from pathlib import Path
import pickle
from typing import Iterable, Literal
from dataclasses import dataclass
import time
from ..proforma.annotation import ProFormaAnnotation
from ..sequence import serialize, digest, mass


@dataclass
class Peptide:
    """Peptide with sequence and optional metadata."""

    sequence: str
    protein_id: str | None = None


class PeptideDatabase:
    def __init__(self, bin_precision: int = 2, verbose: bool = False):
        self.bin_precision = bin_precision
        self.verbose = verbose
        self.db: dict[int, list[str]] = {}
        self._peptide_cache: dict[str, float] = {}  # Cache peptide masses

    def add_peptide(self, peptide: str | ProFormaAnnotation) -> None:
        """Add a single peptide to the database."""
        serialized = serialize([peptide], include_plus=False, precision=5)[0]

        # Skip if already exists
        if serialized in self._peptide_cache:
            return

        # Calculate mass
        neutral_mass = mass([serialized], charge=0, precision=5, ion_type="p")[0]

        # Add to cache
        self._peptide_cache[serialized] = neutral_mass

        # Add to bin
        bin_key: int = int(
            round(neutral_mass, self.bin_precision) * 10**self.bin_precision
        )
        if bin_key not in self.db:
            self.db[bin_key] = []
        self.db[bin_key].append(serialized)

    def add_peptides(self, peptides: Iterable[str | ProFormaAnnotation]) -> None:
        """Add multiple peptides to the database."""
        t0 = time.time()

        serialized_peptides = serialize(list(peptides), include_plus=False, precision=5)
        t1 = time.time()

        if self.verbose:
            print(
                f"Serialized {len(serialized_peptides)} peptides in {t1 - t0:.2f} seconds."
            )

        unique_peptides = [
            p for p in set(serialized_peptides) if p not in self._peptide_cache
        ]

        if not unique_peptides:
            if self.verbose:
                print("No new peptides to add (all duplicates).")
            return

        t2 = time.time()
        neutral_masses: list[float] = mass(
            unique_peptides, charge=0, precision=5, ion_type="p"
        )
        t3 = time.time()

        if self.verbose:
            print(
                f"Calculated masses for {len(unique_peptides)} peptides in {t3 - t2:.2f} seconds."
            )

        for peptide, neutral_mass in zip(unique_peptides, neutral_masses):
            self._peptide_cache[peptide] = neutral_mass
            bin_key: int = int(
                round(neutral_mass, self.bin_precision) * 10**self.bin_precision
            )
            if bin_key not in self.db:
                self.db[bin_key] = []
            self.db[bin_key].append(peptide)

        t4 = time.time()

        if self.verbose:
            print(f"Added peptides to bins in {t4 - t3:.2f} seconds.")
            print(f"Total add_peptides time: {t4 - t0:.2f} seconds.")
            self.print_stats()

    def digest_proteins(
        self,
        proteins: Iterable[str | ProFormaAnnotation],
        cleave_on: str,
        restrict_before: str = "",
        restrict_after: str = "",
        cterminal: bool = True,
        missed_cleavages: int = 1,
        semi: bool = False,
        min_len: int | None = 6,
        max_len: int | None = 50,
    ) -> None:
        """Digest proteins and add resulting peptides to the database."""
        t0 = time.time()

        peptides: list[list[str]] = digest(
            list(proteins),
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )

        t1 = time.time()

        if self.verbose:
            print(f"Digested proteins in {t1 - t0:.2f} seconds.")

        # flatten list of lists
        flat_peptides = [pep for sublist in peptides for pep in sublist]

        if self.verbose:
            print(f"Generated {len(flat_peptides)} peptides from digestion.")

        self.add_peptides(flat_peptides)

    def remove_peptide(self, peptide: str | ProFormaAnnotation) -> bool:
        """Remove a peptide from the database. Returns True if removed, False if not found."""
        serialized = serialize([peptide], include_plus=False, precision=5)[0]

        if serialized not in self._peptide_cache:
            return False

        # Get mass and bin
        neutral_mass = self._peptide_cache[serialized]
        bin_key: int = int(
            round(neutral_mass, self.bin_precision) * 10**self.bin_precision
        )

        # Remove from bin
        if bin_key in self.db and serialized in self.db[bin_key]:
            self.db[bin_key].remove(serialized)
            # Clean up empty bins
            if not self.db[bin_key]:
                del self.db[bin_key]

        # Remove from cache
        del self._peptide_cache[serialized]

        return True

    def remove_peptides(self, peptides: Iterable[str | ProFormaAnnotation]) -> int:
        """Remove multiple peptides. Returns count of peptides removed."""
        count = 0
        for peptide in peptides:
            if self.remove_peptide(peptide):
                count += 1
        return count

    def clear(self) -> None:
        """Remove all peptides from the database."""
        self.db.clear()
        self._peptide_cache.clear()

    @property
    def peptides(self) -> list[str]:
        """Return list of all peptides in the database."""
        return list(self._peptide_cache.keys())

    def query(
        self,
        neutral_mass: float,
        tolerance: float,
        tolerance_type: Literal["ppm", "da"] = "ppm",
    ) -> list[str]:
        """Query peptides by neutral mass."""

        if tolerance_type == "ppm":
            tol_da: float = neutral_mass * tolerance / 1e6
        elif tolerance_type == "da":
            tol_da = tolerance
        else:
            raise ValueError("tolerance_type must be 'ppm' or 'da'")

        min_mass: float = neutral_mass - tol_da
        max_mass: float = neutral_mass + tol_da
        min_bin: int = int(round(min_mass, self.bin_precision) * 10**self.bin_precision)
        max_bin: int = int(round(max_mass, self.bin_precision) * 10**self.bin_precision)

        candidate_peptides: list[str] = []
        for bin_key in range(min_bin, max_bin + 1):
            if bin_key in self.db:
                candidate_peptides.extend(self.db[bin_key])

        final_peptides = [
            pep
            for pep in candidate_peptides
            if min_mass <= self._get_mass(pep) <= max_mass
        ]
        return final_peptides

    def query_mz(
        self,
        precursor_mz: float,
        charge: int,
        tolerance: float,
        tolerance_type: Literal["ppm", "da"] = "ppm",
    ) -> list[str]:
        """Query peptides by precursor m/z with specified charge state."""
        if charge <= 0:
            raise ValueError("Charge must be a positive integer.")

        neutral_mass: float = precursor_mz * charge - charge * 1.007276466812
        return self.query(neutral_mass, tolerance, tolerance_type)

    def _get_mass(self, peptide: str) -> float:
        """Retrieve cached mass for a peptide."""
        if peptide not in self._peptide_cache:
            raise ValueError(f"Peptide {peptide} not found in cache.")
        return self._peptide_cache[peptide]

    def mass_range_query(self, min_mass: float, max_mass: float) -> list[str]:
        """Query peptides within a mass range."""
        min_bin: int = int(round(min_mass, self.bin_precision) * 10**self.bin_precision)
        max_bin: int = int(round(max_mass, self.bin_precision) * 10**self.bin_precision)

        results: list[str] = []
        for bin_key in range(min_bin, max_bin + 1):
            if bin_key in self.db:
                results.extend(self.db[bin_key])

        results = [
            pep for pep in results if min_mass <= self._get_mass(pep) <= max_mass
        ]
        return results

    def get_stats(self) -> dict:
        """Get database statistics."""
        total_bins = len(self.db)
        avg_peptides_per_bin = (
            sum(len(peptides) for peptides in self.db.values()) / total_bins
            if total_bins > 0
            else 0
        )
        avg_peptide_length = (
            sum(len(p) for p in self._peptide_cache.keys()) / len(self._peptide_cache)
            if self._peptide_cache
            else 0
        )
        avg_mass = (
            sum(self._peptide_cache.values()) / len(self._peptide_cache)
            if self._peptide_cache
            else 0
        )

        return {
            "num_peptides": len(self._peptide_cache),
            "num_bins": total_bins,
            "avg_peptides_per_bin": avg_peptides_per_bin,
            "avg_peptide_length": avg_peptide_length,
            "avg_peptide_mass": avg_mass,
            "bin_precision": self.bin_precision,
        }

    def print_stats(self) -> None:
        """Print database statistics."""
        stats = self.get_stats()
        print("\n=== Peptide Database Statistics ===")
        print(f"Peptides: {stats['num_peptides']:,}")
        print(f"Mass bins: {stats['num_bins']:,}")
        print(f"Avg peptides per bin: {stats['avg_peptides_per_bin']:.1f}")
        print(f"Avg peptide length: {stats['avg_peptide_length']:.1f}")
        print(f"Avg peptide mass: {stats['avg_peptide_mass']:.2f} Da")
        print(f"Bin precision: {stats['bin_precision']} decimal places")
        print("====================================\n")

    def __len__(self) -> int:
        """Return the number of unique peptides in the database."""
        return len(self._peptide_cache)

    def __contains__(self, peptide: str) -> bool:
        """Check if a peptide exists in the database."""
        return peptide in self._peptide_cache

    def save(self, filepath: str | Path) -> None:
        """
        Save the database to a file using pickle.

        :param filepath: Path to save the database
        """
        t0 = time.time()
        filepath = Path(filepath)

        # Create parent directories if they don't exist
        filepath.parent.mkdir(parents=True, exist_ok=True)

        # Save database state
        state = {
            "bin_precision": self.bin_precision,
            "db": self.db,
            "_peptide_cache": self._peptide_cache,
        }

        with open(filepath, "wb") as f:
            pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)

        t1 = time.time()

        if self.verbose:
            file_size_mb = filepath.stat().st_size / (1024 * 1024)
            print(
                f"Saved database to {filepath} ({file_size_mb:.2f} MB) in {t1 - t0:.2f} seconds."
            )

    @classmethod
    def load(cls, filepath: str | Path, verbose: bool = False) -> "PeptideDatabase":
        """
        Load a database from a pickle file.

        :param filepath: Path to the saved database
        :param verbose: Whether to print verbose output
        :return: Loaded PeptideDatabase instance
        """
        t0 = time.time()
        filepath = Path(filepath)

        if not filepath.exists():
            raise FileNotFoundError(f"Database file not found: {filepath}")

        with open(filepath, "rb") as f:
            state = pickle.load(f)

        # Create new instance
        db = cls(bin_precision=state["bin_precision"], verbose=verbose)
        db.db = state["db"]
        db._peptide_cache = state["_peptide_cache"]

        t1 = time.time()

        if verbose:
            file_size_mb = filepath.stat().st_size / (1024 * 1024)
            print(
                f"Loaded database from {filepath} ({file_size_mb:.2f} MB) in {t1 - t0:.2f} seconds."
            )
            db.print_stats()

        return db
