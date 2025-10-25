import multiprocessing as mp
import pickle
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Literal

from ..fasta import parse_fasta
from ..proforma.annotation import ProFormaAnnotation
from ..sequence import (
    build_mods,  # Import here to avoid circular imports
    digest,
    mass,
    serialize,
)
from ..sequence.decoy import debruijin_sequence, reverse_sequence, shuffle_sequence
from ..sequence.mod_builder import MOD_BUILDER_INPUT_TYPE, strip_mods

DECOY_STRING = "DECOY_"


@dataclass(frozen=True)
class Protein:
    """Protein with sequence and optional header."""

    sequence: str
    header: str | None = None
    decoy: bool = False

    def __hash__(self):
        return hash(self.header) if self.header is not None else hash(self.sequence)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def accession(self) -> str:
        acc = ""
        if self.header:
            acc = self.header.split()[0].replace(">", "")
            if self.decoy:
                acc = f"{DECOY_STRING}{acc}"
        return acc


@dataclass(frozen=True)
class Peptide:
    """Peptide with reference to its parent protein."""

    sequence: str
    protein_idx: int
    start_pos: int
    end_pos: int

    @property
    def mass(self) -> float:
        """Calculate the monoisotopic mass of the peptide."""
        return mass(self.sequence, charge=0, precision=5, ion_type="p")


class UnifiedDatabase:
    """
    Unified database for proteins and peptides with O(1) lookups.
    - Query by mass -> get peptides (O(1) with binning)
    - Query by peptide sequence -> get source proteins (O(1) with dict lookup)
    """

    def __init__(self, bin_precision: int = 2, verbose: bool = False):
        self.bin_precision = bin_precision
        self.verbose = verbose

        # Protein storage
        self._proteins: list[Protein] = []
        self._sequence_set: set[str] = set()

        # Peptide storage - O(1) lookups
        self._peptide_to_mass: dict[str, float] = {}  # peptide seq -> mass
        self._peptide_to_proteins: dict[
            str, list[Peptide]
        ] = {}  # peptide seq -> list of Peptide objects
        self.mass_bins: dict[int, list[str]] = {}  # mass bin -> list of peptide seqs

        # NEW: Track peptides per protein
        self._protein_to_peptides: dict[
            int, set[str]
        ] = {}  # protein_idx -> set of peptide sequences

    def load_fasta(
        self,
        fasta_path: str,
        decoy: bool = False,
        decoy_strategy: Literal["reverse", "shuffle", "debruijn"] = "reverse",
        static_residues: str = "KR",
        debruijn_kmer: int = 3,
        random_seed: int | None = None,
        method: Literal["sequential", "thread", "process"] = "sequential",
    ) -> None:
        """Load proteins from a FASTA file."""
        t0 = time.time()
        entries = list(parse_fasta(fasta_path))

        # Extract all sequences at once
        sequences = [entry.sequence for entry in entries]

        t1 = time.time()

        if self.verbose:
            print(f"Parsed FASTA in {t1 - t0:.2f} seconds.")

        # Apply decoy strategy to all sequences at once (utilizes multiprocessing)
        if decoy:
            t2 = time.time()
            if decoy_strategy == "reverse":
                sequences = reverse_sequence(sequences)
            elif decoy_strategy == "shuffle":
                sequences = shuffle_sequence(sequences)
            elif decoy_strategy == "debruijn":
                sequences = debruijin_sequence(
                    sequences,
                    k=debruijn_kmer,
                    seed=random_seed,
                    static_residues=static_residues,
                )
            else:
                raise ValueError(
                    "decoy_strategy must be 'reverse', 'shuffle', or 'debruijn'"
                )

            t3 = time.time()

            if self.verbose:
                print(
                    f"Applied {decoy_strategy} decoy strategy to {len(sequences)} sequences in {t3 - t2:.2f} seconds."
                )

        # Create Protein objects
        proteins: list[Protein] = []
        for i, entry in enumerate(entries):
            protein = Protein(sequence=sequences[i], header=entry.header, decoy=decoy)
            proteins.append(protein)

        t4 = time.time()

        if self.verbose:
            print(f"Loading {len(proteins)} proteins from {fasta_path}...")

        self.add_proteins(proteins, method=method)

    def add_proteins(
        self,
        proteins: Iterable[Protein | str | ProFormaAnnotation],
        method: Literal["sequential", "thread", "process"] = "sequential",
    ) -> None:
        """Add multiple proteins."""
        t0 = time.time()

        # Convert to Protein objects if needed
        protein_list: list[Protein] = []
        for p in proteins:
            if isinstance(p, Protein):
                protein_list.append(p)
            else:
                protein_list.append(Protein(sequence=str(p), header=None))

        # Serialize sequences
        sequences = [p.sequence for p in protein_list]
        serialized = serialize(
            sequences, include_plus=False, precision=5, method=method
        )

        t1 = time.time()

        if self.verbose:
            print(f"Serialized {len(serialized)} proteins in {t1 - t0:.2f} seconds.")

        # Filter duplicates
        new_proteins: list[Protein] = []
        for i, serialized_seq in enumerate(serialized):
            if serialized_seq not in self._sequence_set:
                new_proteins.append(
                    Protein(sequence=serialized_seq, header=protein_list[i].header)
                )

        t2 = time.time()

        if self.verbose:
            print(
                f"Filtered duplicates in {t2 - t1:.2f} seconds. {len(new_proteins)} new proteins."
            )

        if not new_proteins:
            return

        self._proteins.extend(new_proteins)
        self._sequence_set.update(p.sequence for p in new_proteins)

        t3 = time.time()

        if self.verbose:
            print(f"Total add_proteins time: {t3 - t0:.2f} seconds.")

    def digest_proteins_slow(
        self,
        cleave_on: str,
        restrict_before: str = "",
        restrict_after: str = "",
        cterminal: bool = True,
        missed_cleavages: int = 1,
        semi: bool = False,
        min_len: int | None = 6,
        max_len: int | None = 50,
        # Modification parameters
        nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
        cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
        internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
        labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
        nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        max_variable_mods: int = 2,
        use_regex: bool = False,
        method: Literal["sequential", "thread", "process"] = "sequential",
    ) -> None:
        """Digest all proteins and add resulting peptides (with optional modifications) to the database."""

        t0 = time.time()

        protein_sequences: list[str] = [p.sequence for p in self._proteins]

        # Digest proteins
        peptides_per_protein: list[list[str]] = digest(
            protein_sequences,
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            method=method,
        )

        t1 = time.time()

        if self.verbose:
            total_peptides = sum(len(peps) for peps in peptides_per_protein)
            print(
                f"Digested {len(self._proteins)} proteins into {total_peptides} peptides in {t1 - t0:.2f} seconds."
            )

        # Check if any modifications are specified
        has_mods = any(
            [
                nterm_static != None,
                cterm_static != None,
                internal_static != None,
                labile_static != None,
                nterm_variable != None,
                cterm_variable != None,
                internal_variable != None,
                labile_variable != None,
            ]
        )

        t4 = 0
        if has_mods:
            # Apply modifications to all peptides
            t2 = time.time()

            # Flatten peptides for batch processing
            all_peptides: list[str] = [
                pep for peps in peptides_per_protein for pep in peps
            ]

            if self.verbose:
                print(f"  Applying modifications to {len(all_peptides)} peptides...")

            # Build modified peptides
            modified_peptides_list = build_mods(
                all_peptides,
                nterm_static=nterm_static,
                cterm_static=cterm_static,
                internal_static=internal_static,
                labile_static=labile_static,
                nterm_variable=nterm_variable,
                cterm_variable=cterm_variable,
                internal_variable=internal_variable,
                labile_variable=labile_variable,
                max_variable_mods=max_variable_mods,
                use_regex=use_regex,
                precision=5,
                include_plus=False,
                method=method,
            )

            t3 = time.time()

            if self.verbose:
                total_modified = sum(len(mods) for mods in modified_peptides_list)
                print(
                    f"  Generated {total_modified} modified peptides in {t3 - t2:.2f} seconds."
                )

            # Reconstruct peptides_per_protein with modifications
            # Map original peptides to their modified versions
            peptide_to_mods = {
                pep: mods for pep, mods in zip(all_peptides, modified_peptides_list)
            }

            # Expand each protein's peptides to include all modified versions
            modified_peptides_per_protein: list[list[str]] = []
            for peps in peptides_per_protein:
                protein_peptides: list[str] = []
                for pep in peps:
                    protein_peptides.extend(peptide_to_mods[pep])
                modified_peptides_per_protein.append(protein_peptides)

            peptides_per_protein = modified_peptides_per_protein

            t4 = time.time()

            if self.verbose:
                total_final = sum(len(peps) for peps in peptides_per_protein)
                print(
                    f"  Reorganized {total_final} modified peptides per protein in {t4 - t3:.2f} seconds."
                )

        # Add peptides with protein references
        self._add_digested_peptides(peptides_per_protein, method=method)

        t5 = time.time()

        if self.verbose:
            print(f"Indexed peptides in {t5 - (t4 if has_mods else t1):.2f} seconds.")
            print(f"Total digest time: {t5 - t0:.2f} seconds.")
            self.print_stats()

    def digest_proteins(
        self,
        cleave_on: str,
        restrict_before: str = "",
        restrict_after: str = "",
        cterminal: bool = True,
        missed_cleavages: int = 1,
        semi: bool = False,
        min_len: int | None = 6,
        max_len: int | None = 50,
        # Modification parameters
        nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
        cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
        internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
        labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
        nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
        max_variable_mods: int = 2,
        use_regex: bool = False,
        n_workers: int | None = None,
    ) -> None:
        """
        Fast digest using pre-initialized workers.
        Avoids function serialization overhead.
        """
        from .worker import DigestWorkerPool

        t0 = time.time()

        has_mods = any(
            [
                nterm_static,
                cterm_static,
                internal_static,
                labile_static,
                nterm_variable,
                cterm_variable,
                internal_variable,
                labile_variable,
            ]
        )

        print(f"has_mods: {has_mods}")

        protein_sequences = [p.sequence for p in self._proteins]

        if self.verbose:
            print(
                f"Digesting {len(protein_sequences)} proteins with {n_workers or mp.cpu_count()} workers..."
            )

        t1 = time.time()

        # Use pre-initialized worker pool
        with DigestWorkerPool(
            n_workers=n_workers or mp.cpu_count(),
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            has_mods=has_mods,
            nterm_static=nterm_static,
            cterm_static=cterm_static,
            internal_static=internal_static,
            labile_static=labile_static,
            nterm_variable=nterm_variable,
            cterm_variable=cterm_variable,
            internal_variable=internal_variable,
            labile_variable=labile_variable,
            max_variable_mods=max_variable_mods,
            use_regex=use_regex,
        ) as pool:
            peptides_mass_per_protein = pool.digest_proteins(protein_sequences)

        t2 = time.time()

        if self.verbose:
            total_peptides = sum(len(peps[0]) for peps in peptides_mass_per_protein)
            print(f"Digested into {total_peptides} peptides in {t2 - t1:.2f} seconds")

        # Now add peptides to database
        # Extract just the peptide sequences for indexing
        peptides_per_protein = [peps[0] for peps in peptides_mass_per_protein]
        peptide_masses_per_protein = [masses[1] for masses in peptides_mass_per_protein]

        # Add to database with pre-calculated masses
        self._add_digested_peptides_with_masses(
            peptides_per_protein, peptide_masses_per_protein
        )

        t3 = time.time()

        if self.verbose:
            print(f"Indexed peptides in {t3 - t2:.2f} seconds")
            print(f"Total digest_proteins_fast time: {t3 - t0:.2f} seconds")
            self.print_stats()

    def _add_digested_peptides_with_masses(
        self,
        peptides_per_protein: list[list[str]],
        masses_per_protein: list[list[float]],
    ) -> None:
        """
        Add digested peptides with pre-calculated masses.
        Much faster than recalculating masses!
        """
        t0 = time.time()

        # Build mapping of peptide -> mass (skip duplicates)
        peptide_to_mass_map: dict[str, float] = {}
        for peptides, masses in zip(peptides_per_protein, masses_per_protein):
            for peptide_seq, mass in zip(peptides, masses):
                if peptide_seq not in self._peptide_to_mass:
                    peptide_to_mass_map[peptide_seq] = mass

        t1 = time.time()

        if self.verbose:
            print(
                f"  Found {len(peptide_to_mass_map)} new unique peptides in {t1 - t0:.2f} seconds"
            )

        # Add to mass cache and bins
        for peptide_seq, peptide_mass in peptide_to_mass_map.items():
            self._peptide_to_mass[peptide_seq] = peptide_mass

            # Add to mass bins
            bin_key = int(
                round(peptide_mass, self.bin_precision) * 10**self.bin_precision
            )
            if bin_key not in self.mass_bins:
                self.mass_bins[bin_key] = []
            self.mass_bins[bin_key].append(peptide_seq)

        t2 = time.time()

        if self.verbose:
            print(f"  Added peptides to mass bins in {t2 - t1:.2f} seconds")

        # Build peptide -> protein mappings
        for protein_idx, peptides in enumerate(peptides_per_protein):
            if protein_idx not in self._protein_to_peptides:
                self._protein_to_peptides[protein_idx] = set()

            protein_seq = self._proteins[protein_idx].sequence

            for peptide_seq in peptides:
                self._protein_to_peptides[protein_idx].add(peptide_seq)
                unmod_seq = strip_mods(peptide_seq)

                # Find all occurrences
                start = 0
                while True:
                    pos = protein_seq.find(unmod_seq, start)
                    if pos == -1:
                        break
                    peptide_obj = Peptide(
                        sequence=peptide_seq,
                        protein_idx=protein_idx,
                        start_pos=pos,
                        end_pos=pos + len(peptide_seq),
                    )

                    if peptide_seq not in self._peptide_to_proteins:
                        self._peptide_to_proteins[peptide_seq] = []
                    self._peptide_to_proteins[peptide_seq].append(peptide_obj)

                    start = pos + 1

        t3 = time.time()

        if self.verbose:
            print(f"  Built peptide->protein mappings in {t3 - t2:.2f} seconds")

    def _add_digested_peptides(
        self,
        peptides_per_protein: list[list[str]],
        method: Literal["sequential", "thread", "process"] = "sequential",
    ) -> None:
        """Add digested peptides with protein references."""
        t0 = time.time()

        # Collect all unique peptide sequences
        unique_peptides: set[str] = set()
        for peptides in peptides_per_protein:
            unique_peptides.update(peptides)

        # Filter out peptides we've already processed
        new_peptides = [p for p in unique_peptides if p not in self._peptide_to_mass]

        t1 = time.time()

        if self.verbose:
            print(
                f"  Found {len(new_peptides)} new unique peptides in {t1 - t0:.2f} seconds."
            )

        if new_peptides:
            # Calculate masses for ALL new peptides at once (better multiprocessing performance)
            t2 = time.time()
            masses = mass(
                new_peptides, charge=0, precision=5, ion_type="p", method=method
            )
            t3 = time.time()

            if self.verbose:
                print(
                    f"  Calculated masses for {len(new_peptides)} peptides in {t3 - t2:.2f} seconds."
                )

            # Add to mass cache and bins
            for peptide_seq, peptide_mass in zip(new_peptides, masses):
                self._peptide_to_mass[peptide_seq] = peptide_mass

                # Add to mass bins
                bin_key = int(
                    round(peptide_mass, self.bin_precision) * 10**self.bin_precision
                )
                if bin_key not in self.mass_bins:
                    self.mass_bins[bin_key] = []
                self.mass_bins[bin_key].append(peptide_seq)

            t4 = time.time()

            if self.verbose:
                print(f"  Added peptides to mass bins in {t4 - t3:.2f} seconds.")

        # Build peptide -> protein mappings
        t5 = time.time()
        for protein_idx, peptides in enumerate(peptides_per_protein):
            # NEW: Initialize peptide set for this protein if not exists
            if protein_idx not in self._protein_to_peptides:
                self._protein_to_peptides[protein_idx] = set()

            protein_seq = self._proteins[protein_idx].sequence
            for peptide_seq in peptides:
                # NEW: Add peptide to protein's set
                self._protein_to_peptides[protein_idx].add(peptide_seq)

                # Find all occurrences of this peptide in the protein
                start = 0
                while True:
                    pos = protein_seq.find(peptide_seq, start)
                    if pos == -1:
                        break

                    peptide_obj = Peptide(
                        sequence=peptide_seq,
                        protein_idx=protein_idx,
                        start_pos=pos,
                        end_pos=pos + len(peptide_seq),
                    )

                    if peptide_seq not in self._peptide_to_proteins:
                        self._peptide_to_proteins[peptide_seq] = []
                    self._peptide_to_proteins[peptide_seq].append(peptide_obj)

                    start = pos + 1

        t6 = time.time()

        if self.verbose:
            print(f"  Built peptide->protein mappings in {t6 - t5:.2f} seconds.")

    def query_mass(
        self,
        neutral_mass: float,
        tolerance: float,
        tolerance_type: Literal["ppm", "da"] = "ppm",
    ) -> list[tuple[Peptide, set[Protein]]]:
        """Query peptides by neutral mass. Returns list of tuples (Peptide, set of Proteins). O(1) average case."""
        if tolerance_type == "ppm":
            tol_da = neutral_mass * tolerance / 1e6
        elif tolerance_type == "da":
            tol_da = tolerance
        else:
            raise ValueError("tolerance_type must be 'ppm' or 'da'")

        min_mass = neutral_mass - tol_da
        max_mass = neutral_mass + tol_da
        min_bin = int(round(min_mass, self.bin_precision) * 10**self.bin_precision)
        max_bin = int(round(max_mass, self.bin_precision) * 10**self.bin_precision)

        candidate_peptides: list[str] = []
        for bin_key in range(min_bin, max_bin + 1):
            if bin_key in self.mass_bins:
                candidate_peptides.extend(self.mass_bins[bin_key])

        matching_peptides: list[str] = [
            pep
            for pep in candidate_peptides
            if min_mass <= self._peptide_to_mass[pep] <= max_mass
        ]

        # Build result with Peptide objects and their associated Proteins
        result: list[tuple[Peptide, set[Protein]]] = []
        for peptide_seq in matching_peptides:
            peptide_objects = self._peptide_to_proteins.get(peptide_seq, [])
            for peptide_obj in peptide_objects:
                proteins = {self._proteins[peptide_obj.protein_idx]}
                result.append((peptide_obj, proteins))

        return result

    def query_mz(
        self,
        precursor_mz: float,
        charge: int,
        tolerance: float,
        tolerance_type: Literal["ppm", "da"] = "ppm",
    ) -> list[str]:
        """Query peptides by precursor m/z. Returns peptide sequences."""
        if charge <= 0:
            raise ValueError("Charge must be a positive integer.")

        neutral_mass = precursor_mz * charge - charge * 1.007276466812
        return self.query_mass(neutral_mass, tolerance, tolerance_type)

    def get_protein_peptides(self, protein_idx: int) -> set[str]:
        """
        Get all peptide sequences for a specific protein.
        Returns set of peptide sequences.
        O(1) lookup.
        """
        return self._protein_to_peptides.get(protein_idx, set())

    def get_protein_peptide_count(self, protein_idx: int) -> int:
        """
        Get the number of unique peptides for a specific protein.
        O(1) lookup.
        """
        return len(self._protein_to_peptides.get(protein_idx, set()))

    def get_all_protein_peptide_counts(self) -> dict[int, int]:
        """
        Get peptide counts for all proteins.
        Returns dict mapping protein_idx -> peptide count.
        """
        return {
            idx: len(peptides) for idx, peptides in self._protein_to_peptides.items()
        }

    def get_peptide_proteins(self, peptide_seq: str) -> list[Peptide]:
        """
        Get Peptide objects for a specific peptide sequence.
        Returns list of Peptide objects.
        O(1) lookup.
        """
        return self._peptide_to_proteins.get(peptide_seq, [])

    def get_proteins_for_peptides(
        self, peptide_seqs: Iterable[str]
    ) -> dict[str, list[Protein]]:
        """
        Get proteins for multiple peptides.
        Returns dict mapping peptide_seq -> list of Protein objects.
        O(1) per peptide.
        """
        result: dict[str, list[Protein]] = {}
        for peptide_seq in peptide_seqs:
            peptide_objects = self.get_peptide_proteins(peptide_seq)
            result[peptide_seq] = [
                self._proteins[p.protein_idx] for p in peptide_objects
            ]
        return result

    def get_stats(self) -> dict[str, float | int]:
        """Get database statistics."""
        avg_protein_length = (
            sum(len(p.sequence) for p in self._proteins) / len(self._proteins)
            if self._proteins
            else 0
        )
        avg_peptide_mass = (
            sum(self._peptide_to_mass.values()) / len(self._peptide_to_mass)
            if self._peptide_to_mass
            else 0
        )

        return {
            "num_proteins": len(self._proteins),
            "num_peptides": len(self._peptide_to_mass),
            "num_mass_bins": len(self.mass_bins),
            "avg_protein_length": avg_protein_length,
            "avg_peptide_mass": avg_peptide_mass,
            "bin_precision": self.bin_precision,
        }

    def print_stats(self) -> None:
        """Print database statistics."""
        stats = self.get_stats()
        print("\n=== Unified Database Statistics ===")
        print(f"Proteins: {stats['num_proteins']:,}")
        print(f"Peptides: {stats['num_peptides']:,}")
        print(f"Mass bins: {stats['num_mass_bins']:,}")
        print(f"Avg protein length: {stats['avg_protein_length']:.1f}")
        print(f"Avg peptide mass: {stats['avg_peptide_mass']:.2f} Da")
        print(f"Bin precision: {stats['bin_precision']}")
        print("===================================\n")

    def clear(self) -> None:
        """Clear all data."""
        self._proteins.clear()
        self._sequence_set.clear()
        self._peptide_to_mass.clear()
        self._peptide_to_proteins.clear()
        self.mass_bins.clear()
        self._protein_to_peptides.clear()  # NEW

    def save(self, filepath: str | Path) -> None:
        """Save the database to a file using pickle."""
        t0 = time.time()
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        state: dict[str, Any] = {
            "bin_precision": self.bin_precision,
            "_proteins": self._proteins,
            "_sequence_set": self._sequence_set,
            "_peptide_to_mass": self._peptide_to_mass,
            "_peptide_to_proteins": self._peptide_to_proteins,
            "mass_bins": self.mass_bins,
            "_protein_to_peptides": self._protein_to_peptides,  # NEW
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
    def load(cls, filepath: str | Path, verbose: bool = False) -> "UnifiedDatabase":
        """Load a database from a pickle file."""
        t0 = time.time()
        filepath = Path(filepath)

        if not filepath.exists():
            raise FileNotFoundError(f"Database file not found: {filepath}")

        with open(filepath, "rb") as f:
            state = pickle.load(f)

        db = cls(bin_precision=state["bin_precision"], verbose=verbose)
        db._proteins = state["_proteins"]
        db._sequence_set = state["_sequence_set"]
        db._peptide_to_mass = state["_peptide_to_mass"]
        db._peptide_to_proteins = state["_peptide_to_proteins"]
        db.mass_bins = state["mass_bins"]
        db._protein_to_peptides = state.get(
            "_protein_to_peptides", {}
        )  # NEW: use .get() for backwards compatibility

        t1 = time.time()

        if verbose:
            file_size_mb = filepath.stat().st_size / (1024 * 1024)
            print(
                f"Loaded database from {filepath} ({file_size_mb:.2f} MB) in {t1 - t0:.2f} seconds."
            )
            db.print_stats()

        return db

    def __len__(self) -> int:
        return len(self._proteins)

    def get_protein_peptide_masses(self, protein_idx: int) -> set[int]:
        """
        Get all unique peptide masses (as integers) for a specific protein.
        Masses are rounded to 5 decimals and converted to int (mass * 100000).
        Returns set of integer masses.
        O(1) lookup.
        """
        peptide_seqs = self._protein_to_peptides.get(protein_idx, set())
        return {
            int(round(self._peptide_to_mass[seq], 5) * 100000) for seq in peptide_seqs
        }

    def get_protein_peptide_masses_float(self, protein_idx: int) -> set[float]:
        """
        Get all unique peptide masses (as floats) for a specific protein.
        Masses are rounded to 5 decimals.
        Returns set of float masses.
        O(1) lookup.
        """
        peptide_seqs = self._protein_to_peptides.get(protein_idx, set())
        return {round(self._peptide_to_mass[seq], 5) for seq in peptide_seqs}

    def get_all_protein_peptide_masses(self) -> dict[int, set[int]]:
        """
        Get unique peptide masses for all proteins.
        Masses are rounded to 5 decimals and converted to int (mass * 100000).
        Returns dict mapping protein_idx -> set of integer masses.
        """
        return {
            idx: {
                int(round(self._peptide_to_mass[seq], 5) * 100000) for seq in peptides
            }
            for idx, peptides in self._protein_to_peptides.items()
        }

    def get_all_protein_peptide_masses_float(self) -> dict[int, set[float]]:
        """
        Get unique peptide masses for all proteins.
        Masses are rounded to 5 decimals.
        Returns dict mapping protein_idx -> set of float masses.
        """
        return {
            idx: {round(self._peptide_to_mass[seq], 5) for seq in peptides}
            for idx, peptides in self._protein_to_peptides.items()
        }
