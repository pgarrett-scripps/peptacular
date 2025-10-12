from pathlib import Path
import pickle
from typing import Iterable
from dataclasses import dataclass
import time
from ..proforma.annotation import ProFormaAnnotation
from ..sequence import serialize
from ..fasta import parse_fasta


@dataclass(frozen=True)
class Protein:
    """Protein with sequence and optional header."""

    sequence: str
    header: str | None = None

    def __hash__(self):
        return hash(self.header) if self.header is not None else hash(self.sequence)


class ProteinDatabase:
    """
    Protein database using k-mer indexing for fast sequence searching.
    """

    def __init__(self, kmer_size: int = 5, kmer_skip: int = 0, verbose: bool = False):
        self.kmer_size = kmer_size
        self.kmer_skip = kmer_skip  # Number of amino acids to skip between k-mers
        self.verbose = verbose
        self.kmer_db: dict[str, list[tuple[int, int]]] = {}
        self._proteins: list[Protein] = []
        self._sequence_set: set[str] = set()  # Track unique sequences

    def load_fasta(self, fasta_path: str) -> None:
        """Load proteins from a FASTA file."""
        t0 = time.time()
        entries = list(parse_fasta(fasta_path))
        proteins = [
            Protein(sequence=entry.sequence, header=entry.header) for entry in entries
        ]
        t1 = time.time()

        if self.verbose:
            print(f"Parsed FASTA in {t1 - t0:.2f} seconds.")
            print(f"Loading {len(proteins)} proteins from {fasta_path}...")

        self.add_proteins(proteins)

    def add_proteins(
        self, proteins: Iterable[Protein | str | ProFormaAnnotation]
    ) -> None:
        """Add multiple proteins."""
        t0 = time.time()

        # Convert to Protein objects if needed
        protein_list = []
        for p in proteins:
            if isinstance(p, Protein):
                protein_list.append(p)
            else:
                # String or ProFormaAnnotation
                protein_list.append(Protein(sequence=str(p), header=None))

        # Serialize sequences
        sequences = [p.sequence for p in protein_list]
        serialized = serialize(sequences, include_plus=False, precision=5)

        t1 = time.time()

        if self.verbose:
            print(f"Serialized {len(serialized)} proteins in {t1 - t0:.2f} seconds.")

        # Filter duplicates and create new Protein objects with serialized sequences
        t2 = time.time()
        new_proteins = []
        for i, serialized_seq in enumerate(serialized):
            if serialized_seq not in self._sequence_set:
                new_proteins.append(
                    Protein(sequence=serialized_seq, header=protein_list[i].header)
                )

        t3 = time.time()

        if self.verbose:
            print(
                f"Filtered duplicates in {t3 - t2:.2f} seconds. {len(new_proteins)} new proteins."
            )

        if not new_proteins:
            return

        t4 = time.time()
        start_idx = len(self._proteins)
        self._proteins.extend(new_proteins)
        self._sequence_set.update(p.sequence for p in new_proteins)
        t5 = time.time()

        if self.verbose:
            print(f"Updated protein lists in {t5 - t4:.2f} seconds.")

        # Build k-mer index
        self._build_kmer_index(new_proteins, start_idx)

        t6 = time.time()

        if self.verbose:
            print(f"Total k-mer indexing time: {t6 - t5:.2f} seconds.")
            print(f"Total add_proteins time: {t6 - t0:.2f} seconds.")
            self.print_stats()

    def _build_kmer_index(self, new_proteins: list[Protein], start_idx: int) -> None:
        """Build k-mer index for new proteins."""
        t0 = time.time()
        kmer_size = self.kmer_size
        kmer_skip = self.kmer_skip

        for protein_idx, protein in enumerate(new_proteins, start=start_idx):
            sequence = protein.sequence
            for i in range(0, len(sequence) - kmer_size + 1, kmer_skip + 1):
                kmer = sequence[i : i + kmer_size]
                if kmer not in self.kmer_db:
                    self.kmer_db[kmer] = []
                self.kmer_db[kmer].append((protein_idx, i))

        t1 = time.time()

        if self.verbose:
            print(f"  Built k-mer index in {t1 - t0:.2f} seconds.")

    def get_min_search_length(self) -> int:
        """
        Get the minimum peptide length that can use k-mer indexing.
        Peptides shorter than this will require exhaustive search.

        The minimum length accounts for the k-mer size and skip:
        - We need at least kmer_size amino acids to form a k-mer
        - We try (kmer_skip + 1) different offsets to find a match
        - Minimum length = kmer_size + kmer_skip

        Returns:
            Minimum searchable length
        """
        return self.kmer_size + self.kmer_skip

    def search(self, query: str) -> list[Protein]:
        """
        Search for proteins containing the query sequence.
        Handles k-mer skip by trying multiple offset windows.
        Returns list of Protein objects.

        Note: Queries shorter than get_min_search_length() will use exhaustive search.
        """
        if len(query) < self.get_min_search_length():
            # Exhaustive search for short queries
            raise NotImplementedError(
                "Exhaustive search for short queries is not implemented."
            )

        # Try different offset windows to account for k-mer skip
        all_candidate_positions: dict[int, list[int]] = {}

        for offset in range(self.kmer_skip + 1):
            if offset >= len(query):
                break

            # Extract first k-mer at this offset
            if offset + self.kmer_size > len(query):
                continue

            first_kmer = query[offset : offset + self.kmer_size]

            if first_kmer in self.kmer_db:
                for protein_idx, pos in self.kmer_db[first_kmer]:
                    if protein_idx not in all_candidate_positions:
                        all_candidate_positions[protein_idx] = []
                    # Adjust position to account for offset
                    all_candidate_positions[protein_idx].append(pos - offset)

        if not all_candidate_positions:
            return []

        # Verify exact match only at candidate positions
        results = []
        for protein_idx, positions in all_candidate_positions.items():
            protein = self._proteins[protein_idx]
            for pos in positions:
                # Check bounds and verify match
                if 0 <= pos <= len(protein.sequence) - len(query):
                    if protein.sequence[pos : pos + len(query)] == query:
                        results.append(protein)

        return results

    def search_multiple(self, queries: Iterable[str]) -> set[Protein]:
        """Search for multiple query sequences and return unique matching proteins."""
        result_set: set[Protein] = set()
        for query in queries:
            matches = self.search(query)
            result_set.update(matches)
        return result_set

    def get_stats(self) -> dict:
        """Get database statistics."""
        total_kmers = sum(len(positions) for positions in self.kmer_db.values())
        avg_protein_length = (
            sum(len(p.sequence) for p in self._proteins) / len(self._proteins)
            if self._proteins
            else 0
        )
        proteins_with_headers = sum(1 for p in self._proteins if p.header is not None)

        return {
            "num_proteins": len(self._proteins),
            "proteins_with_headers": proteins_with_headers,
            "num_unique_kmers": len(self.kmer_db),
            "total_kmer_entries": total_kmers,
            "avg_kmers_per_protein": total_kmers / len(self._proteins)
            if self._proteins
            else 0,
            "avg_protein_length": avg_protein_length,
            "kmer_size": self.kmer_size,
            "kmer_skip": self.kmer_skip,
            "min_search_length": self.get_min_search_length(),
        }

    def print_stats(self) -> None:
        """Print database statistics."""
        stats = self.get_stats()
        print("\n=== Database Statistics ===")
        print(f"Proteins: {stats['num_proteins']:,}")
        print(f"Proteins with headers: {stats['proteins_with_headers']:,}")
        print(f"Unique k-mers: {stats['num_unique_kmers']:,}")
        print(f"Total k-mer entries: {stats['total_kmer_entries']:,}")
        print(f"Avg k-mers per protein: {stats['avg_kmers_per_protein']:.1f}")
        print(f"Avg protein length: {stats['avg_protein_length']:.1f}")
        print(f"K-mer size: {stats['kmer_size']}")
        print(f"K-mer skip: {stats['kmer_skip']}")
        print(
            f"Min search length: {stats['min_search_length']} (shorter queries use exhaustive search)"
        )
        print("===========================\n")

    def clear(self) -> None:
        """Remove all proteins."""
        self.kmer_db.clear()
        self._proteins.clear()
        self._sequence_set.clear()

    def __len__(self) -> int:
        return len(self._proteins)

    def __iter__(self):
        return iter(self._proteins)

    @property
    def proteins(self) -> list[Protein]:
        """Return list of all proteins."""
        return self._proteins.copy()

    @property
    def protein_sequences(self) -> list[str]:
        """Return list of all protein sequences."""
        return [p.sequence for p in self._proteins]

    @property
    def protein_headers(self) -> list[str | None]:
        """Return list of all protein headers."""
        return [p.header for p in self._proteins]

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
            "kmer_size": self.kmer_size,
            "kmer_skip": self.kmer_skip,
            "kmer_db": self.kmer_db,
            "_proteins": self._proteins,
            "_sequence_set": self._sequence_set,
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
    def load(cls, filepath: str | Path, verbose: bool = False) -> "ProteinDatabase":
        """
        Load a database from a pickle file.

        :param filepath: Path to the saved database
        :param verbose: Whether to print verbose output
        :return: Loaded ProteinDatabase instance
        """
        t0 = time.time()
        filepath = Path(filepath)

        if not filepath.exists():
            raise FileNotFoundError(f"Database file not found: {filepath}")

        with open(filepath, "rb") as f:
            state = pickle.load(f)

        # Create new instance
        db = cls(
            kmer_size=state["kmer_size"], kmer_skip=state["kmer_skip"], verbose=verbose
        )
        db.kmer_db = state["kmer_db"]
        db._proteins = state["_proteins"]
        db._sequence_set = state["_sequence_set"]

        t1 = time.time()

        if verbose:
            file_size_mb = filepath.stat().st_size / (1024 * 1024)
            print(
                f"Loaded database from {filepath} ({file_size_mb:.2f} MB) in {t1 - t0:.2f} seconds."
            )
            db.print_stats()

        return db
