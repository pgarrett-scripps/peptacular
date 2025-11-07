import multiprocessing as mp
import pickle
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Literal

from ..fasta import parse_fasta
from ..proforma.annotation import ProFormaAnnotation
from ..sequence import (
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
    header: str
    protein_idx: int

    def __hash__(self):
        return hash((self.sequence, self.header))

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def accession(self) -> str:
        acc = ""
        if self.header:
            acc = self.header.split()[0].replace(">", "")
        return acc


@dataclass(frozen=True)
class Peptide:
    """Peptide with reference to its parent protein."""

    sequence: str
    protein_idx: int
    start_pos: int
    end_pos: int
    mass: float

    def __hash__(self):
        return hash((self.sequence, self.protein_idx, self.start_pos, self.end_pos))


class UnifiedDatabase:
    """
    Unified database for proteins and peptides with O(1) lookups.
    """

    def __init__(self, bin_precision: int = 2):
        self.bin_precision = bin_precision

        # Protein storage with multiple lookup keys
        self._proteins: list[Protein] = []
        self._sequence_to_idx: dict[str, int] = {}  # sequence -> protein_idx
        self._header_to_idx: dict[str, int] = {}    # header -> protein_idx
        self._accession_to_idx: dict[str, int] = {} # accession -> protein_idx

        # Peptide storage
        self._peptide_to_mass: dict[str, float] = {}
        self._peptide_to_objects: dict[str, list[Peptide]] = {}  # Renamed for clarity
        self._peptide_idx: dict[Peptide, int] = {}  # peptide object -> unique index
        self._idx_to_peptide: dict[int, Peptide] = {}  # index -> peptide object
        self.mass_bins: dict[int, list[str]] = {}

        self._protein_to_peptides: dict[int, set[str]] = {}

    def load_fasta(
        self,
        fasta_path: str,
        decoy: bool = False,
        decoy_strategy: Literal["reverse", "shuffle", "debruijn"] = "reverse",
        static_residues: str = "KR",
        debruijn_kmer: int = 3,
        random_seed: int | None = None,
        method: Literal["sequential", "thread", "process"] = "sequential",
        max_proteins: int | None = None,
    ) -> None:
        """Load proteins from a FASTA file."""
        entries = list(parse_fasta(fasta_path))

        if max_proteins is not None:
            entries = entries[:max_proteins]

        # Extract all sequences at once
        sequences = [entry.sequence for entry in entries]

        # Apply decoy strategy to all sequences at once (utilizes multiprocessing)
        if decoy:
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
        decoy_prefix = DECOY_STRING if decoy else ""

        # Create Protein objects
        proteins: list[Protein] = []
        for i, entry in enumerate(entries):
            protein = Protein(sequence=sequences[i], header=f"{decoy_prefix}{entry.header.lstrip(">")}", protein_idx=len(self._proteins) + i)
            proteins.append(protein)

        self._add_proteins(proteins)

    def _add_proteins(
        self,
        proteins: Iterable[Protein],
    ) -> None:
        """Add multiple proteins."""
        new_proteins: list[Protein] = []
        for p in proteins:
            if p.sequence in self._sequence_to_idx:
                raise ValueError(f"Duplicate protein sequence detected: {p.sequence}")
            if p.header in self._header_to_idx:
                raise ValueError(f"Duplicate protein header detected: {p.header}")
            if p.accession in self._accession_to_idx:
                raise ValueError(f"Duplicate protein accession detected: {p.accession}")
            
            new_proteins.append(p)

        if not new_proteins:
            return

        # Add to storage and lookup dicts
        for p in new_proteins:
            idx = p.protein_idx
            self._proteins.append(p)
            self._sequence_to_idx[p.sequence] = idx
            self._header_to_idx[p.header] = idx
            self._accession_to_idx[p.accession] = idx

    def add_proteins(
        self,
        proteins: Iterable[tuple[str | ProFormaAnnotation, str]],
        method: Literal["sequential", "thread", "process"] = "sequential",
    ) -> None:
        """Add multiple proteins."""

        protein_seqs: list[str | ProFormaAnnotation] = []
        headers: list[str] = []
        for p, h in proteins:
            protein_seqs.append(p)
            headers.append(h)

        serialized = serialize(
            protein_seqs, include_plus=False, precision=5, method=method
        )

        # Convert to Protein objects if needed
        protein_list: list[Protein] = []
        for i, (p, h) in enumerate(zip(serialized, headers)):
            protein_list.append(Protein(sequence=p, header=h, protein_idx=len(self._proteins) + i))

        self._add_proteins(protein_list)


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
        condense_peptidoforms: bool = False,
        use_regex: bool = False,
        n_workers: int | None = None,
        invalid_residues: set[str] | None = None,
    ) -> None:
        """
        Fast digest using pre-initialized workers.
        Avoids function serialization overhead.
        """
        from .worker import DigestWorkerPool

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

        protein_sequences = [p.sequence for p in self._proteins]


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
            condense_peptidoforms=condense_peptidoforms,
            invalid_residues=invalid_residues,
            use_static_notation=True,
        ) as pool:
            peptides_mass_per_protein = pool.digest_proteins(protein_sequences)


        # Now add peptides to database
        # Extract just the peptide sequences for indexing
        peptides_per_protein = [peps[0] for peps in peptides_mass_per_protein]
        peptide_masses_per_protein = [masses[1] for masses in peptides_mass_per_protein]

        # Add to database with pre-calculated masses
        self._add_digested_peptides_with_masses(
            peptides_per_protein, peptide_masses_per_protein
        )


    def _add_digested_peptides_with_masses(
        self,
        peptides_per_protein: list[list[str]],
        masses_per_protein: list[list[float]],
    ) -> None:
        """Add digested peptides with pre-calculated masses."""
        
        # OPTIMIZATION 1: Flatten all peptides and strip once with multiprocessing
        # Build index mapping to reconstruct per-protein structure
        all_peptides: list[str] = []
        peptide_indices: list[tuple[int, int]] = []  # (protein_idx, peptide_idx_in_protein)
        
        for protein_idx, peptides in enumerate(peptides_per_protein):
            for peptide_idx, peptide_seq in enumerate(peptides):
                all_peptides.append(peptide_seq)
                peptide_indices.append((protein_idx, peptide_idx))
        
        # Strip all peptides at once (uses multiprocessing internally)
        all_stripped = strip_mods(all_peptides)
        
        # Reconstruct per-protein stripped peptides
        stripped_peptides_per_protein: list[list[str]] = [
            [None] * len(peptides) for peptides in peptides_per_protein
        ]
        
        for idx, (protein_idx, peptide_idx) in enumerate(peptide_indices):
            stripped_peptides_per_protein[protein_idx][peptide_idx] = all_stripped[idx]
        
        # OPTIMIZATION 2: Build peptide_to_mass_map in a single pass
        peptide_to_mass_map: dict[str, float] = {}
        for peptides, masses in zip(peptides_per_protein, masses_per_protein):
            for peptide_seq, mass in zip(peptides, masses):
                if peptide_seq not in self._peptide_to_mass and peptide_seq not in peptide_to_mass_map:
                    peptide_to_mass_map[peptide_seq] = mass
        
        # OPTIMIZATION 3: Batch update mass bins (pre-calculate bin_keys)
        bin_precision_mult = 10 ** self.bin_precision
        
        for peptide_seq, peptide_mass in peptide_to_mass_map.items():
            self._peptide_to_mass[peptide_seq] = peptide_mass
            bin_key = int(round(peptide_mass, self.bin_precision) * bin_precision_mult)
            
            if bin_key not in self.mass_bins:
                self.mass_bins[bin_key] = []
            self.mass_bins[bin_key].append(peptide_seq)
        
        # OPTIMIZATION 4: Process peptides with pre-stripped sequences
        next_peptide_idx = len(self._peptide_idx)
        
        for protein_idx, (peptides, stripped_peptides) in enumerate(
            zip(peptides_per_protein, stripped_peptides_per_protein)
        ):
            if protein_idx not in self._protein_to_peptides:
                self._protein_to_peptides[protein_idx] = set()
            
            # Cache protein sequence lookup
            protein_seq = self._proteins[protein_idx].sequence
            
            # Process peptides for this protein
            for peptide_seq, unmod_seq in zip(peptides, stripped_peptides):
                self._protein_to_peptides[protein_idx].add(peptide_seq)
                
                # Find all occurrences of this peptide
                unmod_len = len(unmod_seq)
                start = 0
                
                while True:
                    pos = protein_seq.find(unmod_seq, start)
                    if pos == -1:
                        break
                    
                    peptide_obj = Peptide(
                        sequence=peptide_seq,
                        protein_idx=protein_idx,
                        start_pos=pos,
                        end_pos=pos + unmod_len,
                        mass=self._peptide_to_mass[peptide_seq],
                    )
                    
                    # Add to peptide object lookups
                    if peptide_obj not in self._peptide_idx:
                        self._peptide_idx[peptide_obj] = next_peptide_idx
                        self._idx_to_peptide[next_peptide_idx] = peptide_obj
                        next_peptide_idx += 1
                    
                    if peptide_seq not in self._peptide_to_objects:
                        self._peptide_to_objects[peptide_seq] = []
                    self._peptide_to_objects[peptide_seq].append(peptide_obj)
                    
                    start = pos + 1

    def get_protein(
        self, identifier: str | int | Protein
    ) -> Protein | None:
        """
        Get protein by any identifier: sequence, header, accession, index, or Protein object.
        Returns Protein object or None if not found.
        """
        if isinstance(identifier, Protein):
            return identifier
        
        if isinstance(identifier, int):
            return self._proteins[identifier] if 0 <= identifier < len(self._proteins) else None
        
        # Try all string-based lookups
        idx = (
            self._sequence_to_idx.get(identifier) or
            self._header_to_idx.get(identifier) or
            self._accession_to_idx.get(identifier)
        )
        
        return self._proteins[idx] if idx is not None else None

    def get_proteins_by_peptide(
        self, identifier: str | Peptide | int
    ) -> list[Protein]:
        """
        Get proteins by peptide sequence, Peptide object, or peptide index.
        Returns list of Protein objects.
        """
        peptide_objects: list[Peptide] = []
        
        if isinstance(identifier, str):
            # Peptide sequence
            peptide_objects = self._peptide_to_objects.get(identifier, [])
        elif isinstance(identifier, Peptide):
            # Peptide object
            peptide_objects = self._peptide_to_objects.get(identifier.sequence, [])
        elif isinstance(identifier, int):
            # Peptide index
            peptide_obj = self._idx_to_peptide.get(identifier)
            if peptide_obj:
                peptide_objects = [peptide_obj]
        
        # Convert to unique proteins
        protein_indices = {p.protein_idx for p in peptide_objects}
        return [self._proteins[idx] for idx in protein_indices]

    def get_peptides_by_protein(self, identifier: str | int | Protein) -> list[Peptide]:
        """
        Get all Peptide objects in the database.
        """
        protein = self.get_protein(identifier)
        if not protein:
            return []

        protein_idx = protein.protein_idx
        peptide_seqs = self._protein_to_peptides.get(protein_idx, set())

        peptides: list[Peptide] = []
        for pep_seq in peptide_seqs:
            peptides.extend(self._peptide_to_objects.get(pep_seq, []))

        return peptides

    def get_proteins_for_peptide(self, peptide: Peptide) -> list[Protein]:
        """
        Convert a Peptide object to its source Protein objects.
        This is useful after query_mass() or query_mz().
        """
        return self.get_proteins_by_peptide(peptide)

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

    def query_mass(
        self,
        neutral_mass: float,
        tolerance: float,
        tolerance_type: Literal["ppm", "da"] = "ppm",
    ) -> list[Peptide]:
        """Query peptides by neutral mass. Returns list of Peptide objects."""
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

        result: list[Peptide] = []
        for peptide_seq in matching_peptides:
            result.extend(self._peptide_to_objects.get(peptide_seq, []))

        return result

    def query_mz(
        self,
        precursor_mz: float,
        charge: int,
        tolerance: float,
        tolerance_type: Literal["ppm", "da"] = "ppm",
    ) -> list[Peptide]:
        """Query peptides by precursor m/z. Returns list of Peptide objects."""
        if charge <= 0:
            raise ValueError("Charge must be a positive integer.")

        neutral_mass = precursor_mz * charge - charge * 1.007276466812
        return self.query_mass(neutral_mass, tolerance, tolerance_type)

    def get_peptide_proteins(self, peptide_seq: str) -> list[Peptide]:
        """Get Peptide objects for a specific peptide sequence. O(1) lookup."""
        return self._peptide_to_objects.get(peptide_seq, [])

    def get_proteins_for_peptides(
        self, peptide_seqs: Iterable[str]
    ) -> dict[str, list[Protein]]:
        """Get proteins for multiple peptides. O(1) per peptide."""
        result: dict[str, list[Protein]] = {}
        for peptide_seq in peptide_seqs:
            result[peptide_seq] = self.get_proteins_by_peptide(peptide_seq)
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


    def save(self, filepath: str | Path) -> None:
        """Save the database to a file using pickle."""
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        state: dict[str, Any] = {
            "bin_precision": self.bin_precision,
            "_proteins": self._proteins,
            "_sequence_to_idx": self._sequence_to_idx,
            "_header_to_idx": self._header_to_idx,
            "_accession_to_idx": self._accession_to_idx,
            "_peptide_to_mass": self._peptide_to_mass,
            "_peptide_to_objects": self._peptide_to_objects,
            "_peptide_idx": self._peptide_idx,
            "_idx_to_peptide": self._idx_to_peptide,
            "mass_bins": self.mass_bins,
            "_protein_to_peptides": self._protein_to_peptides,
        }

        with open(filepath, "wb") as f:
            pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, filepath: str | Path) -> "UnifiedDatabase":
        """Load a database from a pickle file."""
        filepath = Path(filepath)

        if not filepath.exists():
            raise FileNotFoundError(f"Database file not found: {filepath}")

        with open(filepath, "rb") as f:
            state = pickle.load(f)

        db = cls(bin_precision=state["bin_precision"])
        db._proteins = state["_proteins"]
        db._sequence_to_idx = state.get("_sequence_to_idx", {})
        db._header_to_idx = state.get("_header_to_idx", {})
        db._accession_to_idx = state.get("_accession_to_idx", {})
        db._peptide_to_mass = state["_peptide_to_mass"]
        db._peptide_to_objects = state.get("_peptide_to_objects", state.get("_peptide_to_proteins", {}))
        db._peptide_idx = state.get("_peptide_idx", {})
        db._idx_to_peptide = state.get("_idx_to_peptide", {})
        db.mass_bins = state["mass_bins"]
        db._protein_to_peptides = state.get("_protein_to_peptides", {})

        return db

    def __len__(self) -> int:
        return len(self._proteins)
