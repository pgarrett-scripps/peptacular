from functools import cache, cached_property
from typing import Iterable
from .dclass import AminoAcidInfo
from .data import AMINO_ACID_INFOS
from ..elements import ElementInfo


class AALookup:
    def __init__(self, data: dict[str, AminoAcidInfo]):
        self.one_letter_to_info = data
        self.three_letter_to_info = {
            info.three_letter_code.lower(): info for info in data.values()
        }
        self.name_to_info = {info.name.lower(): info for info in data.values()}

    def _query_one_letter(self, code: str) -> AminoAcidInfo | None:
        return self.one_letter_to_info.get(code.upper())

    def _query_three_letter(self, code: str) -> AminoAcidInfo | None:
        return self.three_letter_to_info.get(code.lower())

    def _query_name(self, name: str) -> AminoAcidInfo | None:
        return self.name_to_info.get(name.lower())

    def one_letter(self, code: str) -> AminoAcidInfo:
        val = self._query_one_letter(code)
        if val is not None:
            return val
        raise KeyError(f"Amino acid with one-letter code '{code}' not found.")

    def three_letter(self, code: str) -> AminoAcidInfo:
        val = self._query_three_letter(code)
        if val is not None:
            return val
        raise KeyError(f"Amino acid with three-letter code '{code}' not found.")

    def name(self, name: str) -> AminoAcidInfo:
        val = self._query_name(name)
        if val is not None:
            return val
        raise KeyError(f"Amino acid with name '{name}' not found.")

    @cache
    def __getitem__(self, key: str) -> AminoAcidInfo:
        info = self._query_one_letter(key)
        if info is not None:
            return info

        info = self._query_three_letter(key)
        if info is not None:
            return info

        info = self._query_name(key)
        if info is not None:
            return info

        raise KeyError(
            f"Amino acid '{key}' not found by one-letter code, three-letter code, or name."
        )

    def __contains__(self, key: str) -> bool:
        try:
            _ = self[key]
            return True
        except KeyError:
            return False

    @cached_property
    def ordered_amino_acids(self) -> tuple[AminoAcidInfo, ...]:
        """Get amino acids in order of one-letter codes A-Z"""
        return tuple(
            self.one_letter_to_info[aa] for aa in sorted(self.one_letter_to_info.keys())
        )

    @cached_property
    def ambiguous_amino_acids(self) -> tuple[AminoAcidInfo, ...]:
        """Get ambiguous amino acids (B, J, X, Z)"""
        return tuple(aa for aa in self.ordered_amino_acids if aa.is_ambiguous)

    @cached_property
    def mass_amino_acids(self) -> tuple[AminoAcidInfo, ...]:
        """Get amino acids that have defined masses"""
        return tuple(
            aa
            for aa in self.ordered_amino_acids
            if aa.monoisotopic_mass is not None and aa.average_mass is not None
        )

    @cached_property
    def unambiguous_amino_acids(self) -> tuple[AminoAcidInfo, ...]:
        """Get unambiguous amino acids (all except B, J, X, Z)"""
        return tuple(aa for aa in self.ordered_amino_acids if not aa.is_ambiguous)

    @cached_property
    def mass_unambiguous_amino_acids(self) -> tuple[AminoAcidInfo, ...]:
        """Get unambiguous amino acids that have defined masses"""
        return tuple(
            aa
            for aa in self.unambiguous_amino_acids
            if aa.monoisotopic_mass is not None and aa.average_mass is not None
        )

    @cache
    def is_ambiguous(self, key: str) -> bool:
        """Check if the amino acid identified by key is ambiguous"""
        aa_info = self[key]
        return aa_info.is_ambiguous

    @cache
    def is_mass_ambiguous(self, key: str) -> bool:
        """Check if the amino acid identified by key has mass ambiguity"""
        aa_info = self[key]
        return aa_info.is_mass_ambiguous

    @cache
    def is_unambiguous(self, key: str) -> bool:
        """Check if the amino acid identified by key is unambiguous"""
        aa_info = self[key]
        return not aa_info.is_ambiguous

    @cache
    def is_mass_unambiguous(self, key: str) -> bool:
        """Check if the amino acid identified by key has no mass ambiguity"""
        aa_info = self[key]
        return not aa_info.is_mass_ambiguous

    def mass(self, key: str, monoisotopic: bool = True) -> float:
        """Get the mass of the amino acid identified by key"""
        aa_info = self[key]
        if monoisotopic:
            if aa_info.monoisotopic_mass is None:
                raise ValueError(
                    f"Amino acid '{key}' does not have a defined monoisotopic mass."
                )
            return aa_info.monoisotopic_mass
        else:
            if aa_info.average_mass is None:
                raise ValueError(
                    f"Amino acid '{key}' does not have a defined average mass."
                )
            return aa_info.average_mass

    def composition(self, key: str) -> dict[ElementInfo, int]:
        """Get the elemental composition of the amino acid identified by key"""
        aa_info = self[key]
        if aa_info.composition is None:
            raise ValueError(
                f"Amino acid '{key}' does not have a defined elemental composition."
            )
        return aa_info.composition

    def __iter__(self) -> Iterable[AminoAcidInfo]:
        """Iterator over all amino acids in order of one-letter codes A-Z"""
        for aa in self.ordered_amino_acids:
            yield aa


AA_LOOKUP = AALookup(AMINO_ACID_INFOS)

# prime the cache for all amino acids
for aa in AMINO_ACID_INFOS.keys():
    AA_LOOKUP[aa]

# load cahced properties
_ = AA_LOOKUP.ordered_amino_acids
_ = AA_LOOKUP.ambiguous_amino_acids
_ = AA_LOOKUP.mass_amino_acids
_ = AA_LOOKUP.unambiguous_amino_acids
_ = AA_LOOKUP.mass_unambiguous_amino_acids
