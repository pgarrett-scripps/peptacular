"""
ProForma Parser Module
"""

import copy
import itertools
import random
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from typing import List, Dict, Optional, Generator, Any, Tuple, Union
from typing import Counter as CounterType
import warnings
import regex as re


from ..chem.chem_constants import AVERAGE_AA_MASSES, MONOISOTOPIC_AA_MASSES
from ..chem.chem_util import chem_mass
from ..mass_calc import adjust_mass, adjust_mz, mod_mass

from ..chem.chem_calc import (
    _parse_charge_adducts_comp,
    _parse_mod_delta_mass_only,
    apply_isotope_mods_to_composition,
    estimate_comp,
    mod_comp,
)
from .utils import _validate_single_mod_multiplier, parse_static_mods


from .input_convert import (
    fix_list_of_mods,
    fix_dict_of_mods,
    fix_intervals_input,
)

from ..proforma_dataclasses import (
    Mod,
    Interval,
    are_mods_equal,
    are_intervals_equal,
)
from ..constants import (
    AA_COMPOSITIONS,
    AMINO_ACIDS,
    AMBIGUOUS_AMINO_ACIDS,
    FRAGMENT_ION_BASE_CHARGE_ADDUCTS,
    MASS_AMBIGUOUS_AMINO_ACIDS,
    NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS,
    ModType,
)
from ..errors import AmbiguousAminoAcidError, ProFormaFormatError, UnknownAminoAcidError
from ..types import (
    ACCEPTED_INTERVAL_INPUT,
    ACCEPTED_MOD_INPUT,
    INTERVAL_VALUE,
    ChemComposition,
    ModValue,
)
from ..util import (
    _combine_ambiguity_intervals,
    _construct_ambiguity_intervals,
    _get_mass_shift_interval,
)


class ProFormaAnnotation:

    def __init__(
        self,
        sequence: str,
        isotope_mods: Optional[List[Mod]] = None,
        static_mods: Optional[List[Mod]] = None,
        labile_mods: Optional[List[Mod]] = None,
        unknown_mods: Optional[List[Mod]] = None,
        nterm_mods: Optional[List[Mod]] = None,
        cterm_mods: Optional[List[Mod]] = None,
        internal_mods: Optional[Dict[int, List[Mod]]] = None,
        intervals: Optional[List[Interval]] = None,
        charge: Optional[int] = None,
        charge_adducts: Optional[List[Mod]] = None,
    ) -> None:

        self._sequence = sequence
        self._isotope_mods = isotope_mods if isotope_mods is not None else []
        self._static_mods = static_mods if static_mods is not None else []
        self._labile_mods = labile_mods if labile_mods is not None else []
        self._unknown_mods = unknown_mods if unknown_mods is not None else []
        self._nterm_mods = nterm_mods if nterm_mods is not None else []
        self._cterm_mods = cterm_mods if cterm_mods is not None else []
        self._internal_mods = internal_mods if internal_mods is not None else {}
        self._intervals = intervals if intervals is not None else []
        self._charge = charge
        self._charge_adducts = charge_adducts if charge_adducts is not None else []

    def has_sequence(self) -> bool:
        """
        Check if the annotation has a sequence.

        :return: True if the annotation has a sequence, False otherwise
        :rtype: bool
        """
        return len(self.sequence) > 0

    def has_isotope_mods(self) -> bool:
        """
        Check if the annotation has any isotope modifications.

        :return: True if the annotation has isotope modifications, False otherwise
        :rtype: bool
        """
        return len(self.isotope_mods) > 0

    def has_static_mods(self) -> bool:
        """
        Check if the annotation has any static modifications.

        :return: True if the annotation has static modifications, False otherwise
        :rtype: bool
        """
        return len(self.static_mods) > 0

    def has_labile_mods(self) -> bool:
        """
        Check if the annotation has any labile modifications.

        :return: True if the annotation has labile modifications, False otherwise
        :rtype: bool
        """
        return len(self.labile_mods) > 0

    def has_unknown_mods(self) -> bool:
        """
        Check if the annotation has any unknown modifications.

        :return: True if the annotation has unknown modifications, False otherwise
        :rtype: bool
        """
        return len(self.unknown_mods) > 0

    def has_nterm_mods(self) -> bool:
        """
        Check if the annotation has any N-terminal modifications.

        :return: True if the annotation has N-terminal modifications, False otherwise
        :rtype: bool
        """
        return len(self.nterm_mods) > 0

    def has_cterm_mods(self) -> bool:
        """
        Check if the annotation has any C-terminal modifications.

        :return: True if the annotation has C-terminal modifications, False otherwise
        :rtype: bool
        """
        return len(self.cterm_mods) > 0

    def has_internal_mods(self) -> bool:
        """
        Check if the annotation has any internal modifications.

        :return: True if the annotation has internal modifications, False otherwise
        :rtype: bool
        """
        return len(self.internal_mods) > 0

    def has_intervals(self) -> bool:
        """
        Check if the annotation has any intervals.

        :return: True if the annotation has intervals, False otherwise
        :rtype: bool
        """
        return len(self.intervals) > 0

    def has_charge(self) -> bool:
        """
        Check if the annotation has a charge value.

        :return: True if the annotation has a charge, False otherwise
        :rtype: bool
        """
        return self._charge is not None

    def has_charge_adducts(self) -> bool:
        """
        Check if the annotation has any charge adducts.

        :return: True if the annotation has charge adducts, False otherwise
        :rtype: bool
        """
        return len(self.charge_adducts) > 0

    def has_mods(self) -> bool:
        """
        Check if the annotation has any modifications of any type.

        :return: True if the annotation has any modifications, False otherwise
        :rtype: bool
        """
        return any(
            [
                self.has_isotope_mods(),
                self.has_static_mods(),
                self.has_labile_mods(),
                self.has_unknown_mods(),
                self.has_nterm_mods(),
                self.has_cterm_mods(),
                self.has_internal_mods(),
                self.has_intervals(),
                self.has_charge(),
                self.has_charge_adducts(),
            ]
        )

    @property
    def sequence(self) -> str:
        """
        Returns the sequence.
        """

        if self._sequence is None:
            return ""

        return self._sequence

    @sequence.setter
    def sequence(self, value: Optional[str]) -> None:
        """
        Sets the sequence.
        """
        self._sequence = value

    @property
    def stripped_sequence(self) -> str:
        """
        Returns the sequence without any modifications.
        """
        return self.sequence

    @property
    def isotope_mods(self) -> List[Mod]:
        """
        Returns the isotope modifications.
        """

        if self._isotope_mods is None:
            return []

        return self._isotope_mods

    @isotope_mods.setter
    def isotope_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the isotope modifications.
        """
        if value is None:
            self._isotope_mods = None
        else:
            value = fix_list_of_mods(value)
            self._isotope_mods = copy.deepcopy(value)

    @property
    def static_mods(self) -> List[Mod]:
        """
        Returns the static modifications.
        """

        if self._static_mods is None:
            return []

        return self._static_mods

    @static_mods.setter
    def static_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the static modifications.
        """
        if value is None:
            self._static_mods = None
        else:
            value = fix_list_of_mods(value)
            self._static_mods = copy.deepcopy(value)

    @property
    def labile_mods(self) -> List[Mod]:
        """
        Returns the labile modifications.
        """

        if self._labile_mods is None:
            return []

        return self._labile_mods

    @labile_mods.setter
    def labile_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the labile modifications.
        """
        if value is None:
            self._labile_mods = None
        else:
            value = fix_list_of_mods(value)
            self._labile_mods = copy.deepcopy(value)

    @property
    def unknown_mods(self) -> List[Mod]:
        """
        Returns the unknown modifications.
        """

        if self._unknown_mods is None:
            return []

        return self._unknown_mods

    @unknown_mods.setter
    def unknown_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the unknown modifications.
        """
        if value is None:
            self._unknown_mods = None
        else:
            value = fix_list_of_mods(value)
            self._unknown_mods = copy.deepcopy(value)

    @property
    def nterm_mods(self) -> List[Mod]:
        """
        Returns the N-terminal modifications.
        """

        if self._nterm_mods is None:
            return []

        return self._nterm_mods

    @nterm_mods.setter
    def nterm_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the N-terminal modifications.
        """
        if value is None:
            self._nterm_mods = None
        else:
            value = fix_list_of_mods(value)
            self._nterm_mods = copy.deepcopy(value)

    @property
    def cterm_mods(self) -> List[Mod]:
        """
        Returns the C-terminal modifications.
        """

        if self._cterm_mods is None:
            return []

        return self._cterm_mods

    @cterm_mods.setter
    def cterm_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the C-terminal modifications.
        """
        if value is None:
            self._cterm_mods = None
        else:
            value = fix_list_of_mods(value)
            self._cterm_mods = copy.deepcopy(value)

    @property
    def internal_mods(self) -> Dict[int, List[Mod]]:
        """
        Returns the internal modifications.
        """

        if self._internal_mods is None:
            return {}

        return self._internal_mods

    @internal_mods.setter
    def internal_mods(
        self, value: Optional[Dict[int, Union[List[ModValue], ModValue]]]
    ) -> None:
        if value is None:
            self._internal_mods = None
        else:
            value = fix_dict_of_mods(value)
            self._internal_mods = copy.deepcopy(value)

    @property
    def intervals(self) -> List[Interval]:
        """
        Returns the intervals.
        """
        if self._intervals is None:
            return []

        return self._intervals

    @intervals.setter
    def intervals(
        self, value: Optional[Union[List[INTERVAL_VALUE], INTERVAL_VALUE]]
    ) -> None:
        """
        Sets the intervals.
        """
        if value is None:
            self._intervals = None
        else:
            value = fix_intervals_input(value)
            self._intervals = copy.deepcopy(value)

    @property
    def charge(self) -> Optional[int]:
        """
        Returns the charge.
        """
        return self._charge

    @charge.setter
    def charge(self, value: Optional[int]):
        """
        Sets the charge.
        """
        self._charge = value

    @property
    def charge_adducts(self) -> List[Mod]:
        """
        Returns the charge adducts.
        """
        if self._charge_adducts is None:
            return []

        return self._charge_adducts

    @charge_adducts.setter
    def charge_adducts(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        if value is None:
            self._charge_adducts = None
        else:
            value = fix_list_of_mods(value)
            self._charge_adducts = copy.deepcopy(value)

    def __repr__(self) -> str:
        """
        Only shows non None types
        """
        seq = f"ProFormaAnnotation(sequence={self.sequence}"

        if self.has_isotope_mods():
            seq += f", isotope_mods={self.isotope_mods}"
        if self.has_static_mods():
            seq += f", static_mods={self.static_mods}"
        if self.has_labile_mods():
            seq += f", labile_mods={self.labile_mods}"
        if self.has_unknown_mods():
            seq += f", unknown_mods={self.unknown_mods}"
        if self.has_nterm_mods():
            seq += f", nterm_mods={self.nterm_mods}"
        if self.has_cterm_mods():
            seq += f", cterm_mods={self.cterm_mods}"
        if self.has_internal_mods():
            seq += f", internal_mods={self.internal_mods}"
        if self.has_intervals():
            seq += f", intervals={self.intervals}"
        if self.has_charge():
            seq += f", charge={self.charge}"
        if self.has_charge_adducts():
            seq += f", charge_adducts={self.charge_adducts}"
        seq += ")"

        return seq

    def __len__(self) -> int:
        return len(self.sequence)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, ProFormaAnnotation):
            return NotImplemented

        # Direct sequence comparison
        if self._sequence != other._sequence:
            return False

        if not are_mods_equal(self.labile_mods, other.labile_mods):
            return False

        if not are_mods_equal(self.unknown_mods, other.unknown_mods):
            return False

        if not are_mods_equal(self.nterm_mods, other.nterm_mods):
            return False

        if not are_mods_equal(self.cterm_mods, other.cterm_mods):
            return False

        if not are_mods_equal(self.charge_adducts, other.charge_adducts):
            return False

        if not are_mods_equal(self.isotope_mods, other.isotope_mods):
            return False

        if not are_mods_equal(self.static_mods, other.static_mods):
            return False

        # Compare internal modifications (more complex due to being a dictionary)
        if self.has_internal_mods() or other.has_internal_mods():
            self_keys = (
                set(self.internal_mods.keys()) if self.has_internal_mods() else set()
            )
            other_keys = (
                set(other.internal_mods.keys()) if other.has_internal_mods() else set()
            )

            for k in self_keys.union(other_keys):
                if not are_mods_equal(
                    self.get_internal_mods_by_index(k),
                    other.get_internal_mods_by_index(k),
                ):
                    return False

        # Compare intervals (if applicable and need to ensure a deep comparison)
        if not are_intervals_equal(self.intervals, other.intervals):
            return False

        # Compare charges
        if self.charge != other.charge:
            return False

        return True

    def copy(self, deep: bool = True) -> "ProFormaAnnotation":
        """
        Returns a deep copy of the ProFormaAnnotation instance.
        """

        if not deep:
            return ProFormaAnnotation(
                sequence=self.sequence,
                isotope_mods=self.isotope_mods,
                static_mods=self.static_mods,
                labile_mods=self.labile_mods,
                unknown_mods=self.unknown_mods,
                nterm_mods=self.nterm_mods,
                cterm_mods=self.cterm_mods,
                internal_mods=self.internal_mods,
                intervals=self.intervals,
                charge=self.charge,
                charge_adducts=self.charge_adducts,
            )
        return deepcopy(self)

    def clear_empty_mods(self) -> None:
        """
        Sets the mods to None if they are an empty list
        """

        if not self.has_isotope_mods():
            self._isotope_mods = None
        if not self.has_static_mods():
            self._static_mods = None
        if not self.has_labile_mods():
            self._labile_mods = None
        if not self.has_unknown_mods():
            self._unknown_mods = None
        if not self.has_nterm_mods():
            self._nterm_mods = None
        if not self.has_cterm_mods():
            self._cterm_mods = None
        if not self.has_charge():
            self._charge = None
        if not self.has_charge_adducts():
            self._charge_adducts = None

        if self.has_internal_mods():
            keys_to_remove = []
            for k, mods in self.internal_mods.items():
                if mods is not None and len(mods) == 0:
                    keys_to_remove.append(k)

            for k in keys_to_remove:
                self.internal_mods.pop(k)

        if not self.has_internal_mods():
            self._internal_mods = None

        if self.has_intervals():
            intervals_to_remove = []
            for i, interval in enumerate(self.intervals):
                if interval.mods is not None and len(interval.mods) == 0:
                    interval.mods = None
                if interval.ambiguous is False and interval.mods is None:
                    intervals_to_remove.append(i)

            for i in reversed(intervals_to_remove):
                self.intervals.pop(i)

        if not self.has_intervals():
            self._intervals = None

    def dict(self) -> Dict[str, Any]:
        """
        Returns a dictionary representation of the ProFormaAnnotation instance.
        """
        return {
            "sequence": copy.deepcopy(self.sequence),
            "isotope_mods": copy.deepcopy(self.isotope_mods),
            "static_mods": copy.deepcopy(self.static_mods),
            "labile_mods": copy.deepcopy(self.labile_mods),
            "unknown_mods": copy.deepcopy(self.unknown_mods),
            "nterm_mods": copy.deepcopy(self.nterm_mods),
            "cterm_mods": copy.deepcopy(self.cterm_mods),
            "intervals": copy.deepcopy(self.intervals),
            "internal_mods": copy.deepcopy(self.internal_mods),
            "charge": copy.deepcopy(self.charge),
            "charge_adducts": copy.deepcopy(self.charge_adducts),
        }

    def mod_dict(
        self, mods: Optional[Union[str, List[str]]] = None, condense: bool = True
    ) -> Dict[str, Any]:
        """
        Returns a dictionary representation of the ProFormaAnnotation instance, including only the
        modification-related attributes.
        """
        result = {
            "isotope": self.isotope_mods,
            "static": self.static_mods,
            "labile": self.labile_mods,
            "unknown": self.unknown_mods,
            "nterm": self.nterm_mods,
            "cterm": self.cterm_mods,
            "intervals": self.intervals,
            "charge": self.charge,
            "charge_adducts": self.charge_adducts,
            "internal": self.internal_mods,
        }

        # remove keys not in mods
        if mods is not None:
            if isinstance(mods, str):
                mods = [mods]
            mods = set(mods)
            result = {k: v for k, v in result.items() if k in mods}

        # remove empty lists/dicts/None values if condense is True
        if condense is True:
            for key in list(result.keys()):
                if (
                    result[key] is None
                    or (isinstance(result[key], list) and len(result[key]) == 0)
                    or (isinstance(result[key], dict) and len(result[key]) == 0)
                ):
                    result.pop(key)

        return result

    def add_mod_dict(
        self, mod_dict: Dict, append: bool = False, inplace: bool = True
    ) -> "ProFormaAnnotation":
        """
        Add mods from a dictionary
        """

        if inplace is False:
            new_annotation = deepcopy(self)
            return new_annotation.add_mod_dict(mod_dict, append, inplace=True)

        if "isotope" in mod_dict:
            self.add_isotope_mods(mod_dict["isotope"], append)
        if "static" in mod_dict:
            self.add_static_mods(mod_dict["static"], append)
        if "labile" in mod_dict:
            self.add_labile_mods(mod_dict["labile"], append)
        if "unknown" in mod_dict:
            self.add_unknown_mods(mod_dict["unknown"], append)
        if "nterm" in mod_dict:
            self.add_nterm_mods(mod_dict["nterm"], append)
        if "cterm" in mod_dict:
            self.add_cterm_mods(mod_dict["cterm"], append)
        if "intervals" in mod_dict:
            self.add_intervals(mod_dict["intervals"], append)
        if "charge" in mod_dict:
            self.charge = mod_dict["charge"]
        if "charge_adducts" in mod_dict:
            self.add_charge_adducts(mod_dict["charge_adducts"], append)
        if "internal" in mod_dict:
            self.add_internal_mods(mod_dict["internal"], append)

        internal_mods = {}
        for k, v in mod_dict.items():
            if isinstance(k, int):
                internal_mods[k] = v

        if len(internal_mods) > 0:
            self.add_internal_mods(internal_mods, append)

        return self

    def condense_static_mods(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Condense static mods into internal mods
        """

        if not self.has_static_mods():
            return self

        if inplace is False:
            return self.copy().condense_static_mods(inplace=True)

        static_mod_dict = parse_static_mods(self.pop_static_mods())
        nterm_mod = static_mod_dict.pop("N-Term", None)
        cterm_mod = static_mod_dict.pop("C-Term", None)

        if nterm_mod is not None:
            self.add_nterm_mods(nterm_mod, append=True)

        if cterm_mod is not None:
            self.add_cterm_mods(cterm_mod, append=True)

        # Handle amino acid specific mods
        for aa, mod in static_mod_dict.items():
            indexes = [m.start() for m in re.finditer(aa, self.sequence)]
            for index in indexes:
                self.add_internal_mods({index: mod}, append=True)

        return self

    def contains_sequence_ambiguity(self) -> bool:
        """
        Check if the sequence contains any ambiguity (modifications or intervals).
        """
        return self.has_intervals() or self.has_unknown_mods()

    def contains_residue_ambiguity(self) -> bool:
        """
        Check if the sequence contains any residue ambiguous amino acids.
        """
        return len(self.get_residue_ambiguity_residues()) > 0

    def get_residue_ambiguity_residues(self) -> List[str]:
        """
        Get a list of residue ambiguous amino acids in the sequence.
        """
        return [aa for aa in self.sequence if aa in AMBIGUOUS_AMINO_ACIDS]

    def contains_mass_ambiguity(self) -> bool:
        """
        Check if the sequence contains any mass ambiguous amino acids.
        """
        return len(self.get_mass_ambiguity_residues()) > 0

    def get_mass_ambiguity_residues(self) -> List[str]:
        """
        Check if the sequence contains any mass ambiguous amino acids.
        """
        return [aa for aa in self.sequence if aa in MASS_AMBIGUOUS_AMINO_ACIDS]

    def pop_labile_mods(self) -> List[Mod]:
        """
        Pop all labile mods and return them in a list
        """
        m = self.labile_mods
        self.labile_mods = None
        return m

    def pop_unknown_mods(self) -> List[Mod]:
        """
        Pop all unknown mods and return them in a list
        """
        m = self.unknown_mods
        self.unknown_mods = None
        return m

    def pop_nterm_mods(self) -> List[Mod]:
        """
        Pop all nterm mods and return them in a list
        """
        m = self.nterm_mods
        self.nterm_mods = None
        return m

    def pop_cterm_mods(self) -> List[Mod]:
        """
        Pop all cterm mods and return them in a list
        """
        m = self.cterm_mods
        self.cterm_mods = None
        return m

    def pop_internal_mods(self) -> Dict[int, List[Mod]]:
        """
        Pop all internal mods and return them in a dictionary
        """
        m = self.internal_mods
        self.internal_mods = None
        return m

    def pop_intervals(self) -> List[Interval]:
        """
        Pop all intervals and return them in a list
        """
        m = self.intervals
        self.intervals = None
        return m

    def pop_charge(self) -> Optional[int]:
        """
        Pop the charge and return it
        """
        m = self.charge
        self.charge = None
        return m

    def pop_charge_adducts(self) -> List[Mod]:
        """
        Pop all charge adducts and return them in a list
        """
        m = self.charge_adducts
        self.charge_adducts = None
        return m

    def pop_isotope_mods(self) -> List[Mod]:
        """
        Pop all isotope mods and return them in a list
        """
        m = self.isotope_mods
        self.isotope_mods = None
        return m

    def pop_static_mods(self) -> List[Mod]:
        """
        Pop all static mods and return them in a list
        """
        m = self.static_mods
        self.static_mods = None
        return m

    def pop_mods(
        self, mods: Optional[Union[str, List[str]]] = None, condense: bool = True
    ) -> Dict[str, Any]:
        """
        Pop all mods and return them in a dictionary
        """

        if mods is None:
            mod_types_to_remove = [mod_type.value for mod_type in ModType]
        elif isinstance(mods, str):
            # Single modification type
            mod_types_to_remove = [mods]
        elif isinstance(mods, list):
            # List of modification types
            mod_types_to_remove = mods
        else:
            raise ValueError(
                f"mods parameter must be str, list of str, or None, got {type(mods)}"
            )

        d = {}
        if ModType.ISOTOPE.value in mod_types_to_remove:
            d["isotope"] = self.pop_isotope_mods()
        else:
            d["isotope"] = []

        if ModType.STATIC.value in mod_types_to_remove:
            d["static"] = self.pop_static_mods()
        else:
            d["static"] = []

        if ModType.LABILE.value in mod_types_to_remove:
            d["labile"] = self.pop_labile_mods()
        else:
            d["labile"] = []

        if ModType.UNKNOWN.value in mod_types_to_remove:
            d["unknown"] = self.pop_unknown_mods()
        else:
            d["unknown"] = []

        if ModType.NTERM.value in mod_types_to_remove:
            d["nterm"] = self.pop_nterm_mods()
        else:
            d["nterm"] = []

        if ModType.CTERM.value in mod_types_to_remove:
            d["cterm"] = self.pop_cterm_mods()
        else:
            d["cterm"] = []

        if ModType.CHARGE_ADDUCTS.value in mod_types_to_remove:
            d["charge_adducts"] = self.pop_charge_adducts()
        else:
            d["charge_adducts"] = []

        if ModType.CHARGE.value in mod_types_to_remove:
            d["charge"] = self.pop_charge()
        else:
            d["charge"] = None

        if ModType.INTERNAL.value in mod_types_to_remove:
            d["internal"] = self.pop_internal_mods()
        else:
            d["internal"] = {}

        if ModType.INTERVAL.value in mod_types_to_remove:
            d["intervals"] = self.pop_intervals()
        else:
            d["intervals"] = []

        if condense is True:
            for key in list(d.keys()):
                if (
                    d[key] is None
                    or (isinstance(d[key], list) and len(d[key]) == 0)
                    or (isinstance(d[key], dict) and len(d[key]) == 0)
                ):
                    d.pop(key)

        return d

    def remove_labile_mods(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Remove all labile mods and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with labile mods removed
        :rtype: ProFormaAnnotation
        """
        annotation = self

        if inplace is False:
            annotation = deepcopy(self)

        _ = annotation.pop_labile_mods()

        return annotation

    def remove_unknown_mods(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Remove all unknown mods and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with unknown mods removed
        :rtype: ProFormaAnnotation
        """
        annotation = self

        if inplace is False:
            annotation = deepcopy(self)

        _ = annotation.pop_unknown_mods()

        return annotation

    def remove_nterm_mods(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Remove all N-terminal mods and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with N-terminal mods removed
        :rtype: ProFormaAnnotation
        """
        annotation = self

        if inplace is False:
            annotation = deepcopy(self)

        _ = annotation.pop_nterm_mods()

        return annotation

    def remove_cterm_mods(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Remove all C-terminal mods and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with C-terminal mods removed
        :rtype: ProFormaAnnotation
        """
        annotation = self

        if inplace is False:
            annotation = deepcopy(self)

        _ = annotation.pop_cterm_mods()

        return annotation

    def remove_internal_mods(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Remove all internal mods and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with internal mods removed
        :rtype: ProFormaAnnotation
        """
        annotation = self

        if inplace is False:
            annotation = deepcopy(self)

        _ = annotation.pop_internal_mods()

        return annotation

    def remove_intervals(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Remove all intervals and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with intervals removed
        :rtype: ProFormaAnnotation
        """
        annotation = self

        if inplace is False:
            annotation = deepcopy(self)

        _ = annotation.pop_intervals()

        return annotation

    def remove_charge(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Remove the charge and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with charge removed
        :rtype: ProFormaAnnotation
        """
        annotation = self

        if inplace is False:
            annotation = deepcopy(self)

        _ = annotation.pop_charge()

        return annotation

    def remove_charge_adducts(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Remove all charge adducts and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with charge adducts removed
        :rtype: ProFormaAnnotation
        """
        annotation = self

        if inplace is False:
            annotation = deepcopy(self)

        _ = annotation.pop_charge_adducts()

        return annotation

    def remove_isotope_mods(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Remove all isotope mods and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with isotope mods removed
        :rtype: ProFormaAnnotation
        """
        annotation = self

        if inplace is False:
            annotation = deepcopy(self)

        _ = annotation.pop_isotope_mods()

        return annotation

    def remove_static_mods(self, inplace: bool = True) -> "ProFormaAnnotation":
        """
        Remove all static mods and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with static mods removed
        :rtype: ProFormaAnnotation
        """
        annotation = self

        if inplace is False:
            annotation = deepcopy(self)

        _ = annotation.pop_static_mods()

        return annotation

    def remove_mods(
        self,
        mods: Optional[Union[str, List[str]]] = None,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Remove all modifications and return the modified annotation.

        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The annotation with all modifications removed
        :rtype: ProFormaAnnotation
        """

        if inplace is False:
            return self.copy().remove_mods(inplace=True, mods=mods)

        _ = self.pop_mods(mods, condense=False)
        return self

    def filter_mods(
        self,
        mods: Optional[Union[str, List[str]]] = None,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Filter the modifications in the annotation based on the provided mod types.
        If inplace is False, a new ProFormaAnnotation object is returned with the filtered mods.
        If inplace is True, the current object is modified.
        """
        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().filter_mods(mods, inplace=True)

        if mods is None:
            # If no mods specified, remove all mods
            return self.remove_mods(inplace=True)

        if isinstance(mods, str):
            mods = [mods]

        mod_dict = self.pop_mods(condense=False)
        filtered_mods = {k: v for k, v in mod_dict.items() if k in mods}

        # Re-add the filtered mods to the annotation
        self.add_mod_dict(filtered_mods, append=False)

        return self

    def add_labile_mods(
        self,
        mods: Optional[Union[List[Mod], Mod]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Add labile mods to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().add_labile_mods(mods, append, inplace=True)

        if mods is None:
            if not append:
                self.labile_mods = None
            return self

        if isinstance(mods, Mod):
            mods = [mods]

        if not append:
            self.labile_mods = mods
        else:
            if self.has_labile_mods():
                self.labile_mods.extend(copy.deepcopy(mods))
            else:
                self.labile_mods = mods

        return self

    def add_unknown_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Add unknown mods to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().add_unknown_mods(mods, append, inplace=True)

        if mods is None:
            if not append:
                self.unknown_mods = None
            return self

        mods = fix_list_of_mods(mods)

        if not append:
            self.unknown_mods = mods
        else:
            if self.has_unknown_mods():
                self.unknown_mods.extend(copy.deepcopy(mods))
            else:
                self.unknown_mods = mods

        return self

    def add_nterm_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Add nterm mods to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().add_nterm_mods(mods, append, inplace=True)

        if mods is None:
            if not append:
                self.nterm_mods = None
            return self

        mods = fix_list_of_mods(mods)

        if not append:
            self.nterm_mods = mods
        else:
            if self.has_nterm_mods():
                self.nterm_mods.extend(copy.deepcopy(mods))
            else:
                self.nterm_mods = mods

        return self

    def add_cterm_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Add cterm mods to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().add_cterm_mods(mods, append, inplace=True)

        if mods is None:
            if not append:
                self.cterm_mods = None
            return self

        mods = fix_list_of_mods(mods)

        if not append:
            self.cterm_mods = mods  # Uses the setter to ensure proper copying
        else:
            if self.has_cterm_mods():
                self.cterm_mods.extend(copy.deepcopy(mods))
            else:
                self.cterm_mods = mods

        return self

    def count_internal_mods(self) -> int:
        """
        Count the number of internal mods in the annotation
        """
        if not self.has_internal_mods():
            return 0
        return sum(len(v) for v in self.internal_mods.values())

    def count_modified_residues(self) -> int:
        """
        Count the number of modified residues in the annotation
        """
        if not self.has_internal_mods():
            return 0
        return len(self.internal_mods)

    def has_internal_mods_at_index(self, index: int) -> bool:
        """
        Check if there are internal mods at a specific index
        """
        if not self.has_internal_mods():
            return False
        return index in self.internal_mods

    def get_internal_mods_by_index(self, index: int) -> Union[List[Mod], None]:
        """
        Get internal mods at a specific index
        """
        if not self.has_internal_mods():
            return None

        if index not in self.internal_mods:
            return None

        return self.internal_mods[index]

    def pop_internal_mod(self, index: int, default: Optional[List[Mod]]) -> List[Mod]:
        """
        Pop internal mods at a specific index
        """
        if not self.has_internal_mods():

            if default is not None:
                return default

            raise KeyError(
                f"No internal mods found in the annotation to pop at index {index}."
            )

        if index not in self.internal_mods:

            if default is not None:
                return default

            raise KeyError(
                f"No internal mods found in the annotation to pop at index {index}."
            )

        return self.internal_mods.pop(index)

    def add_internal_mod(
        self,
        index: int,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Add internal mods to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.add_internal_mod(index, mods, append, inplace=True)

        if mods is None:
            if not append:
                if not self.has_internal_mods():
                    self._internal_mods = {}
                self._internal_mods.pop(index, None)
            return self

        mods = fix_list_of_mods(mods)

        if not self.has_internal_mods():
            self._internal_mods = {}

        if not append:
            self._internal_mods[index] = copy.deepcopy(mods)
        else:
            if index in self.internal_mods:
                self._internal_mods[index].extend(copy.deepcopy(mods))
            else:
                self._internal_mods[index] = copy.deepcopy(mods)

        return self

    def add_internal_mods(
        self,
        mods: Optional[Dict[int, Union[List[ModValue], ModValue]]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Add internal mods to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.add_internal_mods(mods, append, inplace=True)

        if mods is None:
            if not append:
                self.internal_mods = None
            return self

        mods = fix_dict_of_mods(mods)

        if not append:
            self.internal_mods = mods
        else:
            if not self.has_internal_mods():
                self.internal_mods = copy.deepcopy(mods)
            else:
                for k, v in mods.items():
                    if k in self.internal_mods:
                        self.internal_mods[k].extend(copy.deepcopy(v))
                    else:
                        self.internal_mods[k] = copy.deepcopy(v)

        return self

    def add_intervals(
        self,
        intervals: Optional[Union[List[INTERVAL_VALUE], INTERVAL_VALUE]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Add intervals to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.add_intervals(intervals, append, inplace=True)

        if intervals is None:
            if append is False:
                self.intervals = None
            return self

        intervals = fix_intervals_input(intervals)

        if append is False:
            self.intervals = intervals  # Uses the setter to ensure proper copying
        else:
            if self.has_intervals():
                self.intervals.extend(copy.deepcopy(intervals))
                self.fix_intervals()
            else:
                self.intervals = intervals

        return self

    def fix_intervals(self) -> None:
        # loop over all intervals and ensure there are not multiple intervals with the same start and end
        # if there are merge them together. taking adding list of mods together and also taking true if any matching
        # intervals have ambiguous=true

        if not self.has_intervals():
            return

        intervals = {}
        for interval in self.intervals:
            key = (interval.start, interval.end)
            if key not in intervals:
                intervals[key] = copy.deepcopy(interval)
                # if mods is None make list
                if intervals[key].mods is None:
                    intervals[key].mods = []
            else:
                intervals[key].mods.extend(copy.deepcopy(interval.mods))
                intervals[key].ambiguous = (
                    intervals[key].ambiguous or interval.ambiguous
                )

        # sort intervals and add back
        intervals = sorted(intervals.values(), key=lambda x: (x.start, x.end))
        self.intervals = intervals

    def add_charge(
        self, charge: Optional[int], inplace: bool = True
    ) -> "ProFormaAnnotation":
        """
        Add charge to the annotation
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.add_charge(charge, inplace=True)

        self.charge = charge
        return self

    def add_charge_adducts(
        self,
        charge_adducts: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Add charge adducts to the annotation
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.add_charge_adducts(charge_adducts, append, inplace=True)

        if charge_adducts is None:
            if not append:
                self.charge_adducts = None
            return self

        charge_adducts = fix_list_of_mods(charge_adducts)

        if not append:
            self.charge_adducts = charge_adducts
        else:
            if self.has_charge_adducts():
                self.charge_adducts.extend(copy.deepcopy(charge_adducts))
            else:
                self.charge_adducts = charge_adducts

        return self

    def add_isotope_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Add isotope mods to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.add_isotope_mods(mods, append, inplace=True)

        if mods is None:
            if not append:
                self.isotope_mods = None
            return self

        mods = fix_list_of_mods(mods)

        if not append:
            self.isotope_mods = mods  # Uses the setter to ensure proper copying
        else:
            if self.has_isotope_mods():
                self.isotope_mods.extend(copy.deepcopy(mods))
            else:
                self.isotope_mods = mods

        return self

    def add_static_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotation":
        """
        Add static mods to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.add_static_mods(mods, append, inplace=True)

        if mods is None:
            if not append:
                self.static_mods = None
            return self

        mods = fix_list_of_mods(mods)

        if not append:
            self.static_mods = mods  # Uses the setter to ensure proper copying
        else:
            if self.has_static_mods():
                self.static_mods.extend(copy.deepcopy(mods))
            else:
                self.static_mods = mods

        return self

    def strip(self, inplace: bool = False) -> "ProFormaAnnotation":
        """
        Remove all modifications from the annotation and return a new annotation with the stripped sequence.
        """

        if inplace is False:
            return self.copy().strip(inplace=True)

        self.isotope_mods = None
        self.static_mods = None
        self.labile_mods = None
        self.unknown_mods = None
        self.nterm_mods = None
        self.cterm_mods = None
        self.internal_mods = None
        self.intervals = None
        self.charge = None
        self.charge_adducts = None

        return self

    def slice(
        self, start: Optional[int], stop: Optional[int], inplace: bool = False
    ) -> "ProFormaAnnotation":
        """
        Slice the annotation sequence and return a new annotation with the sliced sequence and modifications.

        :param start: Start index for slicing (inclusive). If None, defaults to 0.
        :type start: Optional[int]
        :param stop: Stop index for slicing (exclusive). If None, defaults to sequence length.
        :type stop: Optional[int]
        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :return: The sliced annotation
        :rtype: ProFormaAnnotation
        """
        if inplace is False:
            annotation = deepcopy(self)
            return annotation.slice(start, stop, inplace=True)

        # Handle default values and negative indices
        seq_len = len(self.sequence)
        if start is None:
            start = 0
        elif start < 0:
            start = max(0, seq_len + start)
        else:
            start = min(start, seq_len)

        if stop is None:
            stop = seq_len
        elif stop < 0:
            stop = max(0, seq_len + stop)
        else:
            stop = min(stop, seq_len)

        # Ensure start <= stop
        if start > stop:
            start, stop = stop, start

        new_sequence = self.sequence[start:stop]

        # Early return if no modifications exist
        if not self.has_mods():
            self.sequence = new_sequence
            return self

        # Adjust internal modifications
        new_internal_mods = {}
        if self.has_internal_mods():
            for pos, mods in self.internal_mods.items():
                if start <= pos < stop:
                    new_internal_mods[pos - start] = copy.deepcopy(mods)

        # Adjust intervals
        new_intervals = None
        if self.has_intervals():
            new_intervals = []
            for interval in self.intervals:
                # Check if interval overlaps with slice range
                interval_start = interval.start
                interval_end = interval.end if interval.end is not None else seq_len

                if interval_start < stop and interval_end > start:
                    # Calculate new positions relative to slice
                    new_start = max(0, interval_start - start)
                    new_end = min(stop - start, interval_end - start)

                    # Only add if the interval has meaningful content in the slice
                    if new_start < new_end or (
                        new_start == new_end and interval.ambiguous
                    ):
                        new_intervals.append(
                            Interval(
                                start=new_start,
                                end=new_end if new_end < (stop - start) else None,
                                ambiguous=interval.ambiguous,
                                mods=(
                                    copy.deepcopy(interval.mods)
                                    if interval.mods
                                    else None
                                ),
                            )
                        )

        # Handle terminal modifications
        if start > 0:
            self._nterm_mods = None
        if stop < seq_len:
            self._cterm_mods = None

        # Update annotation
        self._sequence = new_sequence
        self._internal_mods = new_internal_mods if new_internal_mods else None
        self._intervals = new_intervals

        return self

    def sliding_windows(
        self, window_size: int
    ) -> Generator["ProFormaAnnotation", None, None]:
        """
        Generate sliding windows of the annotation with a specified size.

        Similar to slicing but creates overlapping windows that slide across the sequence.

        :param window_size: Size of each window/subsequence
        :type window_size: int
        :return: Generator yielding ProFormaAnnotation objects for each window
        :rtype: Generator["ProFormaAnnotation", None, None]

        .. code-block:: python

            >>> annotation = parse("PEPTIDE")
            >>> windows = list(annotation.sliding_windows(3))
            >>> [w.sequence for w in windows]
            ['PEP', 'TID', 'E']

            >>> # Windows preserve modifications
            >>> annotation = parse("PEP[Phospho]TIDE")
            >>> windows = list(annotation.sliding_windows(3))
            >>> [w.serialize() for w in windows]
            ['PEP[Phospho]', 'TID', 'E']

        """
        if not self.sequence:
            return

        if window_size <= 0:
            raise ValueError("Window size must be positive")

        seq_len = len(self.sequence)
        for start in range(0, seq_len, window_size):
            stop = min(start + window_size, seq_len)
            yield self.slice(start, stop, inplace=False)

    def shift(self, n: int, inplace: bool = False) -> "ProFormaAnnotation":
        """
        Shift the annotation by n positions in a cyclic manner.
        """

        if inplace is False:
            annotation = deepcopy(self)
            return annotation.shift(n, inplace=True)

        seq_len = len(self.sequence)
        effective_shift = n % seq_len
        shifted_sequence = (
            self.sequence[effective_shift:] + self.sequence[:effective_shift]
        )

        new_internal_mods = {}
        for mod_index, mods in self.internal_mods.items():
            shifted_index = (mod_index - effective_shift) % seq_len
            new_internal_mods[shifted_index] = copy.deepcopy(mods)

        # Adjust intervals considering the effective shift and sequence length
        new_intervals = []
        for interval in self.intervals:
            new_start = (interval.start - effective_shift) % seq_len
            new_end = (
                (interval.end - effective_shift) % seq_len
                if interval.end is not None
                else None
            )
            # Ensure the start is always less than the end for non-ambiguous intervals
            if new_end is not None and new_start > new_end:
                new_end, new_start = new_start, new_end
            new_intervals.append(
                Interval(
                    start=new_start,
                    end=new_end,
                    ambiguous=interval.ambiguous,
                    mods=copy.deepcopy(interval.mods),
                )
            )

        self._sequence = shifted_sequence
        self._internal_mods = new_internal_mods  # already a copy
        self._intervals = new_intervals  # already a copy
        return self

    def shuffle(
        self, seed: Optional[Any] = None, inplace: bool = False
    ) -> "ProFormaAnnotation":
        """
        Shuffle the annotation sequence and return a new annotation with the shuffled sequence.
        """

        if inplace is False:
            annotation = deepcopy(self)
            return annotation.shuffle(seed, inplace=True)

        if seed is not None:
            random.seed(seed)

        # Convert sequence to a list of characters for shuffling
        sequence_list = list(self.sequence)
        # Track original positions
        original_positions = list(range(len(sequence_list)))
        # Shuffle the sequence list
        combined = list(zip(sequence_list, original_positions))
        random.shuffle(combined)
        shuffled_sequence, shuffled_positions = zip(*combined)

        # Shuffle internal modifications based on new positions
        new_internal_mods = {}
        # Create a mapping from original to new positions
        position_mapping = {
            original: new for new, original in enumerate(shuffled_positions)
        }
        for original_pos, mods in self.internal_mods.items():
            # Map each original position to its new position
            new_pos = position_mapping[original_pos]
            new_internal_mods[new_pos] = copy.deepcopy(mods)

        if len(new_internal_mods) == 0:
            new_internal_mods = None

        new_sequence = "".join(shuffled_sequence)
        self._sequence = new_sequence
        self._internal_mods = new_internal_mods
        return self

    def reverse(
        self, inplace: bool = False, swap_terms: bool = False
    ) -> "ProFormaAnnotation":
        """
        Reverse the annotation sequence and return a new annotation with the reversed sequence.
        """

        if inplace is False:
            annotation = deepcopy(self)
            return annotation.reverse(inplace=True, swap_terms=swap_terms)

        reversed_sequence = self.sequence[::-1]

        # Reverse internal modifications based on new positions
        new_internal_mods = {}
        for original_pos, mods in self.internal_mods.items():
            new_pos = len(self.sequence) - original_pos - 1
            new_internal_mods[new_pos] = copy.deepcopy(mods)

        # reverse intervals
        new_intervals = []
        for interval in self.intervals:
            new_start = len(self.sequence) - interval.start - 1
            new_end = (
                len(self.sequence) - interval.end - 1
                if interval.end is not None
                else None
            )
            # Ensure the start is always less than the end for non-ambiguous intervals
            if new_end is not None and new_start > new_end:
                new_end, new_start = new_start, new_end
            new_intervals.append(
                Interval(
                    start=new_start,
                    end=new_end,
                    ambiguous=interval.ambiguous,
                    mods=copy.deepcopy(interval.mods),
                )
            )

        if swap_terms:
            nterm_mods = self.cterm_mods
            cterm_mods = self.nterm_mods
        else:
            nterm_mods = self.nterm_mods
            cterm_mods = self.cterm_mods

        self._sequence = reversed_sequence  # already a copy
        self._internal_mods = new_internal_mods  # already a copy
        self._nterm_mods = nterm_mods
        self._cterm_mods = cterm_mods
        self._intervals = new_intervals
        return self

    def split(self) -> Generator["ProFormaAnnotation", None, None]:
        """
        Split each amino acid in the sequence into a separate ProFormaAnnotation
        """

        # for labile mods only include on first amino acid
        labile_mods = self.pop_labile_mods()

        for i, _ in enumerate(self.sequence):
            s = self.slice(i, i + 1, inplace=False)
            if i == 0 and labile_mods:
                s.add_labile_mods(labile_mods, append=True, inplace=True)
            yield s

    def count_residues(self) -> CounterType:
        """
        Count the occurrences of each residue in the sequence.
        """
        return Counter([a.serialize() for a in self.split()])

    def sort(self, inplace: bool = False) -> "ProFormaAnnotation":
        """
        Sort the residues in the annotation sequence and return a new annotation with the sorted sequence.
        """

        if inplace is False:
            new_annotation = deepcopy(self)
            return new_annotation.sort(inplace=True)

        # Mapping original positions to their new positions after sorting
        original_to_new_positions = {
            old: new
            for new, old in enumerate(
                sorted(range(len(self.sequence)), key=lambda x: self.sequence[x])
            )
        }

        # Creating new internal mods with adjusted positions
        new_internal_mods = {}
        for pos, mods in self.internal_mods.items():
            new_pos = original_to_new_positions[pos]
            new_internal_mods[new_pos] = copy.deepcopy(mods)

        # Generating sorted sequence
        sorted_sequence = "".join(sorted(self.sequence))

        self.sequence = sorted_sequence
        self.internal_mods = new_internal_mods
        return self

    def serialize(
        self, include_plus: bool = False, precision: Optional[float] = None
    ) -> str:
        """
        Serialize the entire annotation
        """
        return _serialize_annotation(self, include_plus, precision)

    def serialize_start(
        self, include_plus: bool = False, precision: Optional[float] = None
    ) -> str:
        """
        Serialize the start of the annotation
        """
        return _serialize_annotation_start(self, include_plus, precision)

    def serialize_middle(
        self, include_plus: bool = False, precision: Optional[float] = None
    ) -> str:
        """
        Serialize the middle of the annotation
        """
        return _serialize_annotation_middle(self, include_plus, precision)

    def serialize_end(
        self, include_plus: bool = False, precision: Optional[float] = None
    ) -> str:
        """
        Serialize the end of the annotation
        """
        return _serialize_annotation_end(self, include_plus, precision)

    def is_subsequence(
        self, other: "ProFormaAnnotation", ignore_mods: bool = False
    ) -> bool:
        """
        Check if the annotation is a subsequence of another annotation.
        """

        # check if unmodified sequence is a subsequence
        if self.sequence in other.sequence:

            # loop over all starting indexes where the sequence is found
            for start in [
                m.start() for m in re.finditer(self.sequence, other.sequence)
            ]:

                # check if all modifications are also a subsequence
                sliced_annot = other.slice(
                    start, start + len(self.sequence), inplace=False
                )

                if ignore_mods:
                    # if ignore_mods, only check the sequence
                    if sliced_annot.sequence == self.sequence:
                        return True
                else:
                    if sliced_annot == self:
                        return True

        return False

    def find_indices(
        self, other: "ProFormaAnnotation", ignore_mods: bool = False
    ) -> List[int]:
        """
        Find all occurrences of the annotation in another annotation.
        """

        if not isinstance(other, ProFormaAnnotation):
            raise TypeError(f"other must be a ProFormaAnnotation, got {type(other)}")

        if not self.sequence or not other.sequence:
            # if either sequence is empty, return empty list
            return []

        # find all starting indexes where the sequence is found and is a subsequence
        return [
            m.start()
            for m in re.finditer(self.sequence, other.sequence)
            if self.is_subsequence(
                other.slice(m.start(), m.start() + len(self.sequence), inplace=False),
                ignore_mods,
            )
        ]

    def permutations(self, size: Optional[int] = None) -> List["ProFormaAnnotation"]:
        """
        Generate all permutations of the annotation sequence.
        """

        if size is None:
            size = len(self)

        start = self.serialize_start()
        end = self.serialize_end()

        mods = self.pop_mods()
        self._internal_mods = mods.get("internal")

        components = [a.serialize() for a in self.split()]

        return [
            parse(start + "".join(i) + end)
            for i in itertools.permutations(components, size)
        ]

    def product(self, repeat: Optional[int] = None) -> List["ProFormaAnnotation"]:
        """
        Generate the product of the annotation sequence with itself.
        """

        if repeat is None:
            repeat = len(self)

        start = self.serialize_start()
        end = self.serialize_end()

        mods = self.pop_mods()
        self._internal_mods = mods.get("internal")

        components = [a.serialize() for a in self.split()]

        return [
            parse(start + "".join(i) + end)
            for i in itertools.product(components, repeat=repeat)
        ]

    def combinations(self, size: Optional[int] = None) -> List["ProFormaAnnotation"]:
        """
        Generate all combinations of the annotation sequence.
        """

        if size is None:
            size = len(self)

        start = self.serialize_start()
        end = self.serialize_end()

        mods = self.pop_mods()
        self._internal_mods = mods.get("internal")

        components = [a.serialize() for a in self.split()]

        return [
            parse(start + "".join(i) + end)
            for i in itertools.combinations(components, size)
        ]

    def combinations_with_replacement(
        self, size: Optional[int] = None
    ) -> List["ProFormaAnnotation"]:
        """
        Generate all combinations of the annotation sequence with replacement.
        """

        if size is None:
            size = len(self)

        start = self.serialize_start()
        end = self.serialize_end()

        mods = self.pop_mods()
        self._internal_mods = mods.get("internal")

        components = [a.serialize() for a in self.split()]

        return [
            parse(start + "".join(i) + end)
            for i in itertools.combinations_with_replacement(components, size)
        ]

    def annotate_ambiguity(
        self,
        forward_coverage: List[int],
        reverse_coverage: List[int],
        mass_shift: Optional[Any] = None,
        inplace: bool = False,
    ) -> Union["ProFormaAnnotation", None]:
        """
        Generates ambiguity intervals based on the coverage of the sequence.
        See :func:`~peptacular.sequence.sequence_funcs.annotate_ambiguity` for more details.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.annotate_ambiguity(
                forward_coverage, reverse_coverage, mass_shift, inplace=True
            )

        # ensure that annotation does not contain intervals
        if self.has_intervals():
            raise ValueError("Annotation should not contain intervals")

        if len(forward_coverage) != len(reverse_coverage) != len(self):
            raise ValueError(
                f"Coverage length does not match sequence length: {len(forward_coverage)} != {len(reverse_coverage)} != {len(self)}"
            )

        forward_intervals = _construct_ambiguity_intervals(
            forward_coverage, reverse=False
        )
        reverse_intervals = _construct_ambiguity_intervals(
            reverse_coverage, reverse=True
        )
        ambiguity_intervals = _combine_ambiguity_intervals(
            forward_intervals, reverse_intervals
        )

        intervals = [
            Interval(start, end + 1, True, None) for start, end in ambiguity_intervals
        ]

        self.add_intervals(intervals, append=True)

        if mass_shift is not None:
            mass_shift_interval = _get_mass_shift_interval(
                forward_coverage, reverse_coverage
            )
            if mass_shift_interval is None:
                self.add_labile_mods(mass_shift, append=True)
            elif mass_shift_interval[0] == mass_shift_interval[1]:
                # add modification to the sequence
                self.add_internal_mod(mass_shift_interval[0], mass_shift, append=True)
            else:
                mod_interval = Interval(
                    mass_shift_interval[0],
                    mass_shift_interval[1] + 1,
                    False,
                    [Mod(mass_shift, 1)],
                )
                self.add_intervals([mod_interval], append=True)

        return self

    def pop_delta_mass_mods(self, inplace: bool = True) -> float:
        """
        Pop the delta mass modifications from the modifications' dictionary. This leaves only modifications which
        can have their elemental composition calculated.

        :param annotation: The ProForma annotation.
        :type annotation: ProFormaAnnotation

        :return: The delta mass.
        :rtype: float

        .. code-block:: python

            >>> mods = ProFormaAnnotation(sequence='', nterm_mods = [42.0, -20.0])
            >>> mods.pop_delta_mass_mods()
            22.0
            >>> mods
            ProFormaAnnotation(sequence=)

        """
        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.pop_delta_mass_mods(inplace=True)

        delta_mass = 0.0

        # Labile Mods
        if self.has_labile_mods():
            for i in range(len(self.labile_mods) - 1, -1, -1):
                mod = self.labile_mods[i]
                val = _parse_mod_delta_mass_only(mod)
                if val is not None:
                    delta_mass += val
                    self.labile_mods.pop(i)

        # Unknown Mods
        if self.has_unknown_mods():
            for i in range(len(self.unknown_mods) - 1, -1, -1):
                mod = self.unknown_mods[i]
                val = _parse_mod_delta_mass_only(mod)
                if val is not None:
                    delta_mass += val
                    self.unknown_mods.pop(i)

        # NTerm
        if self.has_nterm_mods():
            for i in range(len(self.nterm_mods) - 1, -1, -1):
                mod = self.nterm_mods[i]
                val = _parse_mod_delta_mass_only(mod)
                if val is not None:
                    delta_mass += val
                    self.nterm_mods.pop(i)

        # CTerm
        if self.has_cterm_mods():
            for i in range(len(self.cterm_mods) - 1, -1, -1):
                mod = self.cterm_mods[i]
                val = _parse_mod_delta_mass_only(mod)
                if val is not None:
                    delta_mass += val
                    self.cterm_mods.pop(i)

        # Intervals
        if self.has_intervals():
            for interval in self.intervals:
                if interval.has_mods():
                    for i in range(len(interval.mods) - 1, -1, -1):
                        mod = interval.mods[i]
                        val = _parse_mod_delta_mass_only(mod)
                        if val is not None:
                            delta_mass += val
                            interval.mods.pop(i)

        # Internal mods
        if self.has_internal_mods():
            for k in self.internal_mods:
                for i in range(len(self.internal_mods[k]) - 1, -1, -1):
                    mod = self.internal_mods[k][i]
                    val = _parse_mod_delta_mass_only(mod)
                    if val is not None:
                        delta_mass += val
                        self.internal_mods[k].pop(i)

        self.clear_empty_mods()

        return delta_mass

    def sequence_comp(
        self,
        ion_type: str,
        isotope: int = 0,
        use_isotope_on_mods: bool = False,
    ) -> ChemComposition:

        # If charge is not provided, set it to 0
        charge = 0
        if self.has_charge():
            charge = self.charge

        # if charge_adducts is not provided, set it to None
        charge_adducts = None
        if self.has_charge_adducts():
            charge_adducts = self.charge_adducts[0]

        if charge_adducts is None:
            if ion_type in ("p", "n"):
                charge_adducts = f"{charge}H+"
            else:
                charge_adducts = (
                    f"{charge-1}H+,{FRAGMENT_ION_BASE_CHARGE_ADDUCTS[ion_type]}"
                )

        if ion_type not in ("p", "n"):
            if charge == 0:
                warnings.warn(
                    "Calculating the comp of a fragment ion with charge state 0. Fragment ions should have a "
                    "charge state greater than 0 since the neutral form doesnt exist."
                )

        if "B" in self.sequence:
            raise AmbiguousAminoAcidError(
                "B",
                "Cannot determine the composition of a sequence with an ambiguous amino acid.",
            )

        if "Z" in self.sequence:
            raise AmbiguousAminoAcidError(
                "Z",
                "Cannot determine the composition of a sequence with an ambiguous amino acid.",
            )

        # Get the composition of the base sequence
        sequence_composition = {}
        for aa in self.sequence:
            try:
                aa_comp = AA_COMPOSITIONS[aa]
            except KeyError as err:
                raise UnknownAminoAcidError(aa) from err
            for k, v in aa_comp.items():
                sequence_composition[k] = sequence_composition.get(k, 0) + v

        # Apply the adjustments for the neutral fragment composition based on strictly the ion dissociation points.
        for k, v in NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[ion_type].items():
            sequence_composition[k] = sequence_composition.get(k, 0) + v

        charge_adduct_comp = _parse_charge_adducts_comp(charge_adducts)

        for k, v in charge_adduct_comp.items():
            sequence_composition[k] = sequence_composition.get(k, 0) + v

        mod_composition = {}
        if self.has_unknown_mods():
            for unknown_mod in self.unknown_mods:
                for k, v in mod_comp(unknown_mod).items():
                    mod_composition[k] = mod_composition.get(k, 0) + v

        if self.has_intervals():
            for interval in self.intervals:
                if interval.has_mods():
                    for interval_mod in interval.mods:
                        for k, v in mod_comp(interval_mod).items():
                            mod_composition[k] = mod_composition.get(k, 0) + v

        if self.has_labile_mods() and ion_type == "p":
            for labile_mod in self.labile_mods:
                for k, v in mod_comp(labile_mod).items():
                    mod_composition[k] = mod_composition.get(k, 0) + v

        if self.has_nterm_mods():
            for nterm_mod in self.nterm_mods:
                for k, v in mod_comp(nterm_mod).items():
                    mod_composition[k] = mod_composition.get(k, 0) + v

        if self.has_cterm_mods():
            for cterm_mod in self.cterm_mods:
                for k, v in mod_comp(cterm_mod).items():
                    mod_composition[k] = mod_composition.get(k, 0) + v

        if self.has_internal_mods():
            for _, internal_mods in self.internal_mods.items():
                for internal_mod in internal_mods:
                    for k, v in mod_comp(internal_mod).items():
                        mod_composition[k] = mod_composition.get(k, 0) + v

        if self.has_static_mods():
            static_map = parse_static_mods(self.static_mods)

            n_term_mod = static_map.get("N-Term")
            if n_term_mod is not None:
                for m in n_term_mod:
                    for k, v in mod_comp(m.val).items():
                        mod_composition[k] = mod_composition.get(k, 0) + v

            c_term_mod = static_map.get("C-Term")
            if c_term_mod is not None:
                for m in c_term_mod:
                    for k, v in mod_comp(m.val).items():
                        mod_composition[k] = mod_composition.get(k, 0) + v

            for aa, mod in static_map.items():
                if aa in ["N-Term", "C-Term"]:
                    continue

                aa_count = self.sequence.count(aa)
                for m in mod:
                    for k, v in mod_comp(m.val).items():
                        mod_composition[k] = mod_composition.get(k, 0) + v * aa_count

        mod_composition["n"] = mod_composition.get("n", 0) + isotope

        # Apply isotopic mods
        if self.has_isotope_mods():
            if use_isotope_on_mods:
                sequence_composition = apply_isotope_mods_to_composition(
                    sequence_composition, self.isotope_mods
                )
                mod_composition = apply_isotope_mods_to_composition(
                    mod_composition, self.isotope_mods
                )
            else:
                sequence_composition = apply_isotope_mods_to_composition(
                    sequence_composition, self.isotope_mods
                )

        composition = {}
        for k, v in sequence_composition.items():
            composition[k] = composition.get(k, 0) + v

        for k, v in mod_composition.items():
            composition[k] = composition.get(k, 0) + v

        composition = {k: v for k, v in composition.items() if v != 0}

        return composition

    def mass(
        self,
        ion_type: str = "p",
        monoisotopic: bool = True,
        isotope: int = 0,
        loss: float = 0.0,
        use_isotope_on_mods: bool = False,
        precision: Optional[int] = None,
        inplace: bool = False,
    ) -> float:

        if self.contains_mass_ambiguity():
            raise AmbiguousAminoAcidError(
                aa=",".join(self.get_residue_ambiguity_residues()),
                msg="Cannot determine the mass of a sequence with ambiguous amino acids: {self.sequence}",
            )

        # inplace should not be true
        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.mass(
                ion_type=ion_type,
                monoisotopic=monoisotopic,
                isotope=isotope,
                loss=loss,
                use_isotope_on_mods=use_isotope_on_mods,
                precision=precision,
                inplace=True,
            )
        else:
            warnings.warn(
                "Inplace mass calculation is not recommended. It may lead to unexpected behavior.",
            )

        self.condense_static_mods(inplace=True)

        # more complex mass calculation (based on chem composition)
        if self.has_isotope_mods():
            peptide_composition, delta_mass = self.comp_mass(
                ion_type, isotope, use_isotope_on_mods, True
            )
            return (
                chem_mass(
                    peptide_composition, monoisotopic=monoisotopic, precision=precision
                )
                + delta_mass
            )

        if ion_type not in ("p", "n"):
            if self.charge is None or self.charge == 0:
                warnings.warn(
                    "Calculating the mass of a fragment ion with charge state 0. Fragment ions should have a "
                    "charge state greater than 0 since the neutral mass doesnt exist."
                )

        m = 0.0
        """
        if self.has_static_mods():
            static_map = parse_static_mods(self.static_mods)

            n_term_mod = static_map.get("N-Term")
            if n_term_mod is not None:
                m += sum(mod_mass(m, monoisotopic, precision=None) for m in n_term_mod)

            c_term_mod = static_map.get("C-Term")
            if c_term_mod is not None:
                m += sum(mod_mass(m, monoisotopic, precision=None) for m in c_term_mod)

            for aa, mod in static_map.items():

                if aa in ["N-Term", "C-Term"]:
                    continue

                aa_count = self.sequence.count(aa)
                m += sum(mod_mass(m, monoisotopic, precision=None) for m in mod) * aa_count
        """
        try:
            m += sum(
                MONOISOTOPIC_AA_MASSES[aa] if monoisotopic else AVERAGE_AA_MASSES[aa]
                for aa in self.sequence
            )
        except KeyError as err:
            raise UnknownAminoAcidError(err) from err

        # Apply labile mods
        if self.has_labile_mods() and ion_type == "p":
            for mod in self.labile_mods:
                m += mod_mass(mod)

        # Apply Unknown mods
        if self.has_unknown_mods():
            for mod in self.unknown_mods:
                m += mod_mass(mod)

        # Apply N-term mods
        if self.has_nterm_mods():
            for mod in self.nterm_mods:
                m += mod_mass(mod)

        # Apply intervals
        if self.has_intervals():
            for interval in self.intervals:
                if interval.mods is not None:
                    for mod in interval.mods:
                        m += mod_mass(mod)

        # apply internal mods
        if self.has_internal_mods():
            for _, mods in self.internal_mods.items():
                for mod in mods:
                    m += mod_mass(mod)

        # apply C-term mods
        if self.has_cterm_mods():
            for mod in self.cterm_mods:
                m += mod_mass(mod)

        return adjust_mass(
            base_mass=m,
            charge=self.charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            loss=loss,
            charge_adducts=self.charge_adducts,
            precision=precision,
        )

    def mz(
        self,
        ion_type: str = "p",
        monoisotopic: bool = True,
        isotope: int = 0,
        loss: float = 0.0,
        precision: Optional[int] = None,
    ) -> float:

        m = self.mass(
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            loss=loss,
            precision=precision,
        )

        return adjust_mz(m, self.charge, precision)

    def comp(
        self,
        ion_type: str = "p",
        estimate_delta: bool = False,
        isotope: int = 0,
        use_isotope_on_mods: bool = False,
        inplace: bool = False,
    ) -> ChemComposition:

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.comp(
                ion_type=ion_type,
                estimate_delta=estimate_delta,
                isotope=isotope,
                use_isotope_on_mods=use_isotope_on_mods,
                inplace=True,
            )

        composition, delta_mass = self.comp_mass(
            ion_type=ion_type,
            isotope=isotope,
            use_isotope_on_mods=use_isotope_on_mods,
            inplace=False,
        )

        if delta_mass != 0:

            if estimate_delta is False:
                raise ValueError(
                    f"Non-zero delta mass ({delta_mass}) encountered without estimation enabled for "
                    f"sequence '{self.serialize(include_plus=True, precision=5)}'."
                )

            if use_isotope_on_mods is True:
                delta_mass_comp = estimate_comp(delta_mass, self.isotope_mods)
                warnings.warn(
                    "Applying isotopic modifications to the predicted composition. This will not be accurate."
                )
            else:
                delta_mass_comp = estimate_comp(delta_mass, None)

            # Combine the compositions of the sequence and the delta mass
            for element in delta_mass_comp:
                composition.setdefault(element, 0)
                composition[element] += delta_mass_comp[element]

        return composition

    def comp_mass(
        self,
        ion_type: str = "p",
        isotope: int = 0,
        use_isotope_on_mods: bool = False,
        inplace: bool = False,
    ) -> Tuple[ChemComposition, float]:
        # if inplace is True -> Will condense static mods

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.comp_mass(
                ion_type=ion_type,
                isotope=isotope,
                use_isotope_on_mods=use_isotope_on_mods,
                inplace=True,
            )

        # condenses static mods
        self.condense_static_mods(inplace=True)
        delta_mass = (
            self.pop_delta_mass_mods()
        )  # sum delta mass mods and remove them from the annotation
        peptide_composition = self.sequence_comp(ion_type, isotope, use_isotope_on_mods)

        if use_isotope_on_mods is True and delta_mass != 0:
            warnings.warn(
                "use_isotope_on_mods=True and delta_mass != 0. Cannot apply isotopic modifications to the delta mass."
            )

        return peptide_composition, delta_mass

    def condense_to_mass_mods(
        self, include_plus: bool = False, inplace: bool = True
    ) -> "ProFormaAnnotation":

        if inplace is False:
            # Create a copy of the annotation to modify
            annotation = deepcopy(self)
            return annotation.condense_to_mass_mods(
                include_plus=include_plus, inplace=True
            )

        labile_mods = self.pop_labile_mods()
        isotope_mods = self.isotope_mods

        # fix this so that isotopes get applied
        n_term_mods = self.pop_nterm_mods()
        c_term_mods = self.pop_cterm_mods()
        n_term_mods_mass, c_term_mods_mass = None, None

        if n_term_mods:
            n_term_annot = ProFormaAnnotation(
                sequence="",
                nterm_mods=n_term_mods,
                isotope_mods=isotope_mods,
            )
            n_term_mods_mass = n_term_annot.mass(ion_type="n")

        if c_term_mods:
            c_term_annot = ProFormaAnnotation(
                sequence="",
                cterm_mods=c_term_mods,
                isotope_mods=isotope_mods,
            )
            c_term_mods_mass = c_term_annot.mass(ion_type="n")

        # Split into segments and process each one
        segments = list(self.split())
        stripped_segments = [seg.strip(inplace=False) for seg in segments]

        # Calculate mass differences
        for i, (segment, stripped) in enumerate(zip(segments, stripped_segments)):
            # Calculate mass difference
            mass_diff = segment.mass() - stripped.mass()
            if abs(mass_diff) > 1e-6:  # Only add if difference is significant
                self.add_internal_mod(
                    index=i, mods=mass_diff, append=False, inplace=True
                )

        if labile_mods:
            labile_mods_mass = sum(mod_mass(mod) for mod in labile_mods)
            self.add_labile_mods(mods=labile_mods_mass, append=False, inplace=True)

        if n_term_mods_mass is not None:
            self.add_nterm_mods(mods=n_term_mods_mass, append=False, inplace=True)

        if c_term_mods_mass is not None:
            self.add_cterm_mods(mods=c_term_mods_mass, append=False, inplace=True)

        self.pop_isotope_mods()

        return self


def create_annotation(
    sequence: str,
    isotope_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    static_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    labile_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    unknown_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    nterm_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    cterm_mods: Optional[ACCEPTED_MOD_INPUT] = None,
    internal_mods: Optional[Dict[int, ACCEPTED_MOD_INPUT]] = None,
    intervals: Optional[ACCEPTED_INTERVAL_INPUT] = None,
    charge: Optional[int] = None,
    charge_adducts: Optional[ACCEPTED_MOD_INPUT] = None,
) -> ProFormaAnnotation:
    """
    Create a ProFormaAnnotation from a sequence and modifications

    .. code-block:: python

        >>> create_annotation('PEPTIDE', static_mods=['Carbamidomethyl'])
        ProFormaAnnotation(sequence=PEPTIDE, static_mods=[Mod('Carbamidomethyl', 1)])

    """

    isotope_mods = fix_list_of_mods(isotope_mods) if isotope_mods is not None else None
    static_mods = fix_list_of_mods(static_mods) if static_mods is not None else None
    labile_mods = fix_list_of_mods(labile_mods) if labile_mods is not None else None
    unknown_mods = fix_list_of_mods(unknown_mods) if unknown_mods is not None else None
    nterm_mods = fix_list_of_mods(nterm_mods) if nterm_mods is not None else None
    cterm_mods = fix_list_of_mods(cterm_mods) if cterm_mods is not None else None
    internal_mods = (
        fix_dict_of_mods(internal_mods) if internal_mods is not None else None
    )
    intervals = fix_intervals_input(intervals) if intervals is not None else None
    charge_adducts = (
        fix_list_of_mods(charge_adducts) if charge_adducts is not None else None
    )

    return ProFormaAnnotation(
        sequence=sequence,
        isotope_mods=isotope_mods,
        static_mods=static_mods,
        labile_mods=labile_mods,
        unknown_mods=unknown_mods,
        nterm_mods=nterm_mods,
        cterm_mods=cterm_mods,
        internal_mods=internal_mods,
        intervals=intervals,
        charge=charge,
        charge_adducts=charge_adducts,
    )


@dataclass
class MultiProFormaAnnotation:
    """
    A multi proforma annotation
    """

    annotations: List[ProFormaAnnotation]
    connections: List[bool]

    def serialize(self, include_plus: bool = False) -> str:
        """
        Convert the multi annotation to a proforma string.

        :return: The serialized multi annotation
        :rtype: str
        """
        seq = ""
        for i, annotation in enumerate(self.annotations):
            seq += annotation.serialize(include_plus=include_plus)
            if i != len(self.annotations) - 1:
                connection = self.connections[i]
                if connection is True:
                    seq += r"\\"
                else:
                    seq += r"+"

        return seq


def create_multi_annotation(
    annotations: List[ProFormaAnnotation], connections: List[bool]
) -> MultiProFormaAnnotation:
    """
    Create a MultiProFormaAnnotation from a list of annotations and connections

    :param annotations: The list of annotations
    :type annotations: List[ProFormaAnnotation]

    :raises ValueError: The number of connections should be one less than the number of annotations

    :param connections: The list of connections
    :type connections: List[bool]

    .. code-block:: python

        >>> annotation = create_multi_annotation([create_annotation('PEP'), create_annotation('TIDE')], [False])
        >>> annotation.annotations[0]
        ProFormaAnnotation(sequence=PEP)
        >>> annotation.annotations[1]
        ProFormaAnnotation(sequence=TIDE)
        >>> annotation.connections
        [False]
    """

    if len(annotations) != len(connections) + 1:
        raise ValueError(
            "The number of connections should be one less than the number of annotations"
        )

    return MultiProFormaAnnotation(annotations=annotations, connections=connections)


def _is_unmodified(proforma_sequence: str) -> bool:
    """
    Check if a proforma sequence is unmodified

    :param proforma_sequence: The proforma sequence
    :type proforma_sequence: str

    :return: True if the sequence is unmodified
    :rtype: bool
    """
    return all(c in AMINO_ACIDS for c in proforma_sequence)


class _ProFormaParser:
    """
    A proforma sequence parser
    """

    def __init__(self, proforma_sequence: str):
        self.sequence = proforma_sequence
        self.position = 0
        self.length = len(proforma_sequence)
        self._amino_acids: List[str] = []
        self._isotope_mods = None
        self._static_mods = None
        self._labile_mods = None
        self._unknown_mods = None
        self._nterm_mods = None
        self._cterm_mods = None
        self._internal_mods = None
        self._charge = None
        self._charge_adducts = None
        self._intervals = None
        self._current_connection = None

    def parse(self) -> Generator[Tuple[ProFormaAnnotation, bool], None, None]:
        """
        Parse the proforma sequence, yielding annotations and connections

        :return: A generator of annotations and connections
        :rtype: Generator[(ProFormaAnnotation, bool), None, None]
        """
        while not self._end_of_sequence():
            self._parse_sequence_start()
            self._parse_sequence_middle()
            self._parse_sequence_end()

            yield self._get_result(), self._current_connection

            if not self._end_of_sequence():
                self._reset_sequence()

    @property
    def _unmod_sequence(self) -> str:
        return "".join(self._amino_acids)

    def _get_result(self) -> ProFormaAnnotation:
        return ProFormaAnnotation(
            sequence=self._unmod_sequence,
            isotope_mods=self._isotope_mods,
            static_mods=self._static_mods,
            labile_mods=self._labile_mods,
            unknown_mods=self._unknown_mods,
            nterm_mods=self._nterm_mods,
            cterm_mods=self._cterm_mods,
            internal_mods=self._internal_mods,
            intervals=self._intervals,
            charge=self._charge,
            charge_adducts=self._charge_adducts,
        )

    def _reset_sequence(self) -> None:
        self._amino_acids = []
        self._isotope_mods = None
        self._static_mods = None
        self._labile_mods = None
        self._unknown_mods = None
        self._nterm_mods = None
        self._cterm_mods = None
        self._internal_mods = None
        self._charge = None
        self._charge_adducts = None
        self._intervals = None

    @_validate_single_mod_multiplier
    def _add_static_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._static_mods is None:
            self._static_mods = []
        if isinstance(mod, list):
            self._static_mods.extend(mod)
        else:
            self._static_mods.append(mod)

    @_validate_single_mod_multiplier
    def _add_isotope_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._isotope_mods is None:
            self._isotope_mods = []
        if isinstance(mod, list):
            self._isotope_mods.extend(mod)
        else:
            self._isotope_mods.append(mod)

    def _add_labile_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._labile_mods is None:
            self._labile_mods = []
        if isinstance(mod, list):
            self._labile_mods.extend(mod)
        else:
            self._labile_mods.append(mod)

    def _add_unknown_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._unknown_mods is None:
            self._unknown_mods = []
        if isinstance(mod, list):
            self._unknown_mods.extend(mod)
        else:
            self._unknown_mods.append(mod)

    def _add_nterm_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._nterm_mods is None:
            self._nterm_mods = []
        if isinstance(mod, list):
            self._nterm_mods.extend(mod)
        else:
            self._nterm_mods.append(mod)

    def _add_cterm_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._cterm_mods is None:
            self._cterm_mods = []
        if isinstance(mod, list):
            self._cterm_mods.extend(mod)
        else:
            self._cterm_mods.append(mod)

    def _add_internal_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._internal_mods is None:
            self._internal_mods = {}
        position = len(self._amino_acids) - 1
        if position not in self._internal_mods:
            self._internal_mods[position] = []
        if isinstance(mod, list):
            self._internal_mods[position].extend(mod)
        else:
            self._internal_mods[position].append(mod)

    def _add_interval(self, interval: Union[Interval, List[Interval]]) -> None:
        if self._intervals is None:
            self._intervals = []
        if isinstance(interval, list):
            self._intervals.extend(interval)
        else:
            self._intervals.append(interval)

    @_validate_single_mod_multiplier
    def _add_charge_adducts(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._charge_adducts is None:
            self._charge_adducts = []
        if isinstance(mod, list):
            self._charge_adducts.extend(mod)
        else:
            self._charge_adducts.append(mod)

    def _parse_sequence_start(self) -> None:
        """
        Parse the start of the sequence, up until the first amino acid or interval
        """

        while not self._end_of_sequence():

            cur = self._current()
            if cur in AMINO_ACIDS or cur == "(":  # End of start sequence
                return
            if cur == "[":  # N-term or unknown mods
                mods = self._parse_modifications("[", "]")
                next_char = self._parse_char()
                if next_char == "-":
                    self._add_nterm_mod(mods)
                elif next_char == "?":
                    self._add_unknown_mod(mods)
                else:
                    raise ProFormaFormatError(
                        f"Expected '-' or '?, but got {cur}",
                        self.position,
                        self.sequence,
                    )
            elif cur == "<":  # Global mods
                for mod in self._parse_modifications("<", ">"):
                    if "@" in mod.val:  # Static mod
                        try:
                            self._add_static_mod(mod)
                        except (
                            ValueError
                        ) as err:  # re-raise error with position, and sequence
                            raise ProFormaFormatError(
                                err, self.position, self.sequence
                            ) from err
                    else:  # Isotope mod
                        try:
                            self._add_isotope_mod(mod)
                        except (
                            ValueError
                        ) as err:  # re-raise error with position, and sequence
                            raise ProFormaFormatError(
                                err, self.position, self.sequence
                            ) from err
            elif cur == "{":  # Labile mods
                self._add_labile_mod(self._parse_modification("{", "}"))
            else:
                raise ProFormaFormatError(
                    r"Expected amino acid, '[', '{', or '<' but got " + cur,
                    self.position,
                    self.sequence,
                )

    def _parse_sequence_middle(self) -> None:
        """
        Parse the middle of the sequence, up until the end of the sequence defined by '/' or '+'
        """
        dummy_interval = None
        while not self._end_of_sequence():
            cur = self._current()
            if cur in AMINO_ACIDS:  # Amino acid
                self._amino_acids.append(self._parse_char())
            elif cur == "[":  # mods for the previous amino acid
                self._add_internal_mod(self._parse_modifications("[", "]"))
            elif cur == "-":  # cterm mods (end of sequence)
                self._skip(1)
                self._add_cterm_mod(self._parse_modifications("[", "]"))
                return
            elif cur in ("/", "+"):  # charge ( end of sequence)
                return
            elif cur == "(":  # Interval start
                if dummy_interval is not None:
                    raise ProFormaFormatError(
                        "Overlapping intervals!", self.position, self.sequence
                    )
                dummy_interval = [len(self._amino_acids), None, False, None]
                self._skip(1)
            elif cur == ")":  # Interval end
                if dummy_interval is None:
                    raise ProFormaFormatError(
                        "Interval ended without starting!", self.position, self.sequence
                    )
                dummy_interval[1] = len(self._amino_acids)

                self._skip(1)
                if not self._end_of_sequence() and self._current() == "[":
                    dummy_interval[3] = self._parse_modifications("[", "]")

                self._add_interval(
                    Interval(
                        start=dummy_interval[0],
                        end=dummy_interval[1],
                        ambiguous=dummy_interval[2],
                        mods=dummy_interval[3],
                    )
                )
                dummy_interval = None

            elif cur == "?":  # unknown mods
                if dummy_interval is None:
                    raise ProFormaFormatError(
                        "Unknown mod outside of interval", self.position, self.sequence
                    )

                # interval is ambiguous
                dummy_interval[2] = True
                self._skip(1)
            else:
                raise ProFormaFormatError(
                    f"Expected either '[', '(', '?', '-', or '/' but got: {cur}",
                    self.position,
                    self.sequence,
                )

    def _parse_sequence_end(self) -> None:
        """
        Parse the end of the sequence, up until the end of the sequence or start of the next sequence
        """
        while not self._end_of_sequence():
            cur = self._current()
            if cur == "/":  # charge
                self._skip(1)

                # check for // (crosslink)
                if not self._end_of_sequence() and self._current() == "/":
                    self._skip(1)
                    self._current_connection = True
                    return

                self._charge = self._parse_integer()

                # check for charge adducts
                if not self._end_of_sequence() and self._current() == "[":
                    self._add_charge_adducts(self._parse_modifications("[", "]"))

            elif cur == "+":  # next sequence
                self._skip(1)
                self._current_connection = False
                return
            else:
                raise ProFormaFormatError(
                    f"Invalid sequence: expected '/' or '+' but got {cur}",
                    self.position,
                    self.sequence,
                )

    def _parse_char(self) -> str:
        # Assuming any character not a '[' or ']' is an amino acid for simplicity
        aa = self._current()
        self.position += 1
        return aa

    def _parse_modifications(
        self, opening_bracket="[", closing_bracket="]"
    ) -> List[Mod]:
        """
        Parses modifications from the sequence starting with the current position. The function will continue parsing
        until it reaches the end of the sequence or the there are no more sequential modifications.
        """
        mods = []
        while not self._end_of_sequence():
            if self._current() == opening_bracket:
                mod = self._parse_modification(opening_bracket, closing_bracket)
                mods.append(mod)
            else:
                break

        return mods

    def _parse_modification(self, opening_bracket="[", closing_bracket="]") -> Mod:
        """
        Parses a single modification from the sequence starting with the current position.
        """
        self.position += 1
        start = self.position
        bracket_depth = 1
        while not self._end_of_sequence() and bracket_depth > 0:
            if self.sequence[self.position] == opening_bracket:
                bracket_depth += 1
            elif self.sequence[self.position] == closing_bracket:
                bracket_depth -= 1
            self.position += 1

        if bracket_depth != 0:
            msg = f"Unmatched {opening_bracket} at position {self.position}"
            raise ProFormaFormatError(msg, self.position, self.sequence)

        mod = self.sequence[start : self.position - 1]

        multiplier = 1
        if not self._end_of_sequence() and self._peek() == "^":
            self.position += 1
            multiplier_start = self.position
            while not self._end_of_sequence() and self._peek().isdigit():
                self.position += 1
            multiplier = int(self.sequence[multiplier_start : self.position])

        return Mod(mod, multiplier)

    def _parse_integer(self) -> int:
        start = self.position
        digit_count = 0
        while not self._end_of_sequence():
            if self._peek().isdigit():
                digit_count += 1
                self.position += 1
            elif digit_count == 0 and self._peek() in ["+", "-"]:
                self.position += 1
            else:
                break
        return int(self.sequence[start : self.position])

    def _current(self) -> str:
        return self.sequence[self.position]

    def _peek(self) -> str:
        return self.sequence[self.position] if not self._end_of_sequence() else None

    def _skip(self, n) -> None:
        self.position += n

    def _end_of_sequence(self) -> bool:
        return self.position >= self.length


def parse(sequence: str) -> Union[ProFormaAnnotation, MultiProFormaAnnotation]:
    """
    Parses a ProForma sequence string and returns its corresponding annotation object.

    Note that the function's behavior and the type of object returned depend on the structure of the input sequence.
    Single sequences result in ProFormaAnnotation objects, while multi-sequences result in MultiProFormaAnnotation
    objects.

    :param sequence: The sequence to parse.
    :type sequence: str

    :raises ProFormaFormatError: If the sequence is not valid.

    :return: Either a ProFormaAnnotation or a MultiProFormaAnnotation, based on the input
    :rtype: Union[ProFormaAnnotation, MultiProFormaAnnotation]

    .. python::

        Parsing a simple peptide sequence:
        >>> isinstance(parse('PEPTIDE'), ProFormaAnnotation)
        True

        Parsing a sequence with multiple peptides or complex modifications:
        >>> isinstance(parse('PEPTIDE+PEPTIDE'), MultiProFormaAnnotation)
        True


    """
    if _is_unmodified(sequence) is True:
        return ProFormaAnnotation(sequence=sequence)

    annotations_connections = list(_ProFormaParser(sequence).parse())
    annotations = [annotation for annotation, _ in annotations_connections]
    annotations_connections = [connection for _, connection in annotations_connections]

    if len(annotations) == 1:
        return annotations[0]

    return MultiProFormaAnnotation(annotations, annotations_connections[:-1])


def serialize(
    annotation: Union[ProFormaAnnotation, MultiProFormaAnnotation],
    include_plus: bool = False,
) -> str:
    """
    Serializes a ProForma annotation or multiple ProForma annotations into a single string representation.

    :param annotation: Either a ProFormaAnnotation or a MultiProFormaAnnotation.
    :type annotation: Union[ProFormaAnnotation, MultiProFormaAnnotation]

    :return: A string representation of the ProForma annotation.
    :rtype: str

    . python::

        Serializing a simple ProForma annotation:
        >>> serialize(ProFormaAnnotation(sequence='PEPTIDE'))
        'PEPTIDE'

        >>> pfa1 = ProFormaAnnotation(sequence='PEPTIDE')
        >>> pfa2 = ProFormaAnnotation(sequence='PEPTIDE')

        Serializing a MultiProFormaAnnotation with chimeric connections:
        >>> multi_annotation = MultiProFormaAnnotation([pfa1, pfa2], [False])
        >>> serialize(multi_annotation)
        'PEPTIDE+PEPTIDE'

        Serializing a MultiProFormaAnnotation with crosslink connections:
        >>> multi_annotation = MultiProFormaAnnotation([pfa1, pfa2], [True])
        >>> p = serialize(multi_annotation)
        >>> p == r'PEPTIDE\\\PEPTIDE'
        True

    """

    return annotation.serialize(include_plus)


def _serialize_annotation_start(
    annotation: ProFormaAnnotation,
    include_plus: bool,
    precision: Optional[float] = None,
) -> str:
    comps = []

    # add labile mods
    if annotation.has_labile_mods():
        for mod in annotation.labile_mods:
            comps.append(mod.serialize("{}", include_plus, precision))

    if annotation.has_static_mods():
        for mod in annotation.static_mods:
            comps.append(mod.serialize("<>", include_plus, precision))

    # Add global mods
    if annotation.has_isotope_mods():
        for mod in annotation.isotope_mods:
            comps.append(mod.serialize("<>", include_plus, precision))

    # Unknown mods
    if annotation.has_unknown_mods():
        for mod in annotation.unknown_mods:
            comps.append(mod.serialize("[]", include_plus, precision))
        comps.append("?")

    # N-term mods
    if annotation.has_nterm_mods():
        for mod in annotation.nterm_mods:
            comps.append(mod.serialize("[]", include_plus, precision))
        comps.append("-")

    return "".join(comps)


def _serialize_annotation_middle(
    annotation: ProFormaAnnotation,
    include_plus: bool,
    precision: Optional[float] = None,
) -> str:
    comps = []
    # Sequence
    for i, aa in enumerate(annotation.sequence):

        if annotation.intervals:
            for interval in annotation.intervals:
                if interval.start == i:
                    comps.append("(")
                    if interval.ambiguous:
                        comps.append("?")
                if interval.end == i:
                    comps.append(")")

                    if interval.mods:
                        for mod in interval.mods:
                            comps.append(mod.serialize("[]", include_plus, precision))

        comps.append(aa)

        # Internal mods
        if annotation.internal_mods and i in annotation.internal_mods:
            for mod in annotation.internal_mods[i]:
                comps.append(mod.serialize("[]", include_plus, precision))

    # add end interval
    i = len(annotation.sequence)
    if annotation.intervals:
        for interval in annotation.intervals:
            if interval.end == i:
                comps.append(")")
                if interval.mods:
                    for mod in interval.mods:
                        comps.append(mod.serialize("[]", include_plus, precision))

    return "".join(comps)


def _serialize_annotation_end(
    annotation: ProFormaAnnotation,
    include_plus: bool,
    precision: Optional[float] = None,
) -> str:
    comps = []
    # C-term mods
    if annotation.cterm_mods:
        comps.append("-")
        for mod in annotation.cterm_mods:
            comps.append(mod.serialize("[]", include_plus, precision))

    # Charge
    if annotation.charge:
        comps.append(f"/{annotation.charge}")

    if annotation.charge_adducts:
        for mod in annotation.charge_adducts:
            comps.append(mod.serialize("[]", include_plus, precision))

    return "".join(comps)


def _serialize_annotation(
    annotation: ProFormaAnnotation,
    include_plus: bool,
    precision: Optional[float] = None,
) -> str:
    return (
        _serialize_annotation_start(annotation, include_plus, precision)
        + _serialize_annotation_middle(annotation, include_plus, precision)
        + _serialize_annotation_end(annotation, include_plus, precision)
    )
