import copy
import collections as co
from copy import deepcopy
from typing import *
import warnings
import regex as re

from ..chem.chem_constants import AVERAGE_AA_MASSES, MONOISOTOPIC_AA_MASSES
from ..chem.chem_util import chem_mass
from ..mass_calc import adjust_mass, adjust_mz, mod_mass
from .mod_list import ModList
from ..chem.chem_calc import (
    _parse_charge_adducts_comp,
    _parse_mod_delta_mass_only,
    apply_isotope_mods_to_composition,
    estimate_comp,
    mod_comp,
)
from .utils import parse_static_mods
from ..proforma_dataclasses import (
    ACCEPTED_MOD_INPUT,
    fix_list_of_mods,
    fix_dict_of_mods,
    fix_intervals_input,
    Mod,
    Interval,
    are_mods_equal,
    are_intervals_equal,
)
from ..constants import (
    AA_COMPOSITIONS,
    AMBIGUOUS_AMINO_ACIDS,
    FRAGMENT_ION_BASE_CHARGE_ADDUCTS,
    MASS_AMBIGUOUS_AMINO_ACIDS,
    NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS,
    ModType,
    get_mod_type,
    get_mods,
)
from ..errors import AmbiguousAminoAcidError, UnknownAminoAcidError
from ..proforma_dataclasses import (
    INTERVAL_VALUE,
    ChemComposition,
    ModValue,
)
from ..util import (
    combine_ambiguity_intervals,
    _construct_ambiguity_intervals,
    get_mass_shift_interval,
)


def _setup_mod_list(mods: Optional[ACCEPTED_MOD_INPUT]) -> ModList:
    """
    Helper function to set up a ModList from various input types.
    """
    mod_list: ModList = ModList()
    if mods is not None:
        if isinstance(mods, (str, int, float, Mod)):
            mod_list.append(mods)
        elif isinstance(mods, Iterable):
            mod_list.extend(mods)
        else:
            raise TypeError(
                f"Invalid type for isotope_mods: {type(mods)}. Must be a string, int, float, Mod, or iterable of these."
            )
    return mod_list


class ProFormaAnnotationBase:

    def __init__(
        self,
        sequence: Optional[str] = None,
        isotope_mods: Optional[ACCEPTED_MOD_INPUT] = None,
        static_mods: Optional[ACCEPTED_MOD_INPUT] = None,
        labile_mods: Optional[ACCEPTED_MOD_INPUT] = None,
        unknown_mods: Optional[ACCEPTED_MOD_INPUT] = None,
        nterm_mods: Optional[ACCEPTED_MOD_INPUT] = None,
        cterm_mods: Optional[ACCEPTED_MOD_INPUT] = None,
        internal_mods: Optional[Dict[int, List[Mod]]] = None,
        intervals: Optional[List[Interval]] = None,
        charge: Optional[int] = None,
        charge_adducts: Optional[ACCEPTED_MOD_INPUT] = None,
    ) -> None:

        self._sequence: Optional[str] = sequence if sequence is not None else ""

        self._isotope_mods = _setup_mod_list(isotope_mods)
        self._static_mods = _setup_mod_list(static_mods)
        self._labile_mods = _setup_mod_list(labile_mods)
        self._unknown_mods = _setup_mod_list(unknown_mods)
        self._nterm_mods = _setup_mod_list(nterm_mods)
        self._cterm_mods = _setup_mod_list(cterm_mods)
        self._charge_adducts = _setup_mod_list(charge_adducts)

        self._internal_mods: Optional[Dict[int, List[Mod]]] = (
            fix_dict_of_mods(internal_mods) if internal_mods is not None else {}
        )
        self._intervals: Optional[List[Interval]] = (
            fix_intervals_input(intervals) if intervals is not None else []
        )
        self._charge: Optional[int] = charge

    """
    Magic Methods
    """

    def __repr__(self) -> str:
        """
        Only shows non-None types
        """
        seq = f"ProFormaAnnotationBase(sequence={self.sequence}"

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
        if not isinstance(other, ProFormaAnnotationBase):
            return NotImplemented

        return (
            self.sequence == other.sequence
            and are_mods_equal(self.labile_mods, other.labile_mods)
            and are_mods_equal(self.unknown_mods, other.unknown_mods)
            and are_mods_equal(self.nterm_mods, other.nterm_mods)
            and are_mods_equal(self.cterm_mods, other.cterm_mods)
            and are_mods_equal(self.charge_adducts, other.charge_adducts)
            and are_mods_equal(self.isotope_mods, other.isotope_mods)
            and are_mods_equal(self.static_mods, other.static_mods)
            and self._compare_internal_mods(other)
            and are_intervals_equal(self.intervals, other.intervals)
            and self.charge == other.charge
        )

    def __hash__(self):
        return hash(
            (
                self.sequence,
                tuple(self.isotope_mods),
                tuple(self.static_mods),
                tuple(self.labile_mods),
                tuple(self.unknown_mods),
                tuple(self.nterm_mods),
                tuple(self.cterm_mods),
                frozenset((k, tuple(v)) for k, v in self.internal_mods.items()),
                tuple(self.intervals),
                self.charge,
                tuple(self.charge_adducts),
            )
        )

    """
    Has Methods:
    """

    def has_sequence(self) -> bool:
        return len(self.sequence) > 0

    def has_isotope_mods(self) -> bool:
        return len(self.isotope_mods) > 0

    def has_static_mods(self) -> bool:
        return len(self.static_mods) > 0

    def has_labile_mods(self) -> bool:
        return len(self.labile_mods) > 0

    def has_unknown_mods(self) -> bool:
        return len(self.unknown_mods) > 0

    def has_nterm_mods(self) -> bool:
        return len(self.nterm_mods) > 0

    def has_cterm_mods(self) -> bool:
        return len(self.cterm_mods) > 0

    def has_charge_adducts(self) -> bool:
        return len(self.charge_adducts) > 0

    def has_internal_mods(self) -> bool:
        return len(self.internal_mods) > 0

    def has_intervals(self) -> bool:
        return len(self.intervals) > 0

    def has_unambiguous_intervals(self) -> bool:
        return len(self.unambiguous_intervals) > 0

    def has_ambiguous_intervals(self) -> bool:
        return len(self.amiguous_intervals) > 0

    def has_charge(self) -> bool:
        return self.charge is not None

    def has_mods(self) -> bool:
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

    """
    Properties
    """

    @property
    def sequence(self) -> str:
        if self._sequence is None:
            return ""
        return self._sequence

    @property
    def stripped_sequence(self) -> str:
        return self.sequence

    @property
    def isotope_mods(self) -> List[Mod]:
        return self._isotope_mods.data

    @property
    def static_mods(self) -> List[Mod]:
        return self._static_mods.data

    @property
    def labile_mods(self) -> List[Mod]:
        return self._labile_mods.data

    @property
    def unknown_mods(self) -> List[Mod]:
        return self._unknown_mods.data

    @property
    def nterm_mods(self) -> List[Mod]:
        return self._nterm_mods.data

    @property
    def cterm_mods(self) -> List[Mod]:
        return self._cterm_mods.data

    @property
    def charge_adducts(self) -> List[Mod]:
        return self._charge_adducts.data

    @property
    def internal_mods(self) -> Dict[int, List[Mod]]:

        if self._internal_mods is None:
            return {}

        return self._internal_mods

    @property
    def intervals(self) -> List[Interval]:
        if self._intervals is None:
            return []

        return self._intervals

    @property
    def ambiguous_intervals(self) -> List[Interval]:
        return [interval for interval in self.intervals if interval.ambiguous is True]

    @property
    def unambiguous_intervals(self) -> List[Interval]:
        return [interval for interval in self.intervals if interval.ambiguous is False]

    @property
    def charge(self) -> Optional[int]:
        return self._charge

    """
    Property Setters
    """

    @sequence.setter
    def sequence(self, value: Optional[str]) -> None:
        self._sequence = value

    @isotope_mods.setter
    def isotope_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        self._isotope_mods = ModList(value)

    @static_mods.setter
    def static_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        self._static_mods = ModList(value)

    @labile_mods.setter
    def labile_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        self._labile_mods = ModList(value)

    @unknown_mods.setter
    def unknown_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        self._unknown_mods = ModList(value)

    @nterm_mods.setter
    def nterm_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        self._nterm_mods = ModList(value)

    @cterm_mods.setter
    def cterm_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        self._cterm_mods = ModList(value)

    @charge_adducts.setter
    def charge_adducts(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        self._charge_adducts = ModList(value)

    @internal_mods.setter
    def internal_mods(
        self, value: Optional[Dict[int, Union[List[ModValue], ModValue]]]
    ) -> None:
        if value is None:
            self._internal_mods = None
        else:
            value = fix_dict_of_mods(value)
            self._internal_mods = copy.deepcopy(value)

    @intervals.setter
    def intervals(
        self, value: Optional[Union[List[INTERVAL_VALUE], INTERVAL_VALUE]]
    ) -> None:
        if value is None:
            self._intervals = None
        else:
            value = fix_intervals_input(value)
            self._intervals = copy.deepcopy(value)

    @charge.setter
    def charge(self, value: Optional[int]):
        self._charge = value

    """
    Pop Methods
    """

    def _pop_mod_helper(self, attr_name: str) -> Any:
        """Helper method to pop and clear a modification attribute."""
        value = getattr(self, attr_name)
        setattr(self, attr_name, None)
        return value

    def pop_isotope_mods(self) -> List[Mod]:
        return self._pop_mod_helper("isotope_mods")

    def pop_static_mods(self) -> List[Mod]:
        return self._pop_mod_helper("static_mods")

    def pop_labile_mods(self) -> List[Mod]:
        return self._pop_mod_helper("labile_mods")

    def pop_unknown_mods(self) -> List[Mod]:
        return self._pop_mod_helper("unknown_mods")

    def pop_nterm_mods(self) -> List[Mod]:
        return self._pop_mod_helper("nterm_mods")

    def pop_cterm_mods(self) -> List[Mod]:
        return self._pop_mod_helper("cterm_mods")

    def pop_charge_adducts(self) -> List[Mod]:
        return self._pop_mod_helper("charge_adducts")

    def pop_internal_mods(self) -> Dict[int, List[Mod]]:
        return self._pop_mod_helper("internal_mods")

    def pop_intervals(self) -> List[Interval]:
        return self._pop_mod_helper("intervals")

    def pop_ambiguous_intervals(self) -> List[Interval]:
        m = [interval for interval in self.intervals if interval.ambiguous is True]
        self.intervals = [
            interval for interval in self.intervals if interval.ambiguous is False
        ]
        return m

    def pop_unambiguous_intervals(self) -> List[Interval]:
        m = [interval for interval in self.intervals if interval.ambiguous is False]
        self.intervals = [
            interval for interval in self.intervals if interval.ambiguous is True
        ]
        return m

    def pop_charge(self) -> Optional[int]:
        return self._pop_mod_helper("charge")

    def pop_mod_by_type(self, mod_type: Union[str, ModType]) -> List[Mod]:
        """
        Get the pop method for a given mod type.
        """
        if isinstance(mod_type, str):
            mod_type = get_mod_type(mod_type)
        if not isinstance(mod_type, ModType):
            raise ValueError(
                f"Invalid mod type: {mod_type}. Must be a ModType or a string representation of a ModType."
            )

        mod_type_to_pop_method = {
            ModType.ISOTOPE: self.pop_isotope_mods,
            ModType.STATIC: self.pop_static_mods,
            ModType.LABILE: self.pop_labile_mods,
            ModType.UNKNOWN: self.pop_unknown_mods,
            ModType.NTERM: self.pop_nterm_mods,
            ModType.CTERM: self.pop_cterm_mods,
            ModType.INTERNAL: self.pop_internal_mods,
            ModType.INTERVAL: self.pop_intervals,
            ModType.CHARGE_ADDUCTS: self.pop_charge_adducts,
            ModType.CHARGE: self.pop_charge,
        }

        mods = mod_type_to_pop_method.get(mod_type, None)()
        if mods is None:
            return []

        return mods

    def pop_mods(
        self,
        mods: Optional[Union[str, ModType, List[str], List[ModType]]] = None,
        condense: bool = True,
    ) -> Dict[str, Any]:
        """
        Pop all mods and return them in a dictionary
        """

        mod_types: List[ModType] = get_mods(mods)

        d = {}
        for mod_type in mod_types:
            d[mod_type.value] = self.pop_mod_by_type(mod_type)

        if condense is True:
            for key in list(d.keys()):
                if (
                    d[key] is None
                    or (isinstance(d[key], list) and len(d[key]) == 0)
                    or (isinstance(d[key], dict) and len(d[key]) == 0)
                ):
                    d.pop(key)

        return d

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

    """
    Add Methods
    """

    def _add_mods_helper(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        mod_type: str,
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        """
        Helper method to add modifications to the annotation.
        """
        if inplace is False:
            return self.copy()._add_mods_helper(mods, mod_type, append, inplace=True)

        if mods is None:
            if not append:
                setattr(self, mod_type, None)
            return self

        mods = fix_list_of_mods(mods)

        if not append:
            setattr(self, mod_type, mods)
        else:
            has_mods = getattr(self, f"has_{mod_type}")()
            if has_mods:
                getattr(self, mod_type).extend(mods)
            else:
                setattr(self, mod_type, mods)

        return self

    def add_isotope_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        return self._add_mods_helper(mods, "isotope_mods", append, inplace)

    def add_static_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        return self._add_mods_helper(mods, "static_mods", append, inplace)

    def add_labile_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        return self._add_mods_helper(mods, "labile_mods", append, inplace)

    def add_unknown_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        return self._add_mods_helper(mods, "unknown_mods", append, inplace)

    def add_nterm_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        return self._add_mods_helper(mods, "nterm_mods", append, inplace)

    def add_cterm_mods(
        self,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        return self._add_mods_helper(mods, "cterm_mods", append, inplace)

    def add_charge_adducts(
        self,
        charge_adducts: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        return self._add_mods_helper(charge_adducts, "charge_adducts", append, inplace)

    def add_internal_mod(
        self,
        index: int,
        mods: Optional[Union[List[ModValue], ModValue]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().add_internal_mod(
                index=index, mods=mods, append=append, inplace=True
            )

        if mods is None:
            if not append:
                if not self.has_internal_mods():
                    self.internal_mods = {}
                self.internal_mods.pop(index, None)
            return self

        mods = fix_list_of_mods(mods)

        if not self.has_internal_mods():
            self.internal_mods = {}

        if not append:
            self.internal_mods[index] = mods
        else:
            if index in self.internal_mods:
                self.internal_mods[index].extend(mods)
            else:
                self.internal_mods[index] = mods

        return self

    def add_internal_mods(
        self,
        mods: Optional[Dict[int, Union[List[ModValue], ModValue]]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        """
        Add internal mods to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().add_internal_mods(mods=mods, append=append, inplace=True)

        if mods is None:
            if not append:
                self.internal_mods = None
            return self

        mods = fix_dict_of_mods(mods)

        if not append:
            self.internal_mods = mods
        else:
            if not self.has_internal_mods():
                self.internal_mods = mods
            else:
                for k, v in mods.items():
                    self.add_internal_mod(k, v, append=True, inplace=True)

        return self

    def add_intervals(
        self,
        intervals: Optional[Union[List[INTERVAL_VALUE], INTERVAL_VALUE]],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        """
        Add intervals to the annotation. If not append, existing mods will be replaced.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().add_intervals(intervals, append, inplace=True)

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

    def add_charge(
        self, charge: Optional[int], inplace: bool = True
    ) -> "ProFormaAnnotationBase":
        """
        Add charge to the annotation
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().add_charge(charge, inplace=True)

        self.charge = charge
        return self

    def add_mod_by_type(
        self,
        mod: Any,
        mod_type: Union[str, ModType],
        append: bool = False,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        """
        Add a modification by type.
        """

        if not inplace:
            return self.copy().add_mod_by_type(mod, mod_type, append, True)

        mod_type = get_mod_type(mod_type)

        mod_type_methods = {
            ModType.ISOTOPE: lambda: self.add_isotope_mods(mod, append),
            ModType.STATIC: lambda: self.add_static_mods(mod, append),
            ModType.LABILE: lambda: self.add_labile_mods(mod, append),
            ModType.UNKNOWN: lambda: self.add_unknown_mods(mod, append),
            ModType.NTERM: lambda: self.add_nterm_mods(mod, append),
            ModType.CTERM: lambda: self.add_cterm_mods(mod, append),
            ModType.INTERNAL: lambda: self.add_internal_mods(mod, append),
            ModType.INTERVAL: lambda: self.add_intervals(mod, append),
            ModType.CHARGE: lambda: setattr(self, "charge", mod),
            ModType.CHARGE_ADDUCTS: lambda: self.add_charge_adducts(mod, append),
        }

        mod_type_methods[mod_type]()
        return self

    def add_mod_dict(
        self, mod_dict: Dict, append: bool = False, inplace: bool = True
    ) -> "ProFormaAnnotationBase":
        """
        Add mods from a dictionary
        """

        if inplace is False:
            return self.copy().add_mod_dict(
                mod_dict=mod_dict, append=append, inplace=True
            )

        # Use add_mod_by_type for known mod types
        for mod_type_str, mod_value in mod_dict.items():

            if isinstance(mod_type_str, int):
                # internal mod
                self.add_internal_mod(mod_type_str, mod_value, append, inplace=True)
                continue

            mod_type = get_mod_type(mod_type_str)
            self.add_mod_by_type(mod_value, mod_type, append, inplace=True)

        return self

    def add_mods(
        self, mods: Dict, append: bool = False, inplace: bool = True
    ) -> "ProFormaAnnotationBase":

        if inplace is False:
            return self.copy().add_mods(mods=mods, append=append, inplace=True)

        if not isinstance(mods, dict):
            raise TypeError(f"Expected a dictionary of mods, got {type(mods)}")

        return self.add_mod_dict(mods, append=append, inplace=True)

    """
    Remove Methods
    """

    def remove_isotope_mods(self, inplace: bool = True) -> "ProFormaAnnotationBase":
        return self._remove_mod_helper("pop_isotope_mods", inplace)

    def remove_static_mods(self, inplace: bool = True) -> "ProFormaAnnotationBase":
        return self._remove_mod_helper("pop_static_mods", inplace)

    def remove_nterm_mods(self, inplace: bool = True) -> "ProFormaAnnotationBase":
        return self._remove_mod_helper("pop_nterm_mods", inplace)

    def remove_cterm_mods(self, inplace: bool = True) -> "ProFormaAnnotationBase":
        return self._remove_mod_helper("pop_cterm_mods", inplace)

    def remove_internal_mods(self, inplace: bool = True) -> "ProFormaAnnotationBase":
        return self._remove_mod_helper("pop_internal_mods", inplace)

    def remove_intervals(self, inplace: bool = True) -> "ProFormaAnnotationBase":
        return self._remove_mod_helper("pop_intervals", inplace)

    def remove_charge(self, inplace: bool = True) -> "ProFormaAnnotationBase":
        return self._remove_mod_helper("pop_charge", inplace)

    def remove_charge_adducts(self, inplace: bool = True) -> "ProFormaAnnotationBase":
        return self._remove_mod_helper("pop_charge_adducts", inplace)

    def remove_mods(
        self,
        mods: Optional[Union[str, ModType, List[str], list[ModType]]] = None,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        if inplace is False:
            return self.copy().remove_mods(mods=mods, inplace=True)
        _ = self.pop_mods(mods=mods, condense=False)
        return self

    def _remove_mod_helper(
        self, method_name: str, inplace: bool
    ) -> "ProFormaAnnotationBase":
        if inplace is False:
            return self.copy()._remove_mod_helper(method_name, True)
        pop_method = getattr(self, method_name)
        _ = pop_method()
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

    def filter_mods(
        self,
        mods: Optional[Union[str, ModType, List[str], list[ModType]]] = None,
        inplace: bool = True,
    ) -> "ProFormaAnnotationBase":
        """
        Filter the modifications in the annotation based on the provided mod types.
        If inplace is False, a new ProFormaAnnotationBase object is returned with the filtered mods.
        If inplace is True, the current object is modified.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().filter_mods(mods=mods, inplace=True)

        mods: List[ModType] = get_mods(mods)
        mod_types_to_remove = {mod_type for mod_type in ModType} - set(mods)

        if len(mod_types_to_remove) == 0:
            # If no mods to remove, return the annotation as is
            return self

        for mod_type in mod_types_to_remove:
            # Pop the mods of the specified type
            _ = self.pop_mod_by_type(mod_type)

        return self

    def _compare_internal_mods(self, other: "ProFormaAnnotationBase") -> bool:
        if not (self.has_internal_mods() or other.has_internal_mods()):
            return True

        self_keys = (
            set(self.internal_mods.keys()) if self.has_internal_mods() else set()
        )
        other_keys = (
            set(other.internal_mods.keys()) if other.has_internal_mods() else set()
        )

        return all(
            are_mods_equal(
                self.get_internal_mods_by_index(k), other.get_internal_mods_by_index(k)
            )
            for k in self_keys.union(other_keys)
        )

    def copy(self, deep: bool = True) -> "ProFormaAnnotationBase":
        """
        Returns a deep copy of the ProFormaAnnotationBase instance.
        """

        if not deep:
            return ProFormaAnnotationBase(
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

        # Clear empty internal mods
        if self.has_internal_mods():
            self.internal_mods = {
                k: v for k, v in self.internal_mods.items() if v and len(v) > 0
            } or None

        # Clear empty intervals
        if self.has_intervals():
            # Clean up interval mods
            for interval in self.intervals:
                if interval.mods is not None and len(interval.mods) == 0:
                    interval.mods = None

            # Remove intervals that are non-ambiguous and have no mods
            self.intervals = [
                interval
                for interval in self.intervals
                if interval.ambiguous or interval.mods is not None
            ] or None

    def dict(self) -> Dict[str, Any]:
        """
        Returns a dictionary representation of the ProFormaAnnotationBase instance.
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
        Returns a dictionary representation of the ProFormaAnnotationBase instance, including only the
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

    def condense_static_mods(self, inplace: bool = True) -> "ProFormaAnnotationBase":
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

    """
    Other Methods
    """

    def count_residues(self, include_mods: bool = True) -> Counter:
        """
        Count the occurrences of each residue in the sequence.
        """
        if not include_mods:
            return co.Counter(self.sequence)

        return co.Counter([a.serialize() for a in self.split()])

    def percent_residues(
        self, include_mods: bool = True, precision: Optional[bool] = None
    ) -> Dict[str, float]:
        """
        Calculate the percentage of each residue in the sequence.
        If include_mods is True, counts residues including modifications.
        """
        residue_counts = self.count_residues(include_mods=include_mods)
        total_count = sum(residue_counts.values())

        if total_count == 0:
            return {}

        d = {
            residue: (count / total_count) * 100
            for residue, count in residue_counts.items()
        }

        if precision is not None:
            d = {k: round(v, precision) for k, v in d.items()}

        return d

    def is_subsequence(
        self,
        other: "ProFormaAnnotationBase",
        ignore_mods: bool = False,
        ignore_intervals: bool = True,
    ) -> bool:
        """
        Check if the annotation is a subsequence of another annotation.
        ignore mods -> stripped sequence (if true intervals wont do anythong)
        ignore intervals -> remove intervals
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

                if ignore_intervals:
                    if sliced_annot.remove_intervals(
                        inplace=False
                    ) == self.remove_intervals(inplace=False):
                        return True

                if sliced_annot == self:
                    return True

        return False

    def find_indices(
        self,
        other: "ProFormaAnnotationBase",
        ignore_mods: bool = False,
        ignore_intervals: bool = True,
    ) -> List[int]:
        """
        Find all occurrences of the annotation in another annotation.
        """

        if not isinstance(other, ProFormaAnnotationBase):
            raise TypeError(
                f"other must be a ProFormaAnnotationBase, got {type(other)}"
            )

        if not self.sequence or not other.sequence:
            # if either sequence is empty, return empty list
            return []

        # find all starting indexes where the sequence is found and is a subsequence
        return [
            m.start()
            for m in re.finditer(self.sequence, other.sequence)
            if self.is_subsequence(
                other=other.slice(
                    m.start(), m.start() + len(self.sequence), inplace=False
                ),
                ignore_mods=ignore_mods,
                ignore_intervals=ignore_intervals,
            )
        ]

    def annotate_ambiguity(
        self,
        forward_coverage: List[int],
        reverse_coverage: List[int],
        mass_shift: Optional[Any] = None,
        inplace: bool = False,
    ) -> Union["ProFormaAnnotationBase", None]:
        """
        Generates ambiguity intervals based on the coverage of the sequence.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().annotate_ambiguity(
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
        ambiguity_intervals = combine_ambiguity_intervals(
            forward_intervals, reverse_intervals
        )

        intervals = [
            Interval(start, end + 1, True, None) for start, end in ambiguity_intervals
        ]

        self.add_intervals(intervals, append=True)

        if mass_shift is not None:
            mass_shift_interval = get_mass_shift_interval(
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
        :type annotation: ProFormaAnnotationBase

        :return: The delta mass.
        :rtype: float

        .. code-block:: python

            >>> mods = ProFormaAnnotationBase(sequence='', nterm_mods = [42.0, -20.0])
            >>> mods.pop_delta_mass_mods()
            22.0
            >>> mods
            ProFormaAnnotationBase(sequence=)

        """
        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().pop_delta_mass_mods(inplace=True)

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

    def coverage(
        self,
        annotations: List["ProFormaAnnotationBase"],
        accumulate: bool = False,
        ignore_mods: bool = False,
        ignore_ambiguity: bool = False,
    ) -> List[int]:

        cov_arr = [0] * len(self)

        for sub_annots in annotations:

            peptide_cov = [1] * len(sub_annots)
            if ignore_ambiguity is False:
                for interval in sub_annots.ambiguous_intervals:
                    for i in range(interval.start, interval.end):
                        peptide_cov[i] = 0

            for subsequence_index in sub_annots.find_indices(
                other=self, ignore_mods=ignore_mods
            ):
                start = subsequence_index

                for i, cov in enumerate(peptide_cov):
                    if accumulate:
                        cov_arr[start + i] += cov
                    else:
                        cov_arr[start + i] = cov

        # apply to slef
        if ignore_ambiguity is False:
            # Set coverage to 0 for ambiguous positions
            for interval in self.ambiguous_intervals:
                for i in range(interval.start, interval.end):
                    cov_arr[i] = 0

        return cov_arr

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
            return self.copy().mass(
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
            return self.copy().comp(
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
            return self.copy().comp_mass(
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
    ) -> "ProFormaAnnotationBase":

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().condense_to_mass_mods(
                include_plus=include_plus, inplace=True
            )

        labile_mods = self.pop_labile_mods()
        isotope_mods = self.isotope_mods

        # fix this so that isotopes get applied
        n_term_mods = self.pop_nterm_mods()
        c_term_mods = self.pop_cterm_mods()
        n_term_mods_mass, c_term_mods_mass = None, None

        if n_term_mods:
            n_term_annot = ProFormaAnnotationBase(
                sequence="",
                nterm_mods=n_term_mods,
                isotope_mods=isotope_mods,
            )
            n_term_mods_mass = n_term_annot.mass(ion_type="n")

        if c_term_mods:
            c_term_annot = ProFormaAnnotationBase(
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

    def _serialize_start(
        self,
        include_plus: bool,
        precision: Optional[float] = None,
    ) -> str:
        comps = []

        # add labile mods
        if self.has_labile_mods():
            for mod in self.labile_mods:
                comps.append(mod.serialize("{}", include_plus, precision))

        if self.has_static_mods():
            for mod in self.static_mods:
                comps.append(mod.serialize("<>", include_plus, precision))

        # Add global mods
        if self.has_isotope_mods():
            for mod in self.isotope_mods:
                comps.append(mod.serialize("<>", include_plus, precision))

        # Unknown mods
        if self.has_unknown_mods():
            for mod in self.unknown_mods:
                comps.append(mod.serialize("[]", include_plus, precision))
            comps.append("?")

        # N-term mods
        if self.has_nterm_mods():
            for mod in self.nterm_mods:
                comps.append(mod.serialize("[]", include_plus, precision))
            comps.append("-")

        return "".join(comps)

    def _serialize_middle(
        self,
        include_plus: bool,
        precision: Optional[float] = None,
    ) -> str:
        comps = []
        # Sequence
        for i, aa in enumerate(self.sequence):

            if self.intervals:
                for interval in self.intervals:
                    if interval.start == i:
                        comps.append("(")
                        if interval.ambiguous:
                            comps.append("?")
                    if interval.end == i:
                        comps.append(")")

                        if interval.mods:
                            for mod in interval.mods:
                                comps.append(
                                    mod.serialize("[]", include_plus, precision)
                                )

            comps.append(aa)

            # Internal mods
            if self.internal_mods and i in self.internal_mods:
                for mod in self.internal_mods[i]:
                    comps.append(mod.serialize("[]", include_plus, precision))

        # add end interval
        i = len(self.sequence)
        if self.intervals:
            for interval in self.intervals:
                if interval.end == i:
                    comps.append(")")
                    if interval.mods:
                        for mod in interval.mods:
                            comps.append(mod.serialize("[]", include_plus, precision))

        return "".join(comps)

    def _serialize_end(
        self,
        include_plus: bool,
        precision: Optional[float] = None,
    ) -> str:
        comps = []
        # C-term mods
        if self.cterm_mods:
            comps.append("-")
            for mod in self.cterm_mods:
                comps.append(mod.serialize("[]", include_plus, precision))

        # Charge
        if self.charge:
            comps.append(f"/{self.charge}")

        if self.charge_adducts:
            for mod in self.charge_adducts:
                comps.append(mod.serialize("[]", include_plus, precision))

        return "".join(comps)

    def serialize(
        self,
        include_plus: bool = False,
        precision: Optional[float] = None,
    ) -> str:
        return (
            self._serialize_start(include_plus, precision)
            + self._serialize_middle(include_plus, precision)
            + self._serialize_end(include_plus, precision)
        )
