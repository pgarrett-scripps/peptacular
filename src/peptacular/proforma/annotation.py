from collections.abc import Callable, Generator
import copy
from copy import deepcopy
from typing import Any, Iterable, Self, cast

from .serializer import serialize_annotation
from .parser import ProFormaParser
from .ambiguity import annotate_ambiguity
from .manipulation import (
    condense_static_mods,
    is_subsequence,
    count_residues,
    percent_residues,
    find_indices,
    coverage,
    percent_coverage,
)
from .slicing import (
    generate_sliding_windows,
    reverse_annotation,
    shift_annotation,
    shuffle_annotation,
    slice_annotation,
    sort_annotation,
    split_annotation,
)
from .combinatorics import (
    generate_permutations,
    generate_product,
    generate_combinations,
    generate_combinations_with_replacement,
)

from ..digestion import (
    DigestionMixin,
)

from .mass_comp import (
    comp,
    mass,
    mz
)


from .dclasses import (
    MOD_VALUE_TYPES,
    MODLIST_DATATYPE,
    ACCEPTED_INTERVALLIST_INPUT_TYPES,
    ACCEPTED_MODDICT_INPUT_TYPES,
    ACCEPTED_MODLIST_INPUT_TYPES,
    ACCEPTED_INTERVAL_DATATYPE,
    ModInterval,
    IntervalList,
    ModDict,
    Mod,
    setup_mod_dict,
    setup_mod_list,
    setup_interval_list,
    ModList,
)
from ..constants import (
    AMBIGUOUS_AMINO_ACIDS,
    MASS_AMBIGUOUS_AMINO_ACIDS,
    IonType,
    ModType,
    ModTypeLiteral,
    get_mod_type,
    get_mods,
)

from ..property import SequencePropertyMixin


class ProFormaAnnotation(SequencePropertyMixin, DigestionMixin):

    def __init__(
        self,
        sequence: str | None = None,
        isotope_mods: ACCEPTED_MODLIST_INPUT_TYPES = None,
        static_mods: ACCEPTED_MODLIST_INPUT_TYPES = None,
        labile_mods: ACCEPTED_MODLIST_INPUT_TYPES = None,
        unknown_mods: ACCEPTED_MODLIST_INPUT_TYPES = None,
        nterm_mods: ACCEPTED_MODLIST_INPUT_TYPES = None,
        cterm_mods: ACCEPTED_MODLIST_INPUT_TYPES = None,
        charge_adducts: ACCEPTED_MODLIST_INPUT_TYPES = None,
        internal_mods: ACCEPTED_MODDICT_INPUT_TYPES = None,
        intervals: ACCEPTED_INTERVALLIST_INPUT_TYPES = None,
        charge: int | None = None,
    ) -> None:

        self._sequence: str = sequence if sequence is not None else ""
        self._isotope_mod_list: ModList = setup_mod_list(isotope_mods)
        self._static_mod_list: ModList = setup_mod_list(static_mods)
        self._labile_mod_list: ModList = setup_mod_list(labile_mods)
        self._unknown_mod_list: ModList = setup_mod_list(unknown_mods)
        self._nterm_mod_list: ModList = setup_mod_list(nterm_mods)
        self._cterm_mod_list: ModList = setup_mod_list(cterm_mods)
        self._adduct_mod_list: ModList = setup_mod_list(charge_adducts)
        self._internal_mod_dict: ModDict = setup_mod_dict(internal_mods)
        self._interval_list: IntervalList = setup_interval_list(intervals)
        self._charge: int | None = charge


    @staticmethod
    def parse(sequence: str) -> "ProFormaAnnotation":
        """
        Parse a ProForma sequence string into a ProFormaAnnotation object.

        :param sequence: The ProForma sequence string to parse.
        :type sequence: str

        :return: A ProFormaAnnotation object representing the parsed sequence.
        :rtype: ProFormaAnnotation
        """

        # Implementation goes here
        for prof_parser, connection in ProFormaParser(sequence).parse():

            if connection is not None:
                raise ValueError(f"Unexpected connection value: {connection}")
            
            return ProFormaAnnotation(
                sequence="".join(prof_parser.amino_acids),
                isotope_mods=prof_parser.isotope_mods.copy(),
                static_mods=prof_parser.static_mods.copy(),
                labile_mods=prof_parser.labile_mods.copy(),
                unknown_mods=prof_parser.unknown_mods.copy(),
                nterm_mods=prof_parser.nterm_mods.copy(),
                cterm_mods=prof_parser.cterm_mods.copy(),
                internal_mods=prof_parser.internal_mods.copy(),
                intervals=prof_parser.intervals.copy(),
                charge=prof_parser.charge,
                charge_adducts=prof_parser.charge_adducts.copy(),
            )
        # If no annotations found, raise an error
        raise ValueError(f"Invalid ProForma sequence: {sequence}")
    """
    Magic Methods
    """

    def __repr__(self) -> str:
        """
        Only shows non-None types
        """
        seq = f"ProFormaAnnotationBase(sequence={self.sequence}"

        if self.has_isotope_mods:
            seq += f", {ModType.ISOTOPE}={self.isotope_mods}"
        if self.has_static_mods:
            seq += f", {ModType.STATIC}={self.static_mods}"
        if self.has_labile_mods:
            seq += f", {ModType.LABILE}={self.labile_mods}"
        if self.has_unknown_mods:
            seq += f", {ModType.UNKNOWN}={self.unknown_mods}"
        if self.has_nterm_mods:
            seq += f", {ModType.NTERM}={self.nterm_mods}"
        if self.has_cterm_mods:
            seq += f", {ModType.CTERM}={self.cterm_mods}"
        if self.has_internal_mods:
            seq += f", {ModType.INTERNAL}={self.internal_mods}"
        if self.has_intervals:
            seq += f", {ModType.INTERVAL}={self.intervals}"
        if self.has_charge:
            seq += f", {ModType.CHARGE}={self.charge}"
        if self.has_charge_adducts:
            seq += f", {ModType.CHARGE_ADDUCTS}={self.charge_adducts}"
        seq += ")"

        return seq

    def __len__(self) -> int:
        return len(self._sequence)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ProFormaAnnotation):
            return NotImplemented

        return (
            self.sequence == other.sequence
            and self._labile_mod_list == other._labile_mod_list
            and self._unknown_mod_list == other._unknown_mod_list
            and self._nterm_mod_list == other._nterm_mod_list
            and self._cterm_mod_list == other._cterm_mod_list
            and self._adduct_mod_list == other._adduct_mod_list
            and self._isotope_mod_list == other._isotope_mod_list
            and self._static_mod_list == other._static_mod_list
            and self._internal_mod_dict == other._internal_mod_dict
            and self._interval_list == other._interval_list
            and self._charge == other._charge
        )

    def __hash__(self):
        return hash(
            (
                self.sequence,
                tuple(self._isotope_mod_list),
                tuple(self._static_mod_list),
                tuple(self._labile_mod_list),
                tuple(self._unknown_mod_list),
                tuple(self._nterm_mod_list),
                tuple(self._cterm_mod_list),
                frozenset((k, tuple(v)) for k, v in self._internal_mod_dict.items()),
                tuple(self._interval_list),
                self._charge,
                tuple(self._adduct_mod_list),
            )
        )

    """
    Has Methods:
    """

    @property
    def has_sequence(self) -> bool:
        return len(self._sequence) > 0

    @property
    def has_isotope_mods(self) -> bool:
        return self._isotope_mod_list.has_mods

    @property
    def has_static_mods(self) -> bool:
        return self._static_mod_list.has_mods

    @property
    def has_labile_mods(self) -> bool:
        return self._labile_mod_list.has_mods

    @property
    def has_unknown_mods(self) -> bool:
        return self._unknown_mod_list.has_mods

    @property
    def has_nterm_mods(self) -> bool:
        return self._nterm_mod_list.has_mods

    @property
    def has_cterm_mods(self) -> bool:
        return self._cterm_mod_list.has_mods

    @property
    def has_charge_adducts(self) -> bool:
        return self._adduct_mod_list.has_mods

    @property
    def has_internal_mods(self) -> bool:
        return self._internal_mod_dict.has_mods

    @property
    def has_intervals(self) -> bool:
        return self._interval_list.has_intervals

    @property
    def has_unambiguous_intervals(self) -> bool:
        return self._interval_list.has_unambiguous_intervals

    @property
    def has_ambiguous_intervals(self) -> bool:
        return self._interval_list.has_ambiguous_intervals

    @property
    def has_charge(self) -> bool:
        return self._charge is not None

    @property
    def has_mods(self) -> bool:
        return any(
            [
                self.has_isotope_mods,
                self.has_static_mods,
                self.has_labile_mods,
                self.has_unknown_mods,
                self.has_nterm_mods,
                self.has_cterm_mods,
                self.has_internal_mods,
                self.has_intervals,
                self.has_charge,
                self.has_charge_adducts,
            ]
        )

    """
    Properties
    """

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def stripped_sequence(self) -> str:
        return self.sequence

    @property
    def isotope_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if not self.has_isotope_mods:
            return None
        return self._isotope_mod_list.flatten()

    @property
    def static_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if not self.has_static_mods:
            return None
        return self._static_mod_list.flatten()

    @property
    def labile_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if not self.has_labile_mods:
            return None
        return self._labile_mod_list.flatten()

    @property
    def unknown_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if not self.has_unknown_mods:
            return None
        return self._unknown_mod_list.flatten()

    @property
    def nterm_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if not self.has_nterm_mods:
            return None
        return self._nterm_mod_list.flatten()

    @property
    def cterm_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if not self.has_cterm_mods:
            return None
        return self._cterm_mod_list.flatten()

    @property
    def charge_adducts(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if not self.has_charge_adducts:
            return None
        return self._adduct_mod_list.flatten()

    @property
    def internal_mods(self) -> dict[int, tuple[MOD_VALUE_TYPES, ...]] | None:
        if not self.has_internal_mods:
            return None
        return {k: v.flatten() for k, v in self._internal_mod_dict.items()}

    @property
    def intervals(self) -> tuple[ModInterval, ...] | None:
        if self.has_intervals is False:
            return None
        return self._interval_list.get_mod_intervals()

    @property
    def ambiguous_intervals(self) -> tuple[ModInterval, ...] | None:
        if self._interval_list.has_ambiguous_intervals is False:
            return None
        return self._interval_list.get_ambiguous_intervals().get_mod_intervals()

    @property
    def unambiguous_intervals(self) -> tuple[ModInterval, ...] | None:
        if self._interval_list.has_unambiguous_intervals is False:
            return None
        return self._interval_list.get_unambiguous_intervals().get_mod_intervals()

    @property
    def charge(self) -> int | None:
        return self._charge

    """
    Property Setters
    """

    @sequence.setter
    def sequence(self, value: str | None) -> None:
        if value is None:
            value = ""
        self._sequence = value

    @isotope_mods.setter
    def isotope_mods(self, value: ACCEPTED_MODLIST_INPUT_TYPES) -> None:
        self._isotope_mod_list = setup_mod_list(value)

    @static_mods.setter
    def static_mods(self, value: ACCEPTED_MODLIST_INPUT_TYPES) -> None:
        self._static_mod_list = setup_mod_list(value)

    @labile_mods.setter
    def labile_mods(self, value: ACCEPTED_MODLIST_INPUT_TYPES) -> None:
        self._labile_mod_list = setup_mod_list(value)

    @unknown_mods.setter
    def unknown_mods(self, value: ACCEPTED_MODLIST_INPUT_TYPES) -> None:
        self._unknown_mod_list = setup_mod_list(value)

    @nterm_mods.setter
    def nterm_mods(self, value: ACCEPTED_MODLIST_INPUT_TYPES) -> None:
        self._nterm_mod_list = setup_mod_list(value)

    @cterm_mods.setter
    def cterm_mods(self, value: ACCEPTED_MODLIST_INPUT_TYPES) -> None:
        self._cterm_mod_list = setup_mod_list(value)

    @charge_adducts.setter
    def charge_adducts(self, value: ACCEPTED_MODLIST_INPUT_TYPES) -> None:
        self._adduct_mod_list = setup_mod_list(value)

    @internal_mods.setter
    def internal_mods(self, value: ACCEPTED_MODDICT_INPUT_TYPES) -> None:
        self._internal_mod_dict = setup_mod_dict(value)

    @intervals.setter
    def intervals(self, value: ACCEPTED_INTERVALLIST_INPUT_TYPES) -> None:
        self._interval_list = setup_interval_list(value)

    @charge.setter
    def charge(self, value: int | None):
        self._charge = value

    """
    Pop Methods
    """

    def pop_isotope_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        value = self.isotope_mods
        self.isotope_mods = None
        return value

    def pop_static_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        value = self.static_mods
        self.static_mods = None
        return value

    def pop_labile_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        value = self.labile_mods
        self.labile_mods = None
        return value

    def pop_unknown_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        value = self.unknown_mods
        self.unknown_mods = None
        return value

    def pop_nterm_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        value = self.nterm_mods
        self.nterm_mods = None
        return value

    def pop_cterm_mods(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        value = self.cterm_mods
        self.cterm_mods = None
        return value

    def pop_charge_adducts(self) -> tuple[MOD_VALUE_TYPES, ...] | None:
        value = self.charge_adducts
        self.charge_adducts = None
        return value

    def pop_internal_mods(self) -> dict[int, tuple[MOD_VALUE_TYPES, ...]] | None:
        value = self.internal_mods
        self.internal_mods = None
        return value

    def pop_intervals(self) -> tuple[ModInterval, ...] | None:
        value = self.intervals
        self.intervals = None
        return value

    def pop_ambiguous_intervals(self) -> tuple[ModInterval, ...] | None:
        intervals = self._interval_list.pop_ambiguous_intervals()
        if len(intervals) == 0:
            return None
        return intervals.get_mod_intervals()

    def pop_unambiguous_intervals(self) -> tuple[ModInterval, ...] | None:
        intervals = self._interval_list.pop_unambiguous_intervals()
        if len(intervals) == 0:
            return None
        return intervals.get_mod_intervals()

    def pop_charge(self) -> int | None:
        value = self._charge
        self._charge = None
        return value

    def pop_mod_by_type(self, mod_type: ModTypeLiteral) -> Any:
        """
        Get the pop method for a given mod type.
        """

        mod_type_to_pop_method: dict[str, Callable[[], Any]] = {
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

        if mod_type not in mod_type_to_pop_method:
            raise ValueError(f"Invalid mod type: {mod_type}")

        return mod_type_to_pop_method[mod_type]()

    def pop_mods(
        self,
        mods: ModTypeLiteral | Iterable[ModTypeLiteral] | None = None,
    ) -> dict[str, Any]:
        """
        Pop all mods and return them in a dictionary
        """

        mod_types: list[ModType] = get_mods(mods)

        d: dict[str, Any] = {}
        for mod_type in mod_types:
            d[mod_type.value] = self.pop_mod_by_type(mod_type.value)

        return d

    def pop_internal_mod(self, index: int) -> tuple[Mod, ...] | None:
        """
        Pop internal mods at a specific index
        """

        if index not in self._internal_mod_dict:
            return None

        mod_list = self._internal_mod_dict.pop(index)

        return mod_list.flatten()

    """
    Set Methods - Replace existing modifications
    """

    def set_isotope_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().set_isotope_mods(mods, inplace=True)
        self.isotope_mods = mods
        return self

    def set_static_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().set_static_mods(mods, inplace=True)
        self.static_mods = mods
        return self

    def set_labile_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().set_labile_mods(mods, inplace=True)
        self.labile_mods = mods
        return self

    def set_unknown_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().set_unknown_mods(mods, inplace=True)
        self.unknown_mods = mods
        return self

    def set_nterm_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().set_nterm_mods(mods, inplace=True)
        self.nterm_mods = mods
        return self

    def set_cterm_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().set_cterm_mods(mods, inplace=True)
        self.cterm_mods = mods
        return self

    def set_charge_adducts(
        self, charge_adducts: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().set_charge_adducts(charge_adducts, inplace=True)
        self.charge_adducts = charge_adducts
        return self

    def set_internal_mods(
        self, mods: ACCEPTED_MODDICT_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().set_internal_mods(mods, inplace=True)
        self.internal_mods = mods
        return self

    def set_internal_mods_at_index(
        self, index: int, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().set_internal_mods_at_index(index, mods, inplace=True)
        self._internal_mod_dict[index] = mods
        return self

    def set_intervals(
        self, intervals: ACCEPTED_INTERVALLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().set_intervals(intervals, inplace=True)
        self.intervals = intervals
        return self

    def set_charge(self, charge: int | None, inplace: bool = True) -> Self:
        if not inplace:
            return self.copy().set_charge(charge, inplace=True)
        self.charge = charge
        return self


    def set_mod_by_type(
        self,
        mod: Any,
        mod_type: ModTypeLiteral,
        inplace: bool = True,
    ) -> Self:
        """Set a modification by type, replacing any existing mods of that type"""
        if not inplace:
            return self.copy().set_mod_by_type(mod, mod_type, inplace=True)

        mod_type_enum = get_mod_type(mod_type)

        mod_type_methods: dict[ModType, Callable[[], Self]] = {
            ModType.ISOTOPE: lambda: self.set_isotope_mods(mod, inplace=True),
            ModType.STATIC: lambda: self.set_static_mods(mod, inplace=True),
            ModType.LABILE: lambda: self.set_labile_mods(mod, inplace=True),
            ModType.UNKNOWN: lambda: self.set_unknown_mods(mod, inplace=True),
            ModType.NTERM: lambda: self.set_nterm_mods(mod, inplace=True),
            ModType.CTERM: lambda: self.set_cterm_mods(mod, inplace=True),
            ModType.INTERNAL: lambda: self.set_internal_mods(mod, inplace=True),
            ModType.INTERVAL: lambda: self.set_intervals(mod, inplace=True),
            ModType.CHARGE: lambda: self.set_charge(mod, inplace=True),
            ModType.CHARGE_ADDUCTS: lambda: self.set_charge_adducts(mod, inplace=True),
        }

        mod_type_methods[mod_type_enum]()
        return self

    """
    Append Methods - Add single modifications
    """

    def append_isotope_mod(self, mod: MODLIST_DATATYPE, inplace: bool = True) -> Self:
        if not inplace:
            return self.copy().append_isotope_mod(mod, inplace=True)
        self._isotope_mod_list.append(mod)
        return self

    def append_static_mod(self, mod: MODLIST_DATATYPE, inplace: bool = True) -> Self:
        if not inplace:
            return self.copy().append_static_mod(mod, inplace=True)
        self._static_mod_list.append(mod)
        return self

    def append_labile_mod(self, mod: MODLIST_DATATYPE, inplace: bool = True) -> Self:
        if not inplace:
            return self.copy().append_labile_mod(mod, inplace=True)
        self._labile_mod_list.append(mod)
        return self

    def append_unknown_mod(self, mod: MODLIST_DATATYPE, inplace: bool = True) -> Self:
        if not inplace:
            return self.copy().append_unknown_mod(mod, inplace=True)
        self._unknown_mod_list.append(mod)
        return self

    def append_nterm_mod(self, mod: MODLIST_DATATYPE, inplace: bool = True) -> Self:
        if not inplace:
            return self.copy().append_nterm_mod(mod, inplace=True)
        self._nterm_mod_list.append(mod)
        return self

    def append_cterm_mod(self, mod: MODLIST_DATATYPE, inplace: bool = True) -> Self:
        if not inplace:
            return self.copy().append_cterm_mod(mod, inplace=True)
        self._cterm_mod_list.append(mod)
        return self

    def append_charge_adduct(self, mod: MODLIST_DATATYPE, inplace: bool = True) -> Self:
        if not inplace:
            return self.copy().append_charge_adduct(mod, inplace=True)
        self._adduct_mod_list.append(mod)
        return self

    def append_internal_mod_at_index(
        self, index: int, mod: MODLIST_DATATYPE, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().append_internal_mod_at_index(index, mod, inplace=True)
        self._internal_mod_dict.append_at_key(index, mod)
        return self

    def append_interval(
        self, interval: ACCEPTED_INTERVAL_DATATYPE, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().append_interval(interval, inplace=True)
        self._interval_list.append(interval)
        return self

    """
    Extend Methods - Add multiple modifications
    """

    def extend_isotope_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().extend_isotope_mods(mods, inplace=True)
        if mods is not None:
            self._isotope_mod_list.extend(setup_mod_list(mods))
        return self

    def extend_static_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().extend_static_mods(mods, inplace=True)
        if mods is not None:
            self._static_mod_list.extend(setup_mod_list(mods))
        return self

    def extend_labile_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().extend_labile_mods(mods, inplace=True)
        if mods is not None:
            self._labile_mod_list.extend(setup_mod_list(mods))
        return self

    def extend_unknown_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().extend_unknown_mods(mods, inplace=True)
        if mods is not None:
            self._unknown_mod_list.extend(setup_mod_list(mods))
        return self

    def extend_nterm_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().extend_nterm_mods(mods, inplace=True)
        if mods is not None:
            self._nterm_mod_list.extend(setup_mod_list(mods))
        return self

    def extend_cterm_mods(
        self, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().extend_cterm_mods(mods, inplace=True)
        if mods is not None:
            self._cterm_mod_list.extend(setup_mod_list(mods))
        return self

    def extend_charge_adducts(
        self, charge_adducts: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().extend_charge_adducts(charge_adducts, inplace=True)
        if charge_adducts is not None:
            self._adduct_mod_list.extend(setup_mod_list(charge_adducts))
        return self

    def extend_internal_mods_at_index(
        self, index: int, mods: ACCEPTED_MODLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().extend_internal_mods_at_index(index, mods, inplace=True)
        if mods is not None:
            self._internal_mod_dict.extend_at_key(index, mods)
        return self

    def extend_intervals(
        self, intervals: ACCEPTED_INTERVALLIST_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().extend_intervals(intervals, inplace=True)
        if intervals is not None:
            self._interval_list.extend(setup_interval_list(intervals))
        return self

    """
    Update/Merge Methods - Combine with existing modifications
    """

    def update_internal_mods(
        self, mods: ACCEPTED_MODDICT_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        """Update internal mods dictionary (replaces existing keys)"""
        if not inplace:
            return self.copy().update_internal_mods(mods, inplace=True)
        if mods is not None:
            self._internal_mod_dict.update(mods)
        return self

    def merge_internal_mods(
        self, mods: ACCEPTED_MODDICT_INPUT_TYPES, inplace: bool = True
    ) -> Self:
        """Merge internal mods dictionary (combines existing keys)"""
        if not inplace:
            return self.copy().merge_internal_mods(mods, inplace=True)
        if mods is not None:
            self._internal_mod_dict.merge(mods)
        return self

    """
    Remove Methods
    """

    def remove_isotope_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_isotope_mods(True)
        _ = self.pop_isotope_mods()
        return self

    def remove_static_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_static_mods(True)
        _ = self.pop_static_mods()
        return self

    def remove_nterm_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_nterm_mods(True)
        _ = self.pop_nterm_mods()
        return self

    def remove_cterm_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_cterm_mods(True)
        _ = self.pop_cterm_mods()
        return self

    def remove_internal_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_internal_mods(True)
        _ = self.pop_internal_mods()
        return self

    def remove_intervals(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_intervals(True)
        _ = self.pop_intervals()
        return self

    def remove_charge(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_charge(True)
        _ = self.pop_charge()
        return self

    def remove_charge_adducts(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_charge_adducts(True)
        _ = self.pop_charge_adducts()
        return self

    def remove_mods(
        self,
        mods: ModTypeLiteral | Iterable[ModTypeLiteral] | None = None,
        inplace: bool = True,
    ) -> Self:
        if inplace is False:
            return self.copy().remove_mods(mods=mods, inplace=True)
        _ = self.pop_mods(mods=mods)
        return self

    def strip(self, inplace: bool = False) -> Self:
        return self.remove_mods(None, inplace=inplace)

    def get_isotope_mod_list(self) -> ModList:
        return self._isotope_mod_list

    def get_static_mod_list(self) -> ModList:
        return self._static_mod_list

    def get_labile_mod_list(self) -> ModList:
        return self._labile_mod_list

    def get_unknown_mod_list(self) -> ModList:
        return self._unknown_mod_list

    def get_nterm_mod_list(self) -> ModList:
        return self._nterm_mod_list

    def get_cterm_mod_list(self) -> ModList:

        return self._cterm_mod_list

    def get_charge_adduct_list(self) -> ModList:
        return self._adduct_mod_list

    def get_internal_mod_dict(self) -> ModDict:
        return self._internal_mod_dict

    def get_interval_list(self) -> IntervalList:
        return self._interval_list

    def contains_sequence_ambiguity(self) -> bool:
        return self.has_intervals or self.has_unknown_mods

    def contains_residue_ambiguity(self) -> bool:
        return len(self.get_residue_ambiguity_residues()) > 0

    def get_residue_ambiguity_residues(self) -> tuple[str, ...]:
        return tuple(aa for aa in self.sequence if aa in AMBIGUOUS_AMINO_ACIDS)

    def contains_mass_ambiguity(self) -> bool:
        return len(self.get_mass_ambiguity_residues()) > 0

    def get_mass_ambiguity_residues(self) -> tuple[str, ...]:
        return tuple(aa for aa in self.sequence if aa in MASS_AMBIGUOUS_AMINO_ACIDS)

    def filter_mods(
        self,
        mods: ModTypeLiteral | Iterable[ModTypeLiteral] | None = None,
        inplace: bool = True,
    ) -> Self:
        """
        Filter the modifications in the annotation based on the provided mod types.
        If inplace is False, a new ProFormaAnnotationBase object is returned with the filtered mods.
        If inplace is True, the current object is modified.
        """

        if inplace is False:
            # Create a copy of the annotation to modify
            return self.copy().filter_mods(mods=mods, inplace=True)

        mod_types_to_remove = {mod_type for mod_type in ModType} - set(get_mods(mods))

        if len(mod_types_to_remove) == 0:
            # If no mods to remove, return the annotation as is
            return self

        for mod_type in mod_types_to_remove:
            # Pop the mods of the specified type
            _ = self.pop_mod_by_type(mod_type.value)

        return self

    def copy(self, deep: bool = True) -> Self:
        """
        Returns a deep copy of the ProFormaAnnotationBase instance.
        """

        if not deep:
            return self.__class__(
                sequence=self._sequence,
                isotope_mods=self._isotope_mod_list,
                static_mods=self._static_mod_list,
                labile_mods=self._labile_mod_list,
                unknown_mods=self._unknown_mod_list,
                nterm_mods=self._nterm_mod_list,
                cterm_mods=self._cterm_mod_list,
                internal_mods=self._internal_mod_dict,
                intervals=self._interval_list,
                charge=self._charge,
                charge_adducts=self._adduct_mod_list,
            )
        return deepcopy(self)

    def to_dict(self) -> dict[str, Any]:
        """
        Returns a dictionary representation of the ProFormaAnnotationBase instance.
        """
        return {
            "sequence": self.sequence,
            ModType.ISOTOPE: copy.deepcopy(self.isotope_mods),
            ModType.STATIC: copy.deepcopy(self.static_mods),
            ModType.LABILE: copy.deepcopy(self.labile_mods),
            ModType.UNKNOWN: copy.deepcopy(self.unknown_mods),
            ModType.NTERM: copy.deepcopy(self.nterm_mods),
            ModType.CTERM: copy.deepcopy(self.cterm_mods),
            ModType.INTERVAL: copy.deepcopy(self.intervals),
            ModType.INTERNAL: copy.deepcopy(self.internal_mods),
            ModType.CHARGE: copy.deepcopy(self.charge),
            ModType.CHARGE_ADDUCTS: copy.deepcopy(self.charge_adducts),
        }


    def serialize(
        self,
        include_plus: bool = False,
        precision: int | None = None,
    ) -> str:
        return serialize_annotation(self, include_plus, precision)

    def condense_static_mods(self, inplace: bool = True) -> Self:
        """
        Condense static mods into internal mods.
        """
        return cast(Self, condense_static_mods(self, inplace=inplace))

    def count_residues(self, include_mods: bool = True) -> dict[str, int]:
        """
        Count the occurrences of each residue in the sequence.
        """
        return count_residues(self, include_mods=include_mods)

    def percent_residues(
        self, include_mods: bool = True, precision: int | None = None
    ) -> dict[str, float]:
        """
        Calculate the percentage of each residue in the sequence.
        """
        return percent_residues(self, include_mods=include_mods, precision=precision)

    def is_subsequence(
        self,
        other: Self,
        ignore_mods: bool = False,
        ignore_intervals: bool = True,
    ) -> bool:
        """
        Check if the annotation is a subsequence of another annotation.
        """
        return is_subsequence(
            self, other, ignore_mods=ignore_mods, ignore_intervals=ignore_intervals
        )

    def find_indices(
        self,
        other: Self,
        ignore_mods: bool = False,
        ignore_intervals: bool = True,
    ) -> list[int]:
        """
        Find all occurrences of the annotation in another annotation.
        """
        return find_indices(
            self, other, ignore_mods=ignore_mods, ignore_intervals=ignore_intervals
        )
    
    def coverage(
        self,
        annotations: Iterable[Self],
        accumulate: bool = False,
        ignore_mods: bool = False,
        ignore_ambiguity: bool = False,
    ) -> list[int]:
        return coverage(
            annotation=self,
            annotations=annotations,
            accumulate=accumulate,
            ignore_mods=ignore_mods,
            ignore_ambiguity=ignore_ambiguity,
        )

    def percent_coverage(
        self,
        annotations: Iterable[Self],
        accumulate: bool = False,
        ignore_mods: bool = False,
        ignore_ambiguity: bool = False,
    ) -> float:
        return percent_coverage(
            annotation=self,
            annotations=annotations,
            accumulate=accumulate,
            ignore_mods=ignore_mods,
            ignore_ambiguity=ignore_ambiguity,
        )

    # Add slice and split convenience methods delegating to module-level functions
    def slice(
        self,
        start: int | None,
        stop: int | None,
        inplace: bool = False,
        keep_terms: bool = False,
        keep_labile: bool = True,
    ) -> Self:
        """
        Return a sliced annotation (delegates to slice_annotation).
        """
        return cast(
            Self,
            slice_annotation(
                self,
                start=start,
                stop=stop,
                inplace=inplace,
                keep_terms=keep_terms,
                keep_labile=keep_labile,
            ),
        )

    def split(self) -> list[Self]:
        """
        Split the annotation into single-residue annotations (delegates to split_annotation).
        """
        return [cast(Self, a) for a in split_annotation(self)]

    def annotate_ambiguity(
        self,
        forward_coverage: list[int],
        reverse_coverage: list[int],
        mass_shift: Any | None = None,
        inplace: bool = False,
    ) -> Self:
        return cast(
            Self,
            annotate_ambiguity(
                self, forward_coverage, reverse_coverage, mass_shift, inplace
            ),
        )

    def permutations(self, size: int | None = None) -> Generator[Self, None, None]:
        for item in generate_permutations(self, size):
            yield cast(Self, item)

    def product(self, repeat: int | None = None) -> Generator[Self, None, None]:
        for item in generate_product(self, repeat):
            yield cast(Self, item)

    def combinations(self, r: int | None = None) -> Generator[Self, None, None]:
        for item in generate_combinations(self, r):
            yield cast(Self, item)

    def combinations_with_replacement(
        self, r: int | None = None
    ) -> Generator[Self, None, None]:
        for item in generate_combinations_with_replacement(self, r):
            yield cast(Self, item)

    def shift(
        self,
        n: int,
        inplace: bool = False,
        shift_intervals: bool = True,
        slice_intervals: bool = True,
    ) -> Self:
        return cast(
            Self, shift_annotation(self, n, inplace, shift_intervals, slice_intervals)
        )

    def shuffle(self, seed: Any = None, inplace: bool = False) -> Self:
        return cast(Self, shuffle_annotation(self, seed, inplace))

    def reverse(self, inplace: bool = False, swap_terms: bool = False) -> Self:
        return cast(Self, reverse_annotation(self, inplace, swap_terms))

    def sort(
        self,
        inplace: bool = False,
        key: Callable[[str], Any] | None = None,
        reverse: bool = False,
    ) -> Self:
        return cast(Self, sort_annotation(self, inplace, key, reverse))

    def sliding_windows(
        self,
        window_size: int,
        keep_terms: bool = False,
        slice_intervals: bool = True,
        keep_labile: bool = True,
        reverse: bool = False,
    ) -> Generator[Self, None, None]:
        for window in generate_sliding_windows(
            self, window_size, keep_terms, slice_intervals, keep_labile, reverse
        ):
            yield cast(Self, window)


    def mass(
        self,
        ion_type: str = IonType.PRECURSOR,
        monoisotopic: bool = True,
        isotope: int = 0,
        loss: float = 0.0,
        use_isotope_on_mods: bool = False,
        precision: int | None = None,
        inplace: bool = False,
    ) -> float:
        """
        Calculate the mass of the annotation.
        
        :param ion_type: The type of ion (default: "p" for proton)
        :param monoisotopic: Whether to use monoisotopic masses (default: True)
        :param isotope: Isotope number (default: 0)
        :param loss: Mass loss to apply (default: 0.0)
        :param use_isotope_on_mods: Whether to apply isotopes to modifications (default: False)
        :param precision: Mass precision (default: None)
        :param inplace: Whether to modify in place (default: False)
        :return: The calculated mass
        """
        return mass(
            self,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            loss=loss,
            use_isotope_on_mods=use_isotope_on_mods,
            precision=precision,
            inplace=inplace,
        )

    def mz(
        self,
        ion_type: str = IonType.PRECURSOR,
        monoisotopic: bool = True,
        isotope: int = 0,
        loss: float = 0.0,
        precision: int | None = None,
    ) -> float:
        """
        Calculate the m/z ratio of the annotation.
        
        :param ion_type: The type of ion (default: "p" for proton)
        :param monoisotopic: Whether to use monoisotopic masses (default: True)
        :param isotope: Isotope number (default: 0)
        :param loss: Mass loss to apply (default: 0.0)
        :param precision: Precision for the result (default: None)
        :return: The calculated m/z ratio
        """
        return mz(
            self,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            loss=loss,
            precision=precision,
        )

    def comp(
        self,
        ion_type: str = IonType.PRECURSOR,
        estimate_delta: bool = False,
        isotope: int = 0,
        use_isotope_on_mods: bool = False,
        inplace: bool = False,
    ) -> dict[str, int | float]:
        """
        Calculate the elemental composition of the annotation.
        
        :param ion_type: The type of ion (default: "p" for proton)
        :param estimate_delta: Whether to estimate delta mass composition (default: False)
        :param isotope: Isotope number (default: 0)
        :param use_isotope_on_mods: Whether to apply isotopes to modifications (default: False)
        :param inplace: Whether to modify in place (default: False)
        :return: Dictionary representing the elemental composition
        """
        return comp(
            self,
            ion_type=ion_type,
            estimate_delta=estimate_delta,
            isotope=isotope,
            use_isotope_on_mods=use_isotope_on_mods,
            inplace=inplace,
        )