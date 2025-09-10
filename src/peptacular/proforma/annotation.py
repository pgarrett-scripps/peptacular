from collections.abc import Callable, Generator
import copy
from copy import deepcopy
from typing import Any, Iterable, Mapping, Self, Sequence, cast

from peptacular.proforma.dclasses.interval import Interval

from ..fragment import FragmenterMixin

from ..util import parse_static_mods

from .serializer import serialize_annotation
from .parser import ProFormaParser
from .ambiguity import annotate_ambiguity, condense_ambiguity_to_xnotation
from .manipulation import (
    condense_mods_to_intervals,
    condense_static_mods,
    is_subsequence,
    count_residues,
    modification_coverage,
    percent_residues,
    find_indices,
    coverage,
    percent_coverage,
)
from .slicing import (
    generate_sliding_windows,
    join_annotations,
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
    mz,
    condense_to_mass_mods,
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
    IonTypeLiteral,
    ModType,
    ModTypeLiteral,
    get_mod_type,
    get_mods,
)

from ..property import SequencePropertyMixin

from .mod_builder import build_mods


def _correct_mods(
    mod: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE],
) -> Iterable[MODLIST_DATATYPE]:
    return mod if isinstance(mod, Iterable) and not isinstance(mod, str) else [mod]


class ProFormaAnnotation(SequencePropertyMixin, DigestionMixin, FragmenterMixin):
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
        self._isotope_mod_list: ModList = setup_mod_list(
            isotope_mods, allow_dups=False, stackable=False
        )
        self._static_mod_list: ModList = setup_mod_list(
            static_mods, allow_dups=True, stackable=True
        )
        self._labile_mod_list: ModList = setup_mod_list(
            labile_mods, allow_dups=True, stackable=True
        )
        self._unknown_mod_list: ModList = setup_mod_list(
            unknown_mods, allow_dups=True, stackable=True
        )
        self._nterm_mod_list: ModList = setup_mod_list(
            nterm_mods, allow_dups=True, stackable=True
        )
        self._cterm_mod_list: ModList = setup_mod_list(
            cterm_mods, allow_dups=True, stackable=False
        )
        self._adduct_mod_list: ModList = setup_mod_list(
            charge_adducts, allow_dups=True, stackable=False
        )
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
        annots: list[ProFormaAnnotation] = []
        connections: list[bool | None] = []
        for prof_parser, connection in ProFormaParser(sequence).parse():
            if connection is True:
                raise ValueError(f"Unexpected connection value: {connection}")

            annot = ProFormaAnnotation(
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
            annots.append(annot)
            connections.append(connection)

        if len(annots) > 1:
            raise ValueError(f"Multiple annotations found: {len(annots)}")
        if len(annots) == 0:
            raise ValueError(f"Invalid ProForma sequence: {sequence}")

        # If all checks pass, return the single annotation
        return annots[0]

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
    def has_interval_mods(self) -> bool:
        return any(interval.has_mods for interval in self._interval_list)

    @property
    def has_unambiguous_intervals(self) -> bool:
        return self._interval_list.has_unambiguous_intervals

    @property
    def has_ambiguous_intervals(self) -> bool:
        return self._interval_list.has_ambiguous_intervals

    @property
    def has_charge(self) -> bool:
        return self._charge is not None

    def _has_mods_by_type(self, mod_type: ModType) -> bool:
        match mod_type:
            case ModType.ISOTOPE:
                return self.has_isotope_mods
            case ModType.STATIC:
                return self.has_static_mods
            case ModType.LABILE:
                return self.has_labile_mods
            case ModType.UNKNOWN:
                return self.has_unknown_mods
            case ModType.NTERM:
                return self.has_nterm_mods
            case ModType.CTERM:
                return self.has_cterm_mods
            case ModType.INTERNAL:
                return self.has_internal_mods
            case ModType.INTERVAL:
                return self.has_intervals
            case ModType.CHARGE:
                return self.has_charge
            case ModType.CHARGE_ADDUCTS:
                return self.has_charge_adducts
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

    def has_mods(
        self,
        mod_types: (
            Iterable[ModTypeLiteral | ModType] | ModType | ModTypeLiteral | None
        ) = None,
    ) -> bool:
        mod_enums = get_mods(mod_types)
        return any(self._has_mods_by_type(mod_enum) for mod_enum in mod_enums)

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

    def _get_mod_by_type(self, mod_type: ModType) -> Any:
        match mod_type:
            case ModType.NTERM:
                return self.nterm_mods
            case ModType.CTERM:
                return self.cterm_mods
            case ModType.ISOTOPE:
                return self.isotope_mods
            case ModType.STATIC:
                return self.static_mods
            case ModType.LABILE:
                return self.labile_mods
            case ModType.UNKNOWN:
                return self.unknown_mods
            case ModType.INTERVAL:
                return self.intervals
            case ModType.INTERNAL:
                return self.internal_mods
            case ModType.CHARGE:
                return self.charge
            case ModType.CHARGE_ADDUCTS:
                return self.charge_adducts
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

    @property
    def mods(self) -> dict[str, Any]:
        mod_enums = get_mods(None)
        return {
            mod_enum.value: self._get_mod_by_type(mod_enum) for mod_enum in mod_enums
        }

    def get_mods(
        self,
        mod_types: (
            Iterable[ModType | ModTypeLiteral] | ModType | ModTypeLiteral | None
        ) = None,
    ) -> dict[str, Any]:
        mod_types = get_mods(mod_types)

        return {
            mod_type.value: self._get_mod_by_type(mod_type) for mod_type in mod_types
        }

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

    @mods.setter
    def mods(self, value: Mapping[ModType | ModTypeLiteral, Any] | None) -> None:
        if value is None:
            mod_enums = get_mods(None)
            for mod_enum in mod_enums:
                self._set_mod_by_type(None, mod_enum)
        else:
            for mod_type, mods in value.items():
                mod_enum = get_mod_type(mod_type)
                self._set_mod_by_type(mods, mod_enum)

    """
    Pop Methods
    """

    def pop_isotope_mods(
        self, inplace: bool = True
    ) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if inplace is False:
            return self.copy().pop_isotope_mods(inplace=True)

        value = self.isotope_mods
        self.isotope_mods = None
        return value

    def pop_static_mods(
        self, inplace: bool = True
    ) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if inplace is False:
            return self.copy().pop_static_mods(inplace=True)

        value = self.static_mods
        self.static_mods = None
        return value

    def pop_labile_mods(
        self, inplace: bool = True
    ) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if inplace is False:
            return self.copy().pop_labile_mods(inplace=True)

        value = self.labile_mods
        self.labile_mods = None
        return value

    def pop_unknown_mods(
        self, inplace: bool = True
    ) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if inplace is False:
            return self.copy().pop_unknown_mods(inplace=True)

        value = self.unknown_mods
        self.unknown_mods = None
        return value

    def pop_nterm_mods(
        self, inplace: bool = True
    ) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if inplace is False:
            return self.copy().pop_nterm_mods(inplace=True)

        value = self.nterm_mods
        self.nterm_mods = None
        return value

    def pop_cterm_mods(
        self, inplace: bool = True
    ) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if inplace is False:
            return self.copy().pop_cterm_mods(inplace=True)

        value = self.cterm_mods
        self.cterm_mods = None
        return value

    def pop_charge_adducts(
        self, inplace: bool = True
    ) -> tuple[MOD_VALUE_TYPES, ...] | None:
        if inplace is False:
            return self.copy().pop_charge_adducts(inplace=True)

        value = self.charge_adducts
        self.charge_adducts = None
        return value

    def pop_internal_mods(
        self, inplace: bool = True
    ) -> dict[int, tuple[MOD_VALUE_TYPES, ...]] | None:
        if inplace is False:
            return self.copy().pop_internal_mods(inplace=True)

        value = self.internal_mods
        self.internal_mods = None
        return value

    def pop_intervals(self, inplace: bool = True) -> tuple[ModInterval, ...] | None:
        if inplace is False:
            return self.copy().pop_intervals(inplace=True)

        value = self.intervals
        self.intervals = None
        return value

    def pop_ambiguous_intervals(
        self, inplace: bool = True
    ) -> tuple[ModInterval, ...] | None:
        if inplace is False:
            return self.copy().pop_ambiguous_intervals(inplace=True)

        intervals = self._interval_list.pop_ambiguous_intervals()
        if len(intervals) == 0:
            return None
        return intervals.get_mod_intervals()

    def pop_unambiguous_intervals(
        self, inplace: bool = True
    ) -> tuple[ModInterval, ...] | None:
        if inplace is False:
            return self.copy().pop_unambiguous_intervals(inplace=True)

        intervals = self._interval_list.pop_unambiguous_intervals()
        if len(intervals) == 0:
            return None
        return intervals.get_mod_intervals()

    def pop_charge(self, inplace: bool = True) -> int | None:
        if inplace is False:
            return self.copy().pop_charge(inplace=True)

        value = self._charge
        self._charge = None
        return value

    def _pop_mod_by_type(self, mod_type: ModType) -> Any:
        match mod_type:
            case ModType.ISOTOPE:
                return self.pop_isotope_mods(inplace=True)
            case ModType.STATIC:
                return self.pop_static_mods(inplace=True)
            case ModType.LABILE:
                return self.pop_labile_mods(inplace=True)
            case ModType.UNKNOWN:
                return self.pop_unknown_mods(inplace=True)
            case ModType.NTERM:
                return self.pop_nterm_mods(inplace=True)
            case ModType.CTERM:
                return self.pop_cterm_mods(inplace=True)
            case ModType.INTERNAL:
                return self.pop_internal_mods(inplace=True)
            case ModType.INTERVAL:
                return self.pop_intervals(inplace=True)
            case ModType.CHARGE_ADDUCTS:
                return self.pop_charge_adducts(inplace=True)
            case ModType.CHARGE:
                return self.pop_charge(inplace=True)
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

    def pop_mods(
        self,
        mod_types: (
            ModTypeLiteral | Iterable[ModTypeLiteral | ModType] | ModType | None
        ) = None,
        inplace: bool = True,
    ) -> dict[str, Any]:
        if inplace is False:
            return self.copy().pop_mods(mod_types=mod_types, inplace=True)

        mod_enums: list[ModType] = get_mods(mod_types)

        d: dict[str, Any] = {}
        for mod_enum in mod_enums:
            d[mod_enum.value] = self._pop_mod_by_type(mod_enum)

        return d

    def pop_internal_mod_at_index(self, index: int) -> tuple[Mod, ...] | None:
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

    def _set_mod_by_type(
        self,
        value: Any,
        mod_type: ModType,
    ) -> Self:
        match mod_type:
            case ModType.ISOTOPE:
                self.isotope_mods = value
            case ModType.STATIC:
                self.static_mods = value
            case ModType.LABILE:
                self.labile_mods = value
            case ModType.UNKNOWN:
                self.unknown_mods = value
            case ModType.NTERM:
                self.nterm_mods = value
            case ModType.CTERM:
                self.cterm_mods = value
            case ModType.INTERNAL:
                self.internal_mods = value
            case ModType.INTERVAL:
                self.intervals = value
            case ModType.CHARGE:
                self.charge = value
            case ModType.CHARGE_ADDUCTS:
                self.charge_adducts = value
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

        return self

    def set_mods(
        self,
        mods: Mapping[ModType | ModTypeLiteral | int, Any] | None,
        inplace: bool = True,
    ) -> Self:
        """Set a modification by type, replacing any existing mods of that type"""

        if not inplace:
            return self.copy().set_mods(mods=mods, inplace=True)

        if mods is None:
            self.remove_mods(inplace=True)
            return self

        for mod_type, mod_value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise IndexError(
                        f"Internal modification index out of range: {mod_type}"
                    )
                self.set_internal_mods_at_index(mod_type, mod_value, inplace=True)
                continue

            self._set_mod_by_type(mod_value, ModType(mod_type))

        return self

    """
    Add Methods (For old compatibility)
    """

    def add_isotope_mods(
        self,
        mods: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_isotope_mods(mods, inplace=True, append=append)
        if append:
            self.extend_isotope_mods(_correct_mods(mods))
        else:
            self.set_isotope_mods(_correct_mods(mods))
        return self

    def add_static_mods(
        self,
        mods: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_static_mods(mods, inplace=True, append=append)
        if append:
            self.extend_static_mods(_correct_mods(mods))
        else:
            self.set_static_mods(_correct_mods(mods))
        return self

    def add_labile_mods(
        self,
        mods: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_labile_mods(mods, inplace=True, append=append)
        if append:
            self.extend_labile_mods(_correct_mods(mods))
        else:
            self.set_labile_mods(_correct_mods(mods))
        return self

    def add_unknown_mods(
        self,
        mods: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_unknown_mods(mods, inplace=True, append=append)
        if append:
            self.extend_unknown_mods(_correct_mods(mods))
        else:
            self.set_unknown_mods(_correct_mods(mods))
        return self

    def add_nterm_mods(
        self,
        mods: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_nterm_mods(mods, inplace=True, append=append)
        if append:
            self.extend_nterm_mods(_correct_mods(mods))
        else:
            self.set_nterm_mods(_correct_mods(mods))
        return self

    def add_cterm_mods(
        self,
        mods: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_cterm_mods(mods, inplace=True, append=append)
        if append:
            self.extend_cterm_mods(_correct_mods(mods))
        else:
            self.set_cterm_mods(_correct_mods(mods))
        return self

    def add_charge_adducts(
        self,
        mods: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_charge_adducts(mods, inplace=True, append=append)
        if append:
            self.extend_charge_adducts(_correct_mods(mods))
        else:
            self.set_charge_adducts(_correct_mods(mods))
        return self

    def add_internal_mods_at_index(
        self,
        index: int,
        mods: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_internal_mods_at_index(
                index, mods, inplace=True, append=append
            )
        if append:
            self.extend_internal_mods_at_index(index, _correct_mods(mods))
        else:
            self.set_internal_mods_at_index(index, _correct_mods(mods))
        return self

    def add_intervals(
        self,
        intervals: ACCEPTED_INTERVAL_DATATYPE | Iterable[ACCEPTED_INTERVAL_DATATYPE],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_intervals(intervals, inplace=True, append=append)

        # Handle single interval vs iterable of intervals
        if isinstance(intervals, (ModInterval, Interval)) or (
            isinstance(intervals, tuple) and len(intervals) == 4
        ):
            if append:
                self.append_interval(intervals)  # type: ignore
            else:
                self.set_intervals([intervals])  # type: ignore
        else:
            if append:
                self.extend_intervals(intervals)
            else:
                self.set_intervals(intervals)
        return self

    def add_internal_mods(
        self,
        mods: Mapping[int, MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE]],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_internal_mods(mods, inplace=True, append=append)
        for index, mod in mods.items():
            self.add_internal_mods_at_index(index, mod, append=append)
        return self

    def _add_mods_by_type(
        self, value: Any, mod_type: ModType, inplace: bool = True, append: bool = True
    ) -> Self:
        if not inplace:
            return self.copy()._add_mods_by_type(value, mod_type, inplace=True)

        match mod_type:
            case ModType.ISOTOPE:
                self.add_isotope_mods(value, inplace=inplace, append=append)
            case ModType.STATIC:
                self.add_static_mods(value, inplace=inplace, append=append)
            case ModType.LABILE:
                self.add_labile_mods(value, inplace=inplace, append=append)
            case ModType.UNKNOWN:
                self.add_unknown_mods(value, inplace=inplace, append=append)
            case ModType.NTERM:
                self.add_nterm_mods(value, inplace=inplace, append=append)
            case ModType.CTERM:
                self.add_cterm_mods(value, inplace=inplace, append=append)
            case ModType.INTERNAL:
                self.add_internal_mods(value, inplace=inplace, append=append)
            case ModType.INTERVAL:
                self.add_intervals(value, inplace=inplace, append=append)
            case ModType.CHARGE:
                self.set_charge(value)
            case ModType.CHARGE_ADDUCTS:
                self.add_charge_adducts(value, inplace=inplace, append=append)
            case _:
                raise ValueError(f"Unknown modification type: {mod_type}")

        return self

    def add_mods(
        self,
        mods: Mapping[ModType | ModTypeLiteral | int, Any],
        inplace: bool = True,
        append: bool = True,
    ) -> Self:
        if not inplace:
            return self.copy().add_mods(mods, inplace=True)

        for mod_type, value in mods.items():
            if isinstance(mod_type, int):
                self.add_internal_mods_at_index(
                    mod_type, value, inplace=inplace, append=append
                )
                continue

            self._add_mods_by_type(
                value, ModType(mod_type), inplace=inplace, append=append
            )

        return self

    """
    Append Methods
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

    def _append_by_type(
        self, value: Any, mod_type: ModType, inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy()._append_by_type(value, mod_type, inplace=True)

        match mod_type:
            case ModType.ISOTOPE:
                self.append_isotope_mod(value)
            case ModType.STATIC:
                self.append_static_mod(value)
            case ModType.LABILE:
                self.append_labile_mod(value)
            case ModType.UNKNOWN:
                self.append_unknown_mod(value)
            case ModType.NTERM:
                self.append_nterm_mod(value)
            case ModType.CTERM:
                self.append_cterm_mod(value)
            case ModType.INTERNAL:
                for key, val in value.items():
                    self.append_internal_mod_at_index(key, val)
            case ModType.INTERVAL:
                self.append_interval(value)
            case ModType.CHARGE:
                self.set_charge(value)
            case ModType.CHARGE_ADDUCTS:
                self.extend_charge_adducts(value)
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

        return self

    def append_mods(
        self, mods: Mapping[ModType | ModTypeLiteral | int, Any], inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().append_mods(mods, inplace=True)

        for mod_type, value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise IndexError(
                        f"Internal modification index out of range: {mod_type}"
                    )
                self.append_internal_mod_at_index(mod_type, value, inplace=True)
                continue

            self._append_by_type(value, ModType(mod_type))

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

    def _extend_by_type(self, value: Any, mod_type: ModType) -> Self:
        match mod_type:
            case ModType.ISOTOPE:
                self.extend_isotope_mods(value)
            case ModType.STATIC:
                self.extend_static_mods(value)
            case ModType.LABILE:
                self.extend_labile_mods(value)
            case ModType.UNKNOWN:
                self.extend_unknown_mods(value)
            case ModType.NTERM:
                self.extend_nterm_mods(value)
            case ModType.CTERM:
                self.extend_cterm_mods(value)
            case ModType.INTERNAL:
                for index, mod in value.items():
                    self.extend_internal_mods_at_index(index, mod)
            case ModType.INTERVAL:
                self.extend_intervals(value)
            case ModType.CHARGE_ADDUCTS:
                self.extend_charge_adducts(value)
            case ModType.CHARGE:
                raise NotImplementedError("Extending charge not supported.")
            case _:
                raise NotImplementedError(f"Appending {mod_type} not supported.")

        return self

    def extend_mods(
        self, mods: Mapping[ModType | ModTypeLiteral | int, Any], inplace: bool = True
    ) -> Self:
        if not inplace:
            return self.copy().extend_mods(mods, inplace=True)

        for mod_type, value in mods.items():
            if isinstance(mod_type, int):
                if mod_type < 0 or mod_type >= len(self.sequence):
                    raise IndexError(
                        f"Internal modification index out of range: {mod_type}"
                    )
                self.extend_internal_mods_at_index(mod_type, value, inplace=True)
                continue
            self._extend_by_type(value, ModType(mod_type))

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

    def remove_isotope_mods(self, inplace: bool = False) -> Self:
        if inplace is False:
            return self.copy().remove_isotope_mods(True)
        self.isotope_mods = None
        return self

    def remove_static_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_static_mods(True)
        self.static_mods = None
        return self

    def remove_nterm_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_nterm_mods(True)
        self.nterm_mods = None
        return self

    def remove_cterm_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_cterm_mods(True)
        self.cterm_mods = None
        return self

    def remove_labile_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_labile_mods(True)
        self.labile_mods = None
        return self

    def remove_unknown_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_unknown_mods(True)
        self.unknown_mods = None
        return self

    def remove_internal_mods(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_internal_mods(True)
        self.internal_mods = None
        return self

    def remove_intervals(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_intervals(True)
        self.intervals = None
        return self

    def remove_charge(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_charge(True)
        self.charge = None
        return self

    def remove_charge_adducts(self, inplace: bool = True) -> Self:
        if inplace is False:
            return self.copy().remove_charge_adducts(True)
        self.charge_adducts = None
        return self

    def _remove_mod_by_type(self, mod_type: ModType) -> None:
        match mod_type:
            case ModType.ISOTOPE:
                self.remove_isotope_mods(inplace=True)
            case ModType.STATIC:
                self.remove_static_mods(inplace=True)
            case ModType.LABILE:
                self.remove_labile_mods(inplace=True)
            case ModType.UNKNOWN:
                self.remove_unknown_mods(inplace=True)
            case ModType.NTERM:
                self.remove_nterm_mods(inplace=True)
            case ModType.CTERM:
                self.remove_cterm_mods(inplace=True)
            case ModType.INTERNAL:
                self.remove_internal_mods(inplace=True)
            case ModType.INTERVAL:
                self.remove_intervals(inplace=True)
            case ModType.CHARGE:
                self.remove_charge(inplace=True)
            case ModType.CHARGE_ADDUCTS:
                self.remove_charge_adducts(inplace=True)
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

    def remove_mods(
        self,
        mods: (
            ModTypeLiteral | ModType | Iterable[ModTypeLiteral | ModType] | None
        ) = None,
        inplace: bool = True,
    ) -> Self:
        if inplace is False:
            return self.copy().remove_mods(mods=mods, inplace=True)
        mod_enums = get_mods(mods)
        for mod_enum in mod_enums:
            self._remove_mod_by_type(mod_enum)
        return self

    def strip_mods(self, inplace: bool = False) -> Self:
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

    """
    Sort mods
    """

    def sort_isotope_mods(
        self, key: Callable[[Mod], Any] | None = None, reverse: bool = False
    ) -> Self:
        self._isotope_mod_list.sort(key=key, reverse=reverse)
        return self

    def sort_static_mods(
        self, key: Callable[[Mod], Any] | None = None, reverse: bool = False
    ) -> Self:
        self._static_mod_list.sort(key=key, reverse=reverse)
        return self

    def sort_labile_mods(
        self, key: Callable[[Mod], Any] | None = None, reverse: bool = False
    ) -> Self:
        self._labile_mod_list.sort(key=key, reverse=reverse)
        return self

    def sort_unknown_mods(
        self, key: Callable[[Mod], Any] | None = None, reverse: bool = False
    ) -> Self:
        self._unknown_mod_list.sort(key=key, reverse=reverse)
        return self

    def sort_nterm_mods(
        self, key: Callable[[Mod], Any] | None = None, reverse: bool = False
    ) -> Self:
        self._nterm_mod_list.sort(key=key, reverse=reverse)
        return self

    def sort_cterm_mods(
        self, key: Callable[[Mod], Any] | None = None, reverse: bool = False
    ) -> Self:
        self._cterm_mod_list.sort(key=key, reverse=reverse)
        return self

    def sort_internal_mods(
        self, key: Callable[[Mod], Any] | None = None, reverse: bool = False
    ) -> Self:
        for mod_list in self._internal_mod_dict.values():
            mod_list.sort(key=key, reverse=reverse)
        return self

    def sort_interval_mods(
        self, key: Callable[[Mod], Any] | None = None, reverse: bool = False
    ) -> Self:
        for interval in self._interval_list:
            interval.mods.sort(key=key, reverse=reverse)
        return self

    def sort_charge_adducts(
        self, key: Callable[[Mod], Any] | None = None, reverse: bool = False
    ) -> Self:
        self._adduct_mod_list.sort(key=key, reverse=reverse)
        return self

    def _sort_mods_by_type(
        self,
        mod_type: ModType,
        key: Callable[[Mod], Any] | None = None,
        reverse: bool = False,
    ) -> Self:
        for mod_type in ModType:
            match mod_type:
                case ModType.ISOTOPE:
                    self.sort_isotope_mods(key=key, reverse=reverse)
                case ModType.STATIC:
                    self.sort_static_mods(key=key, reverse=reverse)
                case ModType.LABILE:
                    self.sort_labile_mods(key=key, reverse=reverse)
                case ModType.UNKNOWN:
                    self.sort_unknown_mods(key=key, reverse=reverse)
                case ModType.NTERM:
                    self.sort_nterm_mods(key=key, reverse=reverse)
                case ModType.CTERM:
                    self.sort_cterm_mods(key=key, reverse=reverse)
                case ModType.INTERNAL:
                    self.sort_internal_mods(key=key, reverse=reverse)
                case ModType.INTERVAL:
                    self.sort_interval_mods(key=key, reverse=reverse)
                case ModType.CHARGE_ADDUCTS:
                    self.sort_charge_adducts(key=key, reverse=reverse)
                case ModType.CHARGE:
                    pass
                case _:
                    raise TypeError(f"Unknown mod type: {mod_type}")
        return self

    def sort_mods(
        self,
        mods: (
            ModTypeLiteral | ModType | Iterable[ModTypeLiteral | ModType] | None
        ) = None,
        key: Callable[[Mod], Any] | None = None,
        reverse: bool = False,
        inplace: bool = True,
    ) -> Self:
        if not inplace:
            return self.sort_mods(mods=mods, key=key, reverse=reverse, inplace=True)

        mod_types: list[ModType] = get_mods(mods)
        for mod_type in mod_types:
            self._sort_mods_by_type(mod_type, key=key, reverse=reverse)
        return self

    """
    add (Depriciated)
    """

    """
    Contains
    """

    def contains_isotope_mod(self, mod: MODLIST_DATATYPE) -> bool:
        return mod in self._isotope_mod_list

    def contains_static_mod(self, mod: MODLIST_DATATYPE) -> bool:
        return mod in self._static_mod_list

    def contains_labile_mod(self, mod: MODLIST_DATATYPE) -> bool:
        return mod in self._labile_mod_list

    def contains_unknown_mod(self, mod: MODLIST_DATATYPE) -> bool:
        return mod in self._unknown_mod_list

    def contains_nterm_mod(self, mod: MODLIST_DATATYPE) -> bool:
        return mod in self._nterm_mod_list

    def contains_cterm_mod(self, mod: MODLIST_DATATYPE) -> bool:
        return mod in self._cterm_mod_list

    def contains_internal_mod(self, mod: MODLIST_DATATYPE) -> bool:
        return any(mod in mod_list for mod_list in self._internal_mod_dict.values())

    def contains_interval_mod(self, mod: MODLIST_DATATYPE) -> bool:
        return any(mod in interval.mods for interval in self._interval_list)

    def contains_mod_at_index(self, mod: MODLIST_DATATYPE, index: int) -> bool:
        return mod in self._internal_mod_dict.get(index, ModList())

    def _contains_mod_by_type(self, mod: MODLIST_DATATYPE, mod_type: ModType) -> bool:
        match mod_type:
            case ModType.ISOTOPE:
                return self.contains_isotope_mod(mod)
            case ModType.STATIC:
                return self.contains_static_mod(mod)
            case ModType.LABILE:
                return self.contains_labile_mod(mod)
            case ModType.UNKNOWN:
                return self.contains_unknown_mod(mod)
            case ModType.NTERM:
                return self.contains_nterm_mod(mod)
            case ModType.CTERM:
                return self.contains_cterm_mod(mod)
            case ModType.INTERNAL:
                return self.contains_internal_mod(mod)
            case ModType.INTERVAL:
                return self.contains_interval_mod(mod)
            case ModType.CHARGE:
                return False
            case _:
                raise TypeError(f"Unknown mod type: {mod_type}")

    def contains_mod(
        self,
        mod: MODLIST_DATATYPE,
        mods: (
            ModTypeLiteral | ModType | Iterable[ModTypeLiteral | ModType] | None
        ) = None,
    ) -> bool:
        mod_types: list[ModType] = get_mods(mods)
        return any(self._contains_mod_by_type(mod, mod_type) for mod_type in mod_types)

    """
    Other
    """

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
        mods: (
            ModTypeLiteral | ModType | Iterable[ModTypeLiteral | ModType] | None
        ) = None,
        inplace: bool = True,
    ) -> Self:
        if inplace is False:
            return self.copy().filter_mods(mods=mods, inplace=True)

        mod_types_to_remove = {mod_type for mod_type in ModType} - set(get_mods(mods))

        if len(mod_types_to_remove) == 0:
            # If no mods to remove, return the annotation as is
            return self

        self.pop_mods(mod_types_to_remove)
        return self

    def copy(self, deep: bool = True) -> Self:
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

    def copy_from(self, other: Self, deep: bool = True) -> None:
        self._sequence = other._sequence
        self._isotope_mod_list = other._isotope_mod_list.copy(deep=deep)
        self._static_mod_list = other._static_mod_list.copy(deep=deep)
        self._labile_mod_list = other._labile_mod_list.copy(deep=deep)
        self._unknown_mod_list = other._unknown_mod_list.copy(deep=deep)
        self._nterm_mod_list = other._nterm_mod_list.copy(deep=deep)
        self._cterm_mod_list = other._cterm_mod_list.copy(deep=deep)
        self._internal_mod_dict = other._internal_mod_dict.copy(deep=deep)
        self._interval_list = other._interval_list.copy(deep=deep)
        self._charge = other._charge
        self._adduct_mod_list = other._adduct_mod_list.copy(deep=deep)

    def to_dict(self) -> dict[str, Any]:
        return {
            "sequence": self.sequence,
            ModType.ISOTOPE.value: copy.deepcopy(self.isotope_mods),
            ModType.STATIC.value: copy.deepcopy(self.static_mods),
            ModType.LABILE.value: copy.deepcopy(self.labile_mods),
            ModType.UNKNOWN.value: copy.deepcopy(self.unknown_mods),
            ModType.NTERM.value: copy.deepcopy(self.nterm_mods),
            ModType.CTERM.value: copy.deepcopy(self.cterm_mods),
            ModType.INTERVAL.value: copy.deepcopy(self.intervals),
            ModType.INTERNAL.value: copy.deepcopy(self.internal_mods),
            ModType.CHARGE.value: copy.deepcopy(self.charge),
            ModType.CHARGE_ADDUCTS.value: copy.deepcopy(self.charge_adducts),
        }

    @staticmethod
    def from_dict(data: dict[str, Any]) -> "ProFormaAnnotation":
        return ProFormaAnnotation(
            sequence=data.get("sequence"),
            isotope_mods=data.get(ModType.ISOTOPE.value),
            static_mods=data.get(ModType.STATIC.value),
            labile_mods=data.get(ModType.LABILE.value),
            unknown_mods=data.get(ModType.UNKNOWN.value),
            nterm_mods=data.get(ModType.NTERM.value),
            cterm_mods=data.get(ModType.CTERM.value),
            internal_mods=data.get(ModType.INTERNAL.value),
            intervals=data.get(ModType.INTERVAL.value),
            charge=data.get(ModType.CHARGE.value),
            charge_adducts=data.get(ModType.CHARGE_ADDUCTS.value),
        )

    def serialize(
        self,
        include_plus: bool = False,
        precision: int | None = None,
    ) -> str:
        return serialize_annotation(self, include_plus, precision)

    def condense_static_mods(self, inplace: bool = True) -> Self:
        return cast(Self, condense_static_mods(self, inplace=inplace))

    def get_static_mod_dict(self) -> ModDict:
        static_mod_dict: ModDict = ModDict()
        for static_aa, mods in parse_static_mods(
            self.get_static_mod_list().data
        ).items():
            for i, aa in enumerate(self.sequence):
                if aa == static_aa:
                    static_mod_dict[i].extend(mods)
        return static_mod_dict

    def count_residues(self, include_mods: bool = True) -> dict[str, int]:
        return count_residues(self, include_mods=include_mods)

    def percent_residues(
        self, include_mods: bool = True, precision: int | None = None
    ) -> dict[str, float]:
        return percent_residues(self, include_mods=include_mods, precision=precision)

    def is_subsequence(
        self,
        other: Self,
        ignore_mods: bool = False,
        ignore_intervals: bool = True,
    ) -> bool:
        return is_subsequence(
            self, other, ignore_mods=ignore_mods, ignore_intervals=ignore_intervals
        )

    def find_indices(
        self,
        other: Self,
        ignore_mods: bool = False,
        ignore_intervals: bool = True,
    ) -> list[int]:
        return find_indices(
            self, other, ignore_mods=ignore_mods, ignore_intervals=ignore_intervals
        )

    def condense_mods_to_intervals(self, inplace: bool = True) -> Self:
        return cast(Self, condense_mods_to_intervals(self, inplace=inplace))

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

    def modification_coverage(
        self,
        annotations: Iterable[Self],
        ignore_ambiguity: bool = False,
        accumulate: bool = False,
    ) -> dict[int, int]:
        return modification_coverage(
            annotation=self,
            annotations=annotations,
            ignore_ambiguity=ignore_ambiguity,
            accumulate=accumulate,
        )

    def slice(
        self,
        start: int | None,
        stop: int | None,
        inplace: bool = False,
    ) -> Self:
        return cast(
            Self,
            slice_annotation(
                self,
                start=start,
                stop=stop,
                inplace=inplace,
            ),
        )

    def split(self) -> list[Self]:
        return [cast(Self, a) for a in split_annotation(self)]

    @staticmethod
    def join(annotations: Sequence["ProFormaAnnotation"]) -> "ProFormaAnnotation":
        return join_annotations(annotations)

    def annotate_ambiguity(
        self,
        forward_coverage: list[int],
        reverse_coverage: list[int],
        mass_shift: Any | None = None,
        add_mods_to_intervals: bool = False,
        sort_mods: bool = True,
        inplace: bool = False,
    ) -> Self:
        return cast(
            Self,
            annotate_ambiguity(
                self,
                forward_coverage=forward_coverage,
                reverse_coverage=reverse_coverage,
                mass_shift=mass_shift,
                add_mods_to_intervals=add_mods_to_intervals,
                sort_mods=sort_mods,
                inplace=inplace,
            ),
        )

    def condense_ambiguity_to_xnotation(self, inplace: bool = True) -> Self:
        return cast(Self, condense_ambiguity_to_xnotation(self, inplace=inplace))

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
    ) -> Self:
        return cast(Self, shift_annotation(self, n, inplace))

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
        reverse: bool = False,
    ) -> Generator[Self, None, None]:
        for window in generate_sliding_windows(self, window_size, reverse):
            yield cast(Self, window)

    def mass(
        self,
        ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
        monoisotopic: bool = True,
        isotope: int = 0,
        loss: float = 0.0,
        use_isotope_on_mods: bool = False,
        precision: int | None = None,
    ) -> float:
        return mass(
            self,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            loss=loss,
            use_isotope_on_mods=use_isotope_on_mods,
            precision=precision,
        )

    def mz(
        self,
        ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
        monoisotopic: bool = True,
        isotope: int = 0,
        loss: float = 0.0,
        precision: int | None = None,
        use_isotope_on_mods: bool = False,
    ) -> float:
        return mz(
            self,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            loss=loss,
            precision=precision,
            use_isotope_on_mods=use_isotope_on_mods,
        )

    def comp(
        self,
        ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
        isotope: int = 0,
        use_isotope_on_mods: bool = False,
    ) -> dict[str, int | float]:
        return comp(
            self,
            ion_type=ion_type,
            isotope=isotope,
            use_isotope_on_mods=use_isotope_on_mods,
        )

    def condense_to_delta_mass(
        self,
        include_plus: bool = False,
        inplace: bool = True,
        use_isotope_on_mods: bool = False,
    ) -> Self:
        """
        Condense all modifications to their mass equivalents.
        """
        return cast(
            Self,
            condense_to_mass_mods(
                self,
                include_plus=include_plus,
                inplace=inplace,
                use_isotope_on_mods=use_isotope_on_mods,
            ),
        )

    def build_mods(
        self,
        nterm_static: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
        cterm_static: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
        internal_static: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
        labile_static: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
        nterm_variable: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
        cterm_variable: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
        internal_variable: (
            Mapping[str, Iterable[str | float | int | Mod]] | None
        ) = None,
        labile_variable: Mapping[str, Iterable[str | float | int | Mod]] | None = None,
        max_variable_mods: int = 2,
        use_regex: bool = False,
    ) -> Generator[Self, None, None]:
        """
        Build all modifications from intervals and mass shifts.
        """
        for annot in build_mods(
            self,
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
        ):
            yield cast(Self, annot)
