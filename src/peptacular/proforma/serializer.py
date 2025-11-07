from __future__ import annotations

from typing import TYPE_CHECKING, Any

from .characters import ProformaChar as PC
from .dclasses.intervallist import IntervalList
from .dclasses.modlist import ModList

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


# Cache enum values at module level to avoid repeated dictionary lookups
_MOD_BRACKETS = PC.MOD_START.value + PC.MOD_END.value
_LABILE_BRACKETS = PC.LABILE_MOD_START.value + PC.LABILE_MOD_END.value
_STATIC_BRACKETS = PC.STATIC_MOD_START.value + PC.STATIC_MOD_END.value
_ISOTOPE_BRACKETS = PC.ISOTOPE_MOD_START.value + PC.ISOTOPE_MOD_END.value
_UNKNOWN_BRACKETS = PC.UNKNOWN_MOD_START.value + PC.UNKNOWN_MOD_END.value
_TERM_MOD = PC.TERM_MOD.value
_INTERVAL_START = PC.INTERVAL_START
_INTERVAL_END = PC.INTERVAL_END
_UNKNOWN = PC.UNKNOWN
_CHARGE_SEP = PC.CHARGE_SEP


def _serialize_labile(
    mod_list: ModList,
    comps: list[str],
    include_plus: bool,
    precision: int | None = None,
) -> None:
    """Serialize labile modifications into comps list."""
    val = mod_list.serialize(
        brackets=_LABILE_BRACKETS,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )
    if val:
        comps.append(val)


def _serialize_static(
    mod_list: ModList,
    comps: list[str],
    include_plus: bool,
    precision: int | None = None,
) -> None:
    """Serialize static modifications into comps list."""
    val = mod_list.serialize(
        brackets=_STATIC_BRACKETS,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )
    if val:
        comps.append(val)


def _serialize_isotope(
    mod_list: ModList,
    comps: list[str],
    include_plus: bool,
    precision: int | None = None,
) -> None:
    """Serialize isotope modifications into comps list."""
    val = mod_list.serialize(
        brackets=_ISOTOPE_BRACKETS,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )
    if val:
        comps.append(val)


def _serialize_unknown(
    mod_list: ModList,
    comps: list[str],
    include_plus: bool,
    precision: int | None = None,
) -> None:
    """Serialize unknown modifications into comps list."""
    val = mod_list.serialize(
        brackets=_UNKNOWN_BRACKETS,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )

    if val:
        comps.append(val + _UNKNOWN)


def _serialize_nterm(
    mod_list: ModList,
    comps: list[str],
    include_plus: bool,
    precision: int | None = None,
) -> None:
    """Serialize N-terminal modifications into comps list."""
    val = mod_list.serialize(
        brackets=_MOD_BRACKETS,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )

    if val:
        comps.append(val + _TERM_MOD)


def _serialize_cterm(
    mod_list: ModList,
    comps: list[str],
    include_plus: bool,
    precision: int | None = None,
) -> None:
    """Serialize C-terminal modifications into comps list."""
    val = mod_list.serialize(
        brackets=_MOD_BRACKETS,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )

    if val:
        comps.append(_TERM_MOD + val)


def _serialize_start(
    annotation: ProFormaAnnotation,
    comps: list[str],
    include_plus: bool,
    precision: int | None = None,
) -> None:
    """Serialize the start portion into comps list."""
    # OPTIMIZED: Check has_* before calling getter to avoid lazy initialization
    if annotation.has_labile_mods:
        _serialize_labile(
            annotation.get_labile_mod_list(), comps, include_plus, precision
        )

    if annotation.has_static_mods:
        _serialize_static(
            annotation.get_static_mod_list(), comps, include_plus, precision
        )

    if annotation.has_isotope_mods:
        _serialize_isotope(
            annotation.get_isotope_mod_list(), comps, include_plus, precision
        )

    if annotation.has_unknown_mods:
        _serialize_unknown(
            annotation.get_unknown_mod_list(), comps, include_plus, precision
        )

    if annotation.has_nterm_mods:
        _serialize_nterm(
            annotation.get_nterm_mod_list(), comps, include_plus, precision
        )


def _serialize_interval_start(
    interval: Any,
    comps: list[str],
) -> None:
    """Serialize interval start markers."""
    comps.append(_INTERVAL_START)
    if interval.ambiguous:
        comps.append(_UNKNOWN)


def _serialize_interval_end(
    interval: Any,
    comps: list[str],
    include_plus: bool,
    precision: int | None,
) -> None:
    """Serialize interval end markers and modifications."""
    comps.append(_INTERVAL_END)
    if interval.mods:
        for mod in interval.mods:
            comps.append(mod.serialize(_MOD_BRACKETS, include_plus, precision))


def _serialize_intervals_at_position(
    interval_list: IntervalList | None,
    position: int,
    comps: list[str],
    include_plus: bool,
    precision: int | None,
) -> None:
    """Serialize all intervals that start or end at the given position."""
    if interval_list is None:
        return

    for interval in interval_list.data:
        if interval.start == position:
            _serialize_interval_start(interval, comps)
        if interval.end == position:
            _serialize_interval_end(interval, comps, include_plus, precision)


def _serialize_end_intervals(
    sequence_length: int,
    interval_list: IntervalList | None,
    comps: list[str],
    include_plus: bool,
    precision: int | None,
) -> None:
    """Serialize intervals that end at the sequence terminus."""
    if interval_list is not None:
        for interval in interval_list:
            if interval.end == sequence_length:
                _serialize_interval_end(interval, comps, include_plus, precision)


def _serialize_middle(
    annotation: ProFormaAnnotation,
    comps: list[str],
    include_plus: bool,
    precision: int | None = None,
) -> None:
    """Serialize the middle portion (sequence and internal mods) into comps list."""
    # FAST PATH: No internal mods or intervals - just add the sequence
    if not annotation.has_internal_mods and not annotation.has_intervals:
        comps.append(annotation.stripped_sequence)
        return

    # COMPLEX PATH: Has modifications or intervals
    seq = annotation.sequence
    seq_len = len(seq)

    # Get data structures once
    interval_list = annotation.get_interval_list() if annotation.has_intervals else None
    internal_mods_data = (
        annotation.get_internal_mod_dict().data
        if annotation.has_internal_mods
        else None
    )

    has_intervals = interval_list is not None
    has_internal_mods = internal_mods_data is not None

    # Main sequence loop with inlined internal mod serialization
    for i, aa in enumerate(seq):
        # Process intervals at this position
        if has_intervals:
            _serialize_intervals_at_position(
                interval_list, i, comps, include_plus, precision
            )

        # Add amino acid
        comps.append(aa)

        # INLINED: Internal mods at this position
        if has_internal_mods and i in internal_mods_data:
            for mod in internal_mods_data[i]:
                comps.append(mod.serialize(_MOD_BRACKETS, include_plus, precision))

    # Serialize end intervals
    _serialize_end_intervals(
        seq_len,
        interval_list,
        comps,
        include_plus,
        precision,
    )


def _serialize_charge(
    annotation: ProFormaAnnotation,
    comps: list[str],
    include_plus: bool = False,
    precision: int | None = None,
) -> None:
    """Serialize charge and charge adducts into comps list."""
    if annotation.has_charge:
        comps.append(f"{_CHARGE_SEP}{annotation.charge}")

    # OPTIMIZED: Check before calling getter
    if annotation.has_charge_adducts:
        for mod in annotation.get_charge_adduct_list():
            comps.append(mod.serialize(_MOD_BRACKETS, include_plus, precision))


def serialize_charge(
    annotation: ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """Serialize only the charge and charge adducts of a ProForma annotation."""
    comps: list[str] = []
    _serialize_charge(annotation, comps, include_plus, precision)
    return "".join(comps)


def _serialize_end(
    annotation: ProFormaAnnotation,
    comps: list[str],
    include_plus: bool,
    precision: int | None = None,
) -> None:
    """Serialize the end portion into comps list."""
    # C-term mods
    if annotation.has_cterm_mods:
        _serialize_cterm(
            annotation.get_cterm_mod_list(), comps, include_plus, precision
        )

    _serialize_charge(annotation, comps, include_plus, precision)


def serialize_annotation(
    annotation: ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """Serialize a complete ProForma annotation."""
    # Single list, single join at the end!
    comps: list[str] = []

    # Start
    _serialize_start(annotation, comps, include_plus, precision)

    # Middle (fast path for unmodified sequences)
    _serialize_middle(annotation, comps, include_plus, precision)

    # End
    _serialize_end(annotation, comps, include_plus, precision)

    # Single join operation!
    return "".join(comps)
