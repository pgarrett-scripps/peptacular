from __future__ import annotations

from typing import TYPE_CHECKING, Any

from .dclasses.modlist import ModList

from .characters import ProformaChar as PC

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


def _serialize_labile(
    mod_list: ModList,
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize labile modifications of a ProForma annotation."""
    return mod_list.serialize(
        brackets=PC.LABILE_MOD_START.value + PC.LABILE_MOD_END.value,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )


def _serialize_static(
    mod_list: ModList,
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize static modifications of a ProForma annotation."""
    return mod_list.serialize(
        brackets=PC.STATIC_MOD_START.value + PC.STATIC_MOD_END.value,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )


def _serialize_isotope(
    mod_list: ModList,
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize isotope modifications of a ProForma annotation."""
    return mod_list.serialize(
        brackets=PC.ISOTOPE_MOD_START.value + PC.ISOTOPE_MOD_END.value,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )


def _serialize_unknown(
    mod_list: ModList,
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize unknown modifications of a ProForma annotation."""
    val = mod_list.serialize(
        brackets=PC.UNKNOWN_MOD_START.value + PC.UNKNOWN_MOD_END.value,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )

    if len(val) > 0:
        val += PC.UNKNOWN
    return val


def _serialize_nterm(
    mod_list: ModList,
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize N-terminal modifications of a ProForma annotation."""
    val = mod_list.serialize(
        brackets=PC.MOD_START.value + PC.MOD_END.value,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )

    if len(val) > 0:
        val += PC.TERM_MOD.value
    return val


def _serialize_cterm(
    mod_list: ModList,
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize C-terminal modifications of a ProForma annotation."""
    val = mod_list.serialize(
        brackets=PC.MOD_START.value + PC.MOD_END.value,
        sort=False,
        include_plus=include_plus,
        precision=precision,
    )

    if len(val) > 0:
        val = PC.TERM_MOD.value + val

    return val


def _serialize_start(
    annotation: ProFormaAnnotation,
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize the start portion of a ProForma annotation."""
    comps: list[str] = []

    # OPTIMIZED: Check has_* before calling getter to avoid lazy initialization
    if annotation.has_labile_mods:
        comps.append(
            _serialize_labile(annotation.get_labile_mod_list(), include_plus, precision)
        )

    if annotation.has_static_mods:
        comps.append(
            _serialize_static(annotation.get_static_mod_list(), include_plus, precision)
        )

    if annotation.has_isotope_mods:
        comps.append(
            _serialize_isotope(
                annotation.get_isotope_mod_list(), include_plus, precision
            )
        )

    if annotation.has_unknown_mods:
        comps.append(
            _serialize_unknown(
                annotation.get_unknown_mod_list(), include_plus, precision
            )
        )

    if annotation.has_nterm_mods:
        comps.append(
            _serialize_nterm(annotation.get_nterm_mod_list(), include_plus, precision)
        )

    return "".join(comps)


MOD_BRACKETS = PC.MOD_START + PC.MOD_END


def _serialize_interval_start(
    interval: Any,  # Replace with actual type
    comps: list[str],
) -> None:
    """Serialize interval start markers."""
    comps.append(PC.INTERVAL_START)
    if interval.ambiguous:
        comps.append(PC.UNKNOWN)


def _serialize_interval_end(
    interval: Any,  # Replace with actual type
    comps: list[str],
    include_plus: bool,
    precision: int | None,
) -> None:
    """Serialize interval end markers and modifications."""
    comps.append(PC.INTERVAL_END)
    if interval.mods:
        for mod in interval.mods:
            comps.append(mod.serialize(MOD_BRACKETS, include_plus, precision))


def _serialize_intervals_at_position(
    interval_list: list,  # Replace with actual type
    position: int,
    comps: list[str],
    include_plus: bool,
    precision: int | None,
) -> None:
    """Serialize all intervals that start or end at the given position."""
    for interval in interval_list:
        if interval.start == position:
            _serialize_interval_start(interval, comps)
        if interval.end == position:
            _serialize_interval_end(interval, comps, include_plus, precision)


def _serialize_internal_mods_at_position(
    internal_mods_data: dict[int, Any],  # Replace Any with ModList type
    position: int,
    comps: list[str],
    include_plus: bool,
    precision: int | None,
) -> None:
    """Serialize internal modifications at the given position."""
    if position in internal_mods_data:
        for mod in internal_mods_data[position]:
            comps.append(mod.serialize(MOD_BRACKETS, include_plus, precision))


def _serialize_sequence_with_mods(
    sequence: str,
    interval_list: list | None,  # Replace with actual type
    internal_mods_data: dict[int, Any] | None,  # Replace Any with ModList type
    comps: list[str],
    include_plus: bool,
    precision: int | None,
) -> None:
    """Serialize the sequence with intervals and internal modifications."""
    has_intervals = interval_list is not None
    has_internal_mods = internal_mods_data is not None

    for i, aa in enumerate(sequence):
        # Process intervals at this position
        if has_intervals:
            _serialize_intervals_at_position(
                interval_list, i, comps, include_plus, precision
            )

        # Add amino acid
        comps.append(aa)

        # Add internal mods at this position
        if has_internal_mods:
            _serialize_internal_mods_at_position(
                internal_mods_data, i, comps, include_plus, precision
            )


def _serialize_end_intervals(
    sequence_length: int,
    interval_list: list | None,  # Replace with actual type
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
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize the middle portion (sequence and internal mods) of a ProForma annotation."""
    comps: list[str] = []
    comps_append = comps.append  # Local reference for speed

    # Get data structures once
    interval_list = annotation.get_interval_list() if annotation.has_intervals else None
    internal_mods_data = (
        annotation.get_internal_mod_dict().data
        if annotation.has_internal_mods
        else None
    )

    # OPTIMIZED: Inline the hot path
    has_intervals = interval_list is not None
    has_internal_mods = internal_mods_data is not None

    # Main sequence loop with inlined internal mod serialization
    for i, aa in enumerate(annotation.sequence):
        # Process intervals at this position
        if has_intervals:
            _serialize_intervals_at_position(
                interval_list, i, comps, include_plus, precision
            )

        # Add amino acid
        comps_append(aa)

        # INLINED: Internal mods at this position (hot path - called 137,760 times)
        if has_internal_mods and i in internal_mods_data:
            for mod in internal_mods_data[i]:
                comps_append(mod.serialize(MOD_BRACKETS, include_plus, precision))

    # Serialize end intervals
    _serialize_end_intervals(
        len(annotation.sequence),
        interval_list,
        comps,
        include_plus,
        precision,
    )

    return "".join(comps)


def serialize_charge(
    annotation: ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    comps: list[str] = []

    if annotation.has_charge:
        comps.append(f"{PC.CHARGE_SEP}{annotation.charge}")

    # OPTIMIZED: Check before calling getter
    if annotation.has_charge_adducts:
        for mod in annotation.get_charge_adduct_list():
            comps.append(
                mod.serialize(PC.MOD_START + PC.MOD_END, include_plus, precision)
            )

    return "".join(comps)


def _serialize_end(
    annotation: ProFormaAnnotation,
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize the end portion of a ProForma annotation."""
    comps: list[str] = []

    # C-term mods
    if annotation.has_cterm_mods:
        comps.append(
            _serialize_cterm(annotation.get_cterm_mod_list(), include_plus, precision)
        )

    comps.append(serialize_charge(annotation, include_plus, precision))
    return "".join(comps)


def serialize_annotation(
    annotation: ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """Serialize a complete ProForma annotation."""
    return (
        _serialize_start(annotation, include_plus, precision)
        + _serialize_middle(annotation, include_plus, precision)
        + _serialize_end(annotation, include_plus, precision)
    )
