from __future__ import annotations

from typing import TYPE_CHECKING

from .characters import ProformaChar as PC

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


def _serialize_start(
    annotation: ProFormaAnnotation,
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize the start portion of a ProForma annotation."""
    comps: list[str] = []

    # Add labile mods
    for mod in annotation.get_labile_mod_list():
        comps.append(mod.serialize(PC.LABILE_MOD_START + PC.LABILE_MOD_END, include_plus, precision))

    for mod in annotation.get_static_mod_list():
        comps.append(mod.serialize(PC.STATIC_MOD_START + PC.STATIC_MOD_END, include_plus, precision))

    # Add global mods
    for mod in annotation.get_isotope_mod_list():
        comps.append(mod.serialize(PC.STATIC_MOD_START + PC.STATIC_MOD_END, include_plus, precision))

    # Unknown mods
    if annotation.has_unknown_mods:
        for mod in annotation.get_unknown_mod_list():
            comps.append(mod.serialize(PC.MOD_START + PC.MOD_END, include_plus, precision))
        comps.append(PC.UNKNOWN)

    # N-term mods
    if annotation.has_nterm_mods:
        for mod in annotation.get_nterm_mod_list():
            comps.append(mod.serialize(PC.MOD_START + PC.MOD_END, include_plus, precision))
        comps.append(PC.TERM_MOD)

    return "".join(comps)


def _serialize_middle(
    annotation: ProFormaAnnotation,
    include_plus: bool,
    precision: int | None = None,
) -> str:
    """Serialize the middle portion (sequence and internal mods) of a ProForma annotation."""
    comps: list[str] = []
    
    # Sequence
    for i, aa in enumerate(annotation.sequence):
        for interval in annotation.get_interval_list():
            if interval.start == i:
                comps.append(PC.INTERVAL_START)
                if interval.ambiguous:
                    comps.append(PC.UNKNOWN)
            if interval.end == i:
                comps.append(PC.INTERVAL_END)
                if interval.mods:
                    for mod in interval.mods:
                        comps.append(
                            mod.serialize(PC.MOD_START + PC.MOD_END, include_plus, precision)
                        )

        comps.append(aa)

        # Internal mods
        internal_mods = annotation.get_internal_mod_dict()
        if i in internal_mods:
            for mod in internal_mods[i]:
                comps.append(mod.serialize(PC.MOD_START + PC.MOD_END, include_plus, precision))

    # Add end interval
    i = len(annotation.sequence)
    for interval in annotation.get_interval_list():
        if interval.end == i:
            comps.append(PC.INTERVAL_END)
            if interval.mods:
                for mod in interval.mods:
                    comps.append(mod.serialize(PC.MOD_START + PC.MOD_END, include_plus, precision))

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
        comps.append(PC.TERM_MOD)
        for mod in annotation.get_cterm_mod_list():
            comps.append(mod.serialize(PC.MOD_START + PC.MOD_END, include_plus, precision))

    # Charge
    if annotation.has_charge:
        comps.append(f"{PC.CHARGE_SEP}{annotation.charge}")

    for mod in annotation.get_charge_adduct_list():
        comps.append(mod.serialize(PC.MOD_START + PC.MOD_END, include_plus, precision))

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