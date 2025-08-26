from __future__ import annotations

from typing import TYPE_CHECKING

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

    comps.append(
        _serialize_labile(annotation.get_labile_mod_list(), include_plus, precision)
    )
    comps.append(
        _serialize_static(annotation.get_static_mod_list(), include_plus, precision)
    )
    comps.append(
        _serialize_isotope(annotation.get_isotope_mod_list(), include_plus, precision)
    )
    comps.append(
        _serialize_unknown(annotation.get_unknown_mod_list(), include_plus, precision)
    )
    comps.append(
        _serialize_nterm(annotation.get_nterm_mod_list(), include_plus, precision)
    )

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
                            mod.serialize(
                                PC.MOD_START + PC.MOD_END, include_plus, precision
                            )
                        )

        comps.append(aa)

        # Internal mods
        internal_mods = annotation.get_internal_mod_dict()
        if i in internal_mods:
            for mod in internal_mods[i]:
                comps.append(
                    mod.serialize(PC.MOD_START + PC.MOD_END, include_plus, precision)
                )

    # Add end interval
    i = len(annotation.sequence)
    for interval in annotation.get_interval_list():
        if interval.end == i:
            comps.append(PC.INTERVAL_END)
            if interval.mods:
                for mod in interval.mods:
                    comps.append(
                        mod.serialize(
                            PC.MOD_START + PC.MOD_END, include_plus, precision
                        )
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
    comps.append(
        _serialize_cterm(annotation.get_cterm_mod_list(), include_plus, precision)
    )

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
