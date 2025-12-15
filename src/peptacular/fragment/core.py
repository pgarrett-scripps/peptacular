from __future__ import annotations

import itertools
import re
from collections import defaultdict
from typing import Counter, Generator, Mapping, Sequence

from ..constants import (
    BACKWARD_ION_TYPES,
    C13_NEUTRON_MASS,
    FORWARD_ION_TYPES,
    INTERNAL_ION_TYPES,
    TERMINAL_ION_TYPES,
    IonType,
    IonTypeLiteral,
)
from ..funcs import get_label, get_number
from ..mass_calc import adjust_mass, adjust_mz
from ..spans import (
    build_left_semi_spans,
    build_non_enzymatic_spans,
    build_right_semi_spans,
)
from .types import (
    FRAGMENT_RETURN_TYPING,
    Fragment,
    FragmentableAnnotation,
    FragmentReturnLiteral,
    FragmentReturnType,
)


def get_deltas(
    annotation: FragmentableAnnotation,
    deltas: Mapping[str, Sequence[int | float | str | dict[str, float | int]]],
    max_losses: int,
    use_regex: bool = False,
) -> set[float]:
    """
    Calculate applicable deltas for a given annotation.
    Returns a list to allow multiple deltas of the same type.
    """
    applicable_deltas: list[tuple[int | float | str | dict[str, float | int], int]] = []

    for restr, loss_list in deltas.items():
        # Count how many times each pattern matches
        if use_regex:
            matches = re.findall(restr, annotation.stripped_sequence)
        else:
            matches = [i for i, aa in enumerate(annotation.stripped_sequence) if aa in restr]
            
        match_count = len(matches)

        # Add each loss value multiplied by the number of matches
        for loss_value in loss_list:
            applicable_deltas.append((loss_value, match_count))

    # Convert Counter to list of individual deltas
    individual_deltas: list[float] = []
    for loss_value, count in applicable_deltas.items():
        individual_deltas.extend([loss_value] * count)

    # Generate combinations if max_losses > 1
    all_deltas: set[float] = set()
    if max_losses > 1 and individual_deltas:
        for loss_count in range(1, min(max_losses + 1, len(individual_deltas) + 1)):
            for loss_combination in itertools.combinations(
                individual_deltas, loss_count
            ):
                all_deltas.add(sum(loss_combination))
    else:
        # If max_losses == 1, just add individual deltas
        all_deltas.update(individual_deltas)

    # Always include 0.0 (no loss)
    all_deltas.add(0.0)

    return all_deltas


def _build_fragments(
    annotation: FragmentableAnnotation,
    spans: Sequence[tuple[int, int, int]],
    ion_types: Sequence[IonType],
    charges: Sequence[int],
    deltas: Mapping[str, Sequence[int | float | str | dict[str, float | int]]],
    isotopes: Sequence[int | float | str | dict[str, float | int]],
    monoisotopic: bool,
    return_type: FragmentReturnType,
    mass_components: Sequence[float],
    max_losses: int,
    precision: None | int = None,
    use_isotope_on_mods: bool = False,
    default_isotopic_mass: float = C13_NEUTRON_MASS,
) -> Generator[FRAGMENT_RETURN_TYPING, None, None]:
    """
    Builds fragments for a given sequence.
    """

    if len(annotation.stripped_sequence) == 0:
        return

    annot_sequence = annotation.serialize()

    for span in spans:
        base_mass = adjust_mass(
            sum(mass_components[span[0] : span[1]]),
            charge=0,
            ion_type=IonType.NEUTRAL,
            monoisotopic=monoisotopic,
        )
        base_annotation = annotation.slice(span[0], span[1], inplace=False)
        base_sequence = base_annotation.serialize()
        base_unmod_sequence = base_annotation.stripped_sequence

        applicable_deltas = get_deltas(
            annotation=base_annotation, deltas=deltas, max_losses=max_losses
        )

        for ion_type in ion_types:
            for iso in isotopes:
                for loss in applicable_deltas:
                    for c in charges:
                        fragment_mass = adjust_mass(
                            base_mass=base_mass,
                            charge=c,
                            ion_type=ion_type,
                            monoisotopic=monoisotopic,
                            precision=precision,
                            isotope=iso,
                            loss=loss,
                        )
                        fragment_mz = adjust_mz(
                            base_mass=fragment_mass, charge=c, precision=precision
                        )

                        if return_type == FragmentReturnType.FRAGMENT:
                            yield Fragment(
                                charge=c,
                                ion_type=ion_type,
                                start=span[0],
                                end=span[1],
                                monoisotopic=monoisotopic,
                                isotope=iso,
                                loss=loss,
                                parent_sequence=annot_sequence,
                                mz=fragment_mz,
                                sequence=base_sequence,
                                parent_length=len(annotation),
                                internal=span[0] != 0 and span[1] != len(annotation),
                            )

                        elif return_type == FragmentReturnType.LABEL:
                            number = get_number(
                                ion_type, len(base_unmod_sequence), span[0], span[1]
                            )
                            yield get_label(
                                ion_type, c, number, loss, iso, precision=precision
                            )

                        elif return_type == FragmentReturnType.MASS:
                            yield fragment_mass

                        elif return_type == FragmentReturnType.MZ:
                            yield fragment_mz

                        elif return_type == FragmentReturnType.MASS_LABEL:
                            number = get_number(
                                ion_type, len(base_unmod_sequence), span[0], span[1]
                            )
                            yield (
                                fragment_mass,
                                get_label(
                                    ion_type, c, number, loss, iso, precision=precision
                                ),
                            )

                        elif return_type == FragmentReturnType.MZ_LABEL:
                            number = get_number(
                                ion_type, len(base_unmod_sequence), span[0], span[1]
                            )
                            yield (
                                fragment_mz,
                                get_label(
                                    ion_type, c, number, loss, iso, precision=precision
                                ),
                            )


def _get_internal_fragments(
    annotation: FragmentableAnnotation,
    ion_types: Sequence[IonType],
    charges: Sequence[int],
    monoisotopic: bool,
    isotopes: Sequence[int],
    deltas: Mapping[str, Sequence[float]],
    return_type: FragmentReturnType,
    mass_components: Sequence[float],
    precision: None | int,
    max_losses: int,
) -> Generator[FRAGMENT_RETURN_TYPING, None, None]:
    """
    Build internal fragments for a given sequence.
    """
    spans = list(build_non_enzymatic_spans((0, len(annotation), 0)))
    internal_spans = [
        span for span in spans if span[0] != 0 and span[1] != len(annotation)
    ]
    yield from _build_fragments(
        annotation=annotation,
        spans=internal_spans,
        ion_types=ion_types,
        charges=charges,
        deltas=deltas,
        isotopes=isotopes,
        monoisotopic=monoisotopic,
        return_type=return_type,
        mass_components=mass_components,
        max_losses=max_losses,
        precision=precision,
    )


def _get_immonium_fragments(
    annotation: FragmentableAnnotation,
    charges: Sequence[int],
    monoisotopic: bool,
    isotopes: Sequence[int],
    deltas: Mapping[str, Sequence[float]],
    return_type: FragmentReturnType,
    mass_components: Sequence[float],
    precision: None | int,
    max_losses: int,
) -> Generator[FRAGMENT_RETURN_TYPING, None, None]:
    """
    Build immonium ions for a given sequence.
    """
    spans = [(i, i + 1, 0) for i in range(len(annotation))]
    yield from _build_fragments(
        annotation=annotation,
        spans=spans,
        ion_types=[IonType.IMMONIUM],
        charges=charges,
        deltas=deltas,
        isotopes=isotopes,
        monoisotopic=monoisotopic,
        return_type=return_type,
        mass_components=mass_components,
        max_losses=max_losses,
        precision=precision,
    )


def _get_forward_fragments(
    annotation: FragmentableAnnotation,
    ion_types: Sequence[IonType],
    charges: Sequence[int],
    monoisotopic: bool,
    isotopes: Sequence[int],
    deltas: Mapping[str, Sequence[float]],
    return_type: FragmentReturnType,
    mass_components: Sequence[float],
    precision: None | int,
    max_losses: int,
) -> Generator[FRAGMENT_RETURN_TYPING, None, None]:
    """
    Build forward fragments for a given sequence.
    """
    start_span = (0, len(annotation), 0)
    spans = [start_span] + list(build_left_semi_spans(start_span))
    yield from _build_fragments(
        annotation=annotation,
        spans=spans,
        ion_types=ion_types,
        charges=charges,
        deltas=deltas,
        isotopes=isotopes,
        monoisotopic=monoisotopic,
        return_type=return_type,
        mass_components=mass_components,
        max_losses=max_losses,
        precision=precision,
    )


def _get_backward_fragments(
    annotation: FragmentableAnnotation,
    ion_types: Sequence[IonType],
    charges: Sequence[int],
    monoisotopic: bool,
    isotopes: Sequence[int],
    deltas: Mapping[str, Sequence[float]],
    return_type: FragmentReturnType,
    mass_components: Sequence[float],
    precision: None | int,
    max_losses: int,
) -> Generator[FRAGMENT_RETURN_TYPING, None, None]:
    """
    Build backward fragments for a given sequence.
    """
    start_span = (0, len(annotation), 0)
    spans = [start_span] + list(build_right_semi_spans(start_span))
    yield from _build_fragments(
        annotation=annotation,
        spans=spans,
        ion_types=ion_types,
        charges=charges,
        deltas=deltas,
        isotopes=isotopes,
        monoisotopic=monoisotopic,
        return_type=return_type,
        mass_components=mass_components,
        max_losses=max_losses,
        precision=precision,
    )


def _get_terminal_fragments(
    annotation: FragmentableAnnotation,
    ion_types: Sequence[IonType],
    charges: Sequence[int],
    monoisotopic: bool,
    isotopes: Sequence[int],
    deltas: Mapping[str, Sequence[float]],
    return_type: FragmentReturnType,
    mass_components: Sequence[float],
    precision: None | int,
    max_losses: int,
) -> Generator[FRAGMENT_RETURN_TYPING, None, None]:
    """
    Build terminal fragments for a given sequence.
    """
    forward_ions = [ion for ion in ion_types if ion in FORWARD_ION_TYPES]
    backward_ions = [ion for ion in ion_types if ion in BACKWARD_ION_TYPES]

    yield from _get_forward_fragments(
        annotation=annotation,
        ion_types=forward_ions,
        charges=charges,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        deltas=deltas,
        return_type=return_type,
        mass_components=mass_components,
        precision=precision,
        max_losses=max_losses,
    )

    yield from _get_backward_fragments(
        annotation=annotation,
        ion_types=backward_ions,
        charges=charges,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        deltas=deltas,
        return_type=return_type,
        mass_components=mass_components,
        precision=precision,
        max_losses=max_losses,
    )


def fragment(
    annotation: FragmentableAnnotation,
    ion_types: Sequence[IonTypeLiteral | IonType] | IonTypeLiteral | IonType,
    charges: Sequence[int] | int,
    monoisotopic: bool = True,
    *,
    isotopes: Sequence[int] | int = 0,
    water_loss: bool = False,
    ammonia_loss: bool = False,
    deltas: None | Mapping[str, Sequence[float]] = None,
    max_losses: int = 1,
    precision: None | int = None,
    _mass_components: None | Sequence[float] = None,
    return_type: (
        FragmentReturnLiteral | FragmentReturnType
    ) = FragmentReturnType.FRAGMENT,
) -> Generator[FRAGMENT_RETURN_TYPING, None, None]:
    ion_type_list: list[IonType] = []

    if isinstance(ion_types, IonType):
        ion_type_list = [ion_types]
    elif isinstance(ion_types, str):
        ion_type_list = [IonType(ion_types)]
    elif isinstance(ion_types, Sequence):  # type: ignore
        ion_type_list = [IonType(it) for it in ion_types]
    else:
        raise TypeError("Invalid ion_types format")

    charge_list: list[int] = []
    if not isinstance(charges, Sequence):
        charge_list = [charges]
    else:
        charge_list = list(charges)

    isotope_list: list[int] = []
    if not isinstance(isotopes, Sequence):
        isotope_list = [isotopes]
    else:
        isotope_list = list(isotopes)

    loss_dict: dict[str, list[float]] = defaultdict(list)

    if deltas is not None:
        if isinstance(deltas, Mapping):  # type: ignore
            for k, v in deltas.items():
                loss_dict[k].extend(v)
        else:
            raise TypeError("Invalid deltas format")

    if water_loss:
        loss_dict["[STED]"].append(-18.01056)

    if ammonia_loss:
        loss_dict["[RKNQ]"].append(-17.02655)

    # Convert back to regular dict if needed
    loss_dict = dict(loss_dict)

    if annotation.contains_sequence_ambiguity():
        raise ValueError("Ambiguous sequence")

    terminal_fragment_types = [i for i in ion_type_list if i in TERMINAL_ION_TYPES]
    internal_fragment_types = [i for i in ion_type_list if i in INTERNAL_ION_TYPES]
    immonium = IonType.IMMONIUM in ion_type_list

    if _mass_components is None:
        components = annotation.split()
        _mass_components = [
            c.set_charge(0).mass(ion_type=IonType.NEUTRAL, monoisotopic=monoisotopic)
            for c in components
        ]

    return_type_enum = FragmentReturnType(return_type)

    if terminal_fragment_types:
        yield from _get_terminal_fragments(
            annotation=annotation,
            ion_types=terminal_fragment_types,
            charges=charge_list,
            monoisotopic=monoisotopic,
            isotopes=isotope_list,
            deltas=loss_dict,
            return_type=return_type_enum,
            mass_components=_mass_components,
            precision=precision,
            max_losses=max_losses,
        )

    if internal_fragment_types:
        yield from _get_internal_fragments(
            annotation=annotation,
            ion_types=internal_fragment_types,
            charges=charge_list,
            monoisotopic=monoisotopic,
            isotopes=isotope_list,
            deltas=loss_dict,
            return_type=return_type_enum,
            mass_components=_mass_components,
            precision=precision,
            max_losses=max_losses,
        )

    if immonium:
        yield from _get_immonium_fragments(
            annotation=annotation,
            charges=charge_list,
            monoisotopic=monoisotopic,
            isotopes=isotope_list,
            deltas=loss_dict,
            return_type=return_type_enum,
            mass_components=_mass_components,
            precision=precision,
            max_losses=max_losses,
        )
