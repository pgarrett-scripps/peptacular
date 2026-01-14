from collections.abc import Sequence
from typing import overload

from ..annotation import ProFormaAnnotation
from ..annotation.annotation import (
    CHARGE_TYPE,
    ION_TYPE,
    ISOTOPE_TYPE,
    LOSS_TYPE,
)
from ..annotation.utils import Fragment
from ..constants import ParrallelMethod, ParrallelMethodLiteral
from ..fragment import IonType
from .parrallel import parallel_apply_internal
from .util import get_annotation_input


def _fragment_single(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] | None = None,
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE] | None = None,
    losses: Sequence[LOSS_TYPE] | None = None,
    calculate_composition: bool = False,
    include_sequence: bool = False,
    max_losses: int = 1,
) -> list[Fragment]:
    annotation = get_annotation_input(sequence=sequence, copy=False)

    return annotation.fragment(
        ion_types=ion_types,
        charges=charges,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        losses=losses,
        calculate_composition=calculate_composition,
        include_sequence=include_sequence,
        max_losses=max_losses,
    )


@overload
def fragment(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[ION_TYPE] | ION_TYPE = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] | None = None,
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE] | None = None,
    losses: Sequence[LOSS_TYPE] | None = None,
    calculate_composition: bool = False,
    include_sequence: bool = False,
    max_losses: int = 1,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[Fragment]: ...


@overload
def fragment(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[ION_TYPE] | ION_TYPE = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] | None = None,
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE] | None = None,
    losses: Sequence[LOSS_TYPE] | None = None,
    calculate_composition: bool = False,
    include_sequence: bool = False,
    max_losses: int = 1,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[Fragment]]: ...


def fragment(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] | None = None,
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE] | None = None,
    losses: Sequence[LOSS_TYPE] | None = None,
    calculate_composition: bool = False,
    include_sequence: bool = False,
    max_losses: int = 1,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[Fragment] | list[list[Fragment]]:
    """
    Builds fragment ions from a given input sequence or list of sequences.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _fragment_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            ion_types=ion_types,
            charges=charges,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            losses=losses,
            max_losses=max_losses,
            calculate_composition=calculate_composition,
            include_sequence=include_sequence,
        )
    else:
        return _fragment_single(
            sequence=sequence,
            ion_types=ion_types,
            charges=charges,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            losses=losses,
            max_losses=max_losses,
            calculate_composition=calculate_composition,
            include_sequence=include_sequence,
        )
