from collections.abc import Sequence
from typing import overload

from tacular import IonType

from ..annotation import ProFormaAnnotation
from ..annotation.annotation import (
    CHARGE_TYPE,
    CUSTOM_LOSS_TYPE,
    ION_TYPE,
    ISOTOPE_TYPE,
    LOSS_TYPE,
)
from ..annotation.utils import Fragment
from ..constants import parallelMethod, parallelMethodLiteral
from .parallel import parallel_apply_internal
from .util import get_annotation_input


def _fragment_single(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] = (1,),
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE] = (),
    calculate_composition: bool = False,
    include_sequence: bool = False,
    max_ndeltas: int = 1,
) -> list[Fragment]:
    annotation = get_annotation_input(sequence=sequence, copy=False)

    return annotation.fragment(
        ion_types=ion_types,
        charges=charges,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        deltas=deltas,
        neutral_deltas=neutral_deltas,
        calculate_composition=calculate_composition,
        include_sequence=include_sequence,
        max_ndeltas=max_ndeltas,
    )


@overload
def fragment(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[ION_TYPE] | ION_TYPE = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] = (1,),
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE | None] = (None,),
    calculate_composition: bool = False,
    include_sequence: bool = False,
    max_ndeltas: int = 1,
    n_workers: None = None,
    chunksize: None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[Fragment]: ...


@overload
def fragment(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[ION_TYPE] | ION_TYPE = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] = (1,),
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE | None] = (None,),
    calculate_composition: bool = False,
    include_sequence: bool = False,
    max_ndeltas: int = 1,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[list[Fragment]]: ...


def fragment(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[ION_TYPE] = (IonType.B, IonType.Y),
    charges: Sequence[CHARGE_TYPE] = (1,),
    monoisotopic: bool = True,
    isotopes: Sequence[ISOTOPE_TYPE | None] = (0,),
    deltas: Sequence[CUSTOM_LOSS_TYPE | None] = (None,),
    neutral_deltas: Sequence[LOSS_TYPE] = (),
    calculate_composition: bool = False,
    include_sequence: bool = False,
    max_ndeltas: int = 1,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: parallelMethod | parallelMethodLiteral | None = None,
) -> list[Fragment] | list[list[Fragment]]:
    """
    Builds fragment ions from a given input sequence or list of sequences.
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str) and not isinstance(sequence, ProFormaAnnotation):
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
            deltas=deltas,
            neutral_deltas=neutral_deltas,
            max_ndeltas=max_ndeltas,
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
            deltas=deltas,
            neutral_deltas=neutral_deltas,
            max_ndeltas=max_ndeltas,
            calculate_composition=calculate_composition,
            include_sequence=include_sequence,
        )
