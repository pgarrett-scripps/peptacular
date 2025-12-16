from collections.abc import Sequence
from typing import Counter, overload

from ..constants import (
    ParrallelMethod,
    ParrallelMethodLiteral,
)
from ..fragment import IonType
from ..annotation import ProFormaAnnotation
from ..annotation.annotation import (
    ION_TYPE,
    CHARGE_TYPE,
    ISOTOPE_TYPE,
    LOSS_TYPE,
)
from .parrallel import parallel_apply_internal
from .util import get_annotation_input
from ..elements import ElementInfo


def _mass_single(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
) -> float:
    """Calculate mass for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.mass(
        ion_type=ion_type,
        charge=charge,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        losses=losses,
    )


# Type hints for single sequence
@overload
def mass(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> float: ...


# Type hints for list of sequences
@overload
def mass(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[float]: ...


# Implementation
def mass(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> float | list[float]:
    """
    Calculate the mass of an amino acid 'sequence'.
    """
    # Check if input is a list
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        # Parallel processing for lists
        return parallel_apply_internal(
            _mass_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            charge=charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            losses=losses,
        )
    else:
        # Single sequence processing
        return _mass_single(
            sequence=sequence,
            charge=charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            losses=losses,
        )


# Single sequence implementation for mz
def _mz_single(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
) -> float:
    """Calculate m/z for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.mz(
        ion_type=ion_type,
        charge=charge,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        losses=losses,
    )


# Type hints for mz
@overload
def mz(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> float: ...


@overload
def mz(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[float]: ...


def mz(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    monoisotopic: bool = True,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> float | list[float]:
    """
    Calculate the m/z (mass-to-charge ratio) of an amino acid 'sequence'.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _mz_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            charge=charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            losses=losses,
        )
    else:
        return _mz_single(
            sequence=sequence,
            charge=charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            losses=losses,
        )


def _comp_single(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
) -> Counter[ElementInfo]:
    """Calculate composition for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=True)

    return annotation.comp(
        ion_type=ion_type,
        charge=charge,
        isotopes=isotopes,
        losses=losses,
    )


@overload
def comp(
    sequence: str | ProFormaAnnotation,
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> Counter[ElementInfo]: ...


@overload
def comp(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[Counter[ElementInfo]]: ...


def comp(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_type: ION_TYPE = IonType.PRECURSOR,
    charge: CHARGE_TYPE | None = None,
    isotopes: ISOTOPE_TYPE | None = None,
    losses: dict[LOSS_TYPE, int] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> Counter[ElementInfo] | list[Counter[ElementInfo]]:
    """
    Calculates the elemental composition of a peptide sequence, including modifications.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _comp_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            ion_type=ion_type,
            charge=charge,
            isotopes=isotopes,
            losses=losses,
        )
    else:
        return _comp_single(
            sequence=sequence,
            ion_type=ion_type,
            charge=charge,
            isotopes=isotopes,
            losses=losses,
        )
