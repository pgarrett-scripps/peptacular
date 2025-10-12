"""
fragmentation.py contains functions for fragmenting peptides
"""

from collections.abc import Sequence
from typing import Literal, Mapping, overload

from ..constants import IonTypeLiteral
from ..fragment import FragmentReturnType
from ..proforma.annotation import ProFormaAnnotation
from .basic import get_annotation_input
from .parrallel import parallel_apply_internal


def _get_losses_single(
    sequence: str | ProFormaAnnotation,
    losses: Mapping[str, Sequence[float]],
    max_losses: int,
) -> set[float]:
    """Get losses for a single sequence"""
    annotation = get_annotation_input(sequence, copy=False)
    return annotation.get_losses(losses=losses, max_losses=max_losses)


@overload
def get_losses(
    sequence: str | ProFormaAnnotation,
    losses: Mapping[str, Sequence[float]],
    max_losses: int,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> set[float]: ...


@overload
def get_losses(
    sequence: Sequence[str | ProFormaAnnotation],
    losses: Mapping[str, Sequence[float]],
    max_losses: int,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[set[float]]: ...


def get_losses(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    losses: Mapping[str, Sequence[float]],
    max_losses: int,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> set[float] | list[set[float]]:
    """
    Returns a set of applicable losses for a given sequence or list of sequences.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence to check for losses, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param losses: A mapping of amino acids to their neutral losses.
    :type losses: Mapping[str, Sequence[float]]
    :param max_losses: The maximum number of losses to consider.
    :type max_losses: int
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: A set of applicable losses, or list of sets for multiple sequences.
    :rtype: set[float] | list[set[float]]

    .. code-block:: python

        >>> get_losses('AA', {'A': [-10, -5]}, 1)
        {0.0, -5, -10}

        >>> get_losses('AA', {'A': [-10, -5]}, 2)
        {0.0, -20, -15, -10, -5}

        >>> get_losses('AA', {'A': [-10, -5]}, 3)
        {0.0, -25, -20, -15, -10, -5}

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _get_losses_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            losses=losses,
            max_losses=max_losses,
        )
    else:
        return _get_losses_single(
            sequence=sequence,
            losses=losses,
            max_losses=max_losses,
        )


def _fragment_single(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[IonTypeLiteral] | IonTypeLiteral,
    charges: Sequence[int] | int,
    monoisotopic: bool = True,
    isotopes: Sequence[int] | int = 0,
    water_loss: bool = False,
    ammonia_loss: bool = False,
    losses: Mapping[str, Sequence[float]] | None = None,
    max_losses: int = 1,
    precision: int | None = None,
    return_type: Literal["mz", "mz-label"]
    | Literal[
        FragmentReturnType.MZ, FragmentReturnType.MZ_LABEL
    ] = FragmentReturnType.MZ,
    _mass_components: Sequence[float] | None = None,
) -> list[float] | list[tuple[float, str]]:
    """Fragment a single sequence"""
    return_type_enum: FragmentReturnType = FragmentReturnType(return_type)

    if (
        return_type_enum != FragmentReturnType.MZ
        and return_type_enum != FragmentReturnType.MZ_LABEL
    ):
        raise TypeError(f"return type not supported: {return_type_enum}")

    return list(
        get_annotation_input(sequence=sequence, copy=False).fragment(
            ion_types=ion_types,
            charges=charges,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            water_loss=water_loss,
            ammonia_loss=ammonia_loss,
            losses=losses,
            max_losses=max_losses,
            precision=precision,
            _mass_components=_mass_components,
            return_type=return_type,
        )
    )


@overload
def fragment(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[IonTypeLiteral] | IonTypeLiteral,
    charges: Sequence[int] | int,
    monoisotopic: bool = True,
    *,
    isotopes: Sequence[int] | int = 0,
    water_loss: bool = False,
    ammonia_loss: bool = False,
    losses: Mapping[str, Sequence[float]] | None = None,
    max_losses: int = 1,
    precision: int | None = None,
    return_type: Literal["mz"] | Literal[FragmentReturnType.MZ] = FragmentReturnType.MZ,
    _mass_components: Sequence[float] | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


@overload
def fragment(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[IonTypeLiteral] | IonTypeLiteral,
    charges: Sequence[int] | int,
    monoisotopic: bool = True,
    *,
    isotopes: Sequence[int] | int = 0,
    water_loss: bool = False,
    ammonia_loss: bool = False,
    losses: Mapping[str, Sequence[float]] | None = None,
    max_losses: int = 1,
    precision: int | None = None,
    return_type: Literal["mz"] | Literal[FragmentReturnType.MZ] = FragmentReturnType.MZ,
    _mass_components: Sequence[float] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[list[float]]: ...


@overload
def fragment(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[IonTypeLiteral] | IonTypeLiteral,
    charges: Sequence[int] | int,
    monoisotopic: bool = True,
    *,
    isotopes: Sequence[int] | int = 0,
    water_loss: bool = False,
    ammonia_loss: bool = False,
    losses: Mapping[str, Sequence[float]] | None = None,
    max_losses: int = 1,
    precision: int | None = None,
    return_type: Literal["mz-label"]
    | Literal[FragmentReturnType.MZ_LABEL] = FragmentReturnType.MZ_LABEL,
    _mass_components: Sequence[float] | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[tuple[float, str]]: ...


@overload
def fragment(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[IonTypeLiteral] | IonTypeLiteral,
    charges: Sequence[int] | int,
    monoisotopic: bool = True,
    *,
    isotopes: Sequence[int] | int = 0,
    water_loss: bool = False,
    ammonia_loss: bool = False,
    losses: Mapping[str, Sequence[float]] | None = None,
    max_losses: int = 1,
    precision: int | None = None,
    return_type: Literal["mz-label"]
    | Literal[FragmentReturnType.MZ_LABEL] = FragmentReturnType.MZ_LABEL,
    _mass_components: Sequence[float] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[list[tuple[float, str]]]: ...


def fragment(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_types: Sequence[IonTypeLiteral] | IonTypeLiteral,
    charges: Sequence[int] | int,
    monoisotopic: bool = True,
    *,
    isotopes: Sequence[int] | int = 0,
    water_loss: bool = False,
    ammonia_loss: bool = False,
    losses: Mapping[str, Sequence[float]] | None = None,
    max_losses: int = 1,
    precision: int | None = None,
    return_type: (
        Literal["mz", "mz-label"]
        | Literal[FragmentReturnType.MZ, FragmentReturnType.MZ_LABEL]
    ) = FragmentReturnType.MZ,
    _mass_components: Sequence[float] | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> (
    list[float]
    | list[tuple[float, str]]
    | list[list[float]]
    | list[list[tuple[float, str]]]
):
    """
    Builds fragment ions from a given input sequence or list of sequences.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The amino acid sequence, ProForma annotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param ion_types: Ion types to consider, e.g., ['b', 'y']
    :type ion_types: Sequence[IonTypeLiteral] | IonTypeLiteral
    :param charges: Charge states for the fragments.
    :type charges: Sequence[int] | int
    :param monoisotopic: If True, use monoisotopic masses. If False, use average masses.
    :type monoisotopic: bool
    :param isotopes: Isotope offsets to consider.
    :type isotopes: Sequence[int] | int
    :param water_loss: If True, consider water loss.
    :type water_loss: bool
    :param ammonia_loss: If True, consider ammonia loss.
    :type ammonia_loss: bool
    :param losses: Mapping of neutral losses to consider.
    :type losses: Mapping[str, Sequence[float]] | None
    :param max_losses: The maximum number of losses to consider.
    :type max_losses: int
    :param return_type: The type of data to return ('mz' or 'mz-label').
    :type return_type: Literal["mz", "mz-label"] | Literal[FragmentReturnType.MZ, FragmentReturnType.MZ_LABEL]
    :param precision: Number of decimal places to round values.
    :type precision: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: List of m/z values or (m/z, label) tuples, or list of lists for multiple sequences.
    :rtype: list[float] | list[tuple[float, str]] | list[list[float]] | list[list[tuple[float, str]]]

    .. code-block:: python

        # Single sequence
        >>> fragment('TIDE', 'b', 1, return_type='mz', precision=3)
        [459.209, 330.166, 215.139, 102.055]

        >>> fragment('TIDE', 'b', 1, return_type='mz-label', precision=3)
        [(459.209, '+b4'), (330.166, '+b3'), (215.139, '+b2'), (102.055, '+b1')]

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
            water_loss=water_loss,
            ammonia_loss=ammonia_loss,
            losses=losses,
            max_losses=max_losses,
            precision=precision,
            return_type=return_type,
            _mass_components=_mass_components,
        )
    else:
        return _fragment_single(
            sequence=sequence,
            ion_types=ion_types,
            charges=charges,
            monoisotopic=monoisotopic,
            isotopes=isotopes,
            water_loss=water_loss,
            ammonia_loss=ammonia_loss,
            losses=losses,
            max_losses=max_losses,
            precision=precision,
            return_type=return_type,
            _mass_components=_mass_components,
        )
