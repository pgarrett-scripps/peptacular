"""
fragmentation.py contains functions for fragmenting peptides
"""

from collections.abc import Sequence
from typing import Any, Generator, Literal, Mapping, cast, overload

from ..constants import IonTypeLiteral
from ..fragment import FragmentReturnType
from ..proforma.annotation import ProFormaAnnotation
from .sequence import get_annotation_input


def get_losses(
    sequence: str | ProFormaAnnotation,
    losses: Mapping[str, Sequence[float]],
    max_losses: int,
) -> set[float]:
    """
    Returns a set of applicable losses for a given sequence.

    :param sequence: The sequence to check for losses.
    :type sequence: str
    :param losses: A list of losses to consider.
    :type losses: List[Tuple[str, float]]
    :param max_losses: The maximum number of losses to consider.
    :type max_losses: int

    :return: A set of applicable losses.
    :rtype: Set[float]

    . code-block:: python

        >>> get_losses('AA', {'A': [-10, -5]}, 1)
        {0.0, -5, -10}

        >>> get_losses('AA', {'A': [-10, -5]}, 2)
        {0.0, -20, -15, -10, -5}

        >>> get_losses('AA', {'A': [-10, -5]}, 3)
        {0.0, -25, -20, -15, -10, -5}

    """
    annotation = get_annotation_input(sequence, copy=False)
    return annotation.get_losses(losses=losses, max_losses=max_losses)


@overload
def fragment(  # type: ignore
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[IonTypeLiteral],
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
) -> Generator[float, None, None]: ...


@overload
def fragment(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[IonTypeLiteral],
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
        Literal["mz-label"] | Literal[FragmentReturnType.MZ_LABEL]
    ) = FragmentReturnType.MZ_LABEL,
    _mass_components: Sequence[float] | None = None,
) -> Generator[tuple[float, str], None, None]: ...


@overload
def fragment(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[IonTypeLiteral],
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
) -> Generator[float, None, None]: ...


def fragment(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[IonTypeLiteral],
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
) -> Generator[float | tuple[float, str], None, None]:
    """
    Builds fragment ions from a given input sequence.

    :param sequence: The amino acid sequence or ProForma annotation.
    :param ion_types: Ion types to consider, e.g., ['b', 'y']
    :param charges: Charge states for the fragments.
    :param monoisotopic: If True, use monoisotopic masses. If False, use average masses.
    :param isotopes: Isotope offsets to consider.
    :param water_loss: If True, consider water loss.
    :param ammonia_loss: If True, consider ammonia loss.
    :param losses: Mapping of neutral losses to consider.
    :param max_losses: The maximum number of losses to consider.
    :param return_type: The type of data to return ('mz' or 'mz-label').
    :param precision: Number of decimal places to round values.

    :return: Generator of m/z values or (m/z, label) tuples.

    .. code-block:: python

        >>> list(fragment('TIDE', 'b', 1, return_type='mz', precision=3))
        [459.209, 330.166, 215.139, 102.055]

        >>> list(fragment('TIDE', 'b', 1, return_type='mz-label', precision=3))
        [(459.209, '+b4'), (330.166, '+b3'), (215.139, '+b2'), (102.055, '+b1')]
    """

    return_type_enum: FragmentReturnType = FragmentReturnType(return_type)

    if (
        return_type_enum != FragmentReturnType.MZ
        and return_type_enum != FragmentReturnType.MZ_LABEL
    ):
        raise TypeError(f"return type not supported: {return_type_enum}")

    return cast(
        Generator[Any, None, None],
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
        ),
    )
