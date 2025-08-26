"""
fragmentation.py contains functions for fragmenting peptides
"""


from collections.abc import Sequence
from typing import Any, Generator, Literal, Mapping, cast, overload

from peptacular.constants import IonTypeLiteral
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

        >>> get_losses('AA', [('A', -10), ('A', -5)], 1)
        {0.0, -5, -10}

        >>> get_losses('AA', [('A', -10), ('A', -5)], 2)
        {0.0, -15, -10, -5, -20}

        >>> get_losses('AA', [('A', -10), ('A', -5)], 3)
        {0.0, -15, -10, -25, -5, -20}

    """
    annotation = get_annotation_input(sequence, copy=False)
    return annotation.get_losses(losses=losses, max_losses=max_losses)

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
    losses: Mapping[str, float] | None = None,
    max_losses: int = 1,
    precision: int | None = None,
    return_type: Literal['mz'],
    _mass_components: Sequence[float] | None = None,
) -> Generator[float, None, None]: ...
    
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
    losses: Mapping[str, float] | None = None,
    max_losses: int = 1,
    precision: int | None = None,
    return_type: Literal['mz-label'],
    _mass_components: Sequence[float] | None = None,
) -> Generator[float, None, None]: ...

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
    losses: Mapping[str, float] | None = None,
    max_losses: int = 1,
    precision: int | None = None,
    return_type: Literal['mz','mz-label'] = 'mz'
) -> Generator[float, None, None]: ...

def fragment(
    sequence: str | ProFormaAnnotation,
    ion_types: Sequence[IonTypeLiteral] | IonTypeLiteral,
    charges: Sequence[int] | int,
    monoisotopic: bool = True,
    *,
    isotopes: Sequence[int] | int = 0,
    water_loss: bool = False,
    ammonia_loss: bool = False,
    losses: Mapping[str, float] | None = None,
    max_losses: int = 1,
    precision: int | None = None,
    return_type: Literal['mz','mz-label'] = 'mz',
    _mass_components: Sequence[float] | None = None,
) -> Generator[float | tuple[float, str], None, None]:
    """
    Builds all Fragment objects or a given input 'sequence'.

    :param sequence: The amino acid sequence or ProForma annotation.
    :type sequence: str | ProFormaAnnotation
    :param ion_types: A list of ion types to consider, e.g., ['b', 'y'], or a single ion type, e.g., 'b'.
    :type ion_types: List[str] | str
    :param charges: A list of charge states for the fragments, or a single charge state.
    :type charges: List[int] | int
    :param monoisotopic: If True, use monoisotopic masses. If False, use average masses, default is [True].
    :type monoisotopic: bool
    :param isotopes: A list of isotope offsets to consider, or a single isotope offset, default is [0].
    :type isotopes: List[int] | int
    :param water_loss: If True, consider water loss, default is False.
    :type water_loss: bool
    :param ammonia_loss: If True, consider ammonia loss, default is False.
    :type ammonia_loss: bool
    :param losses: A list of neutral losses to consider, or a single neutral loss, default is [0.0].
    :type losses: List[float] | float
    :param max_losses: The maximum number of losses to consider, default is 1.
    :type max_losses: int
    :param return_type: The type of data to return, either 'fragment', 'mass', or 'mz', default is 'fragment'.
    :type return_type: str
    :param precision: The number of decimal places to round the masses or m/z values to, default is None.
    :type precision: int
    :param _mass_components: The split mass components. Used for more efficient fragmenting.
    :type _mass_components: List[float]

    :return: List of Fragment objects or a list of masses or m/z values.
    :rtype: List[Fragment] | List[float]

    .. code-block:: python

        # By default a fragment object is returned
        >>> len(fragment("[1.0]-P[2.0]E[3.0]-[4.0]", 'y', 1))
        2
        >>> fragment("[1.0]-P[2.0]E[3.0]-[4.0]", 'y', 1)[0].sequence
        '[1.0]-P[2.0]E[3.0]-[4.0]'
        >>> fragment("[1.0]-P[2.0]E[3.0]-[4.0]", 'y', 1)[1].sequence
        'E[3.0]-[4.0]'

        # Charge 1 B-Ions mass
        >>> fragment(sequence='TIDE', ion_types='b', charges=1, return_type='mass', precision=3)
        [459.209, 330.166, 215.139, 102.055]

        # Charge 1 B-Ions mz
        >>> fragment(sequence='TIDE', ion_types='b', charges=1, return_type='mz', precision=3)
        [459.209, 330.166, 215.139, 102.055]

        # Charge 2 B-Ions mz
        >>> fragment(sequence='TIDE', ion_types='b', charges=2, return_type='mz', precision=3)
        [230.108, 165.587, 108.073, 51.531]

        # Charge 1 Y-Ions mass
        >>> fragment(sequence='T[10]IDE', ion_types='y', charges=1, return_type='mass', precision=3)
        [487.219, 376.171, 263.087, 148.06]

        # Get fragment objects by default
        >>> fragments = fragment(sequence='TIDE', ion_types="i", charges=1, return_type='fragment')
        >>> list(map(lambda frag: frag.sequence, fragments))
        ['T', 'I', 'D', 'E']

        # Immonium ions
        >>> fragment(sequence='TIDE', ion_types="i", charges=1, return_type='mz', precision=3)
        [74.06, 86.096, 88.039, 102.055]

        # Internal fragment ions and 'mz-label' return type
        >>> fragment(sequence='TIDE', ion_types="by", charges=1, return_type='mz-label', precision=3)
        [(114.091, '+by1-2'), (229.118, '+by1-3'), (116.034, '+by2-3')]

        # regex loss losses
        >>> fragment(sequence='TIDE', ion_types="b", charges=2, return_type='mz', precision=3, losses=('E', -10))
        [230.108, 225.108, 165.587, 108.073, 51.531]

        # Multiple losses (when max_losses == 1, only one loss is considered at a time)
        >>> l = [('A', -10), ('A', -5)]
        >>> fragment('AA', "b", 1, return_type='mz', precision=3, losses=l)
        [143.082, 138.082, 133.082, 72.044, 67.044, 62.044]
        >>> fragment('AA', "b", 1, return_type='label', precision=3, losses=l)
        ['+b2', '+b2(-5)', '+b2(-10)', '+b1', '+b1(-5)', '+b1(-10)']

        # Multiple losses (when max_losses > 1, all unique combinations of losses are considered)
        >>> l = [('A', -10), ('A', -5)]
        >>> fragment('AA', "b", 1, return_type='mz', precision=3, losses=l, max_losses=2)
        [143.082, 128.082, 133.082, 138.082, 123.082, 72.044, 57.044, 67.044, 62.044]
        >>> fragment('AA', "b", 1, return_type='label', precision=3, losses=l, max_losses=2)
        ['+b2', '+b2(-15)', '+b2(-10)', '+b2(-5)', '+b2(-20)', '+b1', '+b1(-15)', '+b1(-5)', '+b1(-10)']

        # Water and ammonia losses
        >>> fragment('AQE', "b", 1, return_type='mz', precision=3, water_loss=True, ammonia_loss=True)
        [329.146, 311.135, 312.119, 200.103, 183.076, 72.044]
        >>> fragment('AQE', "b", 1, return_type='label', precision=3, water_loss=True, ammonia_loss=True)
        ['+b3', '+b3(-18.01056)', '+b3(-17.02655)', '+b2', '+b2(-17.02655)', '+b1']

    """

    return cast(Generator[Any, None, None], get_annotation_input(sequence=sequence, copy=False).fragment(
        ion_types=ion_types,
        charges=charges,
        monoisotopic=monoisotopic,
        isotopes=isotopes,
        water_loss=water_loss,
        ammonia_loss=ammonia_loss,
        losses=losses,
        max_losses=max_losses,
        return_type=return_type,
        precision=precision,
        _mass_components=_mass_components,
    )   
    )
