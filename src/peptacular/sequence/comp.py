from collections.abc import Sequence
from typing import overload

from ..constants import IonType, IonTypeLiteral, ParrallelMethod, ParrallelMethodLiteral
from ..proforma import ProFormaAnnotation
from .parrallel import parallel_apply_internal
from .util import get_annotation_input, override_annotation_properties


def _comp_single(
    sequence: str | ProFormaAnnotation,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    estimate_delta: bool = False,
    charge: int | None = None,
    isotope: int = 0,
    use_isotope_on_mods: bool = False,
) -> dict[str, int | float]:
    """Calculate composition for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=True)

    override_annotation_properties(
        annotation=annotation,
        charge=charge,
    )

    return annotation.comp(
        ion_type=ion_type,
        isotope=isotope,
        use_isotope_on_mods=use_isotope_on_mods,
    )


@overload
def comp(
    sequence: str | ProFormaAnnotation,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    estimate_delta: bool = False,
    charge: int | None = None,
    isotope: int = 0,
    use_isotope_on_mods: bool = False,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> dict[str, int | float]: ...


@overload
def comp(
    sequence: Sequence[str | ProFormaAnnotation],
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    estimate_delta: bool = False,
    charge: int | None = None,
    isotope: int = 0,
    use_isotope_on_mods: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[dict[str, int | float]]: ...


def comp(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    estimate_delta: bool = False,
    charge: int | None = None,
    isotope: int = 0,
    use_isotope_on_mods: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> dict[str, int | float] | list[dict[str, int | float]]:
    """
    Calculates the elemental composition of a peptide sequence, including modifications,
    and optionally estimates the composition based on the delta mass from modifications.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param ion_type: The type of ion. Default is IonType.PRECURSOR.
    :type ion_type: IonTypeLiteral | IonType
    :param estimate_delta: If True, estimate the composition based on the delta mass from modifications.
    Default is False.
    :type estimate_delta: bool
    :param charge: The charge state of the ion. Default is None.
    :type charge: int | None
    :param isotope: The number of Neutrons to add/subtract from the final mass. Default is 0.
    :type isotope: int
    :param use_isotope_on_mods: If True, use the isotope on the modifications. Default is False.
    :type use_isotope_on_mods: bool
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If delta_mass is nonzero and estimate_delta is False, indicating an unaccounted modification.

    :return: The elemental composition of the peptide sequence, or list of compositions.
    :rtype: dict[str, int | float] | list[dict[str, int | float]]

    .. code-block:: python

        # Single sequence
        >>> comp('PEPTIDE')
        {'C': 34, 'H': 53, 'N': 7, 'O': 15}

        >>> comp('PEPTIDE[1.0]', estimate_delta=True)['C']
        34.04446833455479

        >>> comp('PEPTIDE[1.0]', estimate_delta=False)  # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        ValueError: Non-zero delta mass (1.0) encountered without estimation enabled for sequence 'PEPTIDE[1.00]'.

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN', 'SEQUENCE']
        >>> compositions = comp(sequences)
        >>> len(compositions)
        3

        # With parallel processing parameters
        >>> sequences = ['PEPTIDE'] * 10000
        >>> compositions = comp(sequences, n_workers=8, method='thread')
        >>> len(compositions)
        10000

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
            estimate_delta=estimate_delta,
            charge=charge,
            isotope=isotope,
            use_isotope_on_mods=use_isotope_on_mods,
        )
    else:
        return _comp_single(
            sequence=sequence,
            ion_type=ion_type,
            estimate_delta=estimate_delta,
            charge=charge,
            isotope=isotope,
            use_isotope_on_mods=use_isotope_on_mods,
        )


def condense_to_mass_mods(
    sequence: str | ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """
    Converts all modifications in a sequence to their mass equivalents by calculating
    the mass difference between modified and unmodified segments.

    :param sequence: The sequence or ProFormaAnnotation object to convert.
    :type sequence: Union[str, ProFormaAnnotation]
    :param include_plus: Whether to include the plus sign for positive mass modifications.
    :type include_plus: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The sequence with all modifications converted to mass modifications.
    :rtype: str

    .. code-block:: python

        >>> condense_to_mass_mods('PEP[Phospho]TIDE', include_plus=False, precision=3)
        'PEP[79.966]TIDE'

        >>> condense_to_mass_mods('PEP[Phospho]TIDE', include_plus=True, precision=3)
        'PEP[+79.966]TIDE'

        >>> condense_to_mass_mods('[Acetyl]-PEPTIDE', precision=3)
        '[42.011]-PEPTIDE'

        >>> condense_to_mass_mods('PEPTIDE-[Amidated]', precision=3)
        'PEPTIDE-[-0.984]'

        >>> condense_to_mass_mods('<13C>PEP[Phospho]TIDE', precision=3)
        'P[5.017]E[5.017]P[84.983]T[4.013]I[6.020]D[4.013]E[5.017]'

    """
    annotation = get_annotation_input(sequence=sequence, copy=True)

    return annotation.condense_to_delta_mass(
        include_plus=include_plus,
        inplace=True,
    ).serialize(include_plus=include_plus, precision=precision)
