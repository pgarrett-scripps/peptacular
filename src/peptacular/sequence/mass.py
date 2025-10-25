from collections.abc import Sequence
from typing import overload

from ..constants import IonType, IonTypeLiteral, ParrallelMethod, ParrallelMethodLiteral
from ..proforma import ProFormaAnnotation
from .parrallel import parallel_apply_internal
from .util import get_annotation_input


def _mass_single(
    sequence: str | ProFormaAnnotation,
    charge: int | None = None,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
) -> float:
    """Calculate mass for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.mass(
        ion_type=ion_type,
        charge=charge,
        monoisotopic=monoisotopic,
        isotope=isotope,
        loss=loss,
        precision=precision,
        use_isotope_on_mods=use_isotope_on_mods,
    )


# Type hints for single sequence
@overload
def mass(
    sequence: str | ProFormaAnnotation,
    charge: int | None = None,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> float: ...


# Type hints for list of sequences
@overload
def mass(
    sequence: Sequence[str | ProFormaAnnotation],
    charge: int | None = None,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[float]: ...


# Implementation
def mass(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    charge: int | None = None,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> float | list[float]:
    """
    Calculate the mass of an amino acid 'sequence'.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param charge: The charge state, default is None.
    :type charge: int | None
    :param ion_type: The ion type. Default is IonType.PRECURSOR.
    :type ion_type: IonTypeLiteral | IonType
    :param monoisotopic: If True, use monoisotopic mass else use average mass. Default is True.
    :type monoisotopic: bool
    :param isotope: The number of Neutrons to add/subtract from the final mass. Default is 0.
    :type isotope: int
    :param loss: The loss to add/subtract to the final mass. Default is 0.0.
    :type loss: float
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None
    :param use_isotope_on_mods: Whether to use isotopes on modifications. Default is False.
    :type use_isotope_on_mods: bool
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raise ValueError: If the ion type is not supported.
    :raise UnknownAminoAcidError: If an unknown amino acid is encountered.
    :raise AmbiguousAminoAcidError: If an ambiguous amino acid is encountered.

    :return: The mass of the sequence, or list of masses for list input.
    :rtype: float | list[float]

    .. code-block:: python

        # Single sequence (original behavior)
        >>> mass('PEPTIDE', precision=3)
        799.36

        >>> mass('PEPTIDE', charge=2, precision=3)
        801.375

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN', 'SEQUENCE']
        >>> masses = mass(sequences, precision=3)
        >>> len(masses)
        3

        # With parallel processing parameters
        >>> sequences = ['PEPTIDE'] * 100000
        >>> masses = mass(sequences, charge=2, precision=3, n_workers=16)
        >>> len(masses)
        100000

        # Auto-detects best method (thread if GIL disabled, process otherwise)
        >>> masses = mass(sequences, charge=2, precision=3)

        # Force threading
        >>> masses = mass(sequences, charge=2, precision=3, method='thread')

        # Force multiprocessing
        >>> masses = mass(sequences, charge=2, precision=3, method='process')

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
            isotope=isotope,
            loss=loss,
            precision=precision,
            use_isotope_on_mods=use_isotope_on_mods,
        )
    else:
        # Single sequence processing
        return _mass_single(
            sequence=sequence,
            charge=charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            loss=loss,
            precision=precision,
            use_isotope_on_mods=use_isotope_on_mods,
        )


# Single sequence implementation for mz
def _mz_single(
    sequence: str | ProFormaAnnotation,
    charge: int | None = None,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
) -> float:
    """Calculate m/z for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.mz(
        ion_type=ion_type,
        charge=charge,
        monoisotopic=monoisotopic,
        isotope=isotope,
        loss=loss,
        precision=precision,
        use_isotope_on_mods=use_isotope_on_mods,
    )


# Type hints for mz
@overload
def mz(
    sequence: str | ProFormaAnnotation,
    charge: int | None = None,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> float: ...


@overload
def mz(
    sequence: Sequence[str | ProFormaAnnotation],
    charge: int | None = None,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[float]: ...


def mz(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    charge: int | None = None,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> float | list[float]:
    """
    Calculate the m/z (mass-to-charge ratio) of an amino acid 'sequence'.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param charge: The charge state, default is None.
    :type charge: int | None
    :param ion_type: The ion type. Default is 'p'.
    :type ion_type: str
    :param monoisotopic: If True, use monoisotopic mass else use average mass. Default is True.
    :type monoisotopic: bool
    :param isotope: The number of Neutrons to add/subtract from the final mass. Default is 0.
    :type isotope: int
    :param loss: The loss to add/subtract to the final mass. Default is 0.0.
    :type loss: float
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None
    :param use_isotope_on_mods: Whether to use isotopes on modifications. Default is False.
    :type use_isotope_on_mods: bool
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raise ValueError: If the ion type is not supported.

    :return: The Mass to Charge ratio (m/z) of the sequence, or list of m/z values.
    :rtype: float | list[float]

    .. code-block:: python

        # Single sequence
        >>> mz('PEPTIDE', charge=1, precision=3)
        800.367

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN', 'SEQUENCE']
        >>> mz_values = mz(sequences, charge=1, precision=3)
        >>> len(mz_values)
        3

        # With parallel processing parameters
        >>> sequences = ['PEPTIDE'] * 100000
        >>> mz_values = mz(sequences, charge=2, precision=3, n_workers=16)

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
            isotope=isotope,
            loss=loss,
            precision=precision,
            use_isotope_on_mods=use_isotope_on_mods,
        )
    else:
        return _mz_single(
            sequence=sequence,
            charge=charge,
            ion_type=ion_type,
            monoisotopic=monoisotopic,
            isotope=isotope,
            loss=loss,
            precision=precision,
            use_isotope_on_mods=use_isotope_on_mods,
        )


# Single sequence implementation for condense_to_mass_mods
def _condense_to_mass_mods_single(
    sequence: str | ProFormaAnnotation,
    use_isotope_on_mods: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """Condense to mass mods for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.condense_to_delta_mass(
        use_isotope_on_mods=use_isotope_on_mods,
        include_plus=include_plus,
        inplace=True,
    ).serialize(include_plus=include_plus, precision=precision)


@overload
def condense_to_mass_mods(
    sequence: str | ProFormaAnnotation,
    use_isotope_on_mods: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str: ...


@overload
def condense_to_mass_mods(
    sequence: Sequence[str | ProFormaAnnotation],
    use_isotope_on_mods: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


def condense_to_mass_mods(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    use_isotope_on_mods: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Converts all modifications in a sequence to their mass equivalents by calculating
    the mass difference between modified and unmodified segments.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences to convert.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param use_isotope_on_mods: Whether to use isotopes on modifications. Default is False.
    :type use_isotope_on_mods: bool
    :param include_plus: Whether to include the plus sign for positive mass modifications.
    :type include_plus: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The sequence(s) with all modifications converted to mass modifications.
    :rtype: str | list[str]

    .. code-block:: python

        # Single sequence
        >>> condense_to_mass_mods('PEP[Phospho]TIDE', include_plus=False, precision=3)
        'PEP[79.966]TIDE'

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEP[Phospho]TIDE', 'PE[Acetyl]PTIDE', 'PEPTIDE[Methyl]']
        >>> condensed = condense_to_mass_mods(sequences, precision=3)
        >>> len(condensed)
        3

        # With parallel processing parameters
        >>> sequences = ['PEP[Phospho]TIDE'] * 10000
        >>> condensed = condense_to_mass_mods(sequences, precision=3, n_workers=8)

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _condense_to_mass_mods_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            use_isotope_on_mods=use_isotope_on_mods,
            include_plus=include_plus,
            precision=precision,
        )
    else:
        return _condense_to_mass_mods_single(
            sequence=sequence,
            use_isotope_on_mods=use_isotope_on_mods,
            include_plus=include_plus,
            precision=precision,
        )
