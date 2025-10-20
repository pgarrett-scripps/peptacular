from typing import Iterable, Mapping, overload, Any
from collections.abc import Sequence
from ..proforma.annotation import ProFormaAnnotation
from .util import get_annotation_input
from ..mod import Mod
from .parrallel import parallel_apply_internal
from ..constants import ModType, ModTypeLiteral


MOD_BUILDER_INPUT_TYPE = Mapping[str, Iterable[str | float | int | Mod]]


def _build_mods_single(
    sequence: ProFormaAnnotation | str,
    nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
    nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
    precision: int | None = None,
    include_plus: bool = False,
) -> list[str]:
    """Build modifications for a single sequence"""
    annotation = get_annotation_input(sequence, copy=True)
    return [
        annot.serialize(precision=precision, include_plus=include_plus)
        for annot in annotation.build_mods(
            nterm_static=nterm_static,
            cterm_static=cterm_static,
            internal_static=internal_static,
            labile_static=labile_static,
            nterm_variable=nterm_variable,
            cterm_variable=cterm_variable,
            internal_variable=internal_variable,
            labile_variable=labile_variable,
            max_variable_mods=max_variable_mods,
            use_regex=use_regex,
            inplace=False,
        )
    ]


@overload
def build_mods(
    sequence: ProFormaAnnotation | str,
    nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
    nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
    precision: int | None = None,
    include_plus: bool = False,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def build_mods(
    sequence: Sequence[ProFormaAnnotation | str],
    nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
    nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
    precision: int | None = None,
    include_plus: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def build_mods(
    sequence: ProFormaAnnotation | str | Sequence[ProFormaAnnotation | str],
    nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
    nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
    precision: int | None = None,
    include_plus: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Build modified sequences by applying static and variable modifications to a sequence or list of sequences.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence to modify, ProFormaAnnotation, or list of sequences.
    :type sequence: ProFormaAnnotation | str | list[ProFormaAnnotation | str]
    :param nterm_static: Static N-terminal modifications.
    :type nterm_static: MOD_BUILDER_INPUT_TYPE | None
    :param cterm_static: Static C-terminal modifications.
    :type cterm_static: MOD_BUILDER_INPUT_TYPE | None
    :param internal_static: Static internal modifications.
    :type internal_static: MOD_BUILDER_INPUT_TYPE | None
    :param labile_static: Static labile modifications.
    :type labile_static: MOD_BUILDER_INPUT_TYPE | None
    :param nterm_variable: Variable N-terminal modifications.
    :type nterm_variable: MOD_BUILDER_INPUT_TYPE | None
    :param cterm_variable: Variable C-terminal modifications.
    :type cterm_variable: MOD_BUILDER_INPUT_TYPE | None
    :param internal_variable: Variable internal modifications.
    :type internal_variable: MOD_BUILDER_INPUT_TYPE | None
    :param labile_variable: Variable labile modifications.
    :type labile_variable: MOD_BUILDER_INPUT_TYPE | None
    :param max_variable_mods: Maximum number of variable modifications per sequence.
    :type max_variable_mods: int
    :param use_regex: Whether to use regex for modification matching.
    :type use_regex: bool
    :param precision: Number of decimal places for mass values.
    :type precision: int | None
    :param include_plus: Whether to include plus signs for positive modifications.
    :type include_plus: bool
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: List of modified sequences as strings, or list of lists for multiple input sequences.
    :rtype: list[str] | list[list[str]]

    .. code-block:: python

        # Single sequence
        >>> results = build_mods('PEPTIDE', internal_variable={'P': [79.966]}, max_variable_mods=1)
        >>> len(results) > 0
        True

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN', 'SEQUENCE']
        >>> results = build_mods(sequences, internal_variable={'P': [79.966]}, max_variable_mods=1)
        >>> len(results)
        3


    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _build_mods_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            nterm_static=nterm_static,
            cterm_static=cterm_static,
            internal_static=internal_static,
            labile_static=labile_static,
            nterm_variable=nterm_variable,
            cterm_variable=cterm_variable,
            internal_variable=internal_variable,
            labile_variable=labile_variable,
            max_variable_mods=max_variable_mods,
            use_regex=use_regex,
            precision=precision,
            include_plus=include_plus,
        )
    else:
        return _build_mods_single(
            sequence=sequence,
            nterm_static=nterm_static,
            cterm_static=cterm_static,
            internal_static=internal_static,
            labile_static=labile_static,
            nterm_variable=nterm_variable,
            cterm_variable=cterm_variable,
            internal_variable=internal_variable,
            labile_variable=labile_variable,
            max_variable_mods=max_variable_mods,
            use_regex=use_regex,
            precision=precision,
            include_plus=include_plus,
        )


def get_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | ModTypeLiteral | None = None,
) -> dict[str, Any]:
    """
    Parses a sequence with modifications and returns a dictionary where keys represent the position/type of the
    modifications.

    Internal modifications are mapped to the unmodified sequence by their index in the unmodified sequence.

    Special modifications are mapped to the sequence by the following:
    - N-terminal modifications are mapped to the index 'nterm'
    - C-terminal modifications are mapped to the index 'cterm'
    - Isotopic modifications are mapped to the index 'isotope'
    - Static modifications are mapped to the index 'static'
    - Labile modifications are mapped to the index 'labile'
    - Unknown modifications are mapped to the index 'unknown'
    - Intervals are mapped to the index 'interval'
    - charge state is mapped to the index 'charge'
    - charge adducts are mapped to the index 'charge_adducts'

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A dictionary with the modifications
    :rtype: ModDict

    .. code-block:: python

        # All modifications will be returned as Mod objects which contain the modification value and multiplier
        >>> get_mods('PEP[Phospho]T[1]IDE[-3.14]')['internal']
        {2: ('Phospho',), 3: (1,), 6: (-3.14,)}

        >>> get_mods('PEP[Phospho][1.0]TIDE')['internal']
        {2: ('Phospho', 1.0)}

        # N-terminal modifications are mapped to the index 'nterm'
        >>> get_mods('[Acetyl]-PEPTIDE')['nterm']
        ('Acetyl',)

        # C-terminal modifications are mapped to the index 'cterm'
        >>> get_mods('PEPTIDE-[Amide]')['cterm']
        ('Amide',)

        # Isotopic modifications are mapped to the index 'isotope'
        >>> get_mods('<13C>PEPTIDE')['isotope']
        ('13C',)

        # Static modifications are mapped to the index 'static'
        >>> get_mods('<[+1.234]@P>PEPTIDE')['static']
        ('[+1.234]@P',)

        # Labile modifications are mapped to the index 'labile'
        >>> get_mods('{Glycan:Hex}PEPTIDE')['labile']
        ('Glycan:Hex',)

        # Unknown modifications are mapped to the index 'unknown'
        >>> get_mods('[Phospho]^3?PEPTIDE')['unknown']
        ('Phospho', 'Phospho', 'Phospho')

        # Intervals are mapped to the index 'interval'
        >>> get_mods('PEP(TI)[Phospho]DE')['interval']
        (ModInterval(start=3, end=5, ambiguous=False, mods=('Phospho',)),)

        # Charge state is mapped to the index 'charge'
        >>> get_mods('PEPTIDE/+2')['charge']
        2

        # Charge adducts are mapped to the index 'charge_adducts'
        >>> get_mods('PEPTIDE/+2[+2Na+,-H+]')['charge_adducts']
        ('+2Na+,-H+',)

    """

    return get_annotation_input(sequence, copy=True).get_mods(mods)


def add_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any],
    append: bool = True,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """
    Adds modifications to the given sequence. The modifications can be of type Mod, str, int, or float, and can be
    a single value or a list of values. The modifications will be added to the sequence in the order they are provided.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]
    :param mods: Dictionary representing the modifications to be added to the sequence.
    :type mods: Dict
    :param append: If True, the modifications will be appended to the existing modifications.
                      If False, the existing modifications will be replaced. Defaults to True.
    :type append: bool
    :param include_plus: If True, the modifications will be serialized with a '+' sign for positive values.
    :type include_plus: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The peptide sequence with the specified modifications.
    :rtype: str

    .. code-block:: python

        >>> from peptacular import Mod

        # Add internal modifications to an unmodified peptide
        >>> add_mods('PEPTIDE', {2: [Mod('phospho', 1)]})
        'PEP[phospho]TIDE'

        # Can also add N and C terminal modifications
        >>> add_mods('PEPTIDE', {'nterm': 'Acetyl', 6: 1.234, 'cterm': 'Amide'})
        '[Acetyl]-PEPTIDE[1.234]-[Amide]'

        >>> add_mods('PEPTIDE', {'nterm': 'Acetyl', 6: 1.234, 'cterm': 'Amide'}, include_plus=True)
        '[Acetyl]-PEPTIDE[+1.234]-[Amide]'

        # Can also add isotopic modifications
        >>> add_mods('PEPTIDE', {'isotope': ['13C', '15N']})
        '<13C><15N>PEPTIDE'

        # Can also add static modifications
        >>> add_mods('PEPTIDE', {'static': '[+1.234]@P'})
        '<[+1.234]@P>PEPTIDE'

        # Can also add labile modifications
        >>> add_mods('PEPTIDE', {'labile': 'Glycan:Hex'})
        '{Glycan:Hex}PEPTIDE'

        # Can also add unknown modifications
        >>> add_mods('PEPTIDE', {'unknown': Mod('Phospho', 3)})
        '[Phospho]^3?PEPTIDE'

        # Can also add intervals
        >>> add_mods('PEPTIDE', {'interval': (3, 5, False, 'Phospho')})
        'PEP(TI)[Phospho]DE'

        # Can also add charge state
        >>> add_mods('PEPTIDE', {'charge': 2})
        'PEPTIDE/2'

        # Can also add charge adducts
        >>> add_mods('PEPTIDE', {'charge': 2, 'charge_adducts': '+2Na+,-H+'})
        'PEPTIDE/2[+2Na+,-H+]'

    """

    annotation = get_annotation_input(sequence, copy=True)
    annotation.add_mods(mods, inplace=True, append=append)
    return annotation.serialize(include_plus=include_plus, precision=precision)


def set_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any] | None,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    return (
        get_annotation_input(sequence, copy=True)
        .set_mods(mods, inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


def append_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any],
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    return (
        get_annotation_input(sequence, copy=True)
        .append_mods(mods, inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


def extend_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any],
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    return (
        get_annotation_input(sequence, copy=True)
        .extend_mods(mods, inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


def condense_static_mods(
    sequence: str | ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """
    Condenses static modifications into internal modifications.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The peptide sequence with condensed static modifications.
    :rtype: str

    .. code-block:: python

        # Condenses static modifications to specified internal modifications
        >>> condense_static_mods('<13C><[100]@P>PEPTIDE')
        '<13C>P[100]EP[100]TIDE'

        # If residue is already modified, the static modification will be appended
        >>> condense_static_mods('<13C><[100]@P>P[10]EPTIDE')
        '<13C>P[10][100]EP[100]TIDE'

        # Example for unmodified sequences
        >>> condense_static_mods('PEPTIDE')
        'PEPTIDE'

        # Example for N-Term static modifications
        >>> condense_static_mods('<[Oxidation]@N-Term>PEPTIDE')
        '[Oxidation]-PEPTIDE'

        # Example for C-Term static modifications
        >>> condense_static_mods('<[Oxidation]@C-Term>PEPTIDE')
        'PEPTIDE-[Oxidation]'

    """

    return (
        get_annotation_input(sequence=sequence, copy=True)
        .condense_static_mods(inplace=False)
        .serialize(include_plus=include_plus, precision=precision)
    )


def pop_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
    include_plus: bool = False,
    precision: int | None = None,
) -> tuple[str, dict[str, Any]]:
    """
    Removes all modifications from the given sequence, returning the unmodified sequence and a dictionary of the
    removed modifications.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A tuple containing the unmodified sequence and a dictionary of the removed modifications.
    :rtype: Tuple[str, ModDict]

    .. code-block:: python

        # Simply combines the functionality of strip_modifications and get_modifications
        >>> pop_mods('PEP[phospho]TIDE')[0]
        'PEPTIDE'
        >>> pop_mods('PEP[phospho]TIDE')[1]['internal']
        {2: ('phospho',)}

        # can specify which modifications to pop
        >>> pop_mods('PEP[phospho]TIDE-[+100]', mods=['internal'])
        ('PEPTIDE-[100]', {'internal': {2: ('phospho',)}})

    """
    annotation = get_annotation_input(sequence=sequence, copy=True)
    mod_dict = annotation.pop_mods(mod_types=mods)  # only include keys that have mods
    return (
        annotation.serialize(include_plus=include_plus, precision=precision),
        mod_dict,
    )


def strip_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """
    Strips all modifications from the given sequence, returning the unmodified sequence.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The stripped sequence
    :rtype: str

    .. code-block:: python

        # Removes internal modifications:
        >>> strip_mods('PEP[phospho]TIDE')
        'PEPTIDE'

        # Also removes N and C terminal modifications:
        >>> strip_mods('[Acetyl]-PEPTIDE[1.234]-[Amide]')
        'PEPTIDE'

        # Also remove labile modifications:
        >>> strip_mods('{1.0}[Acetyl]-PEPTIDE[1.234]-[Amide]')
        'PEPTIDE'

        # Also remove isotope notations:
        >>> strip_mods('<C13>[Acetyl]-PEPTIDE[1.234]-[Amide]')
        'PEPTIDE'

        # Using a sequence without modifications will return the same sequence:
        >>> strip_mods('PEPTIDE')
        'PEPTIDE'

        >>> strip_mods('PEP[Formula:[13C]H12]TIDE')
        'PEPTIDE'

        >>> strip_mods('(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K')
        'DQNGTWEMESNENFEGYMK'

        >>> strip_mods('[1][2]^2?[100]^3-PEP[1]^2TIDE')
        'PEPTIDE'

    """

    annotation = get_annotation_input(sequence=sequence, copy=True)

    return annotation.remove_mods(mods=mods, inplace=True).serialize(
        include_plus=include_plus, precision=precision
    )


def filter_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """
    Keeps only the specified modifications in the sequence, removing all others.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]
    :param mods: The modifications to keep. If None, all modifications will be kept.
    :type mods: Optional[Union[str, List[str]]]
    :param include_plus: If True, the modifications will be serialized with a '+' sign for positive values.
    :type include_plus: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The sequence with only the specified modifications kept.
    :rtype: str

    .. code-block:: python

        # Keeps only internal modifications:
        >>> filter_mods('PEP[phospho]TIDE', mods='internal')
        'PEP[phospho]TIDE'

        # Keeps only N and C terminal modifications:
        >>> filter_mods('[Acetyl]-PEPTIDE[1.234]-[Amide]', mods=['nterm', 'cterm'])
        '[Acetyl]-PEPTIDE-[Amide]'

        # Keeps only labile modifications:
        >>> filter_mods('{1.0}[Acetyl]-PEPTIDE[1.234]-[Amide]', mods='labile')
        '{1.0}PEPTIDE'

        # Keeps only isotope notations:
        >>> filter_mods('<C13>[Acetyl]-PEPTIDE[1.234]-[Amide]', mods='isotope')
        '<C13>PEPTIDE'

        # Using a sequence without modifications will return the same sequence:
        >>> filter_mods('PEPTIDE')
        'PEPTIDE'

    """
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .filter_mods(mods=mods, inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )
