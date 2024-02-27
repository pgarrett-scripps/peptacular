"""
The `term.py` module provides utilities for handling terminal modifications in amino acid sequences.
These modifications are represented by square brackets either at the N-terminus or the C-terminus.

Key Features:
    - Extract, strip, add, or condense modifications at the N or C terminus.
    - Support for string, float, and integer type modifications.
    - Index-based representation for easy reference.

Term Modification Notation:
    - N-Terminus: [Acetyl]PEPTIDE or [1.234]PEPTIDE
    - C-Terminus: PEPTIDE[Amide] or PEPTIDE[3.1415]
    - N-Terminus modifications are denoted with a -1 index.
    - C-Terminus modifications use the index based on the length of the unmodified sequence.
"""

from typing import Union, Any, Tuple

from peptacular.util import convert_type


def get_n_term_modification_index(sequence: str) -> int:
    """
    Returns the end index of the N-terminal amino acid in the sequence.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: The end index of the N-terminal amino acid in the sequence.
    :rtype: int

    .. code-block:: python

        # For sequences with N-terminal modifications:
        >>> get_n_term_modification_index('[Acetyl]-PEPTIDE')
        9

        >>> "[Acetyl]-PEPTIDE"[get_n_term_modification_index('[Acetyl]-PEPTIDE'):]
        'PEPTIDE'

        # For sequences without N-terminal modifications:
        >>> get_n_term_modification_index('PEPTIDE')
        0

        >>> "PEPTIDE"[get_n_term_modification_index('PEPTIDE'):]
        'PEPTIDE'

    """

    if sequence.startswith('['):
        return sequence.find(']-') + 2

    return 0


def get_c_term_modification_index(sequence: str) -> int:
    """
    Return the start index of the C-terminal amino acid in the sequence.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: The start index of the C-terminal amino acid in the sequence.
    :rtype: int

    .. code-block:: python

        # For sequences with C-terminal modifications:
        >>> get_c_term_modification_index('PEPTIDE-[Amide]')
        7

        >>> get_c_term_modification_index('PEPTIDE[3.14]-[Amide]')
        13

        >>> get_c_term_modification_index('PEPTIDE[3.14]')
        13

        >>> 'PEPTIDE[Amide]'[:get_c_term_modification_index('PEPTIDE-[Amide]')]
        'PEPTIDE'

        # For sequences without C-terminal modifications:
        >>> get_c_term_modification_index('PEPTIDE')
        7

        >>> 'PEPTIDE'[:get_c_term_modification_index('PEPTIDE')]
        'PEPTIDE'

    """

    if sequence.endswith(']'):
        start = sequence.rfind('-[')
        if start != -1:
            return start

    return len(sequence)


def get_n_term_modification(sequence: str) -> Union[str, float, int, None]:
    """
    Extracts the N-terminal modification from the sequence.

    Given a sequence with potential N-terminal notation, this function will parse and return the N-terminal
    modification, if available. If no modification is present, it returns None.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: N-terminal modification if present, else None.
    :rtype: Union[str, float, int, None]

    .. code-block:: python

        # For string-based modifications:
        >>> get_n_term_modification("[Acetyl]-PEPTIDE")
        'Acetyl'

        # For float-based modifications:
        >>> get_n_term_modification("[3.1415]-PEPTIDE")
        3.1415

        # For int-based modifications:
        >>> get_n_term_modification("[100]-PEPTIDE")
        100

        # When no N-Terminus modification is present:
        >>> get_n_term_modification("PEPTIDE") # returns None

    """

    mod = convert_type(sequence[:get_n_term_modification_index(sequence)][1:-2])

    if mod == '':
        return None

    return mod


def get_c_term_modification(sequence: str) -> Union[str, float, int, None]:
    """
    Extracts the C-terminal modification from the provided sequence.

    This function parses a sequence and retrieves the C-terminal modification if it exists.
    If the sequence lacks a C-terminal modification, the function returns None.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: C-terminal modification if present, else None.
    :rtype: Union[str, float, int, None]

    .. code-block:: python

        # For sequences with string-based modifications:
        >>> get_c_term_modification("PEPTIDE-[Amide]")
        'Amide'

        # For sequences with float-based modifications:
        >>> get_c_term_modification("PEPTIDE-[3.1415]")
        3.1415

        # For sequences with int-based modifications:
        >>> get_c_term_modification("PEPTIDE-[100]")
        100

        # When no C-terminal modification is present:
        >>> get_c_term_modification("PEPTIDE")

    """

    mod = convert_type(sequence[get_c_term_modification_index(sequence):][2:-1])

    if mod == '':
        return None

    return mod


def get_term_modifications(sequence: str) -> Tuple[Union[str, float, int, None], Union[str, float, int, None]]:
    """
    Extracts both N-terminal and C-terminal modifications from the sequence.

    This function parses the input sequence and returns both the N-terminal and C-terminal modifications
    as a tuple. If no modifications are present, the function returns a tuple of None values.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: Tuple of N-terminal and C-terminal modifications.
    :rtype: Tuple[Union[str, float, int, None], Union[str, float, int, None]]

    .. code-block:: python

        # For sequences with different modification types (string, float, int):
        >>> get_term_modifications("[Acetyl]-PEPTIDE-[Amide]")
        ('Acetyl', 'Amide')
        >>> get_term_modifications("[3.1415]-PEPTIDE-[100]")
        (3.1415, 100)
        >>> get_term_modifications("[100]-PEPTIDE-[Acetyl]")
        (100, 'Acetyl')

        # For sequences without any terminal modification:
        >>> get_term_modifications("PEPTIDE")
        (None, None)

    """

    n_term = get_n_term_modification(sequence)
    c_term = get_c_term_modification(sequence)
    return n_term, c_term


def strip_n_term_modification(sequence: str) -> str:
    """
    Removes any N-terminal modification notation from the sequence.

    This function takes a sequence as input and returns the sequence without its N-terminal modification
    notation. If no N-terminal modification exists, the original sequence is returned unchanged.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: Sequence devoid of the N-terminal modification notation.
    :rtype: str

    .. code-block:: python

        # For sequences with different modification types (string, float, int):
        >>> strip_n_term_modification("[Acetyl]-PEPTIDE")
        'PEPTIDE'
        >>> strip_n_term_modification("[3.1415]-PEPTIDE")
        'PEPTIDE'
        >>> strip_n_term_modification("[100]-PEPTIDE")
        'PEPTIDE'

        # If a residue modification is present at the N-terminus, only the terminal notation is removed:
        >>> strip_n_term_modification("[Acetyl]-P[1]EPTIDE")
        'P[1]EPTIDE'

        # For sequences without any N-terminal modification:
        >>> strip_n_term_modification("PEPTIDE")
        'PEPTIDE'

    """

    return sequence[get_n_term_modification_index(sequence):]


def strip_c_term_modification(sequence: str) -> str:
    """
    Removes the C-terminal modification notation from the provided sequence.

    This function parses the input sequence and returns the sequence devoid of
    its C-terminal modification notation. If no C-terminal modification is present in
    the sequence, the original sequence is returned unchanged.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: Sequence without the C-terminal modification notation.
    :rtype: str

    .. code-block:: python

        # For sequences with different modification types (string, float, int):
        >>> strip_c_term_modification("PEPTIDE-[Acetyl]")
        'PEPTIDE'
        >>> strip_c_term_modification("PEPTIDE-[3.1415]")
        'PEPTIDE'
        >>> strip_c_term_modification("PEPTIDE-[100]")
        'PEPTIDE'

        # If a residue modification is present at the C-terminus, only the terminal notation is removed:
        >>> strip_c_term_modification("P[1]EPTIDE-[Amide]")
        'P[1]EPTIDE'

        # For sequences without any C-terminal modification:
        >>> strip_c_term_modification("PEPTIDE")
        'PEPTIDE'

    """

    return sequence[:get_c_term_modification_index(sequence)]


def strip_term_modifications(sequence: str) -> str:
    """
    Removes both N-terminal and C-terminal modification notations from the provided sequence.

    This function parses the input sequence and returns the sequence without its N-terminal and C-terminal
    modification notations. If no terminal modifications are present, the original sequence is returned unchanged.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: Sequence without the terminal modification notations.
    :rtype: str

    .. code-block:: python

        # For sequences with different modification types (string, float, int):
        >>> strip_term_modifications("[Acetyl]-PEPTIDE-[Amide]")
        'PEPTIDE'
        >>> strip_term_modifications("[3.1415]-PEPTIDE-[100]")
        'PEPTIDE'
        >>> strip_term_modifications("[100]-PEPTIDE-[Acetyl]")
        'PEPTIDE'

        # If a residue modification is present at both termini, only the terminal notations are removed:
        >>> strip_term_modifications("[Acetyl]-P[1]EPTIDE-[Amide]")
        'P[1]EPTIDE'

        # For sequences without any terminal modification:
        >>> strip_term_modifications("PEPTIDE")
        'PEPTIDE'

    """

    sequence = strip_n_term_modification(sequence)
    sequence = strip_c_term_modification(sequence)
    return sequence


def add_n_term_modification(sequence: str, mod: Any, overwrite: bool = True) -> str:
    """
    Appends the specified N-terminal modification to the provided sequence.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param mod: Modification to be appended at the N-terminus of the sequence.
    :type mod: Any
    :param overwrite: If True, the new modification will overwrite the existing one. If False, the new modification
    will not be added if an existing modification is present.
    :type overwrite: bool

    :return: Modified sequence with the added N-terminal notation.
    :rtype: str

    .. code-block:: python

        # Adding string-based modifications:
        >>> add_n_term_modification("PEPTIDE", "Acetyl")
        '[Acetyl]-PEPTIDE'

        # Adding float-based modifications:
        >>> add_n_term_modification("PEPTIDE", 3.1415)
        '[3.1415]-PEPTIDE'

        # Adding int-based modifications:
        >>> add_n_term_modification("PEPTIDE", 100)
        '[100]-PEPTIDE'

        # by default, the new modification will overwrite the existing one:
        >>> add_n_term_modification("[Acetyl]-PEPTIDE-[Amide]", "Amide")
        '[Amide]-PEPTIDE-[Amide]'

        # If the overwrite parameter is set to False, the new modification will not be added:
        >>> add_n_term_modification("[Acetyl]-PEPTIDE-[Amide]", "Acetyl", overwrite=False)
        '[Acetyl]-PEPTIDE-[Amide]'

        # No change if the modification is None:
        >>> add_n_term_modification("PEPTIDE", None)
        'PEPTIDE'

    """

    existing_mod = get_n_term_modification(sequence)

    if overwrite:

        if mod is None:
            return strip_n_term_modification(sequence)

        return f"[{str(mod)}]-{strip_n_term_modification(sequence)}"

    else:
        if existing_mod is not None:
            return sequence

        else:
            if mod is None:
                return sequence
            return f"[{str(mod)}]-{sequence}"


def add_c_term_modification(sequence: str, mod: Any, overwrite: bool = True) -> str:
    """
    Adds the specified C-terminal modification notation to the provided sequence.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param mod: Modification to be appended at the C-terminus of the sequence. Can be of type str, int, or float.
    :type mod: Any
    :param overwrite: If True, the new modification will overwrite the existing one. If False, the new modification
    will not be added if an existing modification is present.
    :type overwrite: bool

    :return: Modified sequence with the added C-terminal notation.
    :rtype: str

    .. code-block:: python

        # Adding string-based modifications:
        >>> add_c_term_modification("PEPTIDE", "Amide")
        'PEPTIDE-[Amide]'

        # Adding float-based modifications:
        >>> add_c_term_modification("PEPTIDE", 3.1415)
        'PEPTIDE-[3.1415]'

        # Adding int-based modifications:
        >>> add_c_term_modification("PEPTIDE", 100)
        'PEPTIDE-[100]'

        # by default, the new modification will overwrite the existing one:
        >>> add_c_term_modification("[Acetyl]-PEPTIDE-[Amide]", "Acetyl")
        '[Acetyl]-PEPTIDE-[Acetyl]'

        # If the overwrite parameter is set to False, the new modification will not be added:
        >>> add_c_term_modification("[Acetyl]-PEPTIDE-[Amide]", "Acetyl", overwrite=False)
        '[Acetyl]-PEPTIDE-[Amide]'

        # No change if the modification is None:
        >>> add_c_term_modification("PEPTIDE", None)
        'PEPTIDE'

    """

    existing_mod = get_c_term_modification(sequence)

    if overwrite is True:

        if mod is None:
            return strip_c_term_modification(sequence)

        return f"{strip_c_term_modification(sequence)}-[{str(mod)}]"

    else:
        if existing_mod is not None:
            return sequence

        else:
            if mod is None:
                return sequence

            return f"{sequence}-[{str(mod)}]"


def add_term_modifications(sequence: str, n_term_mod: Any, c_term_mod: Any, overwrite: bool = True) -> str:
    """
    Appends both N-terminal and C-terminal modification notations to the provided sequence.

    If the sequence already contains terminal modifications, the new modifications will be combined with the existing
    ones. Ensure that the types of the modifications are compatible to prevent errors.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param n_term_mod: Modification to be appended at the N-terminus of the sequence.
    :type n_term_mod: Any
    :param c_term_mod: Modification to be appended at the C-terminus of the sequence.
    :type c_term_mod: Any
    :param overwrite: If True, the new modifications will overwrite the existing ones. If False, the new modifications
    will not be added if existing modifications are present.
    :type overwrite: bool

    :return: Modified sequence with the added terminal notations.
    :rtype: str

    .. code-block:: python

        # Adding string-based modifications:
        >>> add_term_modifications("PEPTIDE", "Acetyl", "Amide")
        '[Acetyl]-PEPTIDE-[Amide]'

        # Adding float-based modifications:
        >>> add_term_modifications("PEPTIDE", 3.1415, 3.1415)
        '[3.1415]-PEPTIDE-[3.1415]'

        # Adding int-based modifications:
        >>> add_term_modifications("PEPTIDE", 100, 100)
        '[100]-PEPTIDE-[100]'


        # by default, the new modification will overwrite the existing one:
        >>> add_term_modifications("[Acetyl]-PEPTIDE-[Amide]", "Amide", "Acetyl")
        '[Amide]-PEPTIDE-[Acetyl]'

        # If the overwrite parameter is set to False, the new modification will not be added:
        >>> add_term_modifications("[Acetyl]-PEPTIDE-[Amide]", "Amide", "Acetyl", overwrite=False)
        '[Acetyl]-PEPTIDE-[Amide]'

        # No change if the modifications are None:
        >>> add_term_modifications("PEPTIDE", None, None)
        'PEPTIDE'

    """

    return add_n_term_modification(add_c_term_modification(sequence, c_term_mod, overwrite), n_term_mod, overwrite)


def condense_n_term_modifications(sequence: str) -> str:
    """
    Merges any N-terminal modification with the modification of the first amino acid.

    If the sequence has both an N-terminal modification and a modification on the first amino acid,
    this function will combine them into a single modification on the first amino acid.
    If there's only an N-terminal modification, it will be transferred to the first amino acid.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: Sequence with merged N-terminal and first amino acid modifications.
    :rtype: str

    .. code-block:: python

        # For sequences with just an N-terminal modification:
        >>> condense_n_term_modifications('[Acetyl]-PEPTIDE')
        'P[Acetyl]EPTIDE'

        # Combining N-terminal and first amino acid modifications:
        >>> condense_n_term_modifications('[1.0]-P[2]EPTIDE')
        'P[3.0]EPTIDE'
        >>> condense_n_term_modifications('[Acetyl]-P[Amide]EPTIDE')
        'P[AcetylAmide]EPTIDE'

        # TypeError when incompatible modifications are attempted to be combined:
        >>> condense_n_term_modifications('[Acetyl]-P[1.0]EPTIDE')
        Traceback (most recent call last):
        ...
        TypeError: can only concatenate str (not "float") to str

    """

    n_term_mod = get_n_term_modification(sequence)
    if n_term_mod is None:
        return sequence

    sequence = strip_n_term_modification(sequence)
    if sequence[1] == '[':  # if the first amino acid is modified
        end = sequence.index(']')
        condensed_mod = n_term_mod + convert_type(sequence[2:end])
        sequence = f"{sequence[:1]}[{condensed_mod}]{sequence[end + 1:]}"
    else:
        sequence = f"{sequence[0]}[{n_term_mod}]{sequence[1:]}"

    return sequence


def condense_c_term_modifications(sequence: str) -> str:
    """
    Merges any C-terminal modification with the modification of the last amino acid.

    If the sequence has both a C-terminal modification and a modification on the last amino acid,
    this function will combine them into a single modification on the last amino acid.
    If there's only a C-terminal modification, it will be transferred to the last amino acid.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: Sequence with merged C-terminal and last amino acid modifications.
    :rtype: str

    .. code-block:: python

        # For sequences with just a C-terminal modification:
        >>> condense_c_term_modifications('PEPTIDE-[Amide]')
        'PEPTIDE[Amide]'

        # Combining C-terminal and last amino acid modifications:
        >>> condense_c_term_modifications('PEPTIDE[1.0]-[2.0]')
        'PEPTIDE[3.0]'
        >>> condense_c_term_modifications('PEPTIDE[Acetyl]-[Amide]')
        'PEPTIDE[AcetylAmide]'

        # TypeError when incompatible modifications are attempted to be combined:
        >>> condense_c_term_modifications('PEPTIDE[1.0]-[Amide]')
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand type(s) for +: 'float' and 'str'

    """

    c_term_mod = get_c_term_modification(sequence)
    if c_term_mod is None:
        return sequence

    sequence = strip_c_term_modification(sequence)
    if sequence[-1] == ']':  # if the last amino acid is modified
        start = sequence.rindex('[')
        condensed_mod = convert_type(sequence[start + 1:-1]) + c_term_mod
        sequence = f"{sequence[:start]}[{condensed_mod}]"
    else:
        sequence = f"{sequence}[{c_term_mod}]"

    return sequence


def condense_term_modifications(sequence: str) -> str:
    """
    Merges both the N-terminal and C-terminal modifications with their respective adjacent amino acid modifications.

    If the sequence has terminal modifications (either N-terminal or C-terminal) and adjacent amino acid modifications,
    this function will combine them into single modifications on the respective amino acids.
    If there's only a terminal modification, it will be transferred to the adjacent amino acid.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: Sequence with merged terminal and adjacent amino acid modifications.
    :rtype: str

    .. code-block:: python

        # For sequences with modifications at both terminals:
        >>> condense_term_modifications('[1.0]-P[2.0]EPTIDE[1.0]-[2.0]')
        'P[3.0]EPTIDE[3.0]'

        # Condense_terms is equivalent to the sequential application of condense_n_term and condense_c_term:
        >>> condense_n_term_modifications(condense_c_term_modifications('[1.0]-P[2.0]EPTIDE[1.0]-[2.0]'))
        'P[3.0]EPTIDE[3.0]'

    """

    return condense_n_term_modifications(condense_c_term_modifications(sequence))


def pop_n_term_modification(sequence: str) -> Tuple[str, Union[str, None]]:
    """
    Removes the N-terminal modification from the sequence and returns it separately.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: The sequence without the N-terminal modification, and the removed N-terminal modification.
    :rtype: Tuple[str, Optional[str]]

    .. code-block:: python

        # For sequences with an N-terminal modification:
        >>> pop_n_term_modification('[Acetyl]-PEPTIDE')
        ('PEPTIDE', 'Acetyl')

        # For sequences without an N-terminal modification:
        >>> pop_n_term_modification('PEPTIDE')
        ('PEPTIDE', None)

    """

    n_term_mod = get_n_term_modification(sequence)
    if n_term_mod is None:
        return sequence, None
    sequence = strip_n_term_modification(sequence)
    return sequence, n_term_mod


def pop_c_term_modification(sequence: str) -> Tuple[str, Union[str, None]]:
    """
    Removes the C-terminal modification from the sequence and returns it separately.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: The sequence without the C-terminal modification, and the removed C-terminal modification.
    :rtype: Tuple[str, Optional[str]]

    .. code-block:: python

        # For sequences with a C-terminal modification:
        >>> pop_c_term_modification('PEPTIDE-[Amide]')
        ('PEPTIDE', 'Amide')

        # For sequences without a C-terminal modification:
        >>> pop_c_term_modification('PEPTIDE')
        ('PEPTIDE', None)

    """

    c_term_mod = get_c_term_modification(sequence)
    if c_term_mod is None:
        return sequence, None
    sequence = strip_c_term_modification(sequence)
    return sequence, c_term_mod


def pop_term_modifications(sequence: str) -> Tuple[str, Union[str, None], Union[str, None]]:
    """
    Removes both the N-terminal and C-terminal modifications from the sequence and returns them separately.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: The sequence without the terminal modifications, and the removed N-terminal and C-terminal modifications.
    :rtype: Tuple[str, Optional[str], Optional[str]]

    .. code-block:: python

        # For sequences with both N-terminal and C-terminal modifications:
        >>> pop_term_modifications('[Acetyl]-PEPTIDE-[Amide]')
        ('PEPTIDE', 'Acetyl', 'Amide')

        # For sequences with only an N-terminal modification:
        >>> pop_term_modifications('[Acetyl]-PEPTIDE')
        ('PEPTIDE', 'Acetyl', None)

        # For sequences with only a C-terminal modification:
        >>> pop_term_modifications('PEPTIDE-[Amide]')
        ('PEPTIDE', None, 'Amide')

        # For sequences without terminal modifications:
        >>> pop_term_modifications('PEPTIDE')
        ('PEPTIDE', None, None)

    """

    sequence, n_term_mod = pop_n_term_modification(sequence)
    sequence, c_term_mod = pop_c_term_modification(sequence)
    return sequence, n_term_mod, c_term_mod
