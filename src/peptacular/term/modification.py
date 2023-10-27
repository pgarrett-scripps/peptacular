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

from typing import Union, Any

from peptacular.util import convert_type

# TODO: Fix functions to work with nested brackets

def _get_n_term_modification_index(sequence: str) -> int:
    if sequence.startswith('['):
        nesting_level = 0
        for i, char in enumerate(sequence):
            if char == '[':
                nesting_level += 1
            elif char == ']':
                nesting_level -= 1
                if nesting_level == 0:
                    return i  # Found the end index of the N-terminal modification
        # If the loop ends and the nesting level is not zero, it means brackets are unbalanced.
        return -1
    else:
        # No N-terminal modification is present
        return -1


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
        >>> get_n_term_modification("[Acetyl]PEPTIDE")
        'Acetyl'

        # For float-based modifications:
        >>> get_n_term_modification("[3.1415]PEPTIDE")
        3.1415

        # For int-based modifications:
        >>> get_n_term_modification("[100]PEPTIDE")
        100

        # When no N-Terminus modification is present:
        >>> get_n_term_modification("PEPTIDE") # returns None

        # For string-based modifications (with internal brackets):
        >>> get_n_term_modification("[Acetyl:C(6)H[12]]PEPTIDE")
        'Acetyl:C(6)H[12]'

    """

    i = _get_n_term_modification_index(sequence)

    if i == -1:
        return None

    mod = sequence[1:i]  # Exclude the outer brackets
    return convert_type(mod)


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
        >>> strip_n_term_modification("[Acetyl]PEPTIDE")
        'PEPTIDE'
        >>> strip_n_term_modification("[Acetyl:C(6)H[12]]PEPTIDE")
        'PEPTIDE'
        >>> strip_n_term_modification("[3.1415]PEPTIDE")
        'PEPTIDE'
        >>> strip_n_term_modification("[100]PEPTIDE")
        'PEPTIDE'

        # If a residue modification is present at the N-terminus, only the terminal notation is removed:
        >>> strip_n_term_modification("[Acetyl]P(1)EPTIDE")
        'P(1)EPTIDE'

        # For sequences without any N-terminal modification:
        >>> strip_n_term_modification("PEPTIDE")
        'PEPTIDE'

    """

    i = _get_n_term_modification_index(sequence)

    if i == -1:
        return sequence

    return sequence[i + 1:]


def _get_c_term_modification_index(sequence: str) -> int:
    if sequence.endswith(']'):
        nesting_level = 0
        for i, char in enumerate(reversed(sequence)):
            if char == ']':
                nesting_level += 1
            elif char == '[':
                nesting_level -= 1
                if nesting_level == 0:
                    return len(sequence) - i  # Found the end index of the C-terminal modification
        # If the loop ends and the nesting level is not zero, it means brackets are unbalanced.
        return -1
    else:
        # No C-terminal modification is present
        return -1


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
        >>> get_c_term_modification("PEPTIDE[Amide]")
        'Amide'

        # For sequences with float-based modifications:
        >>> get_c_term_modification("PEPTIDE[3.1415]")
        3.1415

        # For sequences with int-based modifications:
        >>> get_c_term_modification("PEPTIDE[100]")
        100

        # When no C-terminal modification is present:
        >>> get_c_term_modification("PEPTIDE")

        # For string-based modifications (with internal brackets):
        >>> get_c_term_modification("PEPTIDE[Amide:C(6)H[12]]")
        'Amide:C(6)H[12]'

    """

    i = _get_c_term_modification_index(sequence)
    if i != -1:
        mod = sequence[i:-1]  # Exclude the outer brackets
        return convert_type(mod)
    return None


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
        >>> strip_c_term_modification("PEPTIDE[Acetyl]")
        'PEPTIDE'
        >>> strip_c_term_modification("PEPTIDE[Acetyl:C(6)H[12]]")
        'PEPTIDE'
        >>> strip_c_term_modification("PEPTIDE[3.1415]")
        'PEPTIDE'
        >>> strip_c_term_modification("PEPTIDE[100]")
        'PEPTIDE'

        # If a residue modification is present at the C-terminus, only the terminal notation is removed:
        >>> strip_c_term_modification("P(1)EPTIDE[Amide]")
        'P(1)EPTIDE'

        # For sequences without any C-terminal modification:
        >>> strip_c_term_modification("PEPTIDE")
        'PEPTIDE'

    """

    i = _get_c_term_modification_index(sequence)
    if i != -1:
        return sequence[:i - 1]  # Exclude the C-terminal modification notation
    return sequence


def add_n_term_modification(sequence: str, mod: Any) -> str:
    """
    Appends the specified N-terminal modification to the provided sequence.

    If the sequence already contains an N-terminal modification, the new modification will be combined with the
    existing one. Ensure that the types of the modifications are compatible to prevent errors.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param mod: Modification to be appended at the N-terminus of the sequence.
    :type mod: Any

    :return: Modified sequence with the added N-terminal notation.
    :rtype: str

    .. code-block:: python

        # Adding string-based modifications:
        >>> add_n_term_modification("PEPTIDE", "Acetyl")
        '[Acetyl]PEPTIDE'

        # Adding float-based modifications:
        >>> add_n_term_modification("PEPTIDE", 3.1415)
        '[3.1415]PEPTIDE'

        # Adding int-based modifications:
        >>> add_n_term_modification("PEPTIDE", 100)
        '[100]PEPTIDE'

        # Combining new modifications with existing N-terminal modifications:
        >>> add_n_term_modification("[Acetyl]PEPTIDE", "Amide")
        '[AcetylAmide]PEPTIDE'
        >>> add_n_term_modification("[3.1415]PEPTIDE", 3.1415)
        '[6.283]PEPTIDE'
        >>> add_n_term_modification("[100]PEPTIDE", 100)
        '[200]PEPTIDE'

        # In case of incompatible modification types:
        >>> add_n_term_modification("[100]PEPTIDE", '100')
        Traceback (most recent call last):
            ...
        TypeError: unsupported operand type(s) for +: 'int' and 'str'

        # No change if the modification is None:
        >>> add_n_term_modification("PEPTIDE", None)
        'PEPTIDE'

    """

    if mod is None:
        return sequence

    existing_mod = get_n_term_modification(sequence)
    if existing_mod is not None:
        new_mod = existing_mod + mod
        return f"[{str(new_mod)}]{strip_n_term_modification(sequence)}"

    return f"[{str(mod)}]{sequence}"


def add_c_term_modification(sequence: str, mod: Any) -> str:
    """
    Appends the specified C-terminal modification notation to the provided sequence.

    If the sequence already contains a C-terminal modification, the new modification will be
    combined with the existing one. Ensure that the types of the modifications are compatible
    to prevent errors.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param mod: Modification to be appended at the C-terminus of the sequence. Can be of type str, int, or float.
    :type mod: Any

    :return: Modified sequence with the added C-terminal notation.
    :rtype: str

    .. code-block:: python

        # Adding string-based modifications:
        >>> add_c_term_modification("PEPTIDE", "Amide")
        'PEPTIDE[Amide]'

        # Adding float-based modifications:
        >>> add_c_term_modification("PEPTIDE", 3.1415)
        'PEPTIDE[3.1415]'

        # Adding int-based modifications:
        >>> add_c_term_modification("PEPTIDE", 100)
        'PEPTIDE[100]'

        # Combining new modifications with existing C-terminal modifications:
        >>> add_c_term_modification("PEPTIDE[Amide]", "Acetyl")
        'PEPTIDE[AmideAcetyl]'
        >>> add_c_term_modification("PEPTIDE[3.1415]", 3.1415)
        'PEPTIDE[6.283]'
        >>> add_c_term_modification("PEPTIDE[100]", 100)
        'PEPTIDE[200]'

        # In case of incompatible modification types:
        >>> add_c_term_modification("PEPTIDE[100]", '100')
        Traceback (most recent call last):
            ...
        TypeError: unsupported operand type(s) for +: 'int' and 'str'

        # No change if the modification is None:
        >>> add_c_term_modification("PEPTIDE", None)
        'PEPTIDE'

    """

    if mod is None:
        return sequence

    existing_mod = get_c_term_modification(sequence)
    if existing_mod is not None:
        new_mod = existing_mod + mod
        return f"{strip_c_term_modification(sequence)}[{str(new_mod)}]"

    return f"{sequence}[{str(mod)}]"


def strip_term_modifications(sequence: str) -> str:
    """
    Removes both N-terminal and C-terminal modification notations from the sequence.

    If the sequence contains modifications at both terminals, this function will strip both
    notations. If there are no terminal modifications, the sequence is returned unchanged.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: Sequence without any terminal modification notations.
    :rtype: str

    :Examples:

    .. code-block:: python

        # Removing both N-terminal and C-terminal modifications:
        >>> strip_term_modifications("[Acetyl]PEPTIDE[Amide]")
        'PEPTIDE'

        # The result is equivalent to applying strip_n_term_modification followed by strip_c_term_modification:
        >>> strip_n_term_modification(strip_c_term_modification("[Acetyl]PEPTIDE[Amide]"))
        'PEPTIDE'

    """

    return strip_n_term_modification(strip_c_term_modification(sequence))


# TODO: Fix condense functions to handle nested modifications
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
        >>> condense_n_term_modifications('[Acetyl]PEPTIDE')
        'P(Acetyl)EPTIDE'

        # Combining N-terminal and first amino acid modifications:
        >>> condense_n_term_modifications('[1.0]P(2)EPTIDE')
        'P(3.0)EPTIDE'
        >>> condense_n_term_modifications('[Acetyl]P(Amide)EPTIDE')
        'P(AcetylAmide)EPTIDE'

        # TypeError when incompatible modifications are attempted to be combined:
        >>> condense_n_term_modifications('[Acetyl]P(1.0)EPTIDE')
        Traceback (most recent call last):
        ...
        TypeError: can only concatenate str (not "float") to str

    """

    n_term_mod = get_n_term_modification(sequence)
    if n_term_mod is None:
        return sequence

    sequence = strip_n_term_modification(sequence)
    if sequence[1] == '(':  # if the first amino acid is modified
        end = sequence.index(')')
        condensed_mod = n_term_mod + convert_type(sequence[2:end])
        sequence = f"{sequence[:1]}({condensed_mod}){sequence[end + 1:]}"
    else:
        sequence = f"{sequence[0]}({n_term_mod}){sequence[1:]}"

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
        >>> condense_c_term_modifications('PEPTIDE[Amide]')
        'PEPTIDE(Amide)'

        # Combining C-terminal and last amino acid modifications:
        >>> condense_c_term_modifications('PEPTIDE(1.0)[2.0]')
        'PEPTIDE(3.0)'
        >>> condense_c_term_modifications('PEPTIDE(Acetyl)[Amide]')
        'PEPTIDE(AcetylAmide)'

        # TypeError when incompatible modifications are attempted to be combined:
        >>> condense_c_term_modifications('PEPTIDE(1.0)[Amide]')
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand type(s) for +: 'float' and 'str'

    """

    c_term_mod = get_c_term_modification(sequence)
    if c_term_mod is None:
        return sequence

    sequence = strip_c_term_modification(sequence)
    if sequence[-1] == ')':  # if the last amino acid is modified
        start = sequence.rindex('(')
        condensed_mod = convert_type(sequence[start + 1:-1]) + c_term_mod
        sequence = f"{sequence[:start]}({condensed_mod})"
    else:
        sequence = f"{sequence}({c_term_mod})"

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
        >>> condense_term_modifications('[1.0]P(2.0)EPTIDE(1.0)[2.0]')
        'P(3.0)EPTIDE(3.0)'

        # Condense_terms is equivalent to the sequential application of condense_n_term and condense_c_term:
        >>> condense_n_term_modifications(condense_c_term_modifications('[1.0]P(2.0)EPTIDE(1.0)[2.0]'))
        'P(3.0)EPTIDE(3.0)'

    """

    return condense_n_term_modifications(condense_c_term_modifications(sequence))
