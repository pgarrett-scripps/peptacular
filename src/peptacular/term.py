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

import peptacular_bindings

from peptacular.util import convert_type


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

    mod = peptacular_bindings.extract_n_terminal_modification_py(sequence)

    if mod == '':
        return None

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

    return peptacular_bindings.remove_n_terminal_modification_py(sequence)


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

    mod = peptacular_bindings.extract_c_terminal_modification_py(sequence)

    if mod == '':
        return None

    return convert_type(mod)


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

    return peptacular_bindings.remove_c_terminal_modification_py(sequence)


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
