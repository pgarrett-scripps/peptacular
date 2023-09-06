from typing import Tuple


def strip_n_term_residue(sequence: str) -> str:
    """
    Removes the N-terminal amino acid from the sequence, retaining any other modifications in their original positions.

    If the sequence begins with an N-terminal modification, it will be removed alongside the first amino acid.
    If no modifications are present, this operation is equivalent to eliminating the first character of the sequence.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: Sequence with the N-terminal amino acid (and any N-terminal modification) removed.
    :rtype: str

    .. code-block:: python

        # For sequences without modifications, the first character is removed:
        >>> strip_n_term_residue('PEPTIDE')
        'EPTIDE'

        # Sequences with N-terminal modifications have them removed along with the first amino acid:
        >>> strip_n_term_residue('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]')
        'EP(phospho)TIDE[Amide]'

        # For sequences with a single amino acid having a C-terminal modification, the result is an empty string:
        >>> strip_n_term_residue('E(1)[Amide]')
        ''
        >>> strip_n_term_residue('[Acetyl]P(phospho)[Amide]')
        ''
        >>> strip_n_term_residue('E[Amide]')
        ''

    """

    return sequence[_get_n_term_index(sequence):]


def strip_c_term_residue(sequence: str) -> str:
    """
    Removes the C-terminal amino acid from the sequence, retaining any other modifications in their original positions.

    If the sequence ends with a C-terminal modification, it will be removed alongside the last amino acid.
    If no modifications are present, this operation is equivalent to eliminating the last character of the sequence.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: Sequence with the C-terminal amino acid (and any C-terminal modification) removed.
    :rtype: str

    .. code-block:: python

        # For sequences without modifications, the last character is removed:
        >>> strip_c_term_residue('PEPTIDE')
        'PEPTID'

        # Sequences with C-terminal modifications have them removed along with the last amino acid:
        >>> strip_c_term_residue('[Acetyl]P(phospho)EP(phospho)TIDE(1)[Amide]')
        '[Acetyl]P(phospho)EP(phospho)TID'

        # For sequences with a single amino acid having a N-terminal modification, the result is an empty string:
        >>> strip_c_term_residue('[Acetyl]P(phospho)')
        ''
        >>> strip_c_term_residue('[Acetyl]P(phospho)[Amide]')
        ''
        >>> strip_c_term_residue('[Acetyl]P')
        ''

    """

    return sequence[:_get_c_term_index(sequence)]


def get_n_term_residue(sequence: str) -> str:
    """
    Extracts the N-terminal amino acid from the sequence, along with any  modifications.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: The N-terminal amino acid along with any modifications.
    :rtype: str

    .. code-block:: python

        # For sequences without modifications, the first character is returned:
        >>> get_n_term_residue('PEPTIDE')
        'P'

        # Sequences with N-terminal modifications return them along with the last amino acid:
        >>> get_n_term_residue('[Acetyl]P(phospho)EP(phospho)TIDE(phospho)[Amide]')
        '[Acetyl]P(phospho)'

        # For sequences with a single amino acid
        >>> get_n_term_residue('[Acetyl]P(phospho)[Amide]')
        '[Acetyl]P(phospho)[Amide]'
        >>> get_n_term_residue('[Acetyl]P')
        '[Acetyl]P'
        >>> get_n_term_residue('P[Amide]')
        'P[Amide]'
        >>> get_n_term_residue('P(phospho)[Amide]')
        'P(phospho)[Amide]'
        >>> get_n_term_residue('P(phospho)')
        'P(phospho)'

    """

    end_index = _get_n_term_index(sequence)
    return sequence[:end_index]


def get_c_term_residue(sequence: str) -> str:
    """
    Extracts the C-terminal amino acid from the sequence, along with any modifications.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: The C-terminal amino acid along with any modifications.
    :rtype: str

    .. code-block:: python

        # For sequences without modifications, the last character is returned:
        >>> get_c_term_residue('PEPTIDE')
        'E'

        # Sequences with C-terminal modifications return them along with the last amino acid:
        >>> get_c_term_residue('[Acetyl]P(phospho)EP(phospho)TIDE(phospho)[Amide]')
        'E(phospho)[Amide]'

        # For sequences with a single amino acid
        >>> get_c_term_residue('[Acetyl]P(phospho)[Amide]')
        '[Acetyl]P(phospho)[Amide]'
        >>> get_c_term_residue('[Acetyl]P')
        '[Acetyl]P'
        >>> get_c_term_residue('P[Amide]')
        'P[Amide]'
        >>> get_c_term_residue('P(phospho)[Amide]')
        'P(phospho)[Amide]'
        >>> get_c_term_residue('P(phospho)')
        'P(phospho)'

    """

    start_index = _get_c_term_index(sequence)
    return sequence[start_index:]


def pop_c_term_residue(sequence: str) -> Tuple[str, str]:
    """
    Removes the C-terminal amino acid from the sequence, along with any modifications, and returns it as a tuple with
    the remaining sequence.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: A tuple containing the C-terminal amino acid along with any modifications, and the remaining sequence.
    :rtype: tuple[str, str]

    .. code-block:: python

        # For sequences without modifications, the last character is returned:
        >>> pop_c_term_residue('PEPTIDE')
        ('E', 'PEPTID')

        # Sequences with C-terminal modifications return them along with the last amino acid:
        >>> pop_c_term_residue('[Acetyl]P(phospho)EP(phospho)TIDE(phospho)[Amide]')
        ('E(phospho)[Amide]', '[Acetyl]P(phospho)EP(phospho)TID')

    """

    c_term_index = _get_c_term_index(sequence)
    return sequence[c_term_index:], sequence[:c_term_index]


def pop_n_term_residue(sequence: str) -> Tuple[str, str]:
    """
    Removes the N-terminal amino acid from the sequence, along with any modifications, and returns it as a tuple with
    the remaining sequence.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: A tuple containing the N-terminal amino acid along with any modifications, and the remaining sequence.
    :rtype: tuple[str, str]

    .. code-block:: python

        # For sequences without modifications, the first character is returned:
        >>> pop_n_term_residue('PEPTIDE')
        ('P', 'EPTIDE')

        # Sequences with N-terminal modifications return them along with the last amino acid:
        >>> pop_n_term_residue('[Acetyl]P(phospho)EP(phospho)TIDE(phospho)[Amide]')
        ('[Acetyl]P(phospho)', 'EP(phospho)TIDE(phospho)[Amide]')

    """

    n_term_index = _get_n_term_index(sequence)
    return sequence[:n_term_index], sequence[n_term_index:]


def _get_c_term_index(sequence: str) -> int:
    """
    Returns the start index of the C-terminal amino acid in the sequence.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: The start index of the C-terminal amino acid in the sequence.
    :rtype: int

    .. code-block:: python

        >>> _get_c_term_index('PEPTIDE')
        6
        >>> 'PEPTIDE'[:6]
        'PEPTID'

        >>> _get_c_term_index('[Acetyl]P(phospho)EP(phospho)TIDE(phospho)[Amide]')
        32
        >>> '[Acetyl]P(phospho)EP(phospho)TIDE(phospho)[Amide]'[:32]
        '[Acetyl]P(phospho)EP(phospho)TID'

    """

    start_index = len(sequence) - 1
    if sequence and sequence[start_index] == ']':
        start_index = sequence.rindex('[') - 1

    if sequence[start_index] == ')':
        start_index = sequence.rindex('(') - 1

    # Check for N-Term modification
    if start_index != 0 and sequence[start_index - 1] == ']':
        start_index = 0

    return start_index


def _get_n_term_index(sequence: str) -> int:
    """
    Returns the end index of the N-terminal amino acid in the sequence.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str

    :return: The end index of the N-terminal amino acid in the sequence.
    :rtype: int

    .. code-block:: python

        >>> _get_n_term_index('PEPTIDE')
        1
        >>> 'PEPTIDE'[1:]
        'EPTIDE'

        >>> _get_n_term_index('[Acetyl]P(phospho)EP(phospho)TIDE(phospho)[Amide]')
        18
        >>> '[Acetyl]P(phospho)EP(phospho)TIDE(phospho)[Amide]'[18:]
        'EP(phospho)TIDE(phospho)[Amide]'

    """

    end_index = 0
    if sequence and sequence[end_index] == '[':
        end_index = sequence.index(']') + 1

    end_index += 1

    if end_index != len(sequence) and sequence[end_index] == '(':
        end_index = sequence.index(')') + 1

    # Check for N-Term modification
    if end_index != len(sequence) and sequence[end_index] == '[':
        end_index = len(sequence)

    return end_index
