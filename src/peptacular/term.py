from typing import Union, Any


def get_n_term_modification(sequence: str) -> Union[str, None]:
    """
    Parsed and returns the N-terminal modification if present, otherwise None.

    :param sequence: The peptide sequence with potential N-terminal notation.
    :type sequence: str

    :return: The N-terminal modification if present, otherwise None.
    :rtype: Union[str, None]

    :Example:

    .. code-block:: python

        >>> get_n_term_modification("[Acetyl]PEPTIDE[Amide]")
        'Acetyl'

        >>> get_n_term_modification("[Acetyl]P(1)EPTIDE(1)[Amide]")
        'Acetyl'

        >>> get_n_term_modification("PEPTIDE[Amide]") # None
    """

    if sequence.startswith('['):
        # find end notation
        end = sequence.find(']')
        return sequence[1:end]
    return None


def strip_n_term_modification(sequence: str) -> str:
    """
    Strip the N-terminal modification from the given sequence.

    :param sequence: The peptide sequence with potential N-terminal notation.
    :type sequence: str

    :return: The sequence without the N-terminal modification notation.
    :rtype: str

    :Example:

    .. code-block:: python

        >>> strip_n_term_modification("[Acetyl]PEPTIDE[Amide]")
        'PEPTIDE[Amide]'

        >>> strip_n_term_modification("[Acetyl]P(1)EPTIDE(1)[Amide]")
        'P(1)EPTIDE(1)[Amide]'

        >>> strip_n_term_modification("PEPTIDE[Amide]")
        'PEPTIDE[Amide]'
    """

    if sequence.startswith('['):
        # find end notation
        end = sequence.find(']')
        return sequence[end + 1:]
    return sequence


def add_n_term_modification(sequence: str, mod: Any) -> str:
    """
    Add the given N-terminal modification to sequence.

    :param sequence: The peptide sequence to add the N-terminal notation to.
    :type sequence: str

    :param mod: The N-terminal notation to add to the sequence.
    :type mod: str

    :return: The sequence with the N-terminal modification notation added.
    :rtype: str

    :Example:

    .. code-block:: python

        >>> add_n_term_modification("PEPTIDE[Amide]", "Acetyl")
        '[Acetyl]PEPTIDE[Amide]'

        >>> add_n_term_modification("P(1)EPTIDE(1)[Amide]", "Acetyl")
        '[Acetyl]P(1)EPTIDE(1)[Amide]'

        >>> add_n_term_modification("PEPTIDE[Amide]", None)
        'PEPTIDE[Amide]'
    """

    if mod is None:
        return sequence

    return f"[{str(mod)}]{sequence}"


def get_c_term_modification(sequence: str) -> Union[str, None]:
    """
    Parses and returns the C-terminal modification from the sequence.

    :param sequence: The peptide sequence with potential C-terminal notation.
    :type sequence: str

    :return: The C-terminal modification if present, otherwise None.
    :rtype: Union[str, None]

    :Example:

    .. code-block:: python

        >>> get_c_term_modification("[Acetyl]PEPTIDE[Amide]")
        'Amide'

        >>> get_c_term_modification("[Acetyl]P(1)EPTIDE(1)[Amide]")
        'Amide'

        >>> get_c_term_modification("[Acetyl]PEPTIDE") # None

    """

    if sequence.endswith(']'):
        # find start notation from back
        start = sequence.rfind('[')
        return sequence[start + 1:-1]

    return None


def strip_c_term_modification(sequence: str) -> str:
    """
    Strip the C-terminal modification from the given sequence.

    :param sequence: The peptide sequence with potential C-terminal notation.
    :type sequence: str

    :return: The sequence without the C-terminal modification notation.
    :rtype: str

    :Example:

    .. code-block:: python

        >>> strip_c_term_modification("[Acetyl]PEPTIDE[Amide]")
        '[Acetyl]PEPTIDE'

        >>> strip_c_term_modification("[Acetyl]P(1)EPTIDE(1)[Amide]")
        '[Acetyl]PEPTIDE'

        >>> strip_c_term_modification("[Acetyl]PEPTIDE")
        '[Acetyl]PEPTIDE'

    """

    if sequence.endswith(']'):
        # find start notation
        start = sequence.rfind('[')
        return sequence[:start]
    return sequence


def add_c_term_modification(sequence: str, mod: Any) -> str:
    """
    Add the given C-terminal modification to sequence.

    :param sequence: The peptide sequence to add the C-terminal notation to.
    :type sequence: str

    :param mod: The C-terminal notation to add to the sequence.
    :type mod: str

    :return: The sequence with the N-terminal modification notation added.
    :rtype: str

    :Example:

    .. code-block:: python

        >>> add_c_term_modification("[Acetyl]PEPTIDE", "Amide")
        '[Acetyl]PEPTIDE[Amide]'

        >>> add_c_term_modification("[Acetyl]P(1)EPTIDE(1)", "Amide")
        '[Acetyl]P(1)EPTIDE(1)[Amide]'

        >>> add_c_term_modification("[Acetyl]PEPTIDE", None)
        '[Acetyl]PEPTIDE'
    """

    if mod is None:
        return sequence

    return f"{sequence}[{str(mod)}]"


def strip_term_modifications(sequence: str) -> str:
    """
    Strip both N-terminal and C-terminal modifications from the sequence.

    :param sequence: The peptide sequence with potential terminal notations.
    :type sequence: str

    :return: The sequence without the terminal modification notations.
    :rtype: str

    :Example:

    .. code-block:: python

        >>> strip_term_modifications("[Acetyl]PEPTIDE[Amide]")
        'PEPTIDE'

        >>> strip_term_modifications("[Acetyl]P(1)EPTIDE(1)[Amide]")
        'P(1)EPTIDE(1)'

        >>> strip_term_modifications("[Acetyl]PEPTIDE")
        'PEPTIDE'

        >>> strip_term_modifications("PEPTIDE[Amide]")
        'PEPTIDE'

        >>> strip_term_modifications("PEPTIDE")
        'PEPTIDE'
    """

    return strip_n_term_modification(strip_c_term_modification(sequence))


def strip_n_term(sequence: str) -> str:
    """
    Removes the first (N-terminal) amino acid from the sequence, while preserving the position of any modifications.

    :param sequence: The sequence to trim.
    :type sequence: str

    :return: The trimmed sequence
    :rtype: str

    :Example:

    .. code-block:: python

        >>> strip_n_term('PEPTIDE')
        'EPTIDE'

        >>> strip_n_term('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]')
        'EP(phospho)TIDE[Amide]'

        >>> strip_n_term('E(1)[Amide]')
        ''

        >>> strip_n_term('E[Amide]')
        ''

    """
    sequence = strip_n_term_modification(sequence)
    sequence = sequence[1:]
    if sequence and sequence[0] == '(':
        sequence = sequence[sequence.index(')') + 1:]
    if sequence and sequence[0] == '[':
        sequence = sequence[sequence.index(']') + 1:]
    return sequence


def strip_c_term(sequence: str) -> str:
    """
    Removes the last (C-terminal) amino acid from the sequence, while preserving the position of any modifications.

    :param sequence: The sequence to trim.
    :type sequence: str

    :return: The trimmed sequence
    :rtype: str

    :Example:

    .. code-block:: python

        >>> strip_c_term('PEPTIDE')
        'PEPTID'

        >>> strip_c_term('[Acetyl]P(phospho)EP(phospho)TIDE(1)[Amide]')
        '[Acetyl]P(phospho)EP(phospho)TID'

        >>> strip_c_term('[Acetyl]P(phospho)')
        ''

        >>> strip_c_term('[Acetyl]P')
        ''

    """
    sequence = strip_c_term_modification(sequence)

    if sequence and sequence[-1] == ')':
        sequence = sequence[:sequence.rindex('(')]

    sequence = sequence[:-1]

    if sequence and sequence[-1] == ']':
        sequence = sequence[:sequence.rindex('[')]

    return sequence


def condense_n_term(sequence: str) -> str:
    """
    Combines the N-terminal modification with the first amino acid modicication. If the first amino acid is modified,
    the modification is added to the N-terminal modification. Both modifications must be float-like in order to be
    combined.

    :param sequence: The sequence to condense.
    :type sequence: str
    :return: The condensed sequence.
    :rtype: str

    :Example:

    .. code-block:: python

        >>> condense_n_term('[Acetyl]PEPTIDE')
        'P(Acetyl)EPTIDE'

        >>> condense_n_term('[Acetyl]P(1.0)EPTIDE')
        Traceback (most recent call last):
        ValueError: could not convert string to float: 'Acetyl'

        >>> condense_n_term('[1.0]P(2.0)EP(phospho)TIDE')
        'P(3.0)EP(phospho)TIDE'

    """
    n_term_mod = get_n_term_modification(sequence)
    if n_term_mod is None:
        return sequence

    sequence = strip_n_term_modification(sequence)
    if sequence[1] == '(':  # if the first amino acid is modified
        end = sequence.index(')')
        condensed_mod = float(sequence[2:end]) + float(n_term_mod)
        sequence = f"{sequence[:1]}({condensed_mod}){sequence[end + 1:]}"
    else:
        sequence = f"{sequence[0]}({n_term_mod}){sequence[1:]}"

    return sequence


def condense_c_term(sequence: str) -> str:
    """
    Combines the C-terminal modification with the last amino acid. If the last amino acid is modified,
    the modification is added to the C-terminal modification. Both modifications must be float-like in order to be
    combined.

    :param sequence: The sequence to condense.
    :type sequence: str
    :return: The condensed sequence.
    :rtype: str

    :Example:

    .. code-block:: python

        >>> condense_c_term('PEPTIDE[Amide]')
        'PEPTIDE(Amide)'

        >>> condense_c_term('PEPTIDE(1)[Amide]')
        Traceback (most recent call last):
        ValueError: could not convert string to float: 'Amide'

        >>> condense_c_term('PEPTIDE(1)[1.0]')
        'PEPTIDE(2.0)'

    """

    c_term_mod = get_c_term_modification(sequence)
    if c_term_mod is None:
        return sequence

    sequence = strip_c_term_modification(sequence)
    if sequence[-1] == ')':  # if the last amino acid is modified
        start = sequence.rindex('(')
        condensed_mod = float(sequence[start + 1:-1]) + float(c_term_mod)
        sequence = f"{sequence[:start]}({condensed_mod})"
    else:
        sequence = f"{sequence}({c_term_mod})"

    return sequence


def condense_terms(sequence: str) -> str:
    """
    Condenses the term modifications onto the associated term residues, adding modifications together if necessary.
    For term modifications to be combined with the term amino acid modification, both modifications must be float-like.


    :param sequence: modified sequence
    :type sequence: str
    :return: condensed sequence
    :rtype: str

    :Example:

    .. code-block:: python

        >>> condense_terms('[Acetyl]PEPTIDE[Amide]')
        'P(Acetyl)EPTIDE(Amide)'

        >>> condense_terms('[Acetyl]PEPTIDE(1)[Amide]')
        Traceback (most recent call last):
        ValueError: could not convert string to float: 'Amide'

        >>> condense_terms('[Acetyl]PEPTIDE(1)[1.0]')
        'P(Acetyl)EPTIDE(2.0)'
    """

    return condense_n_term(condense_c_term(sequence))
