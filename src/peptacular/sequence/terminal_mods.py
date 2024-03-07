from typing import Tuple, List, Union
from peptacular.sequence.sequence import pop_mods, add_mods
from peptacular.types import ModValue


def pop_nterm_mods(sequence: str) -> Tuple[str, List[ModValue]]:
    """
    Remove the n-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: A Tuple containing the sequence with the n-term modifications removed and the n-term modifications.
    :rtype: Tuple[str, List[ModValue]]

    .. code-block:: python

        >>> pop_nterm_mods('[Acetyl]-PEPTIDE')
        ('PEPTIDE', ['Acetyl'])

        >>> pop_nterm_mods('[Acetyl][1]-PEPTIDE')
        ('PEPTIDE', ['Acetyl', 1])

        >>> pop_nterm_mods('PEPTIDE')
        ('PEPTIDE', [])

        >>> pop_nterm_mods('{Oxidation}<C13>[Acetyl]-PEPTIDE')
        ('{Oxidation}<C13>PEPTIDE', ['Acetyl'])

    """
    stripped_sequence, mods = pop_mods(sequence)
    n_term_mod = mods.pop('n', [])
    sequence = add_mods(stripped_sequence, mods)
    return sequence, n_term_mod


def get_nterm_mods(sequence: str) -> List[ModValue]:
    """
    Get the n-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The n-term modifications.
    :rtype: List[ModValue]

    .. code-block:: python

        >>> get_nterm_mods('[Acetyl]-PEPTIDE')
        ['Acetyl']

        >>> get_nterm_mods('[Acetyl][1]-PEPTIDE')
        ['Acetyl', 1]

        >>> get_nterm_mods('PEPTIDE')
        []

        >>> get_nterm_mods('{Oxidation}<C13>[Acetyl]-PEPTIDE')
        ['Acetyl']

    """
    return pop_nterm_mods(sequence)[1]


def strip_nterm_mods(sequence: str) -> str:
    """
    Remove the n-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with the n-term modifications removed.
    :rtype: str

    .. code-block:: python

        >>> strip_nterm_mods('[Acetyl]-PEPTIDE')
        'PEPTIDE'

        >>> strip_nterm_mods('[Acetyl][1]-PEPTIDE')
        'PEPTIDE'

        >>> strip_nterm_mods('PEPTIDE')
        'PEPTIDE'

        >>> strip_nterm_mods('{Oxidation}<C13>[Acetyl]-PEPTIDE')
        '{Oxidation}<C13>PEPTIDE'

    """
    return pop_nterm_mods(sequence)[0]


def add_nterm_mods(sequence: str, mods: Union[ModValue, List[ModValue]], overwrite: bool = False) -> str:
    """
    Add n-term modifications to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param mods: The n-term modifications to add.
    :type mods: Union[ModValue, List[ModValue]]
    :param overwrite: Whether to overwrite an existing n-term modifications.
    :type overwrite: bool

    :return: The sequence with the n-term modifications added.
    :rtype: str

    .. code-block:: python

        >>> add_nterm_mods('PEPTIDE', ['Acetyl'])
        '[Acetyl]-PEPTIDE'

        >>> add_nterm_mods('PEPTIDE', 'Acetyl')
        '[Acetyl]-PEPTIDE'

        >>> add_nterm_mods('PEPTIDE', ['Acetyl', 1])
        '[Acetyl][1]-PEPTIDE'

        >>> add_nterm_mods('PEPTIDE', [])
        'PEPTIDE'

        >>> add_nterm_mods('{Oxidation}<C13>PEPTIDE', ['Acetyl'])
        '{Oxidation}<C13>[Acetyl]-PEPTIDE'

        >>> add_nterm_mods('[Acetyl]-PEPTIDE', 1)
        '[Acetyl][1]-PEPTIDE'

        >>> add_nterm_mods('[Acetyl]-PEPTIDE', 1, overwrite=True)
        '[1]-PEPTIDE'

    """

    if not isinstance(mods, list):
        mods = [mods]

    return add_mods(sequence, {'n': mods}, overwrite)


def pop_cterm_mods(sequence: str) -> Tuple[str, List[ModValue]]:
    """
    Remove the c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with the c-term modifications and the c-term modifications.
    :rtype: Tuple[str, List[ModValue]]

    .. code-block:: python

        >>> pop_cterm_mods('PEPTIDE-[Amide]')
        ('PEPTIDE', ['Amide'])

        >>> pop_cterm_mods('PEPTIDE-[Amide][1]')
        ('PEPTIDE', ['Amide', 1])

        >>> pop_cterm_mods('PEPTIDE')
        ('PEPTIDE', [])

        >>> pop_cterm_mods('<15C>PEPTIDE[Oxidation]-[Amide]')
        ('<15C>PEPTIDE[Oxidation]', ['Amide'])

    """
    stripped_sequence, mods = pop_mods(sequence)
    c_term_mod = mods.pop('c', [])
    sequence = add_mods(stripped_sequence, mods)
    return sequence, c_term_mod


def get_cterm_mods(sequence: str) -> List[ModValue]:
    """
    Get the c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The c-term modifications.
    :rtype: List[ModValue]

    .. code-block:: python

        >>> get_cterm_mods('PEPTIDE-[Amide]')
        ['Amide']

        >>> get_cterm_mods('PEPTIDE-[Amide][1]')
        ['Amide', 1]

        >>> get_cterm_mods('PEPTIDE')
        []

        >>> get_cterm_mods('<15C>PEPTIDE[Oxidation]-[Amide]')
        ['Amide']

    """

    return pop_cterm_mods(sequence)[1]


def strip_cterm_mods(sequence: str) -> str:
    """
    Remove the c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with the c-term modifications removed.
    :rtype: str

    .. code-block:: python

        >>> strip_cterm_mods('PEPTIDE-[Amide]')
        'PEPTIDE'

        >>> strip_cterm_mods('PEPTIDE-[Amide][1]')
        'PEPTIDE'

        >>> strip_cterm_mods('PEPTIDE')
        'PEPTIDE'

        >>> strip_cterm_mods('<15C>PEPTIDE[Oxidation]-[Amide]')
        '<15C>PEPTIDE[Oxidation]'

    """
    return pop_cterm_mods(sequence)[0]


def add_cterm_mods(sequence: str, mods: Union[ModValue, List[ModValue]], overwrite: bool = False) -> str:
    """
    Add a c-term modification to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param mods: The c-term modification to add.
    :type mods: Union[ModValue, List[ModValue]]
    :param overwrite: If True, overwrite the existing c-term modification. If False, add the modification to the
     existing c-term modification.
    :type overwrite: bool

    :return: The sequence with the c-term modification added.
    :rtype: str

    .. code-block:: python

        >>> add_cterm_mods('PEPTIDE', ['Amide'])
        'PEPTIDE-[Amide]'

        >>> add_cterm_mods('PEPTIDE', 'Amide')
        'PEPTIDE-[Amide]'

        >>> add_cterm_mods('PEPTIDE', ['Amide', 1])
        'PEPTIDE-[Amide][1]'

        >>> add_cterm_mods('PEPTIDE[Oxidation]', ['Amide'])
        'PEPTIDE[Oxidation]-[Amide]'

        >>> add_cterm_mods('PEPTIDE-[Amide]', [1])
        'PEPTIDE-[Amide][1]'

        >>> add_cterm_mods('PEPTIDE-[Amide]', [1], overwrite=True)
        'PEPTIDE-[1]'

    """

    if not isinstance(mods, list):
        mods = [mods]

    return add_mods(sequence, {'c': mods}, overwrite)


def pop_terminal_mods(sequence: str) -> Tuple[str, List[ModValue], List[ModValue]]:
    """
    Remove the n-term and c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: Stripped sequence, n-term modifications, c-term modifications.
    :rtype: Tuple[str, List[ModValue], List[ModValue]]

    .. code-block:: python

        >>> pop_terminal_mods('[Acetyl]-PEPTIDE-[Amide]')
        ('PEPTIDE', ['Acetyl'], ['Amide'])

        >>> pop_terminal_mods('[Acetyl][1]-PEPTIDE-[Amide][1]')
        ('PEPTIDE', ['Acetyl', 1], ['Amide', 1])

        >>> pop_terminal_mods('PEPTIDE')
        ('PEPTIDE', [], [])

        >>> pop_terminal_mods('{Oxidation}<C13>[Acetyl]-PEPTIDE[Oxidation]-[Amide]')
        ('{Oxidation}<C13>PEPTIDE[Oxidation]', ['Acetyl'], ['Amide'])

    """
    stripped_sequence, mods = pop_mods(sequence)
    c_term_mod = mods.pop('c', [])
    n_term_mod = mods.pop('n', [])
    sequence = add_mods(stripped_sequence, mods)
    return sequence, n_term_mod, c_term_mod


def get_terminal_mods(sequence: str) -> Tuple[List[ModValue], List[ModValue]]:
    """
    Get the n-term and c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The n-term modifications, c-term modifications.
    :rtype: Tuple[List[ModValue], List[ModValue]]

    .. code-block:: python

        >>> get_terminal_mods('[Acetyl]-PEPTIDE-[Amide]')
        (['Acetyl'], ['Amide'])

        >>> get_terminal_mods('[Acetyl][1]-PEPTIDE-[Amide][1]')
        (['Acetyl', 1], ['Amide', 1])

        >>> get_terminal_mods('PEPTIDE')
        ([], [])

        >>> get_terminal_mods('{Oxidation}<C13>[Acetyl]-PEPTIDE[Oxidation]-[Amide]')
        (['Acetyl'], ['Amide'])

    """
    return pop_terminal_mods(sequence)[1:]


def strip_terminal_mods(sequence: str) -> str:
    """
    Remove the n-term and c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: Stripped sequence.
    :rtype: str

    .. code-block:: python

        >>> strip_terminal_mods('[Acetyl]-PEPTIDE-[Amide]')
        'PEPTIDE'

        >>> strip_terminal_mods('[Acetyl][1]-PEPTIDE-[Amide][1]')
        'PEPTIDE'

        >>> strip_terminal_mods('PEPTIDE')
        'PEPTIDE'

        >>> strip_terminal_mods('{Oxidation}<C13>[Acetyl]-PEPTIDE[Oxidation]-[Amide]')
        '{Oxidation}<C13>PEPTIDE[Oxidation]'

    """
    return pop_terminal_mods(sequence)[0]


def add_terminal_mods(sequence: str, n_term_mods: Union[ModValue, List[ModValue]],
                      c_term_mods: Union[ModValue, List[ModValue]], overwrite: bool = False) -> str:
    """
    Add n-term and c-term modifications to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param n_term_mods: The n-term modification to add.
    :type n_term_mods: str
    :param c_term_mods: The c-term modification to add.
    :type c_term_mods: str
    :param overwrite: Whether to overwrite existing modifications.
    :type overwrite: bool

    :return: The sequence with the n-term and c-term modifications added.
    :rtype: str

    .. code-block:: python

        >>> add_terminal_mods('PEPTIDE', ['Acetyl'], ['Amide'])
        '[Acetyl]-PEPTIDE-[Amide]'

        >>> add_terminal_mods('PEPTIDE', 'Acetyl', 'Amide')
        '[Acetyl]-PEPTIDE-[Amide]'

        >>> add_terminal_mods('PEPTIDE', ['Acetyl', 1], ['Amide', 1])
        '[Acetyl][1]-PEPTIDE-[Amide][1]'

        >>> add_terminal_mods('PEPTIDE', ['Acetyl', 1], 'Amide')
        '[Acetyl][1]-PEPTIDE-[Amide]'

        >>> add_terminal_mods('PEPTIDE[Oxidation]', ['Acetyl'], ['Amide'])
        '[Acetyl]-PEPTIDE[Oxidation]-[Amide]'

        >>> add_terminal_mods('[Acetyl]-PEPTIDE-[Amide]', [1], [1])
        '[Acetyl][1]-PEPTIDE-[Amide][1]'

        >>> add_terminal_mods('PEPTIDE-[Amide]', [1], [1], overwrite=True)
        '[1]-PEPTIDE-[1]'

    """

    if not isinstance(n_term_mods, list):
        n_term_mods = [n_term_mods]

    if not isinstance(c_term_mods, list):
        c_term_mods = [c_term_mods]

    return add_mods(sequence, {'n': n_term_mods, 'c': c_term_mods}, overwrite)
