from typing import Tuple, List, Union
from peptacular.sequence.sequence import pop_modifications, add_modifications
from peptacular.types import ModValue


def pop_n_term_modifications(sequence: str) -> Tuple[str, List[ModValue]]:
    """
    Remove the n-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: A Tuple containing the sequence with the n-term modifications removed and the n-term modifications.
    :rtype: Tuple[str, List[ModValue]]

    .. code-block:: python

        >>> pop_n_term_modifications('[Acetyl]-PEPTIDE')
        ('PEPTIDE', ['Acetyl'])

        >>> pop_n_term_modifications('[Acetyl][1]-PEPTIDE')
        ('PEPTIDE', ['Acetyl', 1])

        >>> pop_n_term_modifications('PEPTIDE')
        ('PEPTIDE', [])

        >>> pop_n_term_modifications('{Oxidation}<C13>[Acetyl]-PEPTIDE')
        ('{Oxidation}<C13>PEPTIDE', ['Acetyl'])

    """
    stripped_sequence, mods = pop_modifications(sequence)
    n_term_mod = mods.pop('n', [])
    sequence = add_modifications(stripped_sequence, mods)
    return sequence, n_term_mod


def get_n_term_modifications(sequence: str) -> List[ModValue]:
    """
    Get the n-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The n-term modifications.
    :rtype: List[ModValue]

    .. code-block:: python

        >>> get_n_term_modifications('[Acetyl]-PEPTIDE')
        ['Acetyl']

        >>> get_n_term_modifications('[Acetyl][1]-PEPTIDE')
        ['Acetyl', 1]

        >>> get_n_term_modifications('PEPTIDE')
        []

        >>> get_n_term_modifications('{Oxidation}<C13>[Acetyl]-PEPTIDE')
        ['Acetyl']

    """
    return pop_n_term_modifications(sequence)[1]


def strip_n_term_modifications(sequence: str) -> str:
    """
    Remove the n-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with the n-term modifications removed.
    :rtype: str

    .. code-block:: python

        >>> strip_n_term_modifications('[Acetyl]-PEPTIDE')
        'PEPTIDE'

        >>> strip_n_term_modifications('[Acetyl][1]-PEPTIDE')
        'PEPTIDE'

        >>> strip_n_term_modifications('PEPTIDE')
        'PEPTIDE'

        >>> strip_n_term_modifications('{Oxidation}<C13>[Acetyl]-PEPTIDE')
        '{Oxidation}<C13>PEPTIDE'

    """
    return pop_n_term_modifications(sequence)[0]


def add_n_term_modifications(sequence: str, mods: Union[ModValue, List[ModValue]], overwrite: bool = False) -> str:
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

        >>> add_n_term_modifications('PEPTIDE', ['Acetyl'])
        '[Acetyl]-PEPTIDE'

        >>> add_n_term_modifications('PEPTIDE', 'Acetyl')
        '[Acetyl]-PEPTIDE'

        >>> add_n_term_modifications('PEPTIDE', ['Acetyl', 1])
        '[Acetyl][1]-PEPTIDE'

        >>> add_n_term_modifications('PEPTIDE', [])
        'PEPTIDE'

        >>> add_n_term_modifications('{Oxidation}<C13>PEPTIDE', ['Acetyl'])
        '{Oxidation}<C13>[Acetyl]-PEPTIDE'

        >>> add_n_term_modifications('[Acetyl]-PEPTIDE', 1)
        '[Acetyl][1]-PEPTIDE'

        >>> add_n_term_modifications('[Acetyl]-PEPTIDE', 1, overwrite=True)
        '[1]-PEPTIDE'

    """

    if not isinstance(mods, list):
        mods = [mods]

    return add_modifications(sequence, {'n': mods}, overwrite)


def pop_c_term_modifications(sequence: str) -> Tuple[str, List[ModValue]]:
    """
    Remove the c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with the c-term modifications and the c-term modifications.
    :rtype: Tuple[str, List[ModValue]]

    .. code-block:: python

        >>> pop_c_term_modifications('PEPTIDE-[Amide]')
        ('PEPTIDE', ['Amide'])

        >>> pop_c_term_modifications('PEPTIDE-[Amide][1]')
        ('PEPTIDE', ['Amide', 1])

        >>> pop_c_term_modifications('PEPTIDE')
        ('PEPTIDE', [])

        >>> pop_c_term_modifications('<15C>PEPTIDE[Oxidation]-[Amide]')
        ('<15C>PEPTIDE[Oxidation]', ['Amide'])

    """
    stripped_sequence, mods = pop_modifications(sequence)
    c_term_mod = mods.pop('c', [])
    sequence = add_modifications(stripped_sequence, mods)
    return sequence, c_term_mod


def get_c_term_modifications(sequence: str) -> List[ModValue]:
    """
    Get the c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The c-term modifications.
    :rtype: List[ModValue]

    .. code-block:: python

        >>> get_c_term_modifications('PEPTIDE-[Amide]')
        ['Amide']

        >>> get_c_term_modifications('PEPTIDE-[Amide][1]')
        ['Amide', 1]

        >>> get_c_term_modifications('PEPTIDE')
        []

        >>> get_c_term_modifications('<15C>PEPTIDE[Oxidation]-[Amide]')
        ['Amide']

    """

    return pop_c_term_modifications(sequence)[1]


def strip_c_term_modification(sequence: str) -> str:
    """
    Remove the c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with the c-term modifications removed.
    :rtype: str

    .. code-block:: python

        >>> strip_c_term_modification('PEPTIDE-[Amide]')
        'PEPTIDE'

        >>> strip_c_term_modification('PEPTIDE-[Amide][1]')
        'PEPTIDE'

        >>> strip_c_term_modification('PEPTIDE')
        'PEPTIDE'

        >>> strip_c_term_modification('<15C>PEPTIDE[Oxidation]-[Amide]')
        '<15C>PEPTIDE[Oxidation]'

    """
    return pop_c_term_modifications(sequence)[0]


def add_c_term_modifications(sequence: str, mods: Union[ModValue, List[ModValue]], overwrite: bool = False) -> str:
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

        >>> add_c_term_modifications('PEPTIDE', ['Amide'])
        'PEPTIDE-[Amide]'

        >>> add_c_term_modifications('PEPTIDE', 'Amide')
        'PEPTIDE-[Amide]'

        >>> add_c_term_modifications('PEPTIDE', ['Amide', 1])
        'PEPTIDE-[Amide][1]'

        >>> add_c_term_modifications('PEPTIDE[Oxidation]', ['Amide'])
        'PEPTIDE[Oxidation]-[Amide]'

        >>> add_c_term_modifications('PEPTIDE-[Amide]', [1])
        'PEPTIDE-[Amide][1]'

        >>> add_c_term_modifications('PEPTIDE-[Amide]', [1], overwrite=True)
        'PEPTIDE-[1]'

    """

    if not isinstance(mods, list):
        mods = [mods]

    return add_modifications(sequence, {'c': mods}, overwrite)


def pop_terminal_modifications(sequence: str) -> Tuple[str, List[ModValue], List[ModValue]]:
    """
    Remove the n-term and c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: Stripped sequence, n-term modifications, c-term modifications.
    :rtype: Tuple[str, List[ModValue], List[ModValue]]

    .. code-block:: python

        >>> pop_terminal_modifications('[Acetyl]-PEPTIDE-[Amide]')
        ('PEPTIDE', ['Acetyl'], ['Amide'])

        >>> pop_terminal_modifications('[Acetyl][1]-PEPTIDE-[Amide][1]')
        ('PEPTIDE', ['Acetyl', 1], ['Amide', 1])

        >>> pop_terminal_modifications('PEPTIDE')
        ('PEPTIDE', [], [])

        >>> pop_terminal_modifications('{Oxidation}<C13>[Acetyl]-PEPTIDE[Oxidation]-[Amide]')
        ('{Oxidation}<C13>PEPTIDE[Oxidation]', ['Acetyl'], ['Amide'])

    """
    stripped_sequence, mods = pop_modifications(sequence)
    c_term_mod = mods.pop('c', [])
    n_term_mod = mods.pop('n', [])
    sequence = add_modifications(stripped_sequence, mods)
    return sequence, n_term_mod, c_term_mod


def get_terminal_modifications(sequence: str) -> Tuple[List[ModValue], List[ModValue]]:
    """
    Get the n-term and c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The n-term modifications, c-term modifications.
    :rtype: Tuple[List[ModValue], List[ModValue]]

    .. code-block:: python

        >>> get_terminal_modifications('[Acetyl]-PEPTIDE-[Amide]')
        (['Acetyl'], ['Amide'])

        >>> get_terminal_modifications('[Acetyl][1]-PEPTIDE-[Amide][1]')
        (['Acetyl', 1], ['Amide', 1])

        >>> get_terminal_modifications('PEPTIDE')
        ([], [])

        >>> get_terminal_modifications('{Oxidation}<C13>[Acetyl]-PEPTIDE[Oxidation]-[Amide]')
        (['Acetyl'], ['Amide'])

    """
    return pop_terminal_modifications(sequence)[1:]


def strip_terminal_modifications(sequence: str) -> str:
    """
    Remove the n-term and c-term modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: Stripped sequence.
    :rtype: str

    .. code-block:: python

        >>> strip_terminal_modifications('[Acetyl]-PEPTIDE-[Amide]')
        'PEPTIDE'

        >>> strip_terminal_modifications('[Acetyl][1]-PEPTIDE-[Amide][1]')
        'PEPTIDE'

        >>> strip_terminal_modifications('PEPTIDE')
        'PEPTIDE'

        >>> strip_terminal_modifications('{Oxidation}<C13>[Acetyl]-PEPTIDE[Oxidation]-[Amide]')
        '{Oxidation}<C13>PEPTIDE[Oxidation]'

    """
    return pop_terminal_modifications(sequence)[0]


def add_terminal_modifications(sequence: str, n_term_mods: Union[ModValue, List[ModValue]],
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

        >>> add_terminal_modifications('PEPTIDE', ['Acetyl'], ['Amide'])
        '[Acetyl]-PEPTIDE-[Amide]'

        >>> add_terminal_modifications('PEPTIDE', 'Acetyl', 'Amide')
        '[Acetyl]-PEPTIDE-[Amide]'

        >>> add_terminal_modifications('PEPTIDE', ['Acetyl', 1], ['Amide', 1])
        '[Acetyl][1]-PEPTIDE-[Amide][1]'

        >>> add_terminal_modifications('PEPTIDE', ['Acetyl', 1], 'Amide')
        '[Acetyl][1]-PEPTIDE-[Amide]'

        >>> add_terminal_modifications('PEPTIDE[Oxidation]', ['Acetyl'], ['Amide'])
        '[Acetyl]-PEPTIDE[Oxidation]-[Amide]'

        >>> add_terminal_modifications('[Acetyl]-PEPTIDE-[Amide]', [1], [1])
        '[Acetyl][1]-PEPTIDE-[Amide][1]'

        >>> add_terminal_modifications('PEPTIDE-[Amide]', [1], [1], overwrite=True)
        '[1]-PEPTIDE-[1]'

    """

    if not isinstance(n_term_mods, list):
        n_term_mods = [n_term_mods]

    if not isinstance(c_term_mods, list):
        c_term_mods = [c_term_mods]

    return add_modifications(sequence, {'n': n_term_mods, 'c': c_term_mods}, overwrite)
