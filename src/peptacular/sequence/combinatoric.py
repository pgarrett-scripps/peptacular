import itertools
from typing import List, Union

from peptacular.sequence import sequence_length
from peptacular.sequence.sequence import _construct_sequence, split, add_mods, pop_mods


def permutate(sequence: str, size: int = None) -> List[str]:
    """
    Generates all permutations of the input sequence. Terminal sequence are kept in place.

    :param sequence: The sequence to be permuted.
    :type sequence: str
    :param size: The size of the permutations.
    :type size: int

    :return: A list of all permutations of the input sequence.
    :rtype: List[str]

    .. code-block:: python

        >>> permutate('PET')
        ['PET', 'PTE', 'EPT', 'ETP', 'TPE', 'TEP']

        >>> permutate('[3]-PET-[1]')
        ['[3]-PET-[1]', '[3]-PTE-[1]', '[3]-EPT-[1]', '[3]-ETP-[1]', '[3]-TPE-[1]', '[3]-TEP-[1]']

        >>> permutate('PE[3.14]T')
        ['PE[3.14]T', 'PTE[3.14]', 'E[3.14]PT', 'E[3.14]TP', 'TPE[3.14]', 'TE[3.14]P']

        >>> permutate('<13C>PET')
        ['<13C>PET', '<13C>PTE', '<13C>EPT', '<13C>ETP', '<13C>TPE', '<13C>TEP']

    """

    if size is None:
        size = sequence_length(sequence)

    sequence, mods = pop_mods(sequence)

    labile_mods = mods.pop('l', [])
    static_mods = mods.pop('s', [])
    isotope_mods = mods.pop('i', [])
    n_term_mods = mods.pop('n', [])
    c_term_mods = mods.pop('c', [])
    unknown_mods = mods.pop('u', [])

    internal_modified_sequence = _construct_sequence(sequence, mods)

    new_mods = {'l': labile_mods, 's': static_mods, 'i': isotope_mods, 'n': n_term_mods, 'c': c_term_mods,
                'u': unknown_mods}

    components = split(internal_modified_sequence)
    return [add_mods(''.join(p), new_mods) for p in itertools.permutations(components, size)]


def product(sequence: str, size: Union[int, None]) -> List[str]:
    """
    Generates all sartesian products of the input sequence of a given size. Terminal sequence are kept in place.

    :param sequence: The sequence to be combined.
    :type sequence: str
    :param size: The size of the combinations to be generated.
    :type size: int

    :return: A list of all combinations of the input sequence of the given size.
    :rtype: List[str]

    .. code-block:: python

        >>> product('PET', 2)
        ['PP', 'PE', 'PT', 'EP', 'EE', 'ET', 'TP', 'TE', 'TT']

        >>> product('[3]-PET-[1]', 2)[:5]
        ['[3]-PP-[1]', '[3]-PE-[1]', '[3]-PT-[1]', '[3]-EP-[1]', '[3]-EE-[1]']

        >>> product('PE[3.14]T', 2)
        ['PP', 'PE[3.14]', 'PT', 'E[3.14]P', 'E[3.14]E[3.14]', 'E[3.14]T', 'TP', 'TE[3.14]', 'TT']

        >>> product('<13C>PET', 2)
        ['<13C>PP', '<13C>PE', '<13C>PT', '<13C>EP', '<13C>EE', '<13C>ET', '<13C>TP', '<13C>TE', '<13C>TT']

    """

    sequence, mods = pop_mods(sequence)

    labile_mods = mods.pop('l', [])
    static_mods = mods.pop('s', [])
    isotope_mods = mods.pop('i', [])
    n_term_mods = mods.pop('n', [])
    c_term_mods = mods.pop('c', [])
    unknown_mods = mods.pop('u', [])

    internal_modified_sequence = _construct_sequence(sequence, mods)

    new_mods = {'l': labile_mods, 's': static_mods, 'i': isotope_mods, 'n': n_term_mods, 'c': c_term_mods,
                'u': unknown_mods}

    components = split(internal_modified_sequence)
    return [add_mods(''.join(p), new_mods) for p in itertools.product(components, repeat=size)]


def combinations(sequence: str, size: Union[int, None]) -> List[str]:
    """
    Generates all combinations of the input sequence of a given size. Terminal sequence are kept in place.

    :param sequence: The sequence to be combined.
    :type sequence: str

    :param size: The size of the combinations to be generated.
    :type size: int

    :return: A list of all combinations of the input sequence of the given size.
    :rtype: List[str]

    .. code-block:: python

        >>> combinations('PET', 2)
        ['PE', 'PT', 'ET']

        >>> combinations('[3]-PET-[1]', 2)
        ['[3]-PE-[1]', '[3]-PT-[1]', '[3]-ET-[1]']

        >>> combinations('PE[3.14]T', 2)
        ['PE[3.14]', 'PT', 'E[3.14]T']

        >>> combinations('<13C>PET', 2)
        ['<13C>PE', '<13C>PT', '<13C>ET']

    """

    if size is None:
        size = sequence_length(sequence)

    sequence, mods = pop_mods(sequence)

    labile_mods = mods.pop('l', [])
    static_mods = mods.pop('s', [])
    isotope_mods = mods.pop('i', [])
    n_term_mods = mods.pop('n', [])
    c_term_mods = mods.pop('c', [])
    unknown_mods = mods.pop('u', [])

    internal_modified_sequence = _construct_sequence(sequence, mods)

    new_mods = {'l': labile_mods, 's': static_mods, 'i': isotope_mods, 'n': n_term_mods, 'c': c_term_mods,
                'u': unknown_mods}
    components = split(internal_modified_sequence)
    return [add_mods(''.join(c), new_mods) for c in itertools.combinations(components, size)]


def combinations_with_replacement(sequence: str, size: Union[int, None]) -> List[str]:
    """
    Generates all combinations with replacement of the input sequence of a given size. Terminal sequence are kept
    in place.

    :param sequence: The sequence to be combined.
    :type sequence: str

    :param size: The size of the combinations to be generated.
    :type size: int

    :return: A list of all combinations of the input sequence of the given size.
    :rtype: List[str]

    .. code-block:: python

        >>> combinations_with_replacement('PET', 2)
        ['PP', 'PE', 'PT', 'EE', 'ET', 'TT']

        >>> combinations_with_replacement('[3]-PET-[1]', 2)
        ['[3]-PP-[1]', '[3]-PE-[1]', '[3]-PT-[1]', '[3]-EE-[1]', '[3]-ET-[1]', '[3]-TT-[1]']

        >>> combinations_with_replacement('PE[3.14]T', 2)
        ['PP', 'PE[3.14]', 'PT', 'E[3.14]E[3.14]', 'E[3.14]T', 'TT']

        >>> combinations_with_replacement('<13C>PET', 2)
        ['<13C>PP', '<13C>PE', '<13C>PT', '<13C>EE', '<13C>ET', '<13C>TT']

    """

    if size is None:
        size = sequence_length(sequence)

    sequence, mods = pop_mods(sequence)

    labile_mods = mods.pop('l', [])
    static_mods = mods.pop('s', [])
    isotope_mods = mods.pop('i', [])
    n_term_mods = mods.pop('n', [])
    c_term_mods = mods.pop('c', [])
    unknown_mods = mods.pop('u', [])

    internal_modified_sequence = _construct_sequence(sequence, mods)

    new_mods = {'l': labile_mods, 's': static_mods, 'i': isotope_mods, 'n': n_term_mods, 'c': c_term_mods,
                'u': unknown_mods}

    components = split(internal_modified_sequence)
    return [add_mods(''.join(p), new_mods) for p in itertools.combinations_with_replacement(components, size)]