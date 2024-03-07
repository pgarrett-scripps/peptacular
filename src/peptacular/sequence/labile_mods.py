import re
from typing import List, Tuple

from peptacular.types import ModValue
from peptacular.util import convert_type


def pop_labile_mods(sequence) -> Tuple[str, List[ModValue]]:
    """
    Remove labile modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with labile modifications removed and the labile modifications.
    :rtype: Tuple[str, List[ModValue]]

    .. code-block:: python

        >>> pop_labile_mods("{Oxidation}PEPTIDE")
        ('PEPTIDE', ['Oxidation'])

        >>> pop_labile_mods("{Oxidation}{1.0}PEPTIDE")
        ('PEPTIDE', ['Oxidation', 1.0])

        >>> pop_labile_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ('<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]', ['Oxidation'])

        >>> pop_labile_mods("{Oxidation}")
        ('', ['Oxidation'])

        >>> pop_labile_mods("PEPTIDE")
        ('PEPTIDE', [])

    """

    labile_mods = []

    labile_mod_pattern = re.compile(r'{([^}]+)}')
    matches = labile_mod_pattern.finditer(sequence)

    for match in matches:
        labile_mod = convert_type(match.group(1))
        labile_mods.append(labile_mod)
        sequence = sequence.replace(match.group(), '')

    return sequence, labile_mods


def get_labile_mods(sequence: str) -> List[ModValue]:
    """
    Get labile modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The labile modifications found in the sequence.
    :rtype: List[ModValue]

    .. code-block:: python

        >>> get_labile_mods("{Oxidation}PEPTIDE")
        ['Oxidation']

        >>> get_labile_mods("{Oxidation}{1.0}PEPTIDE")
        ['Oxidation', 1.0]

        >>> get_labile_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ['Oxidation']

        >>> get_labile_mods("{Oxidation}")
        ['Oxidation']

        >>> get_labile_mods("PEPTIDE")
        []

    """
    return pop_labile_mods(sequence)[1]


def strip_labile_mods(sequence: str) -> str:
    """
    Remove labile modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with labile modifications removed.
    :rtype: str

    .. code-block:: python

        >>> strip_labile_mods("{Oxidation}PEPTIDE")
        'PEPTIDE'

        >>> strip_labile_mods("{Oxidation}{1.0}PEPTIDE")
        'PEPTIDE'

        >>> strip_labile_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        '<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]'

        >>> strip_labile_mods("{Oxidation}")
        ''

        >>> strip_labile_mods("PEPTIDE")
        'PEPTIDE'
    """

    sequence, _ = pop_labile_mods(sequence)
    return sequence


def add_labile_mods(sequence: str, labile_mods: List[ModValue]) -> str:
    """
    Add labile modifications to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param labile_mods: The labile modifications to add to the sequence.
    :type labile_mods: List[ModValue]

    :return: The sequence with labile modifications added.
    :rtype: str

    .. code-block:: python

        >>> add_labile_mods("PEPTIDE", ['Oxidation'])
        '{Oxidation}PEPTIDE'

        >>> add_labile_mods("PEPTIDE", ['Oxidation', 1.0])
        '{Oxidation}{1.0}PEPTIDE'

        >>> add_labile_mods("<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]", ['Oxidation'])
        '{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]'

        >>> add_labile_mods("", ['Oxidation'])
        '{Oxidation}'

        >>> add_labile_mods("PEPTIDE", [])
        'PEPTIDE'
    """

    for labile_mod in labile_mods[::-1]:
        sequence = f'{{{labile_mod}}}{sequence}'

    return sequence
