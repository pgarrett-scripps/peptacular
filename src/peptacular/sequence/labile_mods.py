import re
from typing import List

from peptacular.types import ModValue
from peptacular.util import convert_type


def pop_labile_modifications(sequence) -> tuple[str, List[ModValue]]:
    """
    Remove labile modifications from the sequence. Labile sequence are denoated by {}

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with labile modifications removed.
    :rtype: str

    .. code-block:: python

        >>> pop_labile_modifications("{Oxidation}PEPTIDE")
        ('PEPTIDE', ['Oxidation'])

        >>> pop_labile_modifications("{Oxidation}{1.0}PEPTIDE")
        ('PEPTIDE', ['Oxidation', 1.0])

        >>> pop_labile_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ('<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]', ['Oxidation'])

    """

    labile_mods = []

    labile_mod_pattern = re.compile(r'\{([^}]+)\}')
    matches = labile_mod_pattern.finditer(sequence)

    for match in matches:
        labile_mod = convert_type(match.group(1))
        labile_mods.append(labile_mod)
        sequence = sequence.replace(match.group(), '')

    return sequence, labile_mods


def get_labile_modifications(sequence: str) -> List[ModValue]:
    """
    Get labile modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The labile modifications found in the sequence.
    :rtype: List[str]

    .. code-block:: python

        >>> get_labile_modifications("{Oxidation}PEPTIDE")
        ['Oxidation']

        >>> get_labile_modifications("{Oxidation}{1.0}PEPTIDE")
        ['Oxidation', 1.0]

        >>> get_labile_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ['Oxidation']

    """
    return pop_labile_modifications(sequence)[1]


def strip_labile_modifications(sequence: str) -> str:
    """
    Remove labile modifications from the sequence. Labile sequence are denoated by {}

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with labile modifications removed.
    :rtype: str

    .. code-block:: python

        >>> strip_labile_modifications("{Oxidation}PEPTIDE")
        'PEPTIDE'

        >>> strip_labile_modifications("{Oxidation}{1.0}PEPTIDE")
        'PEPTIDE'

        >>> strip_labile_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        '<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]'
    """

    sequence, _ = pop_labile_modifications(sequence)
    return sequence


def add_labile_modifications(sequence: str, labile_mods: List[ModValue]) -> str:
    """
    Add labile modifications to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param labile_mods: The labile modifications to add to the sequence.
    :type labile_mods: List[str]

    :return: The sequence with labile modifications added.
    :rtype: str

    .. code-block:: python

        >>> add_labile_modifications("PEPTIDE", ['Oxidation'])
        '{Oxidation}PEPTIDE'

        >>> add_labile_modifications("PEPTIDE", ['Oxidation', 1.0])
        '{Oxidation}{1.0}PEPTIDE'

        >>> add_labile_modifications("<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]", ['Oxidation'])
        '{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]'
    """

    for labile_mod in labile_mods[::-1]:
        sequence = f'{{{labile_mod}}}{sequence}'

    return sequence
