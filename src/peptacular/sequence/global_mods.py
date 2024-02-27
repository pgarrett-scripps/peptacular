import re
from typing import Tuple, List, Dict

from peptacular.types import ModValue
from peptacular.util import convert_type


def pop_isotopic_modifications(sequence: str) -> Tuple[str, List[str]]:
    """
    Remove isotopic modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with isotopic modifications removed.
    :rtype: str

    .. code-block:: python

        >>> pop_isotopic_modifications("<15C>PEPTIDE")
        ('PEPTIDE', ['15C'])

        >>> pop_isotopic_modifications("<15C><13C>PEPTIDE")
        ('PEPTIDE', ['15C', '13C'])

        >>> pop_isotopic_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ('{Oxidation}<[Phospho]@P>PEPT[1.0]IDE-[Amide]', ['15C'])

    """

    global_mods = []

    global_mod_pattern = re.compile(r'\<([^>]+)\>')
    matches = global_mod_pattern.finditer(sequence)

    mods = []

    for match in matches:
        mod = match.group(1)
        global_mods.append(mod)
        if '@' not in mod:
            sequence = sequence.replace(match.group(), '')
            mods.append(mod)

    return sequence, mods


def get_isotopic_modifications(sequence: str) -> List[str]:
    """
    Get isotopic modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The isotopic modifications found in the sequence.
    :rtype: List[str]

    .. code-block:: python

        >>> get_isotopic_modifications("<15C>PEPTIDE")
        ['15C']

        >>> get_isotopic_modifications("<15C><13C>PEPTIDE")
        ['15C', '13C']

        >>> get_isotopic_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ['15C']

    """
    return pop_isotopic_modifications(sequence)[1]


def strip_isotopic_modifications(sequence: str) -> str:
    """
    Remove isotopic modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with isotopic modifications removed.
    :rtype: str

    .. code-block:: python

        >>> strip_isotopic_modifications("<15C>PEPTIDE")
        'PEPTIDE'

        >>> strip_isotopic_modifications("<15C><13C>PEPTIDE")
        'PEPTIDE'

        >>> strip_isotopic_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        '{Oxidation}<[Phospho]@P>PEPT[1.0]IDE-[Amide]'

    """
    return pop_isotopic_modifications(sequence)[0]


def add_isotopic_modifications(sequence: str, mods: List[str], overwrite: bool = False) -> str:
    """
    Add isotopic modifications to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param mods: The isotopic modifications to add.
    :type mods: List[str]

    :return: The sequence with isotopic modifications added.
    :rtype: str

    .. code-block:: python

        >>> add_isotopic_modifications("PEPTIDE", ['15C'])
        '<15C>PEPTIDE'

        >>> add_isotopic_modifications("PEPTIDE", ['15C', '13C'])
        '<15C><13C>PEPTIDE'

        >>> add_isotopic_modifications("{Oxidation}PEPTIDE", ['15C'])
        '<15C>{Oxidation}PEPTIDE'

        >>> add_isotopic_modifications("<13C>PEPTIDE", ['15C'])
        '<15C><13C>PEPTIDE'

        >>> add_isotopic_modifications("<13C>PEPTIDE", ['15C'], overwrite = True)
        '<15C>PEPTIDE'

    """

    if overwrite is True:
        sequence = strip_isotopic_modifications(sequence)

    for mod in mods[::-1]:
        sequence = f'<{mod}>{sequence}'
    return sequence


def pop_static_modifications(sequence: str) -> Tuple[str, List[str]]:
    """
    Remove static modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with static modifications removed.
    :rtype: str

    .. code-block:: python

        >>> pop_static_modifications("<[Phospho]@P>PEPTIDE")
        ('PEPTIDE', ['[Phospho]@P'])

        >>> pop_static_modifications("<[Phospho]@P><[Acetyl]@N>PEPTIDE")
        ('PEPTIDE', ['[Phospho]@P', '[Acetyl]@N'])

        >>> pop_static_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ('{Oxidation}<15C>PEPT[1.0]IDE-[Amide]', ['[Phospho]@P'])

        >>> pop_static_modifications("<[Formula:[13C]H6]@P>PEPTIDE")
        ('PEPTIDE', ['[Formula:[13C]H6]@P'])

    """

    global_mods = []

    global_mod_pattern = re.compile(r'\<([^>]+)\>')
    matches = global_mod_pattern.finditer(sequence)

    mods = []

    for match in matches:
        mod = match.group(1)
        global_mods.append(mod)
        if '@' in mod:
            sequence = sequence.replace(match.group(), '')
            mods.append(mod)

    return sequence, mods


def get_static_modifications(sequence: str) -> List[str]:
    """
    Get static modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The static modifications found in the sequence.
    :rtype: List[str]

    .. code-block:: python

        >>> get_static_modifications("<[Phospho]@P>PEPTIDE")
        ['[Phospho]@P']

        >>> get_static_modifications("<[Phospho]@P><[Acetyl]@N>PEPTIDE")
        ['[Phospho]@P', '[Acetyl]@N']

        >>> get_static_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ['[Phospho]@P']

    """
    return pop_static_modifications(sequence)[1]


def strip_static_modifications(sequence: str) -> str:
    """
    Remove static modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with isotopic modifications removed.
    :rtype: str

    .. code-block:: python

        >>> strip_static_modifications("<[Phospho]@P>PEPTIDE")
        'PEPTIDE'

        >>> strip_static_modifications("<[Phospho]@P><[Acetyl]@N>PEPTIDE")
        'PEPTIDE'

        >>> strip_static_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        '{Oxidation}<15C>PEPT[1.0]IDE-[Amide]'

    """
    return pop_static_modifications(sequence)[0]


def add_static_modifications(sequence: str, mods: List[str], overwrite: bool = False) -> str:
    """
    Add static modifications to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :param mods: The static modifications to add.
    :type mods: List[str]

    :return: The sequence with static modifications added.
    :rtype: str

    .. code-block:: python

        >>> add_static_modifications("PEPTIDE", ['[Phospho]@P'])
        '<[Phospho]@P>PEPTIDE'

        >>> add_static_modifications("PEPTIDE", ['[Phospho]@P', '[Acetyl]@N'])
        '<[Phospho]@P><[Acetyl]@N>PEPTIDE'

        >>> add_static_modifications("{Oxidation}PEPTIDE", ['[Phospho]@P'])
        '<[Phospho]@P>{Oxidation}PEPTIDE'

        >>> add_static_modifications("<[Oxidation]@M>PEPTIDE", ['[Phospho]@P'])
        '<[Phospho]@P><[Oxidation]@M>PEPTIDE'

        >>> add_static_modifications("<[Oxidation]@M>PEPTIDE", ['[Phospho]@P'], overwrite=True)
        '<[Phospho]@P>PEPTIDE'

    """

    if overwrite is True:
        sequence = strip_static_modifications(sequence)

    for mod in mods[::-1]:
        sequence = f'<{mod}>{sequence}'
    return sequence


def pop_global_modifications(sequence: str) -> Tuple[str, List[str], List[str]]:
    """
    Remove isotopic notation from the sequence. Isotopic notation is denoted by <>.

    :param sequence: The amino acid sequence.
    :return: The sequence with isotopic notation removed.

    :return: The sequence
    :rtype: str

    .. code-block:: python

        >>> pop_global_modifications("<13C>PEPTIDE")
        ('PEPTIDE', ['13C'], [])

        >>> pop_global_modifications("<13C><15N>PEPTIDE")
        ('PEPTIDE', ['13C', '15N'], [])

        >>> pop_global_modifications("<13C><15N><[Acetyl]@C>PEPTIDE")
        ('PEPTIDE', ['13C', '15N'], ['[Acetyl]@C'])

        >>> pop_global_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ('{Oxidation}PEPT[1.0]IDE-[Amide]', ['15C'], ['[Phospho]@P'])

    """

    global_mods = []

    global_mod_pattern = re.compile(r'\<([^>]+)\>')
    matches = global_mod_pattern.finditer(sequence)

    for match in matches:
        isotop_mod = match.group(1)
        global_mods.append(isotop_mod)
        sequence = sequence.replace(match.group(), '')

    # separate isotopic and static global sequence
    isotope_mods = []
    static_mods = []
    for mod in global_mods:
        if '@' in mod:
            static_mods.append(mod)
        else:
            isotope_mods.append(mod)

    return sequence, isotope_mods, static_mods


def get_global_modifications(sequence: str) -> Tuple[List[str], List[str]]:
    """
    Get global modifications from the sequence.

    :param sequence: The amino acid sequence.
    :return: The global modifications found in the sequence.

    :return: The isotopic modifications
    :rtype: List[str]

    .. code-block:: python

        >>> get_global_modifications("<13C>PEPTIDE")
        (['13C'], [])

        >>> get_global_modifications("<13C><15N>PEPTIDE")
        (['13C', '15N'], [])

        >>> get_global_modifications("<13C><15N><[Acetyl]@C>PEPTIDE")
        (['13C', '15N'], ['[Acetyl]@C'])

        >>> get_global_modifications("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        (['15C'], ['[Phospho]@P'])

    """
    return pop_global_modifications(sequence)[1:]


def strip_global_modifications(sequence: str) -> str:
    """
    Remove isotopic notation from the sequence. Isotopic notation is denoted by <>.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with isotopic notation removed.
    :rtype: str

    .. code-block:: python

        >>> strip_global_modifications("<13C>PEPTIDE")
        'PEPTIDE'

        >>> strip_global_modifications("<13C><15N>PEPTIDE")
        'PEPTIDE'

        >>> strip_global_modifications("<13C><15N><[Acetyl]@C>PEPTIDE")
        'PEPTIDE'

    """
    return pop_global_modifications(sequence)[0]


def add_global_modifications(sequence: str, isotope_mods: List[str], static_mods: List[str], overwrite: bool = False) -> str:
    """
    Add isotopic notation to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param isotope_mods: The isotopic modifications to add.
    :type isotope_mods: List[str]
    :param static_mods: The static modifications to add.
    :type static_mods: List[str]

    :return: The sequence with isotopic notation added.
    :rtype: str

    .. code-block:: python

        >>> add_global_modifications("PEPTIDE", ['13C'], [])
        '<13C>PEPTIDE'

        >>> add_global_modifications("PEPTIDE", ['13C', '15N'], [])
        '<13C><15N>PEPTIDE'

        >>> add_global_modifications("PEPTIDE", ['13C', '15N'], ['[Acetyl]@C'])
        '<[Acetyl]@C><13C><15N>PEPTIDE'

        >>> add_global_modifications("<[Acetyl]@C><13C>PEPTIDE", ['15N'], ['[Oxidation]@M'])
        '<[Oxidation]@M><15N><[Acetyl]@C><13C>PEPTIDE'

        >>> add_global_modifications("<[Acetyl]@C><13C>PEPTIDE", ['15N'], ['[Oxidation]@M'], overwrite=True)
        '<[Oxidation]@M><15N>PEPTIDE'

    """

    if overwrite is True:
        sequence = strip_global_modifications(sequence)

    for mod in isotope_mods[::-1]:
        sequence = f'<{mod}>{sequence}'
    for mod in static_mods[::-1]:
        sequence = f'<{mod}>{sequence}'
    return sequence


def parse_static_modifications(mods) -> Dict[str, ModValue]:
    """
    Parse static modifications into a dictionary.

    :param mods: List of static modifications.
    :return: Dictionary of static modifications.

    :return: Dictionary of static modifications.
    :rtype: Dict[str, str]

    .. code-block:: python

        >>> parse_static_modifications(['[1]@C', '[3.1415]@S,D'])
        {'C': 1, 'S': 3.1415, 'D': 3.1415}

        >>> parse_static_modifications(['[Acetyl]@C', '[Phospho]@S,D'])
        {'C': 'Acetyl', 'S': 'Phospho', 'D': 'Phospho'}

        >>> parse_static_modifications(['[Acetyl]@N', '[Phospho]@S,D', '[Methyl]@K'])
        {'N': 'Acetyl', 'S': 'Phospho', 'D': 'Phospho', 'K': 'Methyl'}

        >>> parse_static_modifications(['[Acetyl]@N', '[Phospho]@S,D', '[Methyl]@K', '[Methyl]@R'])
        {'N': 'Acetyl', 'S': 'Phospho', 'D': 'Phospho', 'K': 'Methyl', 'R': 'Methyl'}

        >>> parse_static_modifications(['[Formula:[13C]H20]@C'])
        {'C': 'Formula:[13C]H20'}

    """

    static_mod_dict = {}

    for mod in mods:
        mod_info, residues = mod.split('@')
        mod_name = mod_info[1:-1]  # Remove the brackets around the modification name
        for residue in residues.split(','):  # Split on ',' for multiple residues
            static_mod_dict[residue] = convert_type(mod_name)  # Assign the modification to each residue

    return static_mod_dict


def parse_isotope_modifications(mods: List[str]) -> Dict[str, str]:
    """
    Parse isotope modifications into a dictionary.

    :param mods: List of isotope modifications.
    :return: Dictionary of isotope modifications.

    :return: Dictionary of isotope modifications.
    :rtype: Dict[str, str]

    .. code-block:: python

        >>> parse_isotope_modifications(['13C', '15N', 'D'])
        {'C': '13C', 'N': '15N', 'H': 'D'}

        >>> parse_isotope_modifications(['13C', '15N', 'T'])
        {'C': '13C', 'N': '15N', 'H': 'T'}

        >>> parse_isotope_modifications(['13C', '15N', 'D', 'T'])
        {'C': '13C', 'N': '15N', 'H': 'T'}

    """

    isotope_map = {}
    for mod in mods:
        # remove digits
        pattern = r'[0-9]'
        base_aa = re.sub(pattern, '', mod)
        isotope_map[base_aa] = mod

    # If any keys are D or T, then replace them with H
    if 'D' in isotope_map:
        isotope_map['H'] = isotope_map.pop('D')

    if 'T' in isotope_map:
        isotope_map['H'] = isotope_map.pop('T')

    return isotope_map
