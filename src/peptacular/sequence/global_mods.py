import re
from typing import Tuple, List, Dict

from peptacular.types import ModValue
from peptacular.util import convert_type, map_bracket_content_to_index


def pop_isotope_mods(sequence: str) -> Tuple[str, List[str]]:
    """
    Remove isotopic modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with isotopic modifications removed.
    :rtype: Tuple[str, List[str]]

    .. code-block:: python

        >>> pop_isotope_mods("<15C>PEPTIDE")
        ('PEPTIDE', ['15C'])

        >>> pop_isotope_mods("<15C><13C>PEPTIDE")
        ('PEPTIDE', ['15C', '13C'])

        >>> pop_isotope_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ('{Oxidation}<[Phospho]@P>PEPT[1.0]IDE-[Amide]', ['15C'])

        >>> pop_isotope_mods("PEPTIDE")
        ('PEPTIDE', [])

        >>> pop_isotope_mods("<15C>")
        ('', ['15C'])

    """

    global_mods = []

    global_mod_pattern = re.compile(r'<([^>]+)>')
    matches = global_mod_pattern.finditer(sequence)

    mods = []

    for match in matches:
        mod = match.group(1)
        global_mods.append(mod)
        if '@' not in mod:
            sequence = sequence.replace(match.group(), '')
            mods.append(mod)

    return sequence, mods


def get_isotope_mods(sequence: str) -> List[str]:
    """
    Get isotopic modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The isotopic modifications found in the sequence.
    :rtype: List[str]

    .. code-block:: python

        >>> get_isotope_mods("<15C>PEPTIDE")
        ['15C']

        >>> get_isotope_mods("<15C><13C>PEPTIDE")
        ['15C', '13C']

        >>> get_isotope_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ['15C']

        >>> get_isotope_mods("PEPTIDE")
        []

        >>> get_isotope_mods("<15C>")
        ['15C']

    """
    return pop_isotope_mods(sequence)[1]


def strip_isotope_mods(sequence: str) -> str:
    """
    Remove isotopic modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with isotopic modifications removed.
    :rtype: str

    .. code-block:: python

        >>> strip_isotope_mods("<15C>PEPTIDE")
        'PEPTIDE'

        >>> strip_isotope_mods("<15C><13C>PEPTIDE")
        'PEPTIDE'

        >>> strip_isotope_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        '{Oxidation}<[Phospho]@P>PEPT[1.0]IDE-[Amide]'

        >>> strip_isotope_mods("PEPTIDE")
        'PEPTIDE'

        >>> strip_isotope_mods("<15C>")
        ''

    """
    return pop_isotope_mods(sequence)[0]


def add_isotope_mods(sequence: str, mods: List[str], overwrite: bool = False) -> str:
    """
    Add isotopic modifications to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param mods: The isotopic modifications to add.
    :type mods: List[str]
    :param overwrite: Whether to overwrite existing modifications.
    :type overwrite: bool

    :return: The sequence with isotopic modifications added.
    :rtype: str

    .. code-block:: python

        >>> add_isotope_mods("PEPTIDE", ['15C'])
        '<15C>PEPTIDE'

        >>> add_isotope_mods("PEPTIDE", ['15C', '13C'])
        '<15C><13C>PEPTIDE'

        >>> add_isotope_mods("{Oxidation}PEPTIDE", ['15C'])
        '<15C>{Oxidation}PEPTIDE'

        >>> add_isotope_mods("<13C>PEPTIDE", ['15C'])
        '<15C><13C>PEPTIDE'

        >>> add_isotope_mods("<13C>PEPTIDE", ['15C'], overwrite = True)
        '<15C>PEPTIDE'

        >>> add_isotope_mods("PEPTIDE", [])
        'PEPTIDE'

        >>> add_isotope_mods("", ['15C'])
        '<15C>'

    """

    if overwrite is True:
        sequence = strip_isotope_mods(sequence)

    for mod in mods[::-1]:
        sequence = f'<{mod}>{sequence}'
    return sequence


def pop_static_mods(sequence: str) -> Tuple[str, List[str]]:
    """
    Remove static modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with static modifications removed.
    :rtype: str

    .. code-block:: python

        >>> pop_static_mods("<[Phospho]@P>PEPTIDE")
        ('PEPTIDE', ['[Phospho]@P'])

        >>> pop_static_mods("<[Phospho]@P><[Acetyl]@N>PEPTIDE")
        ('PEPTIDE', ['[Phospho]@P', '[Acetyl]@N'])

        >>> pop_static_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ('{Oxidation}<15C>PEPT[1.0]IDE-[Amide]', ['[Phospho]@P'])

        >>> pop_static_mods("<[Formula:[13C]H6]@P>PEPTIDE")
        ('PEPTIDE', ['[Formula:[13C]H6]@P'])

        >>> pop_static_mods("PEPTIDE")
        ('PEPTIDE', [])

        >>> pop_static_mods("<[Phospho]@P>")
        ('', ['[Phospho]@P'])

    """

    global_mods = []

    global_mod_pattern = re.compile(r'<([^>]+)>')
    matches = global_mod_pattern.finditer(sequence)

    mods = []

    for match in matches:
        mod = match.group(1)
        global_mods.append(mod)
        if '@' in mod:
            sequence = sequence.replace(match.group(), '')
            mods.append(mod)

    return sequence, mods


def get_static_mods(sequence: str) -> List[str]:
    """
    Get static modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The static modifications found in the sequence.
    :rtype: List[str]

    .. code-block:: python

        >>> get_static_mods("<[Phospho]@P>PEPTIDE")
        ['[Phospho]@P']

        >>> get_static_mods("<[Phospho]@P><[Acetyl]@N>PEPTIDE")
        ['[Phospho]@P', '[Acetyl]@N']

        >>> get_static_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ['[Phospho]@P']

        >>> get_static_mods("PEPTIDE")
        []

        >>> get_static_mods("<[Phospho]@P>")
        ['[Phospho]@P']

    """
    return pop_static_mods(sequence)[1]


def strip_static_mods(sequence: str) -> str:
    """
    Remove static modifications from the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with isotopic modifications removed.
    :rtype: str

    .. code-block:: python

        >>> strip_static_mods("<[Phospho]@P>PEPTIDE")
        'PEPTIDE'

        >>> strip_static_mods("<[Phospho]@P><[Acetyl]@N>PEPTIDE")
        'PEPTIDE'

        >>> strip_static_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        '{Oxidation}<15C>PEPT[1.0]IDE-[Amide]'

        >>> strip_static_mods("PEPTIDE")
        'PEPTIDE'

        >>> strip_static_mods("<[Phospho]@P>")
        ''

    """
    return pop_static_mods(sequence)[0]


def add_static_mods(sequence: str, mods: List[str], overwrite: bool = False) -> str:
    """
    Add static modifications to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param mods: The static modifications to add.
    :type mods: List[str]
    :param overwrite: Whether to overwrite existing modifications.
    :type overwrite: bool


    :return: The sequence with static modifications added.
    :rtype: str

    .. code-block:: python

        >>> add_static_mods("PEPTIDE", ['[Phospho]@P'])
        '<[Phospho]@P>PEPTIDE'

        >>> add_static_mods("PEPTIDE", ['[Phospho]@P', '[Acetyl]@N'])
        '<[Phospho]@P><[Acetyl]@N>PEPTIDE'

        >>> add_static_mods("{Oxidation}PEPTIDE", ['[Phospho]@P'])
        '<[Phospho]@P>{Oxidation}PEPTIDE'

        >>> add_static_mods("<[Oxidation]@M>PEPTIDE", ['[Phospho]@P'])
        '<[Phospho]@P><[Oxidation]@M>PEPTIDE'

        >>> add_static_mods("<[Oxidation]@M>PEPTIDE", ['[Phospho]@P'], overwrite=True)
        '<[Phospho]@P>PEPTIDE'

        >>> add_static_mods("PEPTIDE", [])
        'PEPTIDE'

        >>> add_static_mods("", ['[Phospho]@P'])
        '<[Phospho]@P>'

    """

    if overwrite is True:
        sequence = strip_static_mods(sequence)

    for mod in mods[::-1]:
        sequence = f'<{mod}>{sequence}'
    return sequence


def pop_global_mods(sequence: str) -> Tuple[str, List[str], List[str]]:
    """
    Remove isotopic notation from the sequence. Isotopic notation is denoted by <>.

    :param sequence: The amino acid sequence.
    :return: The sequence with isotopic notation removed.

    :return: The sequence
    :rtype: str

    .. code-block:: python

        >>> pop_global_mods("<13C>PEPTIDE")
        ('PEPTIDE', ['13C'], [])

        >>> pop_global_mods("<13C><15N>PEPTIDE")
        ('PEPTIDE', ['13C', '15N'], [])

        >>> pop_global_mods("<13C><15N><[Acetyl]@C>PEPTIDE")
        ('PEPTIDE', ['13C', '15N'], ['[Acetyl]@C'])

        >>> pop_global_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        ('{Oxidation}PEPT[1.0]IDE-[Amide]', ['15C'], ['[Phospho]@P'])

        >>> pop_global_mods("PEPTIDE")
        ('PEPTIDE', [], [])

        >>> pop_global_mods("<13C><[Acetyl]@C>")
        ('', ['13C'], ['[Acetyl]@C'])

    """

    mods = []
    pattern = re.compile(r'<([^>]+)>')

    for m in pattern.finditer(sequence):
        isotope_mod = m.group(1)
        mods.append(isotope_mod)
        sequence = sequence.replace(m.group(), '')

    # separate isotopic and static global sequence
    isotope_mods = []
    static_mods = []
    for mod in mods:
        if '@' in mod:
            static_mods.append(mod)
        else:
            isotope_mods.append(mod)

    return sequence, isotope_mods, static_mods


def get_global_mods(sequence: str) -> Tuple[List[str], List[str]]:
    """
    Get global modifications from the sequence.

    :param sequence: The amino acid sequence.
    :return: The global modifications found in the sequence.

    :return: The isotopic modifications
    :rtype: List[str]

    .. code-block:: python

        >>> get_global_mods("<13C>PEPTIDE")
        (['13C'], [])

        >>> get_global_mods("<13C><15N>PEPTIDE")
        (['13C', '15N'], [])

        >>> get_global_mods("<13C><15N><[Acetyl]@C>PEPTIDE")
        (['13C', '15N'], ['[Acetyl]@C'])

        >>> get_global_mods("{Oxidation}<15C><[Phospho]@P>PEPT[1.0]IDE-[Amide]")
        (['15C'], ['[Phospho]@P'])

        >>> get_global_mods("PEPTIDE")
        ([], [])

        >>> get_global_mods("<13C><[Acetyl]@C>")
        (['13C'], ['[Acetyl]@C'])

    """
    return pop_global_mods(sequence)[1:]


def strip_global_mods(sequence: str) -> str:
    """
    Remove isotopic notation from the sequence. Isotopic notation is denoted by <>.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :return: The sequence with isotopic notation removed.
    :rtype: str

    .. code-block:: python

        >>> strip_global_mods("<13C>PEPTIDE")
        'PEPTIDE'

        >>> strip_global_mods("<13C><15N>PEPTIDE")
        'PEPTIDE'

        >>> strip_global_mods("<13C><15N><[Acetyl]@C>PEPTIDE")
        'PEPTIDE'

        >>> strip_global_mods("PEPTIDE")
        'PEPTIDE'

        >>> strip_global_mods("<13C><[Acetyl]@C>")
        ''

    """
    return pop_global_mods(sequence)[0]


def add_global_mods(sequence: str, isotope_mods: List[str], static_mods: List[str],
                    overwrite: bool = False) -> str:
    """
    Add isotopic notation to the sequence.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param isotope_mods: The isotopic modifications to add.
    :type isotope_mods: List[str]
    :param static_mods: The static modifications to add.
    :type static_mods: List[str]
    :param overwrite: Whether to overwrite existing modifications.
    :type overwrite: bool


    :return: The sequence with isotopic notation added.
    :rtype: str

    .. code-block:: python

        >>> add_global_mods("PEPTIDE", ['13C'], [])
        '<13C>PEPTIDE'

        >>> add_global_mods("PEPTIDE", ['13C', '15N'], [])
        '<13C><15N>PEPTIDE'

        >>> add_global_mods("PEPTIDE", ['13C', '15N'], ['[Acetyl]@C'])
        '<[Acetyl]@C><13C><15N>PEPTIDE'

        >>> add_global_mods("<[Acetyl]@C><13C>PEPTIDE", ['15N'], ['[Oxidation]@M'])
        '<[Oxidation]@M><15N><[Acetyl]@C><13C>PEPTIDE'

        >>> add_global_mods("<[Acetyl]@C><13C>PEPTIDE", ['15N'], ['[Oxidation]@M'], overwrite=True)
        '<[Oxidation]@M><15N>PEPTIDE'

        >>> add_global_mods("PEPTIDE", [], [])
        'PEPTIDE'

        >>> add_global_mods("", ['13C'], ['[Oxidation]@M'])
        '<[Oxidation]@M><13C>'

    """

    if overwrite is True:
        sequence = strip_global_mods(sequence)

    for mod in isotope_mods[::-1]:
        sequence = f'<{mod}>{sequence}'
    for mod in static_mods[::-1]:
        sequence = f'<{mod}>{sequence}'
    return sequence


def parse_static_mods(mods: List[str]) -> Dict[str, List[ModValue]]:
    """
    Parse static modifications into a dictionary.

    :param mods: List of static modifications.
    :type mods: List[str]

    :return: Dictionary of static modifications.
    :rtype: Dict[str, str]

    .. code-block:: python

        >>> parse_static_mods(['[1]@C', '[3.1415]@S,D'])
        {'C': [1], 'S': [3.1415], 'D': [3.1415]}

        >>> parse_static_mods(['[1]^2@C', '[3.1415]@S,D'])
        {'C': [1, 1], 'S': [3.1415], 'D': [3.1415]}

        >>> parse_static_mods(['[1][2]@C', '[3.1415]@S,D'])
        {'C': [1, 2], 'S': [3.1415], 'D': [3.1415]}

        >>> parse_static_mods(['[Acetyl]@C', '[Phospho]@S,D'])
        {'C': ['Acetyl'], 'S': ['Phospho'], 'D': ['Phospho']}

        >>> parse_static_mods(['[Formula:[13C]H20]@C'])
        {'C': ['Formula:[13C]H20']}

        >>> parse_static_mods(['[100]@P'])
        {'P': [100]}

        >>> parse_static_mods([])
        {}

    """

    static_mod_dict = {}

    for mod in mods:
        mod_info, residues = mod.split('@')

        _, bracket_mods = map_bracket_content_to_index(mod_info, True)
        for bracket_mod in bracket_mods.values():
            bracket_mod = [convert_type(m) for m in bracket_mod]
            for residue in residues.split(','):  # Split on ',' for multiple residues
                static_mod_dict[residue] = bracket_mod  # Assign the modification to each residue

    return static_mod_dict


def write_static_mods(mods: Dict[str, List[ModValue]]) -> List[str]:
    """
    Write static modifications from a dictionary.

    :param mods: Dictionary of static modifications.
    :type: Dict[str, ModValue]

    :return: List of static modifications.
    :rtype: List[str]

    .. code-block:: python

        >>> write_static_mods({'C': [1], 'S': [3.1415], 'D': [3.1415]})
        ['[1]@C', '[3.1415]@S,D']

        >>> write_static_mods({'C': ['Formula:[13C]H20']})
        ['[Formula:[13C]H20]@C']

        >>> write_static_mods({})
        []

    """

    reverse_mods = {}

    for residue, mod in mods.items():
        mod = tuple(mod)
        if isinstance(mod, str):
            reverse_mods.setdefault(mod, []).append(residue)
        else:
            reverse_mods.setdefault(mod, []).append(residue)

    mod_list = []
    for mod, residues in reverse_mods.items():
        mod_str = ''.join([f'[{m}]' for m in mod])
        if len(residues) == 1:
            mod_list.append(f'{mod_str}@{residues[0]}')
        else:
            mod_list.append(f'{mod_str}@{",".join(residues)}')

    return mod_list


def parse_isotope_mods(mods: List[str]) -> Dict[str, str]:
    """
    Parse isotope modifications into a dictionary.

    :param mods: List of isotope modifications.
    :return: Dictionary of isotope modifications.

    :return: Dictionary of isotope modifications.
    :rtype: Dict[str, str]

    .. code-block:: python

        >>> parse_isotope_mods(['13C', '15N', 'D'])
        {'C': '13C', 'N': '15N', 'H': 'D'}

        >>> parse_isotope_mods(['13C', '15N', 'T'])
        {'C': '13C', 'N': '15N', 'H': 'T'}

        >>> parse_isotope_mods(['13C', '15N', 'H'])
        {'C': '13C', 'N': '15N', 'H': 'H'}

        >>> parse_isotope_mods([])
        {}

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


def write_isotope_mods(mods: Dict[str, str]) -> List[str]:
    """
    Write isotope modifications from a dictionary.

    :param mods: Dictionary of isotope modifications.
    :type: Dict[str, str]

    :return: List of isotope modifications.
    :rtype: List[str]

    .. code-block:: python

        >>> write_isotope_mods({'C': '13C', 'N': '15N', 'H': 'D'})
        ['13C', '15N', 'D']

        >>> write_isotope_mods({'C': '13C', 'N': '15N', 'H': 'T'})
        ['13C', '15N', 'T']

        >>> write_isotope_mods({})
        []

    """

    return list(mods.values())
