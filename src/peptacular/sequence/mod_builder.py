from copy import deepcopy
from typing import Dict, List, Generator

from peptacular.sequence import pop_mods, add_mods
from peptacular.types import ModDict, ModIndex, ModValue
from peptacular.util import get_regex_match_indices


def apply_static_mods(sequence: str,
                      residue_mods: ModDict,
                      nterm_mods: ModDict = None,
                      cterm_mods: ModDict = None,
                      overwrite: bool = True) -> str:
    """
    Add static modifications to an amino acid sequence. If a modification is already present in the sequence,
    it will be replaced by the new one.

    :param sequence: Original amino acid sequence.
    :type sequence: str
    :param residue_mods: Dictionary mapping amino acids to the mass of their modifications.
    :type residue_mods: Dict[str, Any]
    :param nterm_mods: Dictionary mapping the N-terminal amino acid to the mass of its modification.
    :type nterm_mods: Dict[str, Any]
    :param cterm_mods: Dictionary mapping the C-terminal amino acid to the mass of its modification.
    :type cterm_mods: Dict[str, Any]
    :param overwrite: If True, existing modifications will be replaced by the new ones.
    :type overwrite: bool

    :return: Modified amino acid sequence.
    :rtype: str

    .. code-block:: python

        # Applies static modifcation to all matching residues in the sequence:
        >>> apply_static_mods('PEPTIDE', {'P': ['phospho']})
        'P[phospho]EP[phospho]TIDE'

        # Can specify multiple static modifications:
        >>> apply_static_mods('PEPTIDE', {'P': ['phospho'], 'E': ['3.0']})
        'P[phospho]E[3.0]P[phospho]TIDE[3.0]'

        # Works with already modified sequences:
        >>> apply_static_mods('[3.14]-PEPTIDE-[1.234]', {'P': ['phospho']})
        '[3.14]-P[phospho]EP[phospho]TIDE-[1.234]'

        # By default, any overlapping modifications will be replaced:
        >>> apply_static_mods('PEP[1.234]TIDE', {'P': ['phospho']})
        'P[phospho]EP[phospho]TIDE'

        # To not reaplce existing modifications, set override to False:
        >>> apply_static_mods('PEP[1.234]TIDE', {'P': ['phospho']}, overwrite=False)
        'P[phospho]EP[1.234][phospho]TIDE'

        # Can also use regular expressions to match multiple residues:
        >>> apply_static_mods('PEPTIDE', {'(?<=P)E': ['phospho']})
        'PE[phospho]PTIDE'

        >>> apply_static_mods('PEPTIDE', {'PE': ['phospho']})
        'P[phospho]EPTIDE'

        >>> apply_static_mods('<13C>PEPTIDE', {'P': ['phospho']})
        '<13C>P[phospho]EP[phospho]TIDE'

        >>> apply_static_mods('(?PE)PTIDE', {'P': ['phospho']})
        '(?P[phospho]E)P[phospho]TIDE'

        # Term mods
        >>> apply_static_mods('PEPTIDE', {}, nterm_mods={'P': ['acetyl']}, cterm_mods={'P': ['amide']})
        '[acetyl]-PEPTIDE'

        >>> apply_static_mods('PEPTIDE', {}, nterm_mods={'PE': ['acetyl']}, cterm_mods={'P': ['amide']})
        '[acetyl]-PEPTIDE'

        >>> apply_static_mods('PEPTIDE', {}, nterm_mods={'PP': ['acetyl']}, cterm_mods={'(?<=D)E': ['amide']})
        'PEPTIDE-[amide]'

        >>> apply_static_mods('PEPTIDE', {}, nterm_mods={'': ['acetyl']}, cterm_mods={'': ['amide']})
        '[acetyl]-PEPTIDE-[amide]'

    """

    stripped_sequence, original_mods = pop_mods(sequence)

    new_mods = {}
    for regex_str, mod in residue_mods.items():
        for mod_index in get_regex_match_indices(stripped_sequence, regex_str):
            new_mods[mod_index] = mod

    if nterm_mods:
        for regex_str, mod in nterm_mods.items():
            for mod_index in get_regex_match_indices(stripped_sequence, regex_str):
                if mod_index == 0:
                    new_mods['n'] = mod

    if cterm_mods:
        for regex_str, mod in cterm_mods.items():
            for mod_index in get_regex_match_indices(stripped_sequence, regex_str):
                if mod_index == len(stripped_sequence) - 1:
                    new_mods['c'] = mod

    return add_mods(sequence, new_mods, overwrite)


def _apply_variable_mods_rec(mods: Dict[ModIndex, List[List[ModValue]]],
                             sequence: str,
                             index: int,
                             current_mods: ModDict,
                             max_mod_count: int,
                             overwrite: bool) -> Generator[ModDict, None, None]:
    """
    Recursively applies variable modifications to the amino acid sequence.

    :param mods: Dictionary mapping amino acids to the mass of their modifications.
    :type mods: Dict[int, Any]
    :param sequence: the unmodified sequence.
    :type sequence: str
    :param index: Current index on the amino acid sequence.
    :type index: int
    :param current_mods: Dictionary mapping the indices of the amino acids to the mass of their modifications.
    :type current_mods: Dict[int, str]
    :param max_mod_count: Maximum number of modifications allowed on the amino acid sequence.
    :type max_mod_count: int
    :param overwrite: If True, existing modifications will be replaced by the new ones.
    :type overwrite: bool

    :yield: Modified dictionary mapping the indices of the amino acid sequence to the mass of their modifications.
    :rtype: Generator[Dict[int, float], None, None]
    """

    if index == len(sequence) or len(current_mods) == max_mod_count:
        yield current_mods
        return

    original_mod_dict = deepcopy(current_mods)

    curr_mods = mods.get(index, None)

    if curr_mods is not None:
        for curr_mod in curr_mods:
            updated_mod_dict = deepcopy(current_mods)

            if index in current_mods:
                if overwrite is True:
                    updated_mod_dict[index] = curr_mod
                    yield from _apply_variable_mods_rec(mods, sequence, index + 1, updated_mod_dict,
                                                        max_mod_count, overwrite)
            else:
                updated_mod_dict[index] = curr_mod
                yield from _apply_variable_mods_rec(mods, sequence, index + 1, updated_mod_dict,
                                                    max_mod_count, overwrite)
    yield from _apply_variable_mods_rec(mods, sequence, index + 1, original_mod_dict, max_mod_count, overwrite)


def _variable_mods_builder(sequence: str,
                           mod_map: Dict[str, List[List[ModValue]]],
                           max_mods: int,
                           overwrite: bool = False) -> List[str]:
    """
    Apply variable modifications to a sequence.

    :param sequence: The sequence to be modified.
    :type sequence: str
    :param mod_map: Dictionary mapping amino acids / regex to a list of modifications
    :type mod_map: Dict[str, List[List[ModValue]]]
    :param max_mods: Maximum number of modifications allowed on the amino acid sequence.
    :type max_mods: int
    :param overwrite: If True, existing modifications will be replaced by the new ones.
    :type overwrite: bool

    :return: List of all possible modified amino acid sequences.
    :rtype: List[str]
    """

    stripped_sequence, original_mods = pop_mods(sequence)

    new_mod_map: Dict[int, List[List[ModValue]]] = {}
    for regex_str, mods in mod_map.items():
        for mod_index in get_regex_match_indices(stripped_sequence, regex_str):
            for mod in mods:
                new_mod_map.setdefault(mod_index, []).append(mod)

    mods = _apply_variable_mods_rec(new_mod_map, stripped_sequence, 0, original_mods,
                                    max_mods + len(original_mods), overwrite)

    return [add_mods(stripped_sequence, mod) for mod in mods]


def apply_variable_mods(sequence: str,
                        mod_map: Dict[str, List[List[ModValue]]],
                        max_mods: int,
                        nterm_mods: Dict[str, List[List[ModValue]]] = None,
                        cterm_mods: Dict[str, List[List[ModValue]]] = None,
                        overwrite: bool = False) -> List[str]:
    """
    Apply variable modifications to a sequence.

    :param sequence: Original amino acid sequence.
    :type sequence: str
    :param mod_map: Dictionary mapping amino acids  /regex to modifications.
    :type mod_map: Dict[str, List[List[ModValue]]]
    :param max_mods: Maximum number of modifications allowed on the peptide.
    :type max_mods: int
    :param nterm_mods: Dictionary mapping the N-terminal amino acids / regex to modifications.
    :type nterm_mods: Dict[str, List[List[ModValue]]]
    :param cterm_mods: Dictionary mapping the C-terminal amino acids / regex to modifications.
    :type cterm_mods: Dict[str, List[List[ModValue]]]
    :param overwrite: If True, existing modifications will be replaced by the new ones.
    :type overwrite: bool

    :return: List of all possible modified peptide sequences.
    :rtype: List[str]

    .. code-block:: python

        >>> apply_variable_mods('PEPTIDE', {'P': [['phospho']]}, 1)
        ['P[phospho]EPTIDE', 'PEP[phospho]TIDE', 'PEPTIDE']

        # When multiple mods are specified in the same group, they are applied together
        >>> apply_variable_mods('PEPTIDE', {'P': [['phospho', 1]]}, 1)
        ['P[phospho][1]EPTIDE', 'PEP[phospho][1]TIDE', 'PEPTIDE']

        # Can also break up groups to apply modifications separately
        >>> apply_variable_mods('PEPTIDE', {'P': [['phospho'], [1]]}, 1)
        ['P[phospho]EPTIDE', 'P[1]EPTIDE', 'PEP[phospho]TIDE', 'PEP[1]TIDE', 'PEPTIDE']

        # Setting max_mods to 0 will return the original sequence:
        >>> apply_variable_mods('PEPTIDE', {'P': [['phospho']]}, 0)
        ['PEPTIDE']

        # Works with already modified sequences, and will not replace existing modifications:
        >>> apply_variable_mods('[Acetyl]-P[3.14]EPTIDE-[Amide]', {'P': [['phospho']]}, 2)
        ['[Acetyl]-P[3.14]EP[phospho]TIDE-[Amide]', '[Acetyl]-P[3.14]EPTIDE-[Amide]']

        # Can also use regular expressions to match multiple residues:
        >>> apply_variable_mods('PEPTIDE', {'(?<=P)E': [['phospho']]}, 2)
        ['PE[phospho]PTIDE', 'PEPTIDE']

        >>> apply_variable_mods('<13C>PEPTIDE', {'P': [['phospho']]}, 1)
        ['<13C>P[phospho]EPTIDE', '<13C>PEP[phospho]TIDE', '<13C>PEPTIDE']

        >>> apply_variable_mods('(?PE)PTIDE', {'P': [['phospho']]}, 1)
        ['(?P[phospho]E)PTIDE', '(?PE)P[phospho]TIDE', '(?PE)PTIDE']

        >>> apply_variable_mods('PEPTIDE', {'P': [['1']]}, 1, nterm_mods={'P': [['acetyl']]}, cterm_mods={'P': [['amide']]})
        ['[acetyl]-P[1]EPTIDE', '[acetyl]-PEP[1]TIDE', '[acetyl]-PEPTIDE', 'P[1]EPTIDE', 'PEP[1]TIDE', 'PEPTIDE']

    """

    seqs = []
    if nterm_mods:
        for regex_str, mods_list in nterm_mods.items():
            for mods in mods_list:
                nterm_seq = apply_static_mods(sequence, {}, nterm_mods={regex_str: mods}, overwrite=overwrite)
                if nterm_seq != sequence:
                    seqs.extend(_variable_mods_builder(nterm_seq, mod_map, max_mods, overwrite))

    if cterm_mods:
        for regex_str, mods_list in cterm_mods.items():
            for mods in mods_list:
                cterm_seq = apply_static_mods(sequence, {}, cterm_mods={regex_str: mods}, overwrite=overwrite)
                if cterm_seq != sequence:
                    seqs.extend(_variable_mods_builder(cterm_seq, mod_map, max_mods, overwrite))

    seqs.extend(_variable_mods_builder(sequence, mod_map, max_mods, overwrite))

    return seqs
