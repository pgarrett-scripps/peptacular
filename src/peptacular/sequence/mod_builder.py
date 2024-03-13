from __future__ import annotations

from typing import Dict, List, Generator

from peptacular.sequence.sequence import parse_single_sequence
from peptacular.sequence.proforma import ProFormaAnnotation, Mod
from peptacular.types import ModIndex, ModValue
from peptacular.util import get_regex_match_indices
from peptacular.input_parser import fix_list_of_list_of_mods, fix_list_of_mods


def apply_static_mods(sequence: str | ProFormaAnnotation,
                      residue_mods: Dict[str, List[ModValue] | ModValue],
                      nterm_mods: Dict[str, List[ModValue] | ModValue] | None = None,
                      cterm_mods: Dict[str, List[ModValue] | ModValue] | None = None,
                      mode: str = 'skip') -> str | ProFormaAnnotation:
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

        >>> apply_static_mods('PEPTIDE', {'P': ['phospho', 'acetyl']})
        'P[phospho][acetyl]EP[phospho][acetyl]TIDE'

        # Can specify multiple static modifications:
        >>> apply_static_mods('PEPTIDE', {'P': ['phospho'], 'E': [3.0]})
        'P[phospho]E[3.0]P[phospho]TIDE[3.0]'

        # Works with already modified sequences:
        >>> apply_static_mods('[3.14]-PEPTIDE-[1.234]', {'P': ['phospho']})
        '[3.14]-P[phospho]EP[phospho]TIDE-[1.234]'

        # By default, any overlapping modifications will be replaced:
        >>> apply_static_mods('PEP[1.234]TIDE', {'P': ['phospho']}, mode='overwrite')
        'P[phospho]EP[phospho]TIDE'

        # To append new modifications to existing ones, use mode='append':
        >>> apply_static_mods('PEP[1.234]TIDE', {'P': ['phospho']}, mode='append')
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

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
        input_type = str
    else:
        annotation = sequence
        input_type = ProFormaAnnotation

    residue_mods = {k: fix_list_of_mods(v) for k, v in residue_mods.items()}
    nterm_mods = {k: fix_list_of_mods(v) for k, v in nterm_mods.items()} if nterm_mods else {}
    cterm_mods = {k: fix_list_of_mods(v) for k, v in cterm_mods.items()} if cterm_mods else {}

    new_annotation = annotation.copy()

    for regex_str, mods in residue_mods.items():
        for mod_index in get_regex_match_indices(annotation.sequence, regex_str):
            if annotation.has_internal_mod(mod_index):  # mod already present
                if mode == 'overwrite':
                    new_annotation.add_internal_mod(mod_index, mods, False)
                elif mode == 'append':
                    new_annotation.add_internal_mod(mod_index, mods, True)
                elif mode == 'skip':
                    continue
                else:
                    raise ValueError(f'Invalid mode: {mode}')
            else:
                new_annotation.add_internal_mod(mod_index, mods, True)

    for regex_str, mods in nterm_mods.items():
        for mod_index in get_regex_match_indices(annotation.sequence, regex_str):
            if mod_index == 0:
                if annotation.has_nterm_mods():
                    if mode == 'overwrite':
                        new_annotation.add_nterm_mods(mods, False)
                    elif mode == 'append':
                        new_annotation.add_nterm_mods(mods, True)
                    elif mode == 'skip':
                        continue
                    else:
                        raise ValueError(f'Invalid mode: {mode}')
                else:
                    new_annotation.add_nterm_mods(mods, True)

    for regex_str, mods in cterm_mods.items():
        for mod_index in get_regex_match_indices(annotation.sequence, regex_str):
            if mod_index == len(annotation.sequence) - 1:
                if annotation.has_cterm_mods():
                    if mode == 'overwrite':
                        new_annotation.add_cterm_mods(mods, False)
                    elif mode == 'append':
                        new_annotation.add_cterm_mods(mods, True)
                    elif mode == 'skip':
                        continue
                    else:
                        raise ValueError(f'Invalid mode: {mode}')
                else:
                    new_annotation.add_cterm_mods(mods, True)

    if input_type == str:
        return new_annotation.serialize()

    return new_annotation


def _apply_variable_mods_rec(mods: Dict[ModIndex, List[List[Mod]]],
                             annotation: ProFormaAnnotation,
                             index: int,
                             max_mod_count: int,
                             mode: str) -> Generator[ProFormaAnnotation, None, None]:
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

    if index == len(annotation.sequence) or annotation.count_modified_residues() == max_mod_count:
        yield annotation
        return

    original_annotation = annotation.copy()

    curr_list_of_mods = mods.get(index, None)

    if curr_list_of_mods is not None:
        for curr_mods in curr_list_of_mods:
            if annotation.has_internal_mod(index):
                if mode == 'overwrite':
                    updated_annotation = annotation.copy()
                    updated_annotation.add_internal_mod(index, curr_mods, False)  # will overwrite
                    yield from _apply_variable_mods_rec(mods, updated_annotation, index + 1, max_mod_count, mode)
                elif mode == 'append':
                    updated_annotation = annotation.copy()
                    updated_annotation.add_internal_mod(index, curr_mods, True)  # will append
                    yield from _apply_variable_mods_rec(mods, updated_annotation, index + 1, max_mod_count, mode)
                elif mode == 'skip':
                    continue
                else:
                    raise ValueError(f'Invalid mode: {mode}')

            else:
                updated_annotation = annotation.copy()
                updated_annotation.add_internal_mod(index, curr_mods, True)
                yield from _apply_variable_mods_rec(mods, updated_annotation, index + 1, max_mod_count, mode)

    yield from _apply_variable_mods_rec(mods, original_annotation, index + 1, max_mod_count, mode)


def _variable_mods_builder(annotation: ProFormaAnnotation,
                           mod_map: Dict[str, List[List[Mod]]],
                           max_mods: int,
                           mode: str) -> Generator[ProFormaAnnotation, None, None]:
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

    new_mod_map: Dict[int, List[List[Mod]]] = {}
    for regex_str, list_of_list_of_mods in mod_map.items():
        for mod_index in get_regex_match_indices(annotation.sequence, regex_str):
            for list_of_mods in list_of_list_of_mods:
                new_mod_map.setdefault(mod_index, []).append(list_of_mods)

    starting_mod_count = annotation.count_modified_residues()
    annotations = _apply_variable_mods_rec(new_mod_map, annotation, 0, max_mods + starting_mod_count, mode)
    return annotations


def apply_variable_mods(sequence: str | ProFormaAnnotation,
                        mod_map: Dict[str, List[List[ModValue]] | List[ModValue] | ModValue],
                        max_mods: int,
                        nterm_mods: Dict[str, List[List[ModValue]] | List[ModValue] | ModValue] | None = None,
                        cterm_mods: Dict[str, List[List[ModValue]] | List[ModValue] | ModValue] | None = None,
                        mode: str = 'skip') -> List[str] | List[ProFormaAnnotation]:
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
    :param mode: overwrite, append, skip
    :type mode: str


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

        # Setting max_mods to 0 will return the original sequence:
        >>> apply_variable_mods('PEPTIDE[Hello]', {'P': [['phospho']]}, 1)
        ['P[phospho]EPTIDE[Hello]', 'PEP[phospho]TIDE[Hello]', 'PEPTIDE[Hello]']

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

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
        input_type = str
    else:
        annotation = sequence
        input_type = ProFormaAnnotation

    mod_map = {k: fix_list_of_list_of_mods(v) for k, v in mod_map.items()}
    nterm_mods = {k: fix_list_of_list_of_mods(v) for k, v in nterm_mods.items()} if nterm_mods else {}
    cterm_mods = {k: fix_list_of_list_of_mods(v) for k, v in cterm_mods.items()} if cterm_mods else {}

    annotations = []
    if nterm_mods:
        for regex_str, mods_list in nterm_mods.items():
            for mods in mods_list:
                nterm_annot = apply_static_mods(annotation, {}, nterm_mods={regex_str: mods}, mode=mode)
                if nterm_annot != annotation:
                    annotations.extend(_variable_mods_builder(nterm_annot, mod_map, max_mods, mode=mode))

    if cterm_mods:
        for regex_str, mods_list in cterm_mods.items():
            for mods in mods_list:
                cterm_annot = apply_static_mods(annotation, {}, cterm_mods={regex_str: mods}, mode=mode)
                if cterm_annot != annotation:
                    annotations.extend(_variable_mods_builder(cterm_annot, mod_map, max_mods, mode=mode))

    annotations.extend(_variable_mods_builder(annotation, mod_map, max_mods, mode=mode))

    if input_type == str:
        return [a.serialize() for a in annotations]

    return annotations
