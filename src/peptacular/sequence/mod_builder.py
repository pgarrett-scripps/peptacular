"""
Modification builder for amino acid sequences.
"""

from typing import Dict, List, Generator, Union, Optional, Literal

from peptacular.sequence.sequence_funcs import sequence_to_annotation
from peptacular.proforma.proforma_parser import ProFormaAnnotation
from peptacular.proforma.proforma_dataclasses import Mod
from peptacular.util import get_regex_match_indices
from peptacular.proforma.input_convert import fix_list_of_list_of_mods, fix_list_of_mods, ModIndex, ModValue, \
    remove_empty_list_of_list_of_mods

ModMode = Literal['skip', 'append', 'overwrite']
MOD_MODE_VALUES = ['skip', 'append', 'overwrite']

ModBuilderReturnType = Literal["str", "annotation"]
MOD_BUILDER_RETURN_TYPE_VALUES = ['str', 'annotation']

STATIC_MOD_INPUT = Union[List[ModValue], ModValue]
VAR_MOD_INPUT = Union[List[List[ModValue]], List[ModValue], ModValue]


def apply_static_mods(sequence: Union[str, ProFormaAnnotation],
                      internal_mods: Optional[Dict[str, STATIC_MOD_INPUT]],
                      nterm_mods: Optional[Union[STATIC_MOD_INPUT, Dict[str, STATIC_MOD_INPUT]]] = None,
                      cterm_mods: Optional[Union[STATIC_MOD_INPUT, Dict[str, STATIC_MOD_INPUT]]] = None,
                      mode: ModMode = 'skip',
                      return_type: ModBuilderReturnType = 'str') -> Union[str, ProFormaAnnotation]:
    """
    Add static modifications to an amino acid sequence. If a modification is already present in the sequence,
    it will be replaced by the new one.

    :param sequence: Original amino acid sequence.
    :type sequence: str | ProFormaAnnotation
    :param internal_mods: Dictionary mapping amino acids to the mass of their modifications.
    :type internal_mods: Dict[str, ModValue | List[ModValue]]
    :param nterm_mods: Dictionary mapping the N-terminal amino acid to the mass of its modification.
    :type nterm_mods: Dict[str, ModValue | List[ModValue]] | List[ModValue] | ModValue | None
    :param cterm_mods: Dictionary mapping the C-terminal amino acid to the mass of its modification.
    :type cterm_mods: Dict[str, ModValue | List[ModValue]] | List[ModValue] | ModValue | None
    :param mode: overwrite, append, skip
    :type mode: str
    :param return_type: Type of the return value. Either 'str' or 'annotation'.
    :type return_type: ModBuilderReturnType

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

        # Term mods, macth must be either start or end index for N and C term
        >>> apply_static_mods('PEPTIDE', None, nterm_mods={'P': ['acetyl']}, cterm_mods={'P': ['amide']})
        '[acetyl]-PEPTIDE'

        >>> apply_static_mods('PEPTIDE', None, nterm_mods={'PE': 'acetyl'}, cterm_mods={'P': 'amide'})
        '[acetyl]-PEPTIDE'

        >>> apply_static_mods('PEPTIDE', None, nterm_mods={'PP': ['acetyl']}, cterm_mods={'(?<=D)E': ['amide']})
        'PEPTIDE-[amide]'

        # Can also specify a terminla mod to apply without any regex info
        >>> apply_static_mods('PEPTIDE', None, nterm_mods='acetyl', cterm_mods='amide')
        '[acetyl]-PEPTIDE-[amide]'

        # Test empty mods
        >>> apply_static_mods('PEPTIDE', None, nterm_mods=[], cterm_mods='amide')
        'PEPTIDE-[amide]'
        >>> apply_static_mods('PEPTIDE', None, nterm_mods={'P': []}, cterm_mods='amide')
        'PEPTIDE-[amide]'

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    if internal_mods is not None:
        internal_mods = {k: fix_list_of_mods(v) for k, v in internal_mods.items()}
    else:
        internal_mods = {}

    internal_mods = {k: v for k, v in internal_mods.items() if v}

    if isinstance(nterm_mods, Dict) and len(nterm_mods) > 0:
        nterm_mods = {k: fix_list_of_mods(v) for k, v in nterm_mods.items()} if nterm_mods else {}
    elif isinstance(nterm_mods, list) and len(nterm_mods) > 0:
        nterm_mods = {'': fix_list_of_mods(nterm_mods)}
    elif isinstance(nterm_mods, (float, int, str, Mod)):
        nterm_mods = {'': fix_list_of_mods(nterm_mods)}
    else:
        nterm_mods = {}

    nterm_mods = {k: v for k, v in nterm_mods.items() if v}

    if isinstance(cterm_mods, Dict) and len(cterm_mods) > 0:
        cterm_mods = {k: fix_list_of_mods(v) for k, v in cterm_mods.items()} if cterm_mods else {}
    elif isinstance(cterm_mods, list) and len(cterm_mods) > 0:
        cterm_mods = {'': fix_list_of_mods(cterm_mods)}
    elif isinstance(cterm_mods, (float, int, str, Mod)):
        cterm_mods = {'': fix_list_of_mods(cterm_mods)}
    else:
        cterm_mods = {}

    cterm_mods = {k: v for k, v in cterm_mods.items() if v}

    new_annotation = annotation.copy()

    for regex_str, mods in internal_mods.items():
        for mod_index in get_regex_match_indices(annotation.sequence, regex_str, offset=-1):
            if annotation.has_internal_mods_at_index(mod_index):  # mod already present
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
        for mod_index in get_regex_match_indices(annotation.sequence, regex_str, offset=-1):
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
        for mod_index in get_regex_match_indices(annotation.sequence, regex_str, offset=-1):
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

    if return_type == 'str':
        return new_annotation.serialize()

    return new_annotation


def _apply_variable_mods_rec(mods: Dict[ModIndex, List[List[Mod]]],
                             annotation: ProFormaAnnotation,
                             index: int,
                             max_mod_count: int,
                             mode: ModMode) -> Generator[ProFormaAnnotation, None, None]:
    """
    Recursively applies variable modifications to the amino acid sequence.

    :param mods: Dictionary mapping amino acids to the mass of their modifications.
    :type mods: Dict[int, Any]
    :param annotation: The annotation to modify.
    :type annotation: ProFormaAnnotation
    :param index: Current index on the amino acid sequence.
    :type index: int
    :param max_mod_count: Maximum number of modifications allowed on the amino acid sequence.
    :type max_mod_count: int
    :param mode: overwrite, append, skip
    :type mode: str

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
            if annotation.has_internal_mods_at_index(index):
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
                           mode: ModMode) -> Generator[ProFormaAnnotation, None, None]:
    """
    Apply variable modifications to a sequence.

    :param annotation: The annotation to modify.
    :type annotation: ProFormaAnnotation
    :param mod_map: Dictionary mapping amino acids / regex to a list of modifications
    :type mod_map: Dict[str, List[List[ModValue]]]
    :param max_mods: Maximum number of modifications allowed on the amino acid sequence.
    :type max_mods: int
    :param mode: overwrite, append, skip
    :type mode: str

    :return: List of all possible modified amino acid sequences.
    :rtype: List[str]
    """

    new_mod_map: Dict[int, List[List[Mod]]] = {}
    for regex_str, list_of_list_of_mods in mod_map.items():
        for mod_index in get_regex_match_indices(annotation.sequence, regex_str, offset=-1):
            for list_of_mods in list_of_list_of_mods:
                new_mod_map.setdefault(mod_index, []).append(list_of_mods)

    starting_mod_count = annotation.count_modified_residues()
    var_annotations = _apply_variable_mods_rec(new_mod_map, annotation, 0, max_mods + starting_mod_count, mode)
    return var_annotations


def apply_variable_mods(sequence: Union[str, ProFormaAnnotation],
                        internal_mods: Optional[Dict[str, VAR_MOD_INPUT]],
                        max_mods: int,
                        nterm_mods: Optional[Union[Dict[str, VAR_MOD_INPUT], VAR_MOD_INPUT]] = None,
                        cterm_mods: Optional[Union[Dict[str, VAR_MOD_INPUT], VAR_MOD_INPUT]] = None,
                        mode: ModMode = 'skip',
                        return_type: ModBuilderReturnType = 'str') -> Union[List[str], List[ProFormaAnnotation]]:
    """
    Apply variable modifications to a sequence.

    :param sequence: Original amino acid sequence.
    :type sequence: str | ProFormaAnnotation
    :param internal_mods: Dictionary mapping amino acids  /regex to modifications.
    :type internal_mods: Dict[str, List[List[ModValue]] | List[ModValue] | ModValue] | None
    :param max_mods: Maximum number of modifications allowed on the peptide.
    :type max_mods: int
    :param nterm_mods: Dictionary mapping the N-terminal amino acids / regex to modifications.
    :type nterm_mods: Dict[str, List[List[ModValue]] | List[ModValue] | ModValue] | List[List[ModValue]] |
    List[ModValue] | ModValue | None
    :param cterm_mods: Dictionary mapping the C-terminal amino acids / regex to modifications.
    :type cterm_mods: Dict[str, List[List[ModValue]] | List[ModValue] | ModValue] | List[List[ModValue]] |
    List[ModValue] | ModValue | None
    :param mode: overwrite, append, skip
    :type mode: str
    :param return_type: Type of the return value. Either 'str' or 'annotation'.
    :type return_type: ModBuilderReturnType

    :return: List of all possible modified peptide sequences.
    :rtype: List[str]

    .. code-block:: python

        # Dict: Double list
        >>> apply_variable_mods('PEPTIDE', {'P': [['phospho']]}, 1)
        ['P[phospho]EPTIDE', 'PEP[phospho]TIDE', 'PEPTIDE']

        # Dict: Single list
        >>> apply_variable_mods('PEPTIDE', {'P': ['phospho']}, 1)
        ['P[phospho]EPTIDE', 'PEP[phospho]TIDE', 'PEPTIDE']

        # Dict: Single Value
        >>> apply_variable_mods('PEPTIDE', {'P': 'phospho'}, 1)
        ['P[phospho]EPTIDE', 'PEP[phospho]TIDE', 'PEPTIDE']

        # Terms can be Dict or List
        >>> apply_variable_mods('PEPTIDE', None, 1, nterm_mods={'': 'phospho'}, cterm_mods='acetyl')
        ['[phospho]-PEPTIDE', '[phospho]-PEPTIDE-[acetyl]', 'PEPTIDE-[acetyl]', 'PEPTIDE']

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

        # Term mods can also have a regex to match terminus residues
        >>> apply_variable_mods('PEPTIDE', {'P': [['1']]}, 1, nterm_mods={'P': [['A']]}, cterm_mods={'P': [['B']]})
        ['[A]-P[1]EPTIDE', '[A]-PEP[1]TIDE', '[A]-PEPTIDE', 'P[1]EPTIDE', 'PEP[1]TIDE', 'PEPTIDE']

        # Test for empty mods
        >>> apply_variable_mods('PEPTIDE', None, 1, nterm_mods=[])
        ['PEPTIDE']
        >>> apply_variable_mods('PEPTIDE', None, 1, nterm_mods={'P': [[]]})
        ['PEPTIDE']

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    if internal_mods is not None:
        internal_mods = {k: fix_list_of_list_of_mods(v) for k, v in internal_mods.items()}
    else:
        internal_mods = {}

    internal_mods = {k: remove_empty_list_of_list_of_mods(v) for k, v in internal_mods.items()}
    internal_mods = {k: v for k, v in internal_mods.items() if v}

    if isinstance(nterm_mods, Dict) and len(nterm_mods) > 0:
        nterm_mods = {k: fix_list_of_list_of_mods(v) for k, v in nterm_mods.items()} if nterm_mods else {}
    elif isinstance(nterm_mods, list) and len(nterm_mods) > 0:
        nterm_mods = {'': fix_list_of_list_of_mods(nterm_mods)}
    elif isinstance(nterm_mods, (float, int, str, Mod)):
        nterm_mods = {'': fix_list_of_list_of_mods(nterm_mods)}
    else:
        nterm_mods = {}

    nterm_mods = {k: remove_empty_list_of_list_of_mods(v) for k, v in nterm_mods.items()}
    nterm_mods = {k: v for k, v in nterm_mods.items() if v}

    if isinstance(cterm_mods, Dict) and len(cterm_mods) > 0:
        cterm_mods = {k: fix_list_of_list_of_mods(v) for k, v in cterm_mods.items()} if cterm_mods else {}
    elif isinstance(cterm_mods, list) and len(cterm_mods) > 0:
        cterm_mods = {'': fix_list_of_list_of_mods(cterm_mods)}
    elif isinstance(cterm_mods, (float, int, str, Mod)):
        cterm_mods = {'': fix_list_of_list_of_mods(cterm_mods)}
    else:
        cterm_mods = {}

    cterm_mods = {k: remove_empty_list_of_list_of_mods(v) for k, v in cterm_mods.items()}
    cterm_mods = {k: v for k, v in cterm_mods.items() if v}

    n_term_annotations = []
    if nterm_mods:
        for regex_str, mods_list in nterm_mods.items():
            for mods in mods_list:
                nterm_annot = apply_static_mods(annotation, {}, nterm_mods={regex_str: mods},
                                                mode=mode, return_type='annotation')
                if nterm_annot != annotation:
                    n_term_annotations.extend(_variable_mods_builder(nterm_annot, internal_mods, max_mods, mode=mode))

    var_annotations = []
    if cterm_mods:
        for regex_str, mods_list in cterm_mods.items():
            for mods in mods_list:

                for n_term_annot in n_term_annotations:
                    cterm_annot = apply_static_mods(n_term_annot, {}, cterm_mods={regex_str: mods},
                                                    mode=mode, return_type='annotation')
                    if cterm_annot != n_term_annot:
                        var_annotations.extend(_variable_mods_builder(cterm_annot, internal_mods, max_mods, mode=mode))

                cterm_annot = apply_static_mods(annotation, {}, cterm_mods={regex_str: mods},
                                                mode=mode, return_type='annotation')
                if cterm_annot != annotation:
                    var_annotations.extend(_variable_mods_builder(cterm_annot, internal_mods, max_mods, mode=mode))

    var_annotations.extend(_variable_mods_builder(annotation, internal_mods, max_mods, mode=mode))
    var_annotations = n_term_annotations + var_annotations

    if return_type == 'str':
        return [a.serialize() for a in var_annotations]

    return var_annotations
