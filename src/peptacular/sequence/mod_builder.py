"""
Refactored modification builder for amino acid sequences.
"""

from typing import Any, Dict, List, Generator, Union, Optional, Literal
from dataclasses import dataclass

from .util import get_annotation_input
from ..proforma.annot import ProFormaAnnotation
from ..dclasses import Mod
from ..util import get_regex_match_indices
from ..dclasses import (
    fix_list_of_list_of_mods,
    fix_list_of_mods,
    remove_empty_list_of_list_of_mods,
)

from ..dclasses import ACCEPTED_MOD_TYPES
from typing import Dict, List, Generator, Union, Optional, Literal, overload

MOD_MODES = Literal["skip", "append", "overwrite"]
MOD_BUILDER_RETURN_TYPES = Literal["str", "annotation"]

STATIC_MOD_INPUT = Union[List[ACCEPTED_MOD_TYPES], ACCEPTED_MOD_TYPES]
VAR_MOD_INPUT = Union[
    List[List[ACCEPTED_MOD_TYPES]], List[ACCEPTED_MOD_TYPES], ACCEPTED_MOD_TYPES
]


@dataclass
class ModificationSpec:
    """Encapsulates modification specifications for cleaner parameter handling."""

    internal_mods: Dict[str, Union[STATIC_MOD_INPUT, VAR_MOD_INPUT]]
    nterm_mods: Union[
        Dict[str, Union[STATIC_MOD_INPUT, VAR_MOD_INPUT]],
        Union[STATIC_MOD_INPUT, VAR_MOD_INPUT],
        None,
    ]
    cterm_mods: Union[
        Dict[str, Union[STATIC_MOD_INPUT, VAR_MOD_INPUT]],
        Union[STATIC_MOD_INPUT, VAR_MOD_INPUT],
        None,
    ]


class StaticModificationProcessor:
    """Handles static modification processing."""

    def normalize_mods(self, mods: STATIC_MOD_INPUT) -> List[Mod]:
        return fix_list_of_mods(mods)

    def apply_internal_mods(
        self,
        annotation: ProFormaAnnotation,
        mods: Dict[str, List[ACCEPTED_MOD_TYPES]],
        mode: MOD_MODES,
    ) -> ProFormaAnnotation:
        new_annotation = annotation.copy()

        for regex_str, mod_list in mods.items():
            for mod_index in get_regex_match_indices(
                annotation.sequence, regex_str, offset=-1
            ):
                self._apply_mod_at_index(new_annotation, mod_index, mod_list, mode)  # type: ignore

        return new_annotation  # type: ignore

    def _apply_mod_at_index(
        self,
        annotation: ProFormaAnnotation,
        index: int,
        mods: List[ACCEPTED_MOD_TYPES],
        mode: MOD_MODES,
    ):
        """Apply modification at specific index based on mode."""
        if annotation.has_internal_mods_at_index(index):
            if mode == "overwrite":
                annotation.add_internal_mods_at_index(index, mods, False)
            elif mode == "append":
                annotation.add_internal_mods_at_index(index, mods, True)
            elif mode == "skip":
                return
            else:
                raise ValueError(f"Invalid mode: {mode}")
        else:
            annotation.add_internal_mods_at_index(index, mods, True)


class VariableModificationProcessor:
    """Handles variable modification processing."""

    def normalize_mods(self, mods: VAR_MOD_INPUT) -> Optional[List[List[Mod]]]:
        normalized = fix_list_of_list_of_mods(mods)
        return remove_empty_list_of_list_of_mods(normalized)

    def apply_internal_mods(
        self,
        annotation: ProFormaAnnotation,
        mods: Dict[str, List[List[ACCEPTED_MOD_TYPES]]],
        mode: MOD_MODES,
    ) -> Generator[ProFormaAnnotation, None, None]:
        # Convert regex-based mods to index-based mods
        index_mods = self._build_index_mod_map(annotation, mods)

        # Apply variable modifications recursively
        starting_mod_count = annotation.count_modified_residues()
        yield from self._apply_variable_mods_recursive(
            index_mods,
            annotation,
            0,
            float("inf"),
            mode,  # max_mods handled at higher level
        )

    def _build_index_mod_map(
        self,
        annotation: ProFormaAnnotation,
        mods: Dict[str, List[List[ACCEPTED_MOD_TYPES]]],
    ) -> Dict[int, List[List[ACCEPTED_MOD_TYPES]]]:
        """Convert regex-based modifications to index-based modifications."""
        index_mods = {}
        for regex_str, list_of_list_of_mods in mods.items():
            for mod_index in get_regex_match_indices(
                annotation.sequence, regex_str, offset=-1
            ):
                for list_of_mods in list_of_list_of_mods:
                    index_mods.setdefault(mod_index, []).append(list_of_mods)
        return index_mods

    def _apply_variable_mods_recursive(
        self,
        mods: Dict[int, List[List[ACCEPTED_MOD_TYPES]]],
        annotation: ProFormaAnnotation,
        index: int,
        max_mod_count: int,
        mode: MOD_MODES,
    ) -> Generator[ProFormaAnnotation, None, None]:
        """Recursively apply variable modifications."""
        if (
            index == len(annotation.sequence)
            or annotation.count_modified_residues() == max_mod_count
        ):
            yield annotation
            return

        original_annotation = annotation.copy()
        curr_mods = mods.get(index)

        if curr_mods is not None:
            for mod_list in curr_mods:
                if annotation.has_internal_mods_at_index(index):
                    if mode == "skip":
                        continue

                    updated_annotation = annotation.copy()
                    append_mode = mode == "append"
                    updated_annotation.add_internal_mods_at_index(
                        index, mod_list, append_mode
                    )

                    yield from self._apply_variable_mods_recursive(
                        mods, updated_annotation, index + 1, max_mod_count, mode
                    )
                else:
                    updated_annotation = annotation.copy()
                    updated_annotation.add_internal_mods_at_index(index, mod_list, True)
                    yield from self._apply_variable_mods_recursive(
                        mods, updated_annotation, index + 1, max_mod_count, mode
                    )

        yield from self._apply_variable_mods_recursive(
            mods, original_annotation, index + 1, max_mod_count, mode
        )


class TerminalModificationHandler:
    """Handles N-terminal and C-terminal modifications."""

    def __init__(
        self,
        processor: Union[StaticModificationProcessor, VariableModificationProcessor],
    ):
        self.processor = processor

    def normalize_terminal_mods(self, mods) -> Dict[str, Any]:
        """Normalize terminal modifications to consistent format."""
        if isinstance(mods, dict) and len(mods) > 0:
            return {k: self.processor.normalize_mods(v) for k, v in mods.items() if v}
        elif isinstance(mods, list) and len(mods) > 0:
            return {"": self.processor.normalize_mods(mods)}
        elif isinstance(mods, (float, int, str, Mod)):
            return {"": self.processor.normalize_mods(mods)}
        else:
            return {}

    def apply_nterm_mods(
        self, annotation: ProFormaAnnotation, mods: Dict[str, Any], mode: MOD_MODES
    ) -> ProFormaAnnotation:
        """Apply N-terminal modifications."""
        new_annotation = annotation.copy()

        for regex_str, mod_list in mods.items():
            for mod_index in get_regex_match_indices(
                annotation.sequence, regex_str, offset=-1
            ):
                if mod_index == 0:  # N-terminal
                    self._apply_terminal_mod(new_annotation, mod_list, mode, "nterm")

        return new_annotation

    def apply_cterm_mods(
        self, annotation: ProFormaAnnotation, mods: Dict[str, Any], mode: MOD_MODES
    ) -> ProFormaAnnotation:
        """Apply C-terminal modifications."""
        new_annotation = annotation.copy()

        for regex_str, mod_list in mods.items():
            for mod_index in get_regex_match_indices(
                annotation.sequence, regex_str, offset=-1
            ):
                if mod_index == len(annotation.sequence) - 1:  # C-terminal
                    self._apply_terminal_mod(new_annotation, mod_list, mode, "cterm")

        return new_annotation

    def _apply_terminal_mod(
        self, annotation: ProFormaAnnotation, mods: Any, mode: MOD_MODES, term_type: str
    ):
        """Apply terminal modification based on mode."""
        has_existing = (
            annotation.has_nterm_mods()
            if term_type == "nterm"
            else annotation.has_cterm_mods()
        )
        add_method = (
            annotation.add_nterm_mods
            if term_type == "nterm"
            else annotation.add_cterm_mods
        )

        if has_existing:
            if mode == "overwrite":
                add_method(mods, False)
            elif mode == "append":
                add_method(mods, True)
            elif mode == "skip":
                return
            else:
                raise ValueError(f"Invalid mode: {mode}")
        else:
            add_method(mods, True)


# Simplified main functions using the new architecture


@overload
def apply_static_mods(
    sequence: Union[str, ProFormaAnnotation],
    internal_mods: Optional[Dict[str, STATIC_MOD_INPUT]],
    nterm_mods: Optional[Union[STATIC_MOD_INPUT, Dict[str, STATIC_MOD_INPUT]]] = None,
    cterm_mods: Optional[Union[STATIC_MOD_INPUT, Dict[str, STATIC_MOD_INPUT]]] = None,
    mode: MOD_MODES = "skip",
    return_type: Literal["str"] = "str",
) -> str: ...


@overload
def apply_static_mods(
    sequence: Union[str, ProFormaAnnotation],
    internal_mods: Optional[Dict[str, STATIC_MOD_INPUT]],
    nterm_mods: Optional[Union[STATIC_MOD_INPUT, Dict[str, STATIC_MOD_INPUT]]] = None,
    cterm_mods: Optional[Union[STATIC_MOD_INPUT, Dict[str, STATIC_MOD_INPUT]]] = None,
    mode: MOD_MODES = "skip",
    return_type: Literal["annotation"] = "annotation",
) -> ProFormaAnnotation: ...


def apply_static_mods(
    sequence: Union[str, ProFormaAnnotation],
    internal_mods: Optional[Dict[str, STATIC_MOD_INPUT]],
    nterm_mods: Optional[Union[STATIC_MOD_INPUT, Dict[str, STATIC_MOD_INPUT]]] = None,
    cterm_mods: Optional[Union[STATIC_MOD_INPUT, Dict[str, STATIC_MOD_INPUT]]] = None,
    mode: MOD_MODES = "skip",
    return_type: MOD_BUILDER_RETURN_TYPES = "str",
) -> Union[str, ProFormaAnnotation]:
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

    annotation = get_annotation_input(sequence, copy=True)
    processor = StaticModificationProcessor()
    terminal_handler = TerminalModificationHandler(processor)

    # Normalize and filter modifications
    normalized_internal = {}
    if internal_mods:
        normalized_internal = {
            k: processor.normalize_mods(v) for k, v in internal_mods.items() if v
        }

    normalized_nterm = terminal_handler.normalize_terminal_mods(nterm_mods)
    normalized_cterm = terminal_handler.normalize_terminal_mods(cterm_mods)

    # Apply modifications
    result = annotation
    if normalized_internal:
        result = processor.apply_internal_mods(result, normalized_internal, mode)
    if normalized_nterm:
        result = terminal_handler.apply_nterm_mods(result, normalized_nterm, mode)
    if normalized_cterm:
        result = terminal_handler.apply_cterm_mods(result, normalized_cterm, mode)

    return result.serialize() if return_type == "str" else result


def apply_variable_mods(
    sequence: Union[str, ProFormaAnnotation],
    internal_mods: Optional[Dict[str, VAR_MOD_INPUT]],
    max_mods: int,
    nterm_mods: Optional[Union[Dict[str, VAR_MOD_INPUT], VAR_MOD_INPUT]] = None,
    cterm_mods: Optional[Union[Dict[str, VAR_MOD_INPUT], VAR_MOD_INPUT]] = None,
    mode: MOD_MODES = "skip",
    return_type: MOD_BUILDER_RETURN_TYPES = "str",
) -> Union[List[str], List[ProFormaAnnotation]]:
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
        >>> output = set(apply_variable_mods('PEPTIDE', {'P': [['phospho']]}, 1))
        >>> expected = {'P[phospho]EPTIDE', 'PEP[phospho]TIDE', 'PEPTIDE'}
        >>> output == expected
        True

        # Dict: Single list
        >>> output = set(apply_variable_mods('PEPTIDE', {'P': ['phospho']}, 1))
        >>> expected = {'P[phospho]EPTIDE', 'PEP[phospho]TIDE', 'PEPTIDE'}
        >>> output == expected
        True

        # Dict: Single Value
        >>> output = set(apply_variable_mods('PEPTIDE', {'P': 'phospho'}, 1))
        >>> expected = {'P[phospho]EPTIDE', 'PEP[phospho]TIDE', 'PEPTIDE'}
        >>> output == expected
        True

        # Terms can be Dict or List
        >>> output = set(apply_variable_mods('PEPTIDE', None, 1, nterm_mods={'': 'phospho'}, cterm_mods='acetyl'))
        >>> expected = {'[phospho]-PEPTIDE', '[phospho]-PEPTIDE-[acetyl]', 'PEPTIDE-[acetyl]', 'PEPTIDE'}
        >>> output == expected
        True

        # When multiple mods are specified in the same group, they are applied together
        >>> output = set(apply_variable_mods('PEPTIDE', {'P': [['phospho', 1]]}, 1))
        >>> expected = {'P[phospho][1]EPTIDE', 'PEP[phospho][1]TIDE', 'PEPTIDE'}
        >>> output == expected
        True

        # Can also break up groups to apply modifications separately
        >>> output = set(apply_variable_mods('PEPTIDE', {'P': [['phospho'], [1]]}, 1))
        >>> expected = {'P[phospho]EPTIDE', 'P[1]EPTIDE', 'PEP[phospho]TIDE', 'PEP[1]TIDE', 'PEPTIDE'}
        >>> output == expected
        True

        # Setting max_mods to 0 will return the original sequence:
        >>> output = set(apply_variable_mods('PEPTIDE', {'P': [['phospho']]}, 0))
        >>> expected = {'PEPTIDE'}
        >>> output == expected
        True

        # Setting max_mods to 0 will return the original sequence:
        >>> output = set(apply_variable_mods('PEPTIDE[Hello]', {'P': [['phospho']]}, 1))
        >>> expected = {'P[phospho]EPTIDE[Hello]', 'PEP[phospho]TIDE[Hello]', 'PEPTIDE[Hello]'}
        >>> output == expected
        True

        # Works with already modified sequences, and will not replace existing modifications:
        >>> output = set(apply_variable_mods('[Acetyl]-P[3.14]EPTIDE-[Amide]', {'P': [['phospho']]}, 2))
        >>> expected = {'[Acetyl]-P[3.14]EP[phospho]TIDE-[Amide]', '[Acetyl]-P[3.14]EPTIDE-[Amide]'}
        >>> output == expected
        True

        # Can also use regular expressions to match multiple residues:
        >>> output = set(apply_variable_mods('PEPTIDE', {'(?<=P)E': [['phospho']]}, 2))
        >>> expected = {'PE[phospho]PTIDE', 'PEPTIDE'}
        >>> output == expected
        True

        >>> output = set(apply_variable_mods('<13C>PEPTIDE', {'P': [['phospho']]}, 1))
        >>> expected = {'<13C>P[phospho]EPTIDE', '<13C>PEP[phospho]TIDE', '<13C>PEPTIDE'}
        >>> output == expected
        True

        >>> output = set(apply_variable_mods('(?PE)PTIDE', {'P': [['phospho']]}, 1))
        >>> expected = {'(?P[phospho]E)PTIDE', '(?PE)P[phospho]TIDE', '(?PE)PTIDE'}
        >>> output == expected
        True

        # Term mods can also have a regex to match terminus residues
        >>> output = set(apply_variable_mods('PEPTIDE', {'P': [['1']]}, 1, nterm_mods={'P': [['A']]}, cterm_mods={'P': [['B']]}))
        >>> expected = {'[A]-P[1]EPTIDE', '[A]-PEP[1]TIDE', '[A]-PEPTIDE', 'P[1]EPTIDE', 'PEP[1]TIDE', 'PEPTIDE'}
        >>> output == expected
        True

        # Test for empty mods
        >>> output = set(apply_variable_mods('PEPTIDE', None, 1, nterm_mods=[]))
        >>> expected = {'PEPTIDE'}
        >>> output == expected
        True
        >>> output = set(apply_variable_mods('PEPTIDE', None, 1, nterm_mods={'P': [[]]}))
        >>> expected = {'PEPTIDE'}
        >>> output == expected
        True

    """

    annotation = get_annotation_input(sequence, copy=True)
    processor = VariableModificationProcessor()
    terminal_handler = TerminalModificationHandler(processor)

    # Normalize modifications
    normalized_internal = {}
    if internal_mods:
        normalized_internal = {
            k: processor.normalize_mods(v) for k, v in internal_mods.items() if v
        }

    normalized_nterm = terminal_handler.normalize_terminal_mods(nterm_mods)
    normalized_cterm = terminal_handler.normalize_terminal_mods(cterm_mods)

    # Generate all combinations
    all_annotations = set()

    # Base case: no terminal modifications
    base_variants = list(
        processor.apply_internal_mods(annotation, normalized_internal, mode)
    )
    all_annotations.update(base_variants)

    # Add N-terminal modifications
    if normalized_nterm:
        nterm_variants = []
        for regex_str, mods_list in normalized_nterm.items():
            if (
                mods_list
                and hasattr(mods_list[0], "__iter__")
                and not isinstance(mods_list[0], (str, Mod))
            ):
                # Variable mods format
                for mods in mods_list:
                    nterm_annot = apply_static_mods(
                        annotation,
                        {},
                        nterm_mods={regex_str: mods},
                        mode=mode,
                        return_type="annotation",
                    )
                    if nterm_annot != annotation:
                        nterm_variants.extend(
                            processor.apply_internal_mods(
                                nterm_annot, normalized_internal, mode
                            )
                        )
        all_annotations.update(nterm_variants)

    # Add C-terminal modifications (with and without N-terminal)
    if normalized_cterm:
        cterm_variants = []
        for regex_str, mods_list in normalized_cterm.items():
            if (
                mods_list
                and hasattr(mods_list[0], "__iter__")
                and not isinstance(mods_list[0], (str, Mod))
            ):
                # Variable mods format
                for mods in mods_list:
                    # C-term only
                    cterm_annot = apply_static_mods(
                        annotation,
                        {},
                        cterm_mods={regex_str: mods},
                        mode=mode,
                        return_type="annotation",
                    )
                    if cterm_annot != annotation:
                        cterm_variants.extend(
                            processor.apply_internal_mods(
                                cterm_annot, normalized_internal, mode
                            )
                        )

                    # C-term + N-term combinations
                    for nterm_variant in nterm_variants:
                        combined_annot = apply_static_mods(
                            nterm_variant,
                            {},
                            cterm_mods={regex_str: mods},
                            mode=mode,
                            return_type="annotation",
                        )
                        if combined_annot != nterm_variant:
                            cterm_variants.append(combined_annot)

        all_annotations.update(cterm_variants)

    # Filter by max_mods constraint
    starting_mod_count = annotation.count_modified_residues()
    filtered_annotations = [
        annot
        for annot in all_annotations
        if annot.count_modified_residues() <= max_mods + starting_mod_count
    ]

    if return_type == "str":
        return [annot.serialize() for annot in filtered_annotations]

    return filtered_annotations
