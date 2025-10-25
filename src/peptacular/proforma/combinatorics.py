from __future__ import annotations

import itertools
from typing import TYPE_CHECKING, Any, Generator

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


def generate_permutations(
    annotation: ProFormaAnnotation, size: int | None = None
) -> Generator[ProFormaAnnotation, None, None]:
    """
    Generate all permutations of the annotation sequence.

    Args:
        annotation: The ProFormaAnnotation to permute
        size: Number of elements in each permutation (default: sequence length)

    Yields:
        ProFormaAnnotation instances representing each permutation

    Note:
        Non-internal modifications (terminal, labile, etc.) are preserved
        and applied to each permutation.
    """
    if size is None:
        size = len(annotation)

    # Extract mods that should be preserved across permutations
    annotation_copy = annotation.copy()
    outside_mods = _extract_non_internal_mods(annotation_copy)

    # Generate permutations of split amino acids
    split_aas = annotation_copy.split()

    for permutation in itertools.permutations(split_aas, size):
        # Create new annotation from permuted sequence
        combined_sequence = "".join(aa.serialize() for aa in permutation)
        result = annotation.__class__(sequence=combined_sequence)
        _apply_extracted_mods(result, outside_mods)
        yield result


def generate_product(
    annotation: ProFormaAnnotation, repeat: int | None = None
) -> Generator[ProFormaAnnotation, None, None]:
    """
    Generate the Cartesian product of the annotation sequence with itself.

    Args:
        annotation: The ProFormaAnnotation to generate products from
        repeat: Number of repetitions (default: sequence length)

    Yields:
        ProFormaAnnotation instances representing each product combination

    Note:
        Non-internal modifications are preserved and applied to each product.
    """
    if repeat is None:
        repeat = len(annotation)

    annotation_copy = annotation.copy()
    outside_mods = _extract_non_internal_mods(annotation_copy)

    split_aas = annotation_copy.split()

    for product in itertools.product(split_aas, repeat=repeat):
        combined_sequence = "".join(aa.serialize() for aa in product)
        result = annotation.__class__(sequence=combined_sequence)
        _apply_extracted_mods(result, outside_mods)
        yield result


def generate_combinations(
    annotation: ProFormaAnnotation, r: int | None = None
) -> Generator[ProFormaAnnotation, None, None]:
    """
    Generate all combinations of the annotation sequence.

    Args:
        annotation: The ProFormaAnnotation to generate combinations from
        r: Length of each combination (default: sequence length)

    Yields:
        ProFormaAnnotation instances representing each combination

    Note:
        Non-internal modifications are preserved and applied to each combination.
    """
    if r is None:
        r = len(annotation)

    annotation_copy = annotation.copy()
    outside_mods = _extract_non_internal_mods(annotation_copy)

    split_aas = annotation_copy.split()

    for combination in itertools.combinations(split_aas, r=r):
        combined_sequence = "".join(aa.serialize() for aa in combination)
        result = annotation.__class__(sequence=combined_sequence)
        _apply_extracted_mods(result, outside_mods)
        yield result


def generate_combinations_with_replacement(
    annotation: ProFormaAnnotation, r: int | None = None
) -> Generator[ProFormaAnnotation, None, None]:
    """
    Generate all combinations of the annotation sequence with replacement.

    Args:
        annotation: The ProFormaAnnotation to generate combinations from
        r: Length of each combination (default: sequence length)

    Yields:
        ProFormaAnnotation instances representing each combination with replacement

    Note:
        Non-internal modifications are preserved and applied to each combination.
    """
    if r is None:
        r = len(annotation)

    annotation_copy = annotation.copy()
    outside_mods = _extract_non_internal_mods(annotation_copy)

    split_aas = annotation_copy.split()

    for combination in itertools.combinations_with_replacement(split_aas, r=r):
        combined_sequence = "".join(aa.serialize() for aa in combination)
        result = annotation.__class__(sequence=combined_sequence)
        _apply_extracted_mods(result, outside_mods)
        yield result


# Helper functions


def _extract_non_internal_mods(annotation: ProFormaAnnotation) -> dict[str, Any]:
    """Extract non-internal modifications from annotation"""
    from ..constants import ModType

    non_internal_mod_types = [
        ModType.NTERM,
        ModType.CTERM,
        ModType.LABILE,
        ModType.ISOTOPE,
        ModType.STATIC,
        ModType.UNKNOWN,
        ModType.CHARGE,
        ModType.CHARGE_ADDUCTS,
    ]

    return annotation.pop_mods(mod_types=non_internal_mod_types)


def _apply_extracted_mods(
    annotation: ProFormaAnnotation, mods_dict: dict[str, Any]
) -> None:
    """Apply extracted modifications back to annotation"""
    for mod_type_str, mod_value in mods_dict.items():
        if mod_value is not None:
            # Convert string back to ModType if needed
            from ..constants import ModType

            mod_type = ModType(mod_type_str)
            annotation._set_mod_by_type(mod_value, mod_type)  # type: ignore
