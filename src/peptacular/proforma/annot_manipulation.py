from __future__ import annotations
from collections import Counter
import regex as re
from typing import TYPE_CHECKING
from .utils import parse_static_mods

if TYPE_CHECKING:
    from .annot_base import ProFormaAnnotation


def condense_static_mods(
    annotation: ProFormaAnnotation, inplace: bool = True
) -> ProFormaAnnotation:
    """
    Condense static mods into internal mods
    """
    if inplace is False:
        return condense_static_mods(annotation.copy(), inplace=True)

    if not annotation.has_static_mods:
        return annotation

    static_mod_dict = parse_static_mods(annotation.get_static_mod_list().data)
    nterm_mod = static_mod_dict.pop("N-Term", None)
    cterm_mod = static_mod_dict.pop("C-Term", None)

    if nterm_mod is not None:
        annotation.extend_nterm_mods(nterm_mod)

    if cterm_mod is not None:
        annotation.extend_cterm_mods(cterm_mod)
        # Handle amino acid specific mods

    for aa, mods in static_mod_dict.items():
        indexes = [m.start() for m in re.finditer(aa, annotation.sequence)]
        for index in indexes:
            annotation.extend_internal_mods_at_index(index, mods)

    return annotation


def count_residues(
    annotation: ProFormaAnnotation, include_mods: bool = True
) -> dict[str, int]:
    """
    Count the occurrences of each residue in the sequence.

    Args:
        annotation: The ProFormaAnnotation to analyze
        include_mods: Whether to include modifications in the count

    Returns:
        Dictionary mapping residue strings to their counts
    """
    if not include_mods:
        return dict(Counter(annotation.sequence))

    split_annotations = annotation.split()
    return dict(Counter(a.serialize() for a in split_annotations))


def percent_residues(
    annotation: ProFormaAnnotation,
    include_mods: bool = True,
    precision: int | None = None,
) -> dict[str, float]:
    """
    Calculate the percentage of each residue in the sequence.

    Args:
        annotation: The ProFormaAnnotation to analyze
        include_mods: Whether to include modifications in the calculation
        precision: Number of decimal places to round to

    Returns:
        Dictionary mapping residue strings to their percentages
    """
    residue_counts = count_residues(annotation, include_mods=include_mods)
    total_count = sum(residue_counts.values())

    if total_count == 0:
        return {}

    percentages = {
        residue: (count / total_count) * 100
        for residue, count in residue_counts.items()
    }

    if precision is not None:
        percentages = {k: round(v, precision) for k, v in percentages.items()}

    return percentages


def is_subsequence(
    annotation: ProFormaAnnotation,
    other: ProFormaAnnotation,
    ignore_mods: bool = False,
    ignore_intervals: bool = True,
) -> bool:
    """
    Check if the annotation is a subsequence of another annotation.

    Args:
        annotation: The annotation to search for
        other: The annotation to search in
        ignore_mods: If True, only compare sequences (ignore modifications)
        ignore_intervals: If True, ignore interval modifications

    Returns:
        True if annotation is a subsequence of other
    """
    if annotation.sequence not in other.sequence:
        return False

    # Loop over all starting indexes where the sequence is found
    for match in re.finditer(annotation.sequence, other.sequence):
        start = match.start()

        # Check if all modifications are also a subsequence
        sliced_annot = other.slice(
            start, start + len(annotation.sequence), inplace=False
        )

        if ignore_mods and sliced_annot.sequence == annotation.sequence:
            return True

        if ignore_intervals:
            sliced_clean = _remove_intervals_copy(sliced_annot)
            annot_clean = _remove_intervals_copy(annotation)
            if sliced_clean == annot_clean:
                return True

        if sliced_annot == annotation:
            return True

    return False


def find_indices(
    annotation: ProFormaAnnotation,
    other: ProFormaAnnotation,
    ignore_mods: bool = False,
    ignore_intervals: bool = True,
) -> list[int]:
    """
    Find all occurrences of the annotation in another annotation.

    Args:
        annotation: The annotation to search for
        other: The annotation to search in
        ignore_mods: If True, only compare sequences
        ignore_intervals: If True, ignore interval modifications

    Returns:
        List of starting indices where annotation occurs in other

    Raises:
        TypeError: If other is not a ProFormaAnnotation
    """
    if not isinstance(other, type(annotation)):
        raise TypeError(
            f"other must be a {type(annotation).__name__}, got {type(other).__name__}"
        )

    if not annotation.sequence or not other.sequence:
        return []

    indices: list[int] = []
    for match in re.finditer(annotation.sequence, other.sequence):
        start = match.start()
        sliced_other = other.slice(
            start, start + len(annotation.sequence), inplace=False
        )

        if is_subsequence(annotation, sliced_other, ignore_mods, ignore_intervals):
            indices.append(start)

    return indices


# Helper functions


def _remove_intervals_copy(annotation: ProFormaAnnotation) -> ProFormaAnnotation:
    """Create a copy of annotation with intervals removed"""
    copy_annot = annotation.copy()
    copy_annot.intervals = None
    return copy_annot
