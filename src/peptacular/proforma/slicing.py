from __future__ import annotations
from collections.abc import Generator
import random
from typing import TYPE_CHECKING, Any, Callable

from .dclasses import ModDict, IntervalList, Interval

if TYPE_CHECKING:
    from .annotation import ProFormaAnnotation


def slice_annotation(
    annotation: ProFormaAnnotation,
    start: int | None,
    stop: int | None,
    inplace: bool = False,
    keep_terms: bool = False,
    keep_labile: bool = True,
) -> ProFormaAnnotation:
    """
    Slice the annotation sequence and return a new annotation with the sliced sequence and modifications.

    Args:
        annotation: The ProFormaAnnotation to slice
        start: Start index (inclusive)
        stop: Stop index (exclusive)
        inplace: Whether to modify the annotation in place
        keep_terms: Whether to keep terminal modifications when slicing
        keep_labile: Whether to keep labile modifications

    Returns:
        Sliced annotation

    Raises:
        ValueError: If indices are out of bounds or invalid
    """
    if not inplace:
        return slice_annotation(
            annotation.copy(),
            start=start,
            stop=stop,
            inplace=True,
            keep_terms=keep_terms,
            keep_labile=keep_labile,
        )

    # Handle default values and negative indices
    seq_len = len(annotation)
    start = _normalize_start_index(start, seq_len)
    stop = _normalize_stop_index(stop, seq_len)

    _validate_slice_indices(start, stop, seq_len)

    if not keep_labile:
        annotation.labile_mods = None

    new_sequence = annotation.sequence[start:stop]

    # Early return if no modifications exist
    if not annotation.has_mods:
        annotation.sequence = new_sequence
        return annotation

    # Adjust internal modifications
    new_internal_mods = _adjust_internal_mods(
        annotation.get_internal_mod_dict(), start, stop
    )

    # Adjust intervals
    new_intervals = _adjust_intervals(annotation.get_interval_list(), start, stop)

    # Handle terminal modifications
    if start > 0 and not keep_terms:
        annotation.nterm_mods = None
    if stop < seq_len and not keep_terms:
        annotation.cterm_mods = None

    # Update annotation
    annotation.sequence = new_sequence
    annotation.internal_mods = new_internal_mods if new_internal_mods.has_mods else None
    annotation.intervals = new_intervals if new_intervals.has_intervals else None

    return annotation


def split_annotation(annotation: ProFormaAnnotation) -> list[ProFormaAnnotation]:
    """
    Split each amino acid in the sequence into a separate ProFormaAnnotation.

    Args:
        annotation: The ProFormaAnnotation to split

    Returns:
        List of single amino acid annotations

    Raises:
        ValueError: If annotation contains intervals
    """
    if annotation.has_intervals:
        raise ValueError("Cannot split annotation with intervals.")

    # For labile mods only include on first amino acid
    labile_mods = annotation.pop_labile_mods()

    result: list[ProFormaAnnotation] = []
    for i, _ in enumerate(annotation.sequence):
        sliced = slice_annotation(annotation, i, i + 1, inplace=False)
        if i == 0 and labile_mods:
            sliced.extend_labile_mods(labile_mods)
        result.append(sliced)

    return result


def shift_annotation(
    annotation: ProFormaAnnotation,
    n: int,
    inplace: bool = False,
    shift_intervals: bool = True,
    slice_intervals: bool = True,
) -> ProFormaAnnotation:
    """
    Shift the annotation by n positions in a cyclic manner.

    Args:
        annotation: The ProFormaAnnotation to shift
        n: Number of positions to shift (positive shifts right, negative shifts left)
        inplace: If True, modify the annotation in place
        shift_intervals: If True, shift intervals with the sequence
        slice_intervals: If True, slice intervals that get cut off by shift

    Returns:
        The shifted annotation

    Raises:
        ValueError: If intervals would be cut off and slice_intervals=False
    """
    if not inplace:
        return shift_annotation(
            annotation.copy(),
            n=n,
            inplace=True,
            shift_intervals=shift_intervals,
            slice_intervals=slice_intervals,
        )

    seq_len = len(annotation.sequence)
    if seq_len == 0:
        return annotation

    effective_shift = n % seq_len
    if effective_shift == 0:
        return annotation

    # Shift sequence
    shifted_sequence = (
        annotation.sequence[effective_shift:] + annotation.sequence[:effective_shift]
    )

    # Shift internal modifications
    new_internal_mods = {}
    if annotation.has_internal_mods:
        for mod_index, mods in annotation.get_internal_mod_dict().items():
            shifted_index = (mod_index - effective_shift) % seq_len
            new_internal_mods[shifted_index] = mods

    # Handle intervals
    new_intervals: list[Interval] = []
    if annotation.has_intervals and shift_intervals:
        new_intervals = _shift_intervals(
            annotation.get_interval_list().data,
            effective_shift,
            seq_len,
            slice_intervals,
        )

    # Update annotation
    annotation.sequence = shifted_sequence
    annotation.internal_mods = new_internal_mods if new_internal_mods else None
    if shift_intervals:
        annotation.intervals = new_intervals if new_intervals else None

    return annotation


def shuffle_annotation(
    annotation: ProFormaAnnotation, seed: Any = None, inplace: bool = False
) -> ProFormaAnnotation:
    """
    Shuffle the annotation sequence randomly.

    Args:
        annotation: The ProFormaAnnotation to shuffle
        seed: Random seed for reproducible shuffling
        inplace: If True, modify the annotation in place

    Returns:
        The shuffled annotation
    """
    if not inplace:
        return shuffle_annotation(annotation.copy(), seed, inplace=True)

    if seed is not None:
        random.seed(seed)

    seq_len = len(annotation.sequence)
    if seq_len <= 1:
        return annotation

    # Create position mapping through shuffling
    sequence_positions = list(enumerate(annotation.sequence))
    random.shuffle(sequence_positions)

    # Extract shuffled sequence and create position mapping
    shuffled_sequence: list[str] = []
    position_mapping: dict[int, int] = {}

    for new_pos, (original_pos, aa) in enumerate(sequence_positions):
        shuffled_sequence.append(aa)
        position_mapping[original_pos] = new_pos

    # Update sequence
    annotation.sequence = "".join(shuffled_sequence)

    # Update internal modifications with new positions
    if annotation.has_internal_mods:
        new_internal_mods = {}
        for original_pos, mods in annotation.get_internal_mod_dict().items():
            new_pos = position_mapping[original_pos]
            new_internal_mods[new_pos] = mods
        annotation.internal_mods = new_internal_mods if new_internal_mods else None

    return annotation


def reverse_annotation(
    annotation: ProFormaAnnotation, inplace: bool = False, swap_terms: bool = False
) -> ProFormaAnnotation:
    """
    Reverse the annotation sequence.

    Args:
        annotation: The ProFormaAnnotation to reverse
        inplace: If True, modify the annotation in place
        swap_terms: If True, swap N-term and C-term modifications

    Returns:
        The reversed annotation
    """
    if not inplace:
        return reverse_annotation(
            annotation.copy(), inplace=True, swap_terms=swap_terms
        )

    # Reverse sequence
    annotation.sequence = annotation.sequence[::-1]

    # Reverse internal modifications
    if annotation.has_internal_mods:
        new_internal_mods = {}
        for pos, mods in annotation.get_internal_mod_dict().items():
            new_pos = len(annotation.sequence) - 1 - pos
            new_internal_mods[new_pos] = mods
        annotation.internal_mods = new_internal_mods

    # Swap terminal modifications if requested
    if swap_terms:
        annotation.nterm_mods, annotation.cterm_mods = (
            annotation.cterm_mods,
            annotation.nterm_mods,
        )

    return annotation


def sort_annotation(
    annotation: ProFormaAnnotation,
    inplace: bool = False,
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
) -> ProFormaAnnotation:
    """
    Sort the residues in the annotation sequence.

    Args:
        annotation: The ProFormaAnnotation to sort
        inplace: If True, modify the annotation in place
        key: Function to use for sorting key (default: alphabetical)
        reverse: If True, sort in descending order

    Returns:
        The sorted annotation
    """
    if not inplace:
        return sort_annotation(
            annotation.copy(), inplace=True, key=key, reverse=reverse
        )

    seq_len = len(annotation.sequence)
    if seq_len <= 1:
        return annotation

    # Remove intervals to avoid complications during sorting
    filtered_annotation = annotation.copy()
    filtered_annotation.intervals = None

    sorted_components = split_annotation(filtered_annotation)

    # Create list of (component, original_position) tuples
    sequence_positions = list(enumerate(sorted_components))

    # Sort based on key function or alphabetically
    if key is None:
        sequence_positions.sort(key=lambda x: x[1].serialize(), reverse=reverse)
    else:
        sequence_positions.sort(key=lambda x: key(x[1].serialize()), reverse=reverse)

    # Extract sorted sequence and create position mapping
    sorted_sequence: list[str] = []
    position_mapping: dict[int, int] = {}

    for new_pos, (original_pos, comp) in enumerate(sequence_positions):
        sorted_sequence.append(comp.sequence)
        position_mapping[original_pos] = new_pos

    # Update sequence
    annotation.sequence = "".join(sorted_sequence)

    # Update internal modifications with new positions
    if annotation.has_internal_mods:
        new_internal_mods = {}
        for original_pos, mods in annotation.get_internal_mod_dict().items():
            new_pos = position_mapping[original_pos]
            new_internal_mods[new_pos] = mods
        annotation.internal_mods = new_internal_mods if new_internal_mods else None

    return annotation


def generate_sliding_windows(
    annotation: ProFormaAnnotation,
    window_size: int,
    keep_terms: bool = False,
    slice_intervals: bool = True,
    keep_labile: bool = True,
    reverse: bool = False,
) -> Generator[ProFormaAnnotation, None, None]:
    """
    Generate sliding windows of the annotation with a specified size.

    Args:
        annotation: The ProFormaAnnotation to generate windows from
        window_size: Size of each window
        keep_terms: Whether to keep terminal modifications in windows
        slice_intervals: Whether to slice intervals at window boundaries
        keep_labile: Whether to keep labile modifications in windows
        reverse: If True, generate windows from right to left

    Yields:
        ProFormaAnnotation instances representing each window

    Raises:
        ValueError: If window_size is invalid or annotation has no sequence
    """
    if not annotation.has_sequence:
        raise ValueError("Annotation must have a sequence to create sliding windows")

    if window_size <= 0:
        raise ValueError("Window size must be positive")

    seq_len = len(annotation.sequence)
    if window_size > seq_len:
        raise ValueError(
            f"Window size {window_size} cannot be greater than sequence length {seq_len}."
        )

    if reverse:
        # Generate windows from right to left
        for start in range(seq_len - window_size, -1, -1):
            stop = start + window_size
            yield slice_annotation(
                annotation,
                start=start,
                stop=stop,
                inplace=False,
                keep_terms=keep_terms,
                keep_labile=keep_labile,
            )
    else:
        # Generate windows from left to right
        for start in range(seq_len - window_size + 1):
            stop = start + window_size
            yield slice_annotation(
                annotation,
                start=start,
                stop=stop,
                inplace=False,
                keep_terms=keep_terms,
                keep_labile=keep_labile,
            )


# Helper functions


def _shift_intervals(
    intervals: list[Interval], effective_shift: int, seq_len: int, slice_intervals: bool
) -> list[Interval]:
    """Helper function to shift intervals during cyclic shift"""
    new_intervals: list[Interval] = []

    for interval in intervals:
        # Shift intervals with the sequence
        new_start = (interval.start - effective_shift) % seq_len
        new_end = (interval.end - effective_shift) % seq_len

        if new_end == 0:
            new_end = seq_len

        # Check if interval wraps around (gets cut off by the shift)
        if new_start >= new_end:
            if slice_intervals:
                # Split the interval that wraps around
                # First part: from new_start to end of sequence
                new_intervals.append(
                    Interval(
                        start=new_start,
                        end=seq_len,
                        ambiguous=interval.ambiguous,
                        mods=interval.mods,
                    )
                )
                # Second part: from start of sequence to new_end
                if new_end > 0:
                    # Create modification copies for the second part
                    new_interval_mods = []
                    if interval.has_mods:
                        # Simple copy without complex labeling logic
                        new_interval_mods = [mod.copy() for mod in interval.mods]

                    new_intervals.append(
                        Interval(
                            start=0,
                            end=new_end,
                            ambiguous=interval.ambiguous,
                            mods=new_interval_mods,
                        )
                    )
            else:
                raise ValueError(
                    f"Interval [{interval.start}:{interval.end}] would be cut off by shift. "
                    f"Set slice_intervals=True to allow interval slicing."
                )
        else:
            # Normal case - interval doesn't wrap around
            new_intervals.append(
                Interval(
                    start=new_start,
                    end=new_end,
                    ambiguous=interval.ambiguous,
                    mods=interval.mods,
                )
            )

    return new_intervals


def _normalize_start_index(start: int | None, seq_len: int) -> int:
    """Normalize start index handling None and negative values"""
    if start is None:
        return 0
    elif start < 0:
        start = seq_len + start
        if start < 0:
            raise ValueError("Start index is out of bounds for the sequence length.")
        return start
    elif start > seq_len:
        raise ValueError(
            f"Start index exceeds the sequence length. "
            f"Sequence length is {seq_len}, but start is {start}."
        )
    return start


def _normalize_stop_index(stop: int | None, seq_len: int) -> int:
    """Normalize stop index handling None and negative values"""
    if stop is None:
        return seq_len
    elif stop < 0:
        stop = seq_len + stop
        if stop < 0:
            raise ValueError("Stop index is out of bounds for the sequence length.")
        return stop
    elif stop > seq_len:
        raise ValueError(
            f"Stop index exceeds the sequence length. "
            f"Sequence length is {seq_len}, but stop is {stop}."
        )
    return stop


def _validate_slice_indices(start: int, stop: int, seq_len: int) -> None:
    """Validate that slice indices are valid"""
    if start > stop:
        raise ValueError(
            f"Start index {start} cannot be greater than stop index {stop}."
        )


def _adjust_internal_mods(internal_mods: ModDict, start: int, stop: int) -> ModDict:
    """Adjust internal modifications for slicing"""
    new_internal_mods = ModDict()
    for pos, mods in internal_mods.items():
        if start <= pos < stop:
            new_internal_mods[pos - start] = mods
    return new_internal_mods


def _adjust_intervals(intervals: IntervalList, start: int, stop: int) -> IntervalList:
    """Adjust intervals for slicing"""
    new_intervals = IntervalList()
    for interval in intervals:
        # Check if interval overlaps with the slice
        if interval.start < stop and interval.end > start:
            # Calculate new positions relative to slice
            new_interval_start = max(0, interval.start - start)
            new_interval_end = min(stop - start, interval.end - start)

            # Only add if the interval has meaningful content in the slice
            if new_interval_start < new_interval_end:
                new_intervals.append(
                    Interval(
                        start=new_interval_start,
                        end=new_interval_end,
                        ambiguous=interval.ambiguous,
                        mods=interval.mods,
                    )
                )
    return new_intervals
