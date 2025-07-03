import itertools
import random
from typing import *

from ..constants import ModType
from ..proforma_dataclasses import Interval
from ..utils2 import ModLabler
from .annot_base import ProFormaAnnotationBase


class ProFormaAnnotationManipulation(ProFormaAnnotationBase):
    """
    Manipulation Methods
    """

    def strip(self, inplace: bool = False) -> "ProFormaAnnotationBase":
        """
        Remove all modifications from the annotation and return a new annotation with the stripped sequence.
        """

        if inplace is False:
            return self.copy().strip(inplace=True)

        self.isotope_mods = None
        self.static_mods = None
        self.labile_mods = None
        self.unknown_mods = None
        self.nterm_mods = None
        self.cterm_mods = None
        self.internal_mods = None
        self.intervals = None
        self.charge = None
        self.charge_adducts = None

        return self

    def slice(
        self,
        start: Optional[int],
        stop: Optional[int],
        inplace: bool = False,
        keep_terms: bool = False,
        slice_intervals: bool = True,
        keep_labile: bool = True,
    ) -> "ProFormaAnnotationBase":
        """
        Slice the annotation sequence and return a new annotation with the sliced sequence and modifications.
        """
        if inplace is False:
            return self.copy().slice(
                start=start,
                stop=stop,
                inplace=True,
                keep_terms=keep_terms,
                slice_intervals=slice_intervals,
                keep_labile=keep_labile,
            )

        # Handle default values and negative indices
        seq_len = len(self.sequence)
        if start is None:
            start = 0
        elif start < 0:
            start = seq_len + start
            if start < 0:
                raise ValueError(
                    "Start index is out of bounds for the sequence length."
                )
        elif start > seq_len:
            # Ensure start does not exceed sequence length
            raise ValueError(
                "Start index exceeds the sequence length. "
                f"Sequence length is {seq_len}, but start is {start}."
            )

        if stop is None:
            stop = seq_len
        elif stop < 0:
            stop = seq_len + stop
            if stop < 0:
                raise ValueError("Stop index is out of bounds for the sequence length.")
        elif stop > seq_len:
            # Ensure stop does not exceed sequence length
            raise ValueError(
                "Stop index exceeds the sequence length. "
                f"Sequence length is {seq_len}, but stop is {stop}."
            )

        # Ensure start <= stop
        if start > stop:
            raise ValueError(
                f"Start index {start} cannot be greater than stop index {stop}."
            )

        if not keep_labile:
            self.labile_mods = None

        new_sequence = self.sequence[start:stop]

        # Early return if no modifications exist
        if not self.has_mods():
            self.sequence = new_sequence
            return self

        # Adjust internal modifications
        new_internal_mods = {}
        if self.has_internal_mods():
            for pos, mods in self.internal_mods.items():
                if start <= pos < stop:
                    new_internal_mods[pos - start] = mods

        # Adjust intervals
        new_intervals = []
        for interval in self.intervals:
            # Check if interval overlaps with the slice
            if interval.start < stop and interval.end > start:
                # Calculate new positions relative to slice
                new_interval_start = interval.start - start
                new_interval_end = interval.end - start

                if slice_intervals == False:
                    if new_interval_start < 0:
                        raise ValueError(
                            f"Interval start {new_interval_start} is out of bounds for the slice."
                        )
                    if new_interval_end > stop - start:
                        raise ValueError(
                            f"Interval end {new_interval_end} is out of bounds for the slice."
                        )

                new_interval_start = max(0, new_interval_start)
                new_interval_end = min(stop - start, new_interval_end)

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

        # Handle terminal modifications
        if start > 0 and keep_terms is False:
            self.nterm_mods = None
        if stop < seq_len and keep_terms is False:
            self.cterm_mods = None

        # Update annotation
        self.sequence = new_sequence
        self.internal_mods = new_internal_mods if new_internal_mods else None
        self.intervals = new_intervals

        return self

    def shift(
        self,
        n: int,
        inplace: bool = False,
        shift_intervals: bool = True,
        slice_intervals: bool = True,
    ) -> "ProFormaAnnotationBase":
        """
        Shift the annotation by n positions in a cyclic manner.

        :param n: Number of positions to shift (positive shifts right, negative shifts left)
        :type n: int
        :param inplace: If True, modify the current annotation. If False, create a copy.
        :type inplace: bool
        :param shift_intervals: If True, shift intervals with the sequence. If False, keep intervals in original positions.
        :type shift_intervals: bool
        :param slice_intervals: If True, slice intervals that get cut off by shift. If False, raise error when this happens.
        :type slice_intervals: bool
        :return: The shifted annotation
        :rtype: ProFormaAnnotationBase
        """

        if inplace is False:
            return self.copy().shift(
                n=n,
                inplace=True,
                shift_intervals=shift_intervals,
                slice_intervals=slice_intervals,
            )

        seq_len = len(self.sequence)
        if seq_len == 0:
            return self  # No shift needed for empty sequence

        effective_shift = n % seq_len  # Normalize shift to within sequence length
        if effective_shift == 0:
            return self

        shifted_sequence = (
            self.sequence[effective_shift:] + self.sequence[:effective_shift]
        )

        new_internal_mods = {}
        for mod_index, mods in self.internal_mods.items():
            shifted_index = (mod_index - effective_shift) % seq_len
            new_internal_mods[shifted_index] = mods

        labler = ModLabler()

        # Handle intervals based on shift_intervals and slice_intervals parameters
        new_intervals = []
        for i, interval in enumerate(self.intervals):
            if shift_intervals:
                # Shift intervals with the sequence
                new_start = (interval.start - effective_shift) % seq_len
                new_end = (interval.end - effective_shift) % seq_len

                if new_end == 0:
                    new_end = len(self.sequence)

                # Check if interval wraps around (gets cut off by the shift)
                if new_start >= new_end:
                    if slice_intervals:
                        # Split the interval that wraps around
                        # First part: from new_start to end of sequence

                        new_interval = interval.copy()
                        new_interval.start = new_start
                        new_interval.end = seq_len

                        new_intervals.append(new_interval)
                        # Second part: from start of sequence to new_end
                        if new_end > 0:

                            new_interval_mods = []
                            if new_intervals[-1].has_mods():

                                mod_labels = []
                                for mod in new_intervals[-1].mods:
                                    mod.tag = labler.label
                                    mod_labels.append(labler.label)
                                    labler.increment_num()

                                labler.increment_letter()

                                for i, mod in enumerate(interval.mods):
                                    new_mod = mod.copy()
                                    new_mod.tag = mod_labels[i]
                                    new_mod.dislpay_val = False
                                    new_interval_mods.append(new_mod)

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
                            f"Interval [{interval.start}:{interval.end}] would be cut off by shift of {n} positions. "
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

        if shift_intervals is True:
            self.intervals = new_intervals

        self.sequence = shifted_sequence  # already a copy
        self.internal_mods = new_internal_mods  # already a copy

        return self

    def shuffle(
        self, seed: Optional[Any] = None, inplace: bool = False
    ) -> "ProFormaAnnotationBase":
        """
        Shuffle the annotation sequence and return a new annotation with the shuffled sequence.
        """

        if inplace is False:
            return self.copy().shuffle(seed, inplace=True)

        if seed is not None:
            random.seed(seed)

        seq_len = len(self.sequence)
        if seq_len <= 1:
            return self

        # Create list of (amino_acid, original_position) tuples
        sequence_positions = list(enumerate(self.sequence))
        random.shuffle(sequence_positions)

        # Extract shuffled sequence and create position mapping
        shuffled_sequence = []
        position_mapping = {}

        for new_pos, (original_pos, aa) in enumerate(sequence_positions):
            shuffled_sequence.append(aa)
            position_mapping[original_pos] = new_pos

        # Update sequence
        self.sequence = "".join(shuffled_sequence)

        # Update internal modifications with new positions
        if self.has_internal_mods():
            new_internal_mods = {}
            for original_pos, mods in self.internal_mods.items():
                new_pos = position_mapping[original_pos]
                new_internal_mods[new_pos] = mods
            self.internal_mods = new_internal_mods if new_internal_mods else None

        return self

    def reverse(
        self, inplace: bool = False, swap_terms: bool = False
    ) -> "ProFormaAnnotationBase":
        """
        Reverse the annotation sequence and return a new annotation with the reversed sequence.
        """

        if inplace is False:
            return self.copy().reverse(inplace=True, swap_terms=swap_terms)

        self.sequence = self.sequence[::-1]

        # Reverse internal modifications
        if self.has_internal_mods():
            new_internal_mods = {}
            for pos, mods in self.internal_mods.items():
                new_pos = len(self.sequence) - 1 - pos
                new_internal_mods[new_pos] = mods
            self.internal_mods = new_internal_mods

        if swap_terms:
            # Swap terminal modifications
            self.nterm_mods, self.cterm_mods = self.cterm_mods, self.nterm_mods

        return self

    def sort(
        self,
        inplace: bool = False,
        key: Optional[callable] = None,
        reverse: bool = False,
    ) -> "ProFormaAnnotationBase":
        """
        Sort the residues in the annotation sequence and return a new annotation with the sorted sequence.
        """

        if inplace is False:
            return self.copy().sort(inplace=True, key=key, reverse=reverse)

        seq_len = len(self.sequence)
        if seq_len <= 1:
            return self

        filtered_annotation = self.filter_mods(mods="interval", inplace=False)
        sorted_comps = list(filtered_annotation.split())

        # Create list of (amino_acid, original_position) tuples
        sequence_positions = list(enumerate(sorted_comps))

        # Sort based on the key function or alphabetically
        if key is None:
            sequence_positions.sort(
                key=lambda x: x[1].serialize()
            )  # Sort by amino acid / mod
        else:
            sequence_positions.sort(
                key=lambda x: key(x[1].serialize())
            )  # Sort by custom key

        if reverse is True:
            sequence_positions.reverse()

        # Extract sorted sequence and create position mapping
        sorted_sequence = []
        position_mapping = {}

        for new_pos, (original_pos, comp) in enumerate(sequence_positions):
            sorted_sequence.append(comp.sequence)
            position_mapping[original_pos] = new_pos

        # Update sequence
        self.sequence = "".join(sorted_sequence)

        # Update internal modifications with new positions
        if self.has_internal_mods():
            new_internal_mods = {}
            for original_pos, mods in self.internal_mods.items():
                new_pos = position_mapping[original_pos]
                new_internal_mods[new_pos] = mods
            self.internal_mods = new_internal_mods if new_internal_mods else None

        return self

    def sliding_windows(
        self,
        window_size: int,
        keep_terms: bool = False,
        slice_intervals: bool = True,
        keep_labile: bool = True,
        reverse: bool = False,
    ) -> Generator["ProFormaAnnotationBase", None, None]:
        """
        Generate sliding windows of the annotation with a specified size.
        """
        if not self.has_sequence():
            raise ValueError(
                "Annotation must have a sequence to create sliding windows"
            )

        if window_size <= 0:
            raise ValueError("Window size must be positive")

        seq_len = len(self.sequence)
        if window_size > seq_len:
            raise ValueError(
                f"Window size {window_size} cannot be greater than sequence length {seq_len}."
            )

        if reverse:
            # Generate windows from right to left
            for start in range(seq_len - window_size, -1, -1):
                stop = start + window_size
                yield self.slice(
                    start=start,
                    stop=stop,
                    inplace=False,
                    keep_terms=keep_terms,
                    slice_intervals=slice_intervals,
                    keep_labile=keep_labile,
                )
        else:
            # Generate windows from left to right
            for start in range(seq_len - window_size + 1):
                stop = start + window_size
                yield self.slice(
                    start=start,
                    stop=stop,
                    inplace=False,
                    keep_terms=keep_terms,
                    slice_intervals=slice_intervals,
                    keep_labile=keep_labile,
                )

    def split(self) -> Generator["ProFormaAnnotationBase", None, None]:
        """
        Split each amino acid in the sequence into a separate ProFormaAnnotationBase
        """

        # for labile mods only include on first amino acid
        labile_mods = self.pop_labile_mods()

        if self.has_intervals():
            raise ValueError(
                "Cannot split annotation with intervals. Please slice the annotation first."
            )

        for i, _ in enumerate(self.sequence):
            s = self.slice(i, i + 1, inplace=False)
            if i == 0 and labile_mods:
                s.add_labile_mods(labile_mods, append=True, inplace=True)
            yield s

    def permutations(
        self, size: Optional[int] = None
    ) -> Generator["ProFormaAnnotationBase", None, None]:
        """
        Generate all permutations of the annotation sequence.
        """

        if size is None:
            size = len(self)

        outside_mods = self.pop_mods(
            mods=[
                ModType.NTERM,
                ModType.CTERM,
                ModType.LABILE,
                ModType.ISOTOPE,
                ModType.STATIC,
                ModType.UNKNOWN,
                ModType.CHARGE,
                ModType.CHARGE_ADDUCTS,
            ]
        )

        for i in itertools.permutations(self.split(), size):
            yield ProFormaAnnotationBase(
                sequence="".join(a.serialize() for a in i)
            ).add_mods(outside_mods)

    def product(
        self, repeat: Optional[int] = None
    ) -> Generator["ProFormaAnnotationBase", None, None]:
        """
        Generate the product of the annotation sequence with itself.
        """

        if repeat is None:
            repeat = len(self)

        outside_mods = self.pop_mods(
            mods=[
                ModType.NTERM,
                ModType.CTERM,
                ModType.LABILE,
                ModType.ISOTOPE,
                ModType.STATIC,
                ModType.UNKNOWN,
                ModType.CHARGE,
                ModType.CHARGE_ADDUCTS,
            ]
        )

        for i in itertools.product(self.split(), repeat=repeat):
            yield ProFormaAnnotationBase(
                sequence="".join(a.serialize() for a in i)
            ).add_mods(outside_mods)

    def combinations(
        self, r: Optional[int] = None
    ) -> Generator["ProFormaAnnotationBase", None, None]:
        """
        Generate all combinations of the annotation sequence.
        """

        if r is None:
            r = len(self)

        outside_mods = self.pop_mods(
            mods=[
                ModType.NTERM,
                ModType.CTERM,
                ModType.LABILE,
                ModType.ISOTOPE,
                ModType.STATIC,
                ModType.UNKNOWN,
                ModType.CHARGE,
                ModType.CHARGE_ADDUCTS,
            ]
        )

        for i in itertools.combinations(self.split(), r=r):
            yield ProFormaAnnotationBase(
                sequence="".join(a.serialize() for a in i)
            ).add_mods(outside_mods)

    def combinations_with_replacement(
        self, r: Optional[int] = None
    ) -> Generator["ProFormaAnnotationBase", None, None]:
        """
        Generate all combinations of the annotation sequence with replacement.
        """

        if r is None:
            r = len(self)

        outside_mods = self.pop_mods(
            mods=[
                ModType.NTERM,
                ModType.CTERM,
                ModType.LABILE,
                ModType.ISOTOPE,
                ModType.STATIC,
                ModType.UNKNOWN,
                ModType.CHARGE,
                ModType.CHARGE_ADDUCTS,
            ]
        )

        for i in itertools.combinations_with_replacement(self.split(), r=r):
            yield ProFormaAnnotationBase(
                sequence="".join(a.serialize() for a in i)
            ).add_mods(outside_mods)
