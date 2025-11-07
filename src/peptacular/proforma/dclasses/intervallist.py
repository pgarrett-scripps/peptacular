from __future__ import annotations

import copy
import warnings
from collections.abc import Iterable, Iterator, MutableSequence
from typing import Any, Self, SupportsIndex

from ..dclasses.modlist import MODLIST_DATATYPE, ModList
from .interval import Interval, ModInterval, setup_interval

INTERVALLIST_DATATYPE = (
    Interval
    | ModInterval
    | tuple[
        int, int, bool, ModList | Iterable[MODLIST_DATATYPE] | MODLIST_DATATYPE | None
    ]
)


class IntervalList(MutableSequence[Interval]):
    """
    A list of Interval objects that automatically merges intervals with identical spans.

    Accepts Interval instances, ModInterval instances, tuples (start, end, ambiguous[, mods]),
    or iterables thereof. Internally stores everything as Interval instances with no duplicate spans.
    """

    __slots__ = ("_data",)

    def __init__(self, data: Iterable[INTERVALLIST_DATATYPE] | None = None) -> None:
        self._data: list[Interval] = []

        if data is not None:
            self.extend(data)

    # Required MutableSequence abstract methods
    def __len__(self) -> int:
        return len(self._data)

    def __getitem__(self, index: int | slice) -> Interval | list[Interval]:
        return self._data[index]

    def __setitem__(
        self, index: int | slice, value: Interval | Iterable[Interval]
    ) -> None:
        self._data[index] = value  # type: ignore

    def __delitem__(self, index: int | slice) -> None:
        del self._data[index]

    def __iter__(self) -> Iterator[Interval]:
        return iter(self._data)

    def _normalize_input(self, item: INTERVALLIST_DATATYPE) -> Interval:
        """Convert input to Interval instance"""
        return setup_interval(item)

    def _find_same_span_index(self, interval: Interval) -> int | None:
        """Find index of existing interval with same span"""
        for i, existing in enumerate(self._data):
            if existing.start == interval.start and existing.end == interval.end:
                return i
        return None

    def _merge_or_append(self, item: Interval) -> None:
        """Merge with existing interval or append new one"""
        idx = self._find_same_span_index(item)
        if idx is not None:
            # Merge into existing interval
            self._data[idx].merge(item)
        else:
            # Append new interval
            self._data.append(item)

    def append(self, item: INTERVALLIST_DATATYPE) -> None:
        """Add an interval, merging with existing interval if same span"""
        interval = self._normalize_input(item)
        self._merge_or_append(interval)

    def extend(self, other: Iterable[INTERVALLIST_DATATYPE]) -> None:
        """Extend with intervals, merging duplicates"""
        if isinstance(other, Iterable):  # type: ignore
            for item in other:
                self.append(item)
        else:
            raise TypeError(f"Cannot extend IntervalList with {type(other)}")

    def insert(self, i: int, item: INTERVALLIST_DATATYPE) -> None:
        """Insert behaves like append (merges duplicates, ignores index)"""
        warnings.warn(
            "IntervalList.insert() does not support index and behaves like append().",
            UserWarning,
            stacklevel=2,
        )
        self.append(item)

    def remove(self, item: INTERVALLIST_DATATYPE) -> None:
        """Remove an interval by exact match"""
        interval = self._normalize_input(item)

        # Find exact match
        for i, existing in enumerate(self._data):
            if existing == interval:
                del self._data[i]
                return

        raise ValueError(f"{item} is not in list")

    def discard(self, item: INTERVALLIST_DATATYPE) -> None:
        """Remove interval without raising error if not found"""
        try:
            self.remove(item)
        except ValueError:
            pass

    def count(self, item: INTERVALLIST_DATATYPE) -> int:
        """Count occurrences of an interval (0 or 1 due to merging)"""
        interval = self._normalize_input(item)
        return 1 if interval in self._data else 0

    def index(
        self,
        item: INTERVALLIST_DATATYPE,
        start: SupportsIndex = 0,
        stop: SupportsIndex | None = None,
    ) -> int:
        """Find index of interval"""
        interval = self._normalize_input(item)

        length = len(self._data)
        start_idx = start if isinstance(start, int) else 0
        stop_idx = stop if isinstance(stop, int) else length

        for i in range(start_idx, min(stop_idx, length)):
            if self._data[i] == interval:
                return i
        raise ValueError(f"{item} is not in list")

    def __contains__(self, item: object) -> bool:
        """Check if interval exists in list"""
        try:
            if isinstance(item, (Interval, ModInterval)):  # type: ignore
                interval = self._normalize_input(item)
                return interval in self._data
        except (TypeError, ValueError):
            pass
        return False

    def __add__(self, other: Iterable[INTERVALLIST_DATATYPE] | Self) -> IntervalList:
        """Return new IntervalList combining both lists"""
        if isinstance(other, IntervalList):
            other_intervallist = other
        elif isinstance(other, Iterable):  # type: ignore
            other_intervallist = IntervalList(other)
        else:
            raise TypeError(f"Cannot add IntervalList with {type(other)}")

        result = IntervalList()
        result.extend(self._data)
        result.extend(other_intervallist._data)
        return result

    def __iadd__(self, other: Iterable[INTERVALLIST_DATATYPE] | Self) -> IntervalList:
        """In-place addition"""
        if isinstance(other, IntervalList):
            self.extend(other._data)
        elif isinstance(other, Iterable):  # type: ignore
            self.extend(other)
        else:
            raise TypeError(f"Cannot add IntervalList with {type(other)}")
        return self

    def __eq__(self, other: Any) -> bool:
        """Equality comparison (order-independent due to merging)"""
        if not isinstance(other, IntervalList):
            try:
                other = IntervalList(other)
            except (TypeError, ValueError):
                return False

        if len(self._data) != len(other._data):
            return False

        # Since spans are unique due to merging, we can compare sets
        return set(self._data) == set(other._data)

    def copy(self, deep: bool = True) -> IntervalList:
        result = IntervalList()
        # Need to copy each Interval (which contains ModLists)
        result._data = [interval.copy(deep=False) for interval in self._data]
        return result

    def clear(self) -> None:
        """Remove all items from the IntervalList"""
        self._data.clear()

    @property
    def has_ambiguous_intervals(self) -> bool:
        """Check if there are any ambiguous intervals"""
        return any(interval.ambiguous for interval in self._data)

    def get_ambiguous_intervals(self) -> IntervalList:
        """Return a list of intervals that are ambiguous"""
        return IntervalList([interval for interval in self._data if interval.ambiguous])

    def pop_ambiguous_intervals(self) -> IntervalList:
        """Remove and return all ambiguous intervals"""
        ambiguous = self.get_ambiguous_intervals()
        self._data = [interval for interval in self._data if not interval.ambiguous]
        return ambiguous

    @property
    def has_unambiguous_intervals(self) -> bool:
        """Check if there are any unambiguous intervals"""
        return any(not interval.ambiguous for interval in self._data)

    def get_unambiguous_intervals(self) -> IntervalList:
        """Return a list of intervals that are unambiguous"""
        return IntervalList(
            [interval for interval in self._data if not interval.ambiguous]
        )

    def pop_unambiguous_intervals(self) -> IntervalList:
        """Remove and return all unambiguous intervals"""
        unambiguous = self.get_unambiguous_intervals()
        self._data = [interval for interval in self._data if interval.ambiguous]
        return unambiguous

    @property
    def has_intervals(self) -> bool:
        """Check if list has any intervals"""
        return len(self._data) > 0

    @property
    def is_empty(self) -> bool:
        """Check if list is empty"""
        return len(self._data) == 0

    @property
    def data(self) -> list[Interval]:
        """Compatibility property for existing code that accesses .data"""
        return self._data

    @data.setter
    def data(self, value: list[Interval]) -> None:
        """Setter for compatibility with existing code that sets .data"""
        self._data = value

    def __repr__(self) -> str:
        return f"IntervalList({self._data})"

    # Additional convenience methods

    def get_intervals_at_position(self, position: int) -> IntervalList:
        """Get all intervals that contain the given position"""
        return IntervalList(
            [
                interval
                for interval in self._data
                if interval.start <= position < interval.end
            ]
        )

    def get_overlapping_intervals(self, start: int, end: int) -> IntervalList:
        """Get all intervals that overlap with the given range"""
        return IntervalList(
            [
                interval
                for interval in self._data
                if interval.start < end and interval.end > start
            ]
        )

    def total_span(self) -> tuple[int, int] | None:
        """Get the total span covered by all intervals"""
        if not self._data:
            return None

        min_start = min(interval.start for interval in self._data)
        max_end = max(interval.end for interval in self._data)
        return (min_start, max_end)

    def get_mod_intervals(self) -> tuple[ModInterval, ...]:
        """Return a list of intervals that have modifications"""
        return tuple(interval.get_immutable_interval() for interval in self._data)


ACCEPTED_INTERVALLIST_INPUT_TYPES = (
    INTERVALLIST_DATATYPE | Iterable[INTERVALLIST_DATATYPE] | IntervalList | None
)


def setup_interval_list(data: ACCEPTED_INTERVALLIST_INPUT_TYPES) -> IntervalList:
    """Helper function to set up an IntervalList from various input types"""

    if isinstance(data, IntervalList):
        return data

    if data is None:
        return IntervalList()

    # Single interval
    if isinstance(data, (Interval, ModInterval)):
        return IntervalList([data])

    if isinstance(data, tuple) and len(data) == 4 and isinstance(data[0], int):
        return IntervalList([data])  # type: ignore

    # Iterable of intervals
    if isinstance(data, Iterable):  # type: ignore
        return IntervalList(data)  # type: ignore

    raise TypeError(
        f"Invalid type for interval data: {type(data)}. "
        f"Must be an Interval, ModInterval, or iterable of these."
    )


def populate_interval_list(
    interval_list: IntervalList, data: ACCEPTED_INTERVALLIST_INPUT_TYPES
) -> IntervalList:
    """Populates and returns an IntervalList from various input types."""

    if data is None:
        return interval_list

    if isinstance(data, IntervalList):
        interval_list += data
    elif isinstance(data, Iterable):  # type: ignore
        interval_list += data  # type: ignore
    elif isinstance(data, (Interval, ModInterval)):
        interval_list.append(data)
    elif isinstance(data, tuple) and len(data) == 4 and isinstance(data[0], int):
        interval_list.append(data)  # type: ignore
    else:
        raise TypeError(f"Cannot populate IntervalList with {type(data)}")

    return interval_list
