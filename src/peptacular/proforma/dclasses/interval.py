from __future__ import annotations

import copy
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any, Self

from ...mod import MOD_VALUE_TYPES
from .modlist import MODLIST_DATATYPE, ModList, setup_mod_list


@dataclass(slots=True)
class ModInterval:
    """An immutable sequence interval for return type"""

    start: int
    end: int
    ambiguous: bool
    mods: tuple[MOD_VALUE_TYPES, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError("Interval start must be non-negative")
        if self.end < self.start:
            raise ValueError("Interval end must be >= start")
        # if not ambiguous and no mods, raise error
        if not self.ambiguous and len(self.mods) == 0:
            raise ValueError("Interval must be ambiguous or have mods")


@dataclass(slots=True)
class Interval:
    """A sequence interval with modlist"""

    start: int
    end: int
    ambiguous: bool
    mods: ModList | Iterable[MODLIST_DATATYPE] | MODLIST_DATATYPE | None = None

    def __post_init__(self) -> None:
        self.mods: ModList = (
            setup_mod_list(self.mods, allow_dups=True, stackable=False)
            if not isinstance(self.mods, ModList)
            else self.mods
        )

        if self.start < 0:
            raise ValueError("Interval start must be non-negative")
        if self.end < self.start:
            raise ValueError("Interval end must be >= start")
        # if not ambiguous and no mods, raise error

    @property
    def has_mods(self) -> bool:
        """Check if the interval has modifications"""
        return len(self.mods) > 0

    def to_dict(self) -> dict[str, Any]:
        """Convert the interval to a dictionary"""
        return {
            "start": self.start,
            "end": self.end,
            "ambiguous": self.ambiguous,
            "mods": copy.deepcopy(self.mods),
        }

    def __repr__(self) -> str:
        return f"Interval({self.start}, {self.end}, {self.ambiguous}, {self.mods})"

    def __eq__(self, other: Any) -> bool:
        """Equality comparison for Interval objects"""
        if not isinstance(other, Interval):
            try:
                other_interval = setup_interval(other)
            except Exception:
                return False
        else:
            other_interval = other

        return (
            self.start == other_interval.start
            and self.end == other_interval.end
            and self.ambiguous == other_interval.ambiguous
            and self.mods == other_interval.mods
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.start,
                self.end,
                self.ambiguous,
                self.mods.flatten(sort=True),
            )
        )

    def copy(self, deep: bool = True) -> Self:
        """Create a copy of the interval"""
        if deep:
            return copy.deepcopy(self)
        return copy.copy(self)

    def merge(self, other: Interval | ModInterval) -> None:
        """
        Merge this Interval with another Interval.
        Merge only allowed when start and end match.
        Mods are combined and ambiguous is the logical OR of the two flags.
        """
        if not isinstance(other, ModInterval):
            other = setup_interval(other)
        else:
            raise TypeError(f"Cannot merge Interval with {type(other).__name__}")

        if self.start != other.start or self.end != other.end:
            raise ValueError(
                f"Cannot merge intervals with different spans: "
                f"({self.start},{self.end}) != ({other.start},{other.end})"
            )

        self.ambiguous = self.ambiguous or other.ambiguous
        self.mods += other.mods

    def get_immutable_interval(self) -> ModInterval:
        """Get an immutable representation of the interval"""
        return ModInterval(
            start=self.start,
            end=self.end,
            ambiguous=self.ambiguous,
            mods=self.mods.flatten(),
        )


ACCEPTED_INTERVAL_DATATYPE = (
    Interval
    | ModInterval
    | tuple[
        int, int, bool, ModList | Iterable[MODLIST_DATATYPE] | MODLIST_DATATYPE | None
    ]
)


# Type definitions for input handling
def setup_interval(interval: ACCEPTED_INTERVAL_DATATYPE) -> Interval:
    """Convert input to Interval instance"""
    if isinstance(interval, Interval):
        return interval
    elif isinstance(interval, ModInterval):  # type: ignore
        return Interval(
            start=interval.start,
            end=interval.end,
            ambiguous=interval.ambiguous,
            mods=interval.mods,
        )
    elif isinstance(interval, tuple):  # type: ignore
        if len(interval) == 4:
            try:
                start = int(interval[0])
                end = int(interval[1])
                ambiguous = bool(interval[2])
                mods = setup_mod_list(interval[3])
            except (ValueError, TypeError) as e:
                raise ValueError(f"Invalid tuple values for Interval: {e}")

            return Interval(
                start=start,
                end=end,
                ambiguous=ambiguous,
                mods=mods,
            )
        else:
            raise ValueError(f"Invalid tuple length for Interval: {len(interval)}")

    else:
        raise TypeError(f"Cannot convert {type(interval).__name__} to Interval")
