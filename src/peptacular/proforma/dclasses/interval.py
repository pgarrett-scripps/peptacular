from __future__ import annotations
import copy
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any, Self

from ...mod import MOD_VALUE_TYPES
from .modlist import ModList, setup_mod_list, MODLIST_DATATYPE


@dataclass
class ModInterval:
    """An immutable sequence interval for return type"""

    start: int
    end: int
    ambiguous: bool
    mods: tuple[MOD_VALUE_TYPES, ...] = field(default_factory=tuple)


class Interval:
    """A sequence interval with modlist"""

    def __init__(
        self,
        start: int,
        end: int,
        ambiguous: bool,
        mods: ModList | Iterable[MODLIST_DATATYPE] | MODLIST_DATATYPE | None = None,
    ) -> None:
        self.start: int = start
        self.end: int = end
        self.ambiguous: bool = ambiguous
        self.mods: ModList = (
            setup_mod_list(mods) if not isinstance(mods, ModList) else mods
        )

    @property
    def has_mods(self) -> bool:
        """Check if the interval has modifications"""
        return len(self.mods) > 0

    def dict(self) -> dict[str, Any]:
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


ACCEPTED_INTERVAL_DATATYPE = Interval | ModInterval


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
    else:
        raise TypeError(f"Cannot convert {type(interval).__name__} to Interval")
