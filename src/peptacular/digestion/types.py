from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Generator, Iterable, Literal, Protocol, Self, TypeAlias, Union


@dataclass
class EnzymeConfig:
    """Configuration for a single enzyme in a digestion process."""
    
    regex: Iterable[str] | str
    missed_cleavages: int = 0
    semi_enzymatic: bool = False
    complete_digestion: bool = True


class DigestReturnType(StrEnum):
    STR = "str"
    ANNOTATION = "annotation"
    SPAN = "span"
    STR_SPAN = "str-span"
    ANNOTATION_SPAN = "annotation-span"



class DigestProtocol(Protocol):
    """Protocol defining the interface for objects that can be digested."""

    @property
    def stripped_sequence(self) -> str:
        """The sequence without modifications (for calculations)."""
        ...

    def slice(
        self, 
        start: int | None, 
        stop: int | None, 
        inplace: bool = False, 
        keep_terms: bool = False, 
        keep_labile: bool = True
    ) -> Self: ...

    def serialize(self) -> str:
        """Serialize the annotation to a ProForma string."""
        ...

    def __len__(self) -> int:
        """Return the length of the sequence."""
        ...


DigestResult: TypeAlias = Union[
    str,
    DigestProtocol,
    tuple[int, int, int],
    tuple[str, tuple[int, int, int]],
    tuple[DigestProtocol, tuple[int, int, int]]
]

DigestGenerator: TypeAlias = Generator[DigestResult, None, None]

ReturnTypeLiteral = Literal[
   "str", "span", "annotation", "str-span", "annotation-span"
]