from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Generator, Literal, Protocol, Self, TypeAlias, Union

from ..spans import Span


@dataclass
class EnzymeConfig:
    """Configuration for a single enzyme in a digestion process."""

    enzyme_regex: str
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
    ) -> Self: ...

    def serialize(self) -> str:
        """Serialize the annotation to a ProForma string."""
        ...

    def __len__(self) -> int:
        """Return the length of the sequence."""
        ...

    def has_mods(self) -> bool:
        """Return True if the sequence has modifications."""
        ...


DigestResult: TypeAlias = Union[
    str,
    DigestProtocol,
    Span,
    tuple[str, Span],
    tuple[DigestProtocol, Span],
]

DigestGenerator: TypeAlias = Generator[DigestResult, None, None]

DigestReturnTypeLiterals = Literal[
    "str", "span", "annotation", "str-span", "annotation-span"
]
DigestReturnTypeStrLiteral = Literal["str"]
DigestReturnTypeSpanLiteral = Literal["span"]
DigestReturnTypeAnnotationLiteral = Literal["annotation"]
DigestReturnTypeStrSpanLiteral = Literal["str-span"]
DigestReturnTypeAnnotationSpanLiteral = Literal["annotation-span"]
