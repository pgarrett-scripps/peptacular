"""Type definitions and enums for protein property calculations."""

from __future__ import annotations
from typing import Protocol, Union, Literal
from enum import Enum


class MissingAAHandling(str, Enum):
    """Strategy for handling missing amino acid values"""

    ZERO = "zero"
    AVG = "avg"
    MIN = "min"
    MAX = "max"
    MEDIAN = "median"
    ERROR = "error"
    SKIP = "skip"


class AggregationMethod(str, Enum):
    """Strategy for aggregating amino acid values"""

    SUM = "sum"
    AVG = "avg"


class SequenceProtocol(Protocol):
    """Protocol defining the interface for objects with sequences"""

    @property
    def stripped_sequence(self) -> str:
        """The sequence without modifications (for property calculations)"""
        ...


# Type aliases
WEIGHTING_SCHEMES = Union[
    list[float],
    Literal[
        "uniform",
        "linear",
        "exponential",
        "gaussian",
        "sigmoid",
        "cosine",
        "sinusoidal",
    ],
]