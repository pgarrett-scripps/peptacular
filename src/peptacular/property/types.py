"""Type definitions and enums for protein property calculations."""

from __future__ import annotations
from enum import StrEnum
from typing import Protocol


class MissingAAHandling(StrEnum):
    """Strategy for handling missing amino acid values"""

    ZERO = "zero"
    AVG = "avg"
    MIN = "min"
    MAX = "max"
    MEDIAN = "median"
    ERROR = "error"
    SKIP = "skip"


class AggregationMethod(StrEnum):
    """Strategy for aggregating amino acid values"""

    SUM = "sum"
    AVG = "avg"


class SequenceProtocol(Protocol):
    """Protocol defining the interface for objects with sequences"""

    @property
    def stripped_sequence(self) -> str:
        """The sequence without modifications (for property calculations)"""
        ...


class WeightingMethods(StrEnum):
    """Enumeration of possible weighting methods for amino acid properties"""

    UNIFORM = "uniform"
    LINEAR = "linear"
    EXPONENTIAL = "exponential"
    GAUSSIAN = "gaussian"
    SIGMOID = "sigmoid"
    COSINE = "cosine"
    SINUSOIDAL = "sinusoidal"
