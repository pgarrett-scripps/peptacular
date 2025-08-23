"""Protein property calculation module."""

from .types import (
    MissingAAHandling,
    AggregationMethod,
    SequenceProtocol,
    WEIGHTING_SCHEMES,
)

from .mixin import SequencePropertyMixin

__all__ = [
    "MissingAAHandling",
    "AggregationMethod",
    "SequenceProtocol",
    "WEIGHTING_SCHEMES",
    "SequencePropertyMixin",
]