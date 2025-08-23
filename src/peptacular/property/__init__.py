"""Protein property calculation module."""

from .types import (
    SequenceProtocol,
)

from .mixin import SequencePropertyMixin

__all__ = [
    "SequenceProtocol",
    "SequencePropertyMixin",
]
