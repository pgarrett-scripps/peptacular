"""Protein property calculation module."""

from .mixin import SequencePropertyMixin
from .types import (
    SequenceProtocol,
)

__all__ = [
    "SequenceProtocol",
    "SequencePropertyMixin",
]
