from __future__ import annotations

from .constants import PROTEASES, PROTEASES_COMPILED
from .types import DigestReturnType, EnzymeConfig, DigestProtocol
from .mixin import DigestionMixin

__all__ = [
    "PROTEASES",
    "PROTEASES_COMPILED",
    "DigestReturnType",
    "EnzymeConfig",
    "DigestProtocol",
    "DigestionMixin",
]
