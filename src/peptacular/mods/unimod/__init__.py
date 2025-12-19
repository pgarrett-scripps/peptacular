"""
UNIMOD module. Enum and Literals are not supported due to the scale of UNIMOD.
"""

from .lookup import UNIMOD_LOOKUP, UnimodLookup
from .data import UnimodInfo

__all__ = [
    "UNIMOD_LOOKUP",
    "UnimodLookup",
    "UnimodInfo",
]
