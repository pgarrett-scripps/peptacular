"""
UNIMOD module. Enum and Literals are not supported due to the scale of UNIMOD.
"""

from .data import UnimodInfo
from .lookup import UNIMOD_LOOKUP, UnimodLookup

__all__ = [
    "UNIMOD_LOOKUP",
    "UnimodLookup",
    "UnimodInfo",
]
