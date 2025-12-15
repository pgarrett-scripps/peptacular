from typing import overload

from peptacular.constants import ParrallelMethod, ParrallelMethodLiteral
from .parser import ProFormaParser, Interval
from .annotation import ProFormaAnnotation
from .serializer import serialize_annotation
from .mod import Mod, Mods
from .utils import Fragment


__all__ = [
    "ProFormaParser",
    "Interval",
    "ProFormaAnnotation",
    "serialize_annotation",
    "Mod",
    "Mods",
    "Fragment",
]
