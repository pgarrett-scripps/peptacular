"""
The annotation module contains classes and functions for handling peptide annotations.

All methods/properties are available on the ProFormaAnnotation class. The one exception is the AnnotationProperties,
which is accessible via the annot.prop attribute.
"""

from .parser import Interval
from .annotation import ProFormaAnnotation, AnnotationProperties
from .mod import Mod, Mods
from .utils import Fragment


__all__ = [
    "Interval",
    "ProFormaAnnotation",
    "AnnotationProperties",
    "Mod",
    "Mods",
    "Fragment",
]
