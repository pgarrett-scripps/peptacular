"""Protein property calculation module."""

from .prop import AnnotationProperties
from .data import (
    SecondaryStructureMethod,
    SecondaryStructureType,
    PropertyScale,
    HydrophobicityScale,
    SecondaryStructureScale,
    SurfaceAccessibilityScale,
    ChargeScale,
    PolarityScale,
    HPLCScale,
    BetaStrandScale,
    PhysicalPropertyScale,
    CompositionScale,
    PROPERTY_SCALES,
    HYDROPHOBICITY_SCALES,
    SURFACE_ACCESSIBILITY_SCALES,
    HPLC_SCALES,
    HYDROPHILICITY_SCALES,
    FLIXIBILITY_SCALES,
    POLARITY_SCALES,
    COMPOSITION_SCALES,
    PHYSICAL_PROPERTY_SCALES,
)
from .types import (
    MissingAAHandling,
    MissingAAHandlingLiteral,
    AggregationMethod,
    AggregationMethodLiteral,
    WeightingMethods,
    WeightingMethodsLiteral,
)

__all__ = [
    # main class
    "AnnotationProperties",
    # data
    "SecondaryStructureMethod",
    "SecondaryStructureType",
    "PropertyScale",
    "HydrophobicityScale",
    "SecondaryStructureScale",
    "SurfaceAccessibilityScale",
    "ChargeScale",
    "PolarityScale",
    "HPLCScale",
    "BetaStrandScale",
    "PhysicalPropertyScale",
    "CompositionScale",
    "PROPERTY_SCALES",
    "HYDROPHOBICITY_SCALES",
    "SURFACE_ACCESSIBILITY_SCALES",
    "HPLC_SCALES",
    "HYDROPHILICITY_SCALES",
    "FLIXIBILITY_SCALES",
    "POLARITY_SCALES",
    "COMPOSITION_SCALES",
    "PHYSICAL_PROPERTY_SCALES",
    # types
    "MissingAAHandling",
    "MissingAAHandlingLiteral",
    "AggregationMethod",
    "AggregationMethodLiteral",
    "WeightingMethods",
    "WeightingMethodsLiteral",
]
