from typing import Union, Dict, Tuple, List

from .proforma_dataclasses import Interval, Mod


# Span Type
Span = Tuple[int, int, int]

# Chem Composition Type
ChemComposition = Dict[str, Union[int, float]]


ModIndex = Union[int, str]
ModValue = Union[str, int, float, Mod]
ModDictValue = Union[List[Mod], List[Interval], int]
ModDict = Dict[ModIndex, ModDictValue]

ACCEPTED_MOD_INPUT = Union[List[ModValue], ModValue]
INTERVAL_VALUE = Union[Tuple[int, int, bool, Union[ACCEPTED_MOD_INPUT, None]], Interval]
ACCEPTED_INTERVAL_INPUT = Union[List[INTERVAL_VALUE], INTERVAL_VALUE]
