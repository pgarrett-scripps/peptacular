from typing import Union, Dict, Tuple, List

from .dataclasses import Mod, Interval

ModIndex = Union[int, str]

# mod type
ModValue = Union[str, int, float, Mod]

# mod dict value type
ModDictValue = Union[ModValue, List[ModValue]]

# ModDict type
ModDict = Dict[ModIndex, ModDictValue]

# Span Type
Span = Tuple[int, int, int]

# Chem Composition Type
Chem_Composition = Dict[str, Union[int, float]]

# Fragments
ChargeType = Union[List[int], int]
IsotopeType = Union[List[int], int]
LossType = Union[List[float], float]
IonTypeType = Union[List[str], str]

ACCEPTED_MOD_INPUT = Union[List[ModValue], ModValue]
INTERVAL_TYPE = Union[Tuple[int, int, bool, ACCEPTED_MOD_INPUT], Interval]

