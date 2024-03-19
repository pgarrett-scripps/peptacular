from typing import Union, Dict, Tuple, List

from peptacular.proforma_dataclasses import Mod, Interval

ModIndex = Union[int, str]

# mod type
ModValue = Union[str, int, float, Mod]

# mod dict value type (Mod dicts can have Mods, Intervals, and charge state
ModDictValue = Union[List[Mod], List[Interval], int]

# ModDict type
ModDict = Dict[ModIndex, ModDictValue]

# Span Type
Span = Tuple[int, int, int]

# Chem Composition Type
ChemComposition = Dict[str, Union[int, float]]

# Fragments
ChargeType = Union[List[int], int]
IsotopeType = Union[List[int], int]
LossType = Union[List[float], float]
IonTypeType = Union[List[str], str]

ACCEPTED_MOD_INPUT = Union[List[ModValue], ModValue]
IntervalValue = Union[Tuple[int, int, bool, Union[ACCEPTED_MOD_INPUT, None]], Interval]
ACCEPTED_INTERVAL_INPUT = Union[List[IntervalValue], IntervalValue]

