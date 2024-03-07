# mod index type
from typing import Union, Dict, Tuple, List

ModIndex = Union[int, str]

# mod type
ModValue = Union[str, int, float]

# mod dict value type
ModDictValue = Union[ModValue, List[ModValue]]

# ModDict type
ModDict = Dict[ModIndex, ModDictValue]

# Span Type
Span = Tuple[int, int, int]

# Chem Composition Type
Chem_Composition = Dict[str, Union[int, float]]