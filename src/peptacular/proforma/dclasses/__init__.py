from ...mod import Mod, MOD_VALUE_TYPES, setup_mod
from .interval import Interval, ModInterval, setup_interval, ACCEPTED_INTERVAL_DATATYPE
from .modlist import (
    ModList,
    setup_mod_list,
    MODLIST_DATATYPE,
    ACCEPTED_MODLIST_INPUT_TYPES,
)
from .intervallist import (
    IntervalList,
    setup_interval_list,
    INTERVALLIST_DATATYPE,
    ACCEPTED_INTERVALLIST_INPUT_TYPES,
)
from .moddict import ModDict, setup_mod_dict, ACCEPTED_MODDICT_INPUT_TYPES

SPAN_TYPE = tuple[int, int, int]

# Chem Composition Type
CHEM_COMPOSITION_TYPE = dict[str, int | float]

__all__ = [
    "Mod",
    "MOD_VALUE_TYPES",
    "setup_mod",
    "Interval",
    "ModInterval",
    "setup_interval",
    "ACCEPTED_INTERVAL_DATATYPE",
    "ModList",
    "setup_mod_list",
    "ACCEPTED_MODLIST_INPUT_TYPES",
    "MODLIST_DATATYPE",
    "IntervalList",
    "setup_interval_list",
    "ACCEPTED_INTERVALLIST_INPUT_TYPES",
    "INTERVALLIST_DATATYPE",
    "ModDict",
    "setup_mod_dict",
    "ACCEPTED_MODDICT_INPUT_TYPES",
    "SPAN_TYPE",
    "CHEM_COMPOSITION_TYPE",
]
