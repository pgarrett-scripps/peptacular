from ...mod import MOD_VALUE_TYPES, Mod, setup_mod
from .interval import ACCEPTED_INTERVAL_DATATYPE, Interval, ModInterval, setup_interval
from .intervallist import (
    ACCEPTED_INTERVALLIST_INPUT_TYPES,
    INTERVALLIST_DATATYPE,
    IntervalList,
    setup_interval_list,
    populate_interval_list
)
from .moddict import ACCEPTED_MODDICT_INPUT_TYPES, ModDict, setup_mod_dict, populate_mod_dict
from .modlist import (
    ACCEPTED_MODLIST_INPUT_TYPES,
    MODLIST_DATATYPE,
    ModList,
    setup_mod_list,
    populate_mod_list,
)

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
    "populate_mod_list",
    "populate_interval_list",
    "populate_mod_dict",
]
