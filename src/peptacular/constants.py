"""
Constants used throughout the peptacular package.
"""

import os
from typing import Iterable, Literal
from enum import StrEnum

import regex as re

from .element_setup import (
    get_element_info,
    map_atomic_symbol_to_average_mass,
    map_atomic_number_to_comp_neutron_offset,
    get_isotopic_atomic_masses,
    map_hill_order,
    map_atomic_number_to_comp,
    map_atomic_number_to_symbol,
)


def _merge_dicts(d1: dict[str, int], d2: dict[str, int]) -> dict[str, int]:
    """
    Merge two dictionaries. And remove any keys with value 0.
    """
    d: dict[str, int] = {}
    for k, v in d1.items():
        d[k] = v
    for k, v in d2.items():
        d[k] = d.get(k, 0) + v

    # remove any keys with value 0
    d = {k: v for k, v in d.items() if v != 0}

    return d


# Partical masses
PROTON_MASS: float = 1.00727646688
ELECTRON_MASS: float = 0.00054857990946
NEUTRON_MASS: float = 1.00866491597
C13_NEUTRON_MASS: float = 1.003350
PEPTIDE_AVERAGINE_NEUTRON_MASS: float = 1.002856

_dir_name = os.path.dirname(__file__)
_element_path = os.path.join(_dir_name, "data", "chem.txt")

# Get element infos
_infos = get_element_info(_element_path)

# calculate average atomic mass for each element
AVERAGE_ATOMIC_MASSES: dict[str, float] = map_atomic_symbol_to_average_mass(_infos)
ISOTOPIC_ATOMIC_MASSES: dict[str, float] = get_isotopic_atomic_masses(_infos)
ATOMIC_NUMBER_TO_SYMBOL: dict[int, str] = map_atomic_number_to_symbol(_infos)
ATOMIC_SYMBOL_TO_NUMBER: dict[str, int] = {
    v: k for k, v in ATOMIC_NUMBER_TO_SYMBOL.items()
}
ATOMIC_SYMBOL_TO_ISOTOPE_NEUTRON_OFFSETS_AND_ABUNDANCES: dict[
    str, list[tuple[int, float]]
] = map_atomic_number_to_comp_neutron_offset(_infos)
ATOMIC_SYMBOL_TO_ISOTOPE_MASSES_AND_ABUNDANCES: dict[str, list[tuple[float, float]]] = (
    map_atomic_number_to_comp(_infos)
)
HILL_ORDER: dict[str, int] = map_hill_order(_infos)

# Neutral Mass
NTERM_COMPOSITION = {"H": 1}
CTERM_COMPOSITION = {"O": 1, "H": 1}


class IonType(StrEnum):
    PRECURSOR = "p"
    NEUTRAL = "n"
    A = "a"
    B = "b"
    C = "c"
    X = "x"
    Y = "y"
    Z = "z"
    AX = "ax"
    AY = "ay"
    AZ = "az"
    BX = "bx"
    BY = "by"
    BZ = "bz"
    CX = "cx"
    CY = "cy"
    CZ = "cz"
    IMMONIUM = "i"


IonTypeLiteral = Literal[
    "p",
    "n",
    "a",
    "b",
    "c",
    "x",
    "y",
    "z",
    "ax",
    "ay",
    "az",
    "bx",
    "by",
    "bz",
    "cx",
    "cy",
    "cz",
    "i",
]


# Additions to get the +1 ion
FRAGMENT_ION_COMPOSITIONS: dict[IonType | str, dict[str, int]] = {
    IonType.PRECURSOR: {"H": 1, "e": -1},  # Proton
    IonType.NEUTRAL: {},  # Nothing
    IonType.A: {"e": -1},  # -electron
    IonType.B: {"e": -1},  # -electron
    IonType.C: {"H": 2, "e": -1},  # Steals proton and hydrogen from neutral fragment
    IonType.X: {"e": -1},  # -electron
    IonType.Y: {"H": 2, "e": -1},  # Steals proton and hydrogen from neutral fragment
    IonType.Z: {"e": -1},
    IonType.AX: {"e": -1},
    IonType.AY: {"H": 1, "e": -1},
    IonType.AZ: {"e": -1},
    IonType.BX: {"e": -1},
    IonType.BY: {"H": 1, "e": -1},
    IonType.BZ: {"e": -1},
    IonType.CX: {"H": 1, "e": -1},
    IonType.CY: {"H": 3, "e": -1},
    IonType.CZ: {"H": 1, "e": -1},
    IonType.IMMONIUM: {"H": 1, "e": -1},
}

FRAGMENT_ION_BASE_CHARGE_ADDUCTS: dict[IonType | str, str] = {
    IonType.PRECURSOR: "+H+",
    IonType.NEUTRAL: "",
    IonType.A: "-e-",
    IonType.B: "-e-",
    IonType.C: "+2H+,+e-",
    IonType.X: "-e-",
    IonType.Y: "+2H+,+e-",
    IonType.Z: "-e-",
    IonType.AX: "-e-",
    IonType.AY: "+H+",
    IonType.AZ: "-e-",
    IonType.BX: "-e-",
    IonType.BY: "+H+",
    IonType.BZ: "-e-",
    IonType.CX: "+H+",
    IonType.CY: "+3H+,+e-",
    IonType.CZ: "+H+",
    IonType.IMMONIUM: "+H+",
}

# According to the dissociation sites
NEUTRAL_FRAGMENT_START_COMPOSITIONS: dict[IonType | str, dict[str, int]] = {
    IonType.PRECURSOR: NTERM_COMPOSITION,
    IonType.NEUTRAL: NTERM_COMPOSITION,
    IonType.A: NTERM_COMPOSITION,
    IonType.B: NTERM_COMPOSITION,
    IonType.C: NTERM_COMPOSITION,
    IonType.X: {"O": 1, "C": 1},
    IonType.Y: {},
    IonType.Z: {"N": -1, "H": -1},
    IonType.IMMONIUM: {},
    IonType.NEUTRAL: {},
}

# According to the dissociation sites
NEUTRAL_FRAGMENT_END_COMPOSITIONS: dict[IonType | str, dict[str, int]] = {
    IonType.PRECURSOR: CTERM_COMPOSITION,
    IonType.A: {"O": -1, "C": -1},
    IonType.B: {},
    IonType.C: {"N": 1, "H": 1},
    IonType.X: CTERM_COMPOSITION,
    IonType.Y: CTERM_COMPOSITION,
    IonType.Z: CTERM_COMPOSITION,
    IonType.IMMONIUM: {"O": -1, "C": -1},
    IonType.NEUTRAL: {},
}

NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS: dict[str | str, dict[str, int]] = {
    IonType.PRECURSOR: _merge_dicts(
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.PRECURSOR],
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.PRECURSOR],
    ),
    IonType.NEUTRAL: _merge_dicts(
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.NEUTRAL],
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.NEUTRAL],
    ),
    IonType.A: _merge_dicts(
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.A],
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.A],
    ),
    IonType.B: _merge_dicts(
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.B],
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.B],
    ),
    IonType.C: _merge_dicts(
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.C],
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.C],
    ),
    IonType.X: _merge_dicts(
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.X],
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.X],
    ),
    IonType.Y: _merge_dicts(
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.Y],
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.Y],
    ),
    IonType.Z: _merge_dicts(
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.Z],
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.Z],
    ),
    IonType.AX: _merge_dicts(
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.A],
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.X],
    ),
    IonType.AY: _merge_dicts(
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.A],
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.Y],
    ),
    IonType.AZ: _merge_dicts(
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.A],
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.Z],
    ),
    IonType.BX: _merge_dicts(
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.B],
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.X],
    ),
    IonType.BY: _merge_dicts(
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.B],
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.Y],
    ),
    IonType.BZ: _merge_dicts(
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.B],
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.Z],
    ),
    IonType.CX: _merge_dicts(
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.C],
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.X],
    ),
    IonType.CY: _merge_dicts(
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.C],
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.Y],
    ),
    IonType.CZ: _merge_dicts(
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.C],
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.Z],
    ),
    IonType.IMMONIUM: _merge_dicts(
        NEUTRAL_FRAGMENT_START_COMPOSITIONS[IonType.IMMONIUM],
        NEUTRAL_FRAGMENT_END_COMPOSITIONS[IonType.IMMONIUM],
    ),
}

FRAGMENT_ION_COMPOSITION_ADJUSTMENTS: dict[IonType | str, dict[str, int]] = {
    IonType.PRECURSOR: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.PRECURSOR],
        FRAGMENT_ION_COMPOSITIONS[IonType.PRECURSOR],
    ),
    IonType.NEUTRAL: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.NEUTRAL],
        FRAGMENT_ION_COMPOSITIONS[IonType.NEUTRAL],
    ),
    IonType.A: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.A],
        FRAGMENT_ION_COMPOSITIONS[IonType.A],
    ),
    IonType.B: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.B],
        FRAGMENT_ION_COMPOSITIONS[IonType.B],
    ),
    IonType.C: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.C],
        FRAGMENT_ION_COMPOSITIONS[IonType.C],
    ),
    IonType.X: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.X],
        FRAGMENT_ION_COMPOSITIONS[IonType.X],
    ),
    IonType.Y: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.Y],
        FRAGMENT_ION_COMPOSITIONS[IonType.Y],
    ),
    IonType.Z: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.Z],
        FRAGMENT_ION_COMPOSITIONS[IonType.Z],
    ),
    IonType.AX: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.AX],
        FRAGMENT_ION_COMPOSITIONS[IonType.AX],
    ),
    IonType.AY: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.AY],
        FRAGMENT_ION_COMPOSITIONS[IonType.AY],
    ),
    IonType.AZ: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.AZ],
        FRAGMENT_ION_COMPOSITIONS[IonType.AZ],
    ),
    IonType.BX: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.BX],
        FRAGMENT_ION_COMPOSITIONS[IonType.BX],
    ),
    IonType.BY: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.BY],
        FRAGMENT_ION_COMPOSITIONS[IonType.BY],
    ),
    IonType.BZ: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.BZ],
        FRAGMENT_ION_COMPOSITIONS[IonType.BZ],
    ),
    IonType.CX: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.CX],
        FRAGMENT_ION_COMPOSITIONS[IonType.CX],
    ),
    IonType.CY: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.CY],
        FRAGMENT_ION_COMPOSITIONS[IonType.CY],
    ),
    IonType.CZ: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.CZ],
        FRAGMENT_ION_COMPOSITIONS[IonType.CZ],
    ),
    IonType.IMMONIUM: _merge_dicts(
        NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[IonType.IMMONIUM],
        FRAGMENT_ION_COMPOSITIONS[IonType.IMMONIUM],
    ),
}

AA_COMPOSITIONS: dict[str, dict[str, int]] = {
    "G": {"C": 2, "H": 3, "N": 1, "O": 1},  # Glycine
    "A": {"C": 3, "H": 5, "N": 1, "O": 1},  # Alanine
    "S": {"C": 3, "H": 5, "N": 1, "O": 2},  # Serine
    "P": {"C": 5, "H": 7, "N": 1, "O": 1},  # Proline
    "V": {"C": 5, "H": 9, "N": 1, "O": 1},  # Valine
    "T": {"C": 4, "H": 7, "N": 1, "O": 2},  # Threonine
    "C": {"C": 3, "H": 5, "N": 1, "O": 1, "S": 1},  # Cysteine
    "I": {"C": 6, "H": 11, "N": 1, "O": 1},  # Isoleucine
    "L": {"C": 6, "H": 11, "N": 1, "O": 1},  # Leucine
    "J": {"C": 6, "H": 11, "N": 1, "O": 1},  # Leucine or Isoleucine
    "N": {"C": 4, "H": 6, "N": 2, "O": 2},  # Asparagine
    "D": {"C": 4, "H": 5, "N": 1, "O": 3},  # Aspartic acid
    "Q": {"C": 5, "H": 8, "N": 2, "O": 2},  # Glutamine
    "K": {"C": 6, "H": 12, "N": 2, "O": 1},  # Lysine
    "E": {"C": 5, "H": 7, "N": 1, "O": 3},  # Glutamic acid
    "M": {"C": 5, "H": 9, "N": 1, "O": 1, "S": 1},  # Methionine
    "H": {"C": 6, "H": 7, "N": 3, "O": 1},  # Histidine
    "F": {"C": 9, "H": 9, "N": 1, "O": 1},  # Phenylalanine
    "R": {"C": 6, "H": 12, "N": 4, "O": 1},  # Arginine
    "Y": {"C": 9, "H": 9, "N": 1, "O": 2},  # Tyrosine
    "W": {"C": 11, "H": 10, "N": 2, "O": 1},  # Tryptophan
    "U": {"C": 3, "H": 5, "N": 1, "O": 1, "Se": 1},  # Selenocysteine
    "O": {"C": 12, "H": 19, "N": 3, "O": 2},  # Pyrrolysine
    "X": {},  # Unknown amino acid (treated as 0 mass)
}

AMINO_ACIDS: set[str] = set(AA_COMPOSITIONS.keys()) | {"B", "Z"}
ORDERED_AMINO_ACIDS = sorted(list(AMINO_ACIDS))
AMBIGUOUS_AMINO_ACIDS: set[str] = {"J", "B", "Z", "X"}
MASS_AMBIGUOUS_AMINO_ACIDS: set[str] = {"B", "Z"}
UNAMBIGUOUS_AMINO_ACIDS: set[str] = AMINO_ACIDS - AMBIGUOUS_AMINO_ACIDS

# 3 letter codes
AA_TO_THREE_LETTER_CODE: dict[str, str] = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "E": "Glu",
    "Q": "Gln",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "B": "Asx",
    "Z": "Glx",
    "X": "Xaa",
    "U": "Sec",
    "O": "Pyl",
    "J": "Xle",
}

# Reverse 3 letter codes
THREE_LETTER_CODE_TO_AA: dict[str, str] = {
    v: k for k, v in AA_TO_THREE_LETTER_CODE.items()
}

# Amino acid to name
AA_TO_NAME: dict[str, str] = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic acid",
    "C": "Cysteine",
    "E": "Glutamic acid",
    "Q": "Glutamine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
    "B": "Asparagine or Aspartic acid",
    "Z": "Glutamine or Glutamic acid",
    "X": "Any amino acid",
    "U": "Selenocysteine",
    "O": "Pyrrolysine",
    "J": "Leucine or Isoleucine",
}

# Reverse amino acid to name
NAME_TO_AA: dict[str, str] = {v: k for k, v in AA_TO_NAME.items()}

FORWARD_ION_TYPES: set[str] = {IonType.A, IonType.B, IonType.C}
BACKWARD_ION_TYPES: set[str] = {IonType.X, IonType.Y, IonType.Z}
INTERNAL_ION_TYPES: set[str] = {
    IonType.AX,
    IonType.AY,
    IonType.AZ,
    IonType.BX,
    IonType.BY,
    IonType.BZ,
    IonType.CX,
    IonType.CY,
    IonType.CZ,
}
TERMINAL_ION_TYPES: set[str] = FORWARD_ION_TYPES | BACKWARD_ION_TYPES
IMMONIUM_ION_TYPES: set[str] = {IonType.IMMONIUM}
VALID_ION_TYPES: set[str] = TERMINAL_ION_TYPES | INTERNAL_ION_TYPES | IMMONIUM_ION_TYPES
AVERAGINE_RATIOS: dict[str, float] = {
    "C": 4.9384,
    "H": 7.7583,
    "N": 1.3577,
    "O": 1.4773,
    "S": 0.0417,
}


# Compiling regex patterns used in your module
ISOTOPE_COMPONENT_PATTERN = re.compile(r"([0-9]*)([A-Za-z]+)(-?\d*\.?\d*)")
CONDENSED_CHEM_FORMULA_PATTERN = re.compile(r"([A-Z][a-z]*|e|p|n)(-?\d*\.?\d*)")
ADDUCT_PATTERN = re.compile(r"([+-]?)(\d*)?([A-Za-z]{1,2}\d*\+?-?)")
ISOTOPE_NUM_PATTERN = re.compile(r"[0-9]")


# str enum
class ModType(StrEnum):
    NTERM = "nterm"
    CTERM = "cterm"
    ISOTOPE = "isotope"
    STATIC = "static"
    LABILE = "labile"
    UNKNOWN = "unknown"
    INTERVAL = "interval"
    INTERNAL = "internal"
    CHARGE = "charge"
    CHARGE_ADDUCTS = "charge_adducts"


ModTypeLiteral = Literal[
    "nterm",
    "cterm",
    "isotope",
    "static",
    "labile",
    "unknown",
    "interval",
    "internal",
    "charge",
    "charge_adducts",
]


def get_mod_type(mod: ModTypeLiteral | ModType) -> ModType:
    # return ModType Enum for the given mod string
    if isinstance(mod, ModType):
        return mod

    if isinstance(mod, str):  # type: ignore
        for mod_type in ModType:
            if mod_type.value == mod:
                return mod_type
    else:
        raise ValueError(f"mod must be a string or ModType, got {type(mod)}")

    raise ValueError(f"Unknown mod type: {mod}")


def get_mods(
    mods: Iterable[ModTypeLiteral | ModType] | ModType | ModTypeLiteral | None,
) -> list[ModType]:
    """
    Get the list of modification types from the input.

    :param mods: Modification types as a ModType, iterable of ModTypes, or None.
    :type mods: None | ModType | Iterable[ModType]
    :return: List of ModType Enum values.
    :rtype: list[ModType]
    :raises ValueError: If mods is not None, ModType, or iterable of ModTypes.
    """

    if mods is None:
        return [mod_type for mod_type in ModType]
    elif isinstance(mods, (str, ModType)):
        # Single modification type
        return [get_mod_type(mods)]
    elif isinstance(mods, Iterable):  # type: ignore
        # List of modification types
        return [get_mod_type(mod) for mod in mods]

    raise ValueError(
        f"mods parameter must be str, list of str, or None, got {type(mods)}"
    )
