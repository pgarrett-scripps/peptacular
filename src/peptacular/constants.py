"""
Constants used throughout the peptacular package.
"""

import os
import re
from enum import StrEnum
from typing import Iterable, Literal, Final

from .element_setup import (
    get_element_info,
    get_isotopic_atomic_masses,
    map_atomic_number_to_comp,
    map_atomic_number_to_comp_neutron_offset,
    map_atomic_number_to_symbol,
    map_atomic_symbol_to_average_mass,
    map_hill_order,
)


# Partical masses
PROTON_MASS: Final[float] = 1.00727646688
ELECTRON_MASS: Final[float] = 0.00054857990946
NEUTRON_MASS: Final[float] = 1.00866491597
C13_NEUTRON_MASS: Final[float] = 1.003350
PEPTIDE_AVERAGINE_NEUTRON_MASS: Final[float] = 1.002856

_dir_name = os.path.dirname(__file__)
_element_path = os.path.join(_dir_name, "data", "chem.txt")

# Get element infos
_infos = get_element_info(_element_path)

# calculate average atomic mass for each element
AVERAGE_ATOMIC_MASSES: Final[dict[str, float]] = map_atomic_symbol_to_average_mass(_infos)
ISOTOPIC_ATOMIC_MASSES: Final[dict[str, float]] = get_isotopic_atomic_masses(_infos)
ATOMIC_NUMBER_TO_SYMBOL: Final[dict[int, str]] = map_atomic_number_to_symbol(_infos)
ATOMIC_SYMBOL_TO_NUMBER: Final[dict[str, int]] = {
    v: k for k, v in ATOMIC_NUMBER_TO_SYMBOL.items()
}
ATOMIC_SYMBOL_TO_ISOTOPE_NEUTRON_OFFSETS_AND_ABUNDANCES: Final[dict[
    str, list[tuple[int, float]]
]] = map_atomic_number_to_comp_neutron_offset(_infos)
ATOMIC_SYMBOL_TO_ISOTOPE_MASSES_AND_ABUNDANCES: Final[dict[str, list[tuple[float, float]]]] = (
    map_atomic_number_to_comp(_infos)
)
HILL_ORDER: Final[dict[str, int]] = map_hill_order(_infos)

# Neutral Mass
NTERM_COMPOSITION: Final[dict[str, int]] = {"H": 1}
CTERM_COMPOSITION: Final[dict[str, int]] = {"O": 1, "H": 1}


class IonType(StrEnum):
    PRECURSOR = "p"
    NEUTRAL = "n"
    A = "a"
    B = "b"
    C = "c"
    X = "x"
    Y = "y"
    Z = "z"
    Z_RADICAL = "z."     # z•+ radical ion (explicit radical notation)
    Z_PLUS_H = "z+H"     # z+1 ion (hydrogen abstraction product)
    C_MINUS_H = "c-H"    # c-H radical ion (complement to z+1)

    # Internal fragment ions
    AX = "ax"
    AY = "ay"
    AZ = "az"
    BX = "bx"
    BY = "by"
    BZ = "bz"
    CX = "cx"
    CY = "cy"
    CZ = "cz"
    IMMONIUM = "i" # single amino acid 'by' fragment ion
    
    # Satellite ions - d series (from a ions, partial side chain loss)
    D = "d"                      # Base d ion: C2H3N
    D_VALINE = "d-valine"        # Valine d+CH3
    DA_THREONINE = "da-threonine"  # Threonine da (heavier) +OH
    DA_ISOLEUCINE = "da-isoleucine"  # Isoleucine da (heavier) +C2H5
    DB_THREONINE = "db-threonine"  # Threonine db (lighter) +CH3
    DB_ISOLEUCINE = "db-isoleucine"  # Isoleucine db (lighter) +CH3
    
    # Satellite ions - v series (from y ions, complete side chain loss)
    V = "v"                      # Base v ion: C2H2NO
    
    # Satellite ions - w series (from z ions, partial side chain loss)
    W = "w"                      # Base w ion: C3H3O
    W_VALINE = "w-valine"        # Valine w+CH3
    WA_THREONINE = "wa-threonine"  # Threonine wa (heavier) +OH
    WA_ISOLEUCINE = "wa-isoleucine"  # Isoleucine wa (heavier) +C2H5
    WB_THREONINE = "wb-threonine"  # Threonine wb (lighter) +CH3
    WB_ISOLEUCINE = "wb-isoleucine"  # Isoleucine wb (lighter) +CH3


IonTypeLiteral = Literal[
    "p",
    "n",
    "a",
    "b",
    "c",
    "x",
    "y",
    "z",
    "z.",
    "z+H",
    "c-H",
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
    "d",
    "v",
    "w",
    "d-valine",
    "da-threonine",
    "da-isoleucine",
    "db-threonine",
    "db-isoleucine",
    "w-valine",
    "wa-threonine",
    "wa-isoleucine",
    "wb-threonine",
    "wb-isoleucine",
]


NEUTRAL_FRAGMENT_ION_COMPOSITIONS: Final[dict[IonType | str, dict[str, int]]] = {

    # Terminal fragments
    IonType.PRECURSOR: {"O": 1, "H": 2},
    IonType.NEUTRAL: {},
    IonType.A: {"C": -1, "O": -1},
    IonType.B: {},
    IonType.C: {"N": 1, "H": 3},
    IonType.X: {"C": 1, "O": 2},
    IonType.Y: {"O": 1, "H": 2},
    IonType.Z: {"N": -1, "O": 1, "H": -1},
    IonType.Z_RADICAL: {"N": -1, "O": 1},
    IonType.Z_PLUS_H: {"N": -1, "O": 1, "H": 1},
    IonType.C_MINUS_H: {"N": 1, "H": 2},
    IonType.IMMONIUM: {"O": -1, "C": -1},
    
    
    #The d ion is partial loss of the side chain from an a ion
    IonType.D: {"C": 2, "N": 1, "H": 3},
    IonType.D_VALINE: {"C": 3, "N": 1, "H": 6}, # Valine d+CH3
    IonType.DA_THREONINE: {"C": 2, "N": 1, "H": 4, "O": 1}, # Threonine d+OH,
    IonType.DA_ISOLEUCINE: {"C": 4, "N": 1, "H": 8}, # Isoleucine d+C2H5
    IonType.DB_THREONINE: {"C": 3, "N": 1, "H": 6}, # Threonine d+CH3
    IonType.DB_ISOLEUCINE: {"C": 3, "N": 1, "H": 6}, # Isoleucine d+CH3

    # The v ion is complete loss of the side chain of the y ion.
    IonType.V: {"C": 2, "N": 1, "H": 2, "O": 1},

    # The w ion is partial loss of the side chain from the z ion
    IonType.W: {"C": 3, "H": 3, "O": 1},
    IonType.W_VALINE: {"C": 4, "H": 6, "O": 1}, # Valine w+CH3
    IonType.WA_THREONINE: {"C": 3, "H": 4, "O": 2}, # Threonine w+OH,
    IonType.WA_ISOLEUCINE: {"C": 5, "H": 8, "O": 1}, # Isoleucine w+C2H5
    IonType.WB_THREONINE: {"C": 4, "H": 6, "O": 1}, # Threonine w+CH3
    IonType.WB_ISOLEUCINE: {"C": 4, "H": 6, "O": 1}, # Isoleucine w+CH3

    # Internal fragments
    IonType.BY: {"O": -1, "C": -1}, # Same as immonium
    IonType.AX: {"O": -1, "C": -1}, # Same as immonium
    IonType.CZ: {"O": -1, "C": -1}, # Same as immonium
    IonType.AY: {"O": -2, "C": -2}, # immonium-CO
    IonType.AZ: {"O": -2, "C": -2, "N": -1, "H": -1}, # immonium-CHNO
    IonType.BX: {}, # immonium+CO
    IonType.BZ: {"O": -1, "C": -1, "H": -1, "N": -1}, # immonium-NH
    IonType.CX: {"N": 1}, # immonium+CHNO
    IonType.CY: {"O": -1, "C": -1, "N": 1, "H": -1}, # immonium+NH
}

AA_COMPOSITIONS: Final[dict[str, dict[str, int]]] = {
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

AMINO_ACIDS: Final[frozenset[str]] = frozenset(AA_COMPOSITIONS.keys()) | {"B", "Z"}
ORDERED_AMINO_ACIDS: Final[list[str]] = sorted(list(AMINO_ACIDS))
AMBIGUOUS_AMINO_ACIDS: Final[frozenset[str]] = frozenset({"J", "B", "Z", "X"})
MASS_AMBIGUOUS_AMINO_ACIDS: Final[frozenset[str]] = frozenset({"B", "Z"})
UNAMBIGUOUS_AMINO_ACIDS: Final[frozenset[str]] = AMINO_ACIDS - AMBIGUOUS_AMINO_ACIDS

# 3 letter codes
AA_TO_THREE_LETTER_CODE: Final[dict[str, str]] = {
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
THREE_LETTER_CODE_TO_AA: Final[dict[str, str]] = {
    v: k for k, v in AA_TO_THREE_LETTER_CODE.items()
}

# Amino acid to name
AA_TO_NAME: Final[dict[str, str]] = {
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
NAME_TO_AA: Final[dict[str, str]] = {v: k for k, v in AA_TO_NAME.items()}

FORWARD_ION_TYPES: Final[frozenset[str]] = frozenset([IonType.A, IonType.B, IonType.C])
BACKWARD_ION_TYPES: Final[frozenset[str]] = frozenset([IonType.X, IonType.Y, IonType.Z])
INTERNAL_ION_TYPES: Final[frozenset[str]] = frozenset(
    [
        IonType.AX,
        IonType.AY,
        IonType.AZ,
        IonType.BX,
        IonType.BY,
        IonType.BZ,
        IonType.CX,
        IonType.CY,
        IonType.CZ,
    ]
)
TERMINAL_ION_TYPES: Final[frozenset[str]] = FORWARD_ION_TYPES | BACKWARD_ION_TYPES
IMMONIUM_ION_TYPES: Final[frozenset[str]] = frozenset([IonType.IMMONIUM])
VALID_ION_TYPES: Final[frozenset[str]] = (
    TERMINAL_ION_TYPES | INTERNAL_ION_TYPES | IMMONIUM_ION_TYPES
)
AVERAGINE_RATIOS: Final[dict[str, float]] = {
    "C": 4.9384,
    "H": 7.7583,
    "N": 1.3577,
    "O": 1.4773,
    "S": 0.0417,
}


# Compiling regex patterns used in your module
ISOTOPE_COMPONENT_PATTERN: Final[re.Pattern[str]] = re.compile(r"([0-9]*)([A-Za-z]+)(-?\d*\.?\d*)")
CONDENSED_CHEM_FORMULA_PATTERN: Final[re.Pattern[str]] = re.compile(r"([A-Z][a-z]*|e|p|n)(-?\d*\.?\d*)")
ADDUCT_PATTERN: Final[re.Pattern[str]] = re.compile(r"([+-]?)(\d*)?([A-Za-z]{1,2}\d*\+?-?)")
ISOTOPE_NUM_PATTERN: Final[re.Pattern[str]] = re.compile(r"[0-9]")


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


class ParrallelMethod(StrEnum):
    PROCESS = "process"
    THREAD = "thread"
    SEQUENTIAL = "sequential"


ParrallelMethodLiteral = Literal["process", "thread", "sequential"]
