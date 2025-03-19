"""
Constants used throughout the peptacular package.
"""

import os
from typing import Dict, List, Tuple, Set

from regex import regex

from peptacular.element_setup import get_element_info, map_atomic_symbol_to_average_mass, \
    map_atomic_number_to_comp_neutron_offset, get_isotopic_atomic_masses, map_hill_order, map_atomic_number_to_comp, \
    map_atomic_number_to_symbol
from peptacular.util import merge_dicts

# Partical masses
PROTON_MASS = 1.00727646688
NEUTRON_MASS = 1.00866491597
ELECTRON_MASS = 0.00054857990946

_dir_name = os.path.dirname(__file__)
_element_path = os.path.join(_dir_name, "data", "chem.txt")

# Get element infos
_infos = get_element_info(_element_path)

# calculate average atomic mass for each element
AVERAGE_ATOMIC_MASSES: Dict[str, float] = map_atomic_symbol_to_average_mass(_infos)
ISOTOPIC_ATOMIC_MASSES: Dict[str, float] = get_isotopic_atomic_masses(_infos)
ATOMIC_NUMBER_TO_SYMBOL: Dict[int, str] = map_atomic_number_to_symbol(_infos)
ATOMIC_SYMBOL_TO_NUMBER: Dict[str, int] = {v: k for k, v in ATOMIC_NUMBER_TO_SYMBOL.items()}
ATOMIC_SYMBOL_TO_ISOTOPE_NEUTRON_OFFSETS_AND_ABUNDANCES: Dict[str, List[Tuple[int, float]]] = \
    map_atomic_number_to_comp_neutron_offset(_infos)
ATOMIC_SYMBOL_TO_ISOTOPE_MASSES_AND_ABUNDANCES: Dict[str, List[Tuple[float, float]]] = \
    map_atomic_number_to_comp(_infos)
HILL_ORDER: Dict[str, int] = map_hill_order(_infos)

# Neutral Mass
NTERM_COMPOSITION = {"H": 1}
CTERM_COMPOSITION = {"O": 1, "H": 1}

# Additions to get the +1 ion
FRAGMENT_ION_COMPOSITIONS: Dict[str, Dict[str, int]] = {
    "p": {'H': 1, 'e': -1},  # Proton
    "n": {},  # Nothing
    "a": {'e': -1},  # -electron
    "b": {'e': -1},  # -electron
    "c": {'H': 2, 'e': -1},  # Steals proton and hydrogen from neutral fragment
    "x": {'e': -1},  # -electron
    "y": {'H': 2, 'e': -1},  # Steals proton and hydrogen from neutral fragment
    "z": {'e': -1},
    "ax": {'e': -1},
    "ay": {'H': 1, 'e': -1},
    "az": {'e': -1},
    "bx": {'e': -1},
    "by": {'H': 1, 'e': -1},
    "bz": {'e': -1},
    "cx": {'H': 1, 'e': -1},
    "cy": {'H': 3, 'e': -1},
    "cz": {'H': 1, 'e': -1},
    'i': {'H': 1, 'e': -1}
}

FRAGMENT_ION_BASE_CHARGE_ADDUCTS: Dict[str, str] = {
    "p": '+H+',
    "n": '',
    "a": '-e-',
    "b": '-e-',
    "c": '+2H+,+e-',
    "x": '-e-',
    "y": '+2H+,+e-',
    "z": '-e-',
    "ax": '-e-',
    "ay": '+H+',
    "az": '-e-',
    "bx": '-e-',
    "by": '+H+',
    "bz": '-e-',
    "cx": '+H+',
    "cy": '+3H+,+e-',
    "cz": '+H+',
    'i': '+H+'
}

# According to the dissociation sites
NEUTRAL_FRAGMENT_START_COMPOSITIONS: Dict[str, Dict[str, int]] = {
    "p": NTERM_COMPOSITION,
    "a": NTERM_COMPOSITION,
    "b": NTERM_COMPOSITION,
    "c": NTERM_COMPOSITION,
    "x": {"O": 1, "C": 1},
    "y": {},
    "z": {"N": -1, "H": -1},
    "i": {},
    'n': {}
}

# According to the dissociation sites
NEUTRAL_FRAGMENT_END_COMPOSITIONS: Dict[str, Dict[str, int]] = {
    "p": CTERM_COMPOSITION,
    "a": {"O": -1, "C": -1},
    "b": {},
    "c": {"N": 1, "H": 1},
    "x": CTERM_COMPOSITION,
    "y": CTERM_COMPOSITION,
    "z": CTERM_COMPOSITION,
    "i": {"O": -1, "C": -1},
    'n': {}
}

NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS: Dict[str, Dict[str, int]] = {
    "p": merge_dicts(NEUTRAL_FRAGMENT_START_COMPOSITIONS['p'], NEUTRAL_FRAGMENT_END_COMPOSITIONS['p']),
    "n": merge_dicts(NEUTRAL_FRAGMENT_START_COMPOSITIONS['n'], NEUTRAL_FRAGMENT_END_COMPOSITIONS['n']),
    "a": merge_dicts(NEUTRAL_FRAGMENT_START_COMPOSITIONS['a'], NEUTRAL_FRAGMENT_END_COMPOSITIONS['a']),
    "b": merge_dicts(NEUTRAL_FRAGMENT_START_COMPOSITIONS['b'], NEUTRAL_FRAGMENT_END_COMPOSITIONS['b']),
    "c": merge_dicts(NEUTRAL_FRAGMENT_START_COMPOSITIONS['c'], NEUTRAL_FRAGMENT_END_COMPOSITIONS['c']),
    "x": merge_dicts(NEUTRAL_FRAGMENT_START_COMPOSITIONS['x'], NEUTRAL_FRAGMENT_END_COMPOSITIONS['x']),
    "y": merge_dicts(NEUTRAL_FRAGMENT_START_COMPOSITIONS['y'], NEUTRAL_FRAGMENT_END_COMPOSITIONS['y']),
    "z": merge_dicts(NEUTRAL_FRAGMENT_START_COMPOSITIONS['z'], NEUTRAL_FRAGMENT_END_COMPOSITIONS['z']),
    "ax": merge_dicts(NEUTRAL_FRAGMENT_END_COMPOSITIONS['a'], NEUTRAL_FRAGMENT_START_COMPOSITIONS['x']),
    "ay": merge_dicts(NEUTRAL_FRAGMENT_END_COMPOSITIONS['a'], NEUTRAL_FRAGMENT_START_COMPOSITIONS['y']),
    "az": merge_dicts(NEUTRAL_FRAGMENT_END_COMPOSITIONS['a'], NEUTRAL_FRAGMENT_START_COMPOSITIONS['z']),
    "bx": merge_dicts(NEUTRAL_FRAGMENT_END_COMPOSITIONS['b'], NEUTRAL_FRAGMENT_START_COMPOSITIONS['x']),
    "by": merge_dicts(NEUTRAL_FRAGMENT_END_COMPOSITIONS['b'], NEUTRAL_FRAGMENT_START_COMPOSITIONS['y']),
    "bz": merge_dicts(NEUTRAL_FRAGMENT_END_COMPOSITIONS['b'], NEUTRAL_FRAGMENT_START_COMPOSITIONS['z']),
    "cx": merge_dicts(NEUTRAL_FRAGMENT_END_COMPOSITIONS['c'], NEUTRAL_FRAGMENT_START_COMPOSITIONS['x']),
    "cy": merge_dicts(NEUTRAL_FRAGMENT_END_COMPOSITIONS['c'], NEUTRAL_FRAGMENT_START_COMPOSITIONS['y']),
    "cz": merge_dicts(NEUTRAL_FRAGMENT_END_COMPOSITIONS['c'], NEUTRAL_FRAGMENT_START_COMPOSITIONS['z']),
    "i": merge_dicts(NEUTRAL_FRAGMENT_START_COMPOSITIONS['i'], NEUTRAL_FRAGMENT_END_COMPOSITIONS['i']),
}

FRAGMENT_ION_COMPOSITION_ADJUSTMENTS: Dict[str, Dict[str, int]] = {
    "p": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['p'], FRAGMENT_ION_COMPOSITIONS['p']),
    "n": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['n'], FRAGMENT_ION_COMPOSITIONS['n']),
    "a": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['a'], FRAGMENT_ION_COMPOSITIONS['a']),
    "b": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['b'], FRAGMENT_ION_COMPOSITIONS['b']),
    "c": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['c'], FRAGMENT_ION_COMPOSITIONS['c']),
    "x": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['x'], FRAGMENT_ION_COMPOSITIONS['x']),
    "y": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['y'], FRAGMENT_ION_COMPOSITIONS['y']),
    "z": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['z'], FRAGMENT_ION_COMPOSITIONS['z']),
    "ax": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['ax'], FRAGMENT_ION_COMPOSITIONS['ax']),
    "ay": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['ay'], FRAGMENT_ION_COMPOSITIONS['ay']),
    "az": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['az'], FRAGMENT_ION_COMPOSITIONS['az']),
    "bx": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['bx'], FRAGMENT_ION_COMPOSITIONS['bx']),
    "by": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['by'], FRAGMENT_ION_COMPOSITIONS['by']),
    "bz": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['bz'], FRAGMENT_ION_COMPOSITIONS['bz']),
    "cx": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['cx'], FRAGMENT_ION_COMPOSITIONS['cx']),
    "cy": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['cy'], FRAGMENT_ION_COMPOSITIONS['cy']),
    "cz": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['cz'], FRAGMENT_ION_COMPOSITIONS['cz']),
    "i": merge_dicts(NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS['i'], FRAGMENT_ION_COMPOSITIONS['i']),
}

AA_COMPOSITIONS: Dict[str, Dict[str, int]] = {
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
    "X": {},  # Unknown amino acid
}

AMINO_ACIDS: Set[str] = set(AA_COMPOSITIONS.keys()) | {'B', 'Z'}
ORDERED_AMINO_ACIDS = sorted(list(AMINO_ACIDS))
AMBIGUOUS_AMINO_ACIDS: Set[str] = {'J', 'B', 'Z', 'X'}
MASS_AMBIGUOUS_AMINO_ACIDS: Set[str] = {'B', 'Z'}

# 3 letter codes
AA_TO_THREE_LETTER_CODE: Dict[str, str] = {
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
THREE_LETTER_CODE_TO_AA: Dict[str, str] = {v: k for k, v in AA_TO_THREE_LETTER_CODE.items()}

# Amino acid to name
AA_TO_NAME: Dict[str, str] = {
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
NAME_TO_AA: Dict[str, str] = {v: k for k, v in AA_TO_NAME.items()}

FORWARD_ION_TYPES: Set[str] = {'a', 'b', 'c'}
BACKWARD_ION_TYPES: Set[str] = {'x', 'y', 'z'}
INTERNAL_ION_TYPES: Set[str] = {'ax', 'ay', 'az', 'bx', 'by', 'bz', 'cx', 'cy', 'cz'}
TERMINAL_ION_TYPES: Set[str] = FORWARD_ION_TYPES | BACKWARD_ION_TYPES
IMMONIUM_ION_TYPES: Set[str] = {'i'}
VALID_ION_TYPES: Set[str] = TERMINAL_ION_TYPES | INTERNAL_ION_TYPES | IMMONIUM_ION_TYPES
AVERAGINE_RATIOS: Dict[str, float] = {'C': 4.9384, 'H': 7.7583, 'N': 1.3577, 'O': 1.4773, 'S': 0.0417}

PROTEASES: Dict[str, str] = \
    {
        'arg-c': '(?<=R)',
        'asp-n': '(?=D)',
        'chymotrypsin': '(?<=[FWYL])(?!P)',
        'chymotrypsin/P': '(?<=[FWYL])',
        'promega-chymotrypsin-high-specificity': '(?<=[YFW])',
        'promega-chymotrypsin-low-specificity': '(?<=[YFWLM])',
        'glu-c': '(?<=E)',
        'lys-c': '(?<=K)',
        'lys-n': '(?=K)',
        'proteinase k': '(?<=[AEFILTVWY])',
        'trypsin': '(?<=[KR])(?=[^P])',
        'trypsin/P': '(?<=[KR])',
        'proalanase': '(?<=[PA])',
        'elastase': r'(?<=[AGSVLI])',
        'pepsin': r'(?<=[FLWY])',
        'thermolysin': r'(?<=[LFIAVM])',
        'proalanase-low-specificity': '(?<=[PASG])',
        'non-specific': '()',
        'no-cleave': '_',
    }

PROTEASES_COMPILED = {}
for key, value in PROTEASES.items():
    try:
        PROTEASES_COMPILED[key] = regex.compile(value)
    except Exception as e:
        print(f"Error compiling regex for {key}: {e}")
        raise e

# Compiling regex patterns used in your module
ISOTOPE_COMPONENT_PATTERN = regex.compile(r'([0-9]*)([A-Za-z]+)(-?\d*\.?\d*)')
CONDENSED_CHEM_FORMULA_PATTERN = regex.compile(r'([A-Z][a-z]*|e|p|n)(-?\d*\.?\d*)')
ADDUCT_PATTERN = regex.compile(r"([+-]?)(\d*)?([A-Za-z]{1,2}\d*\+?-?)")
ISOTOPE_NUM_PATTERN = regex.compile(r'[0-9]')
