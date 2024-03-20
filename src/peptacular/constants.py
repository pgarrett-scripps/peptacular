import os
from copy import deepcopy
from typing import Dict, List, Tuple, Set

from peptacular.chem_setup import get_element_info, map_atomic_symbol_to_average_mass, \
    map_atomic_number_to_comp_neutron_offset, get_isotopic_atomic_masses, map_hill_order, map_atomic_number_to_comp, \
    map_atomic_number_to_symbol

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
ATOMIC_SYMBOL_TO_ISOTOPE_NEUTRON_OFFSETS_AND_ABUNDANCES: Dict[str, List[Tuple[int, float]]] = \
    map_atomic_number_to_comp_neutron_offset(_infos)
ATOMIC_SYMBOL_TO_ISOTOPE_MASSES_AND_ABUNDANCES: Dict[str, List[Tuple[float, float]]] = \
    map_atomic_number_to_comp(_infos)
HILL_ORDER: Dict[str, int] = map_hill_order(_infos)

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
    "W": {"C": 11, "H": 10, "N": 2},  # Tryptophan
    "U": {"C": 3, "H": 5, "N": 1, "O": 1, "Se": 1},  # Selenocysteine
    "O": {"C": 12, "H": 19, "N": 3, "O": 2},  # Pyrrolysine
    "X": {},  # Unknown amino acid
}


ION_TYPE_START_COMPOSITIONS: Dict[str, Dict[str, int]] = {
    "a": {"H": 1},
    "b": {"H": 1},
    "c": {"H": 1},
    "x": {"O": 1, "C": 1},
    "y": {"H": 2},
    "z": {"H": -1, "N": -1},
    "p": {'H': 2},
}

ION_TYPE_END_COMPOSITIONS: Dict[str, Dict[str, int]] = {
    "a": {"C": -1, "O": -1},
    "b": {},
    "c": {"H": 3, "N": 1},
    "x": {"O": 1, "H": 1},
    "y": {"O": 1, "H": 1},
    "z": {"O": 1, "H": 1},
    "p": {"O": 1, "H": 1},
}


def _create_ion_adjustments(atomic_masses, start_comps, end_comps):
    adjustments = {}
    for ion_type in start_comps:
        start_mass = sum([atomic_masses[k] * v for k, v in start_comps[ion_type].items()])
        end_mass = sum([atomic_masses[k] * v for k, v in end_comps[ion_type].items()])
        adjustments[ion_type] = start_mass + end_mass - atomic_masses['H']

    for forward_ion in 'abc':
        for backward_ion in 'xyz':
            start_mass = sum([atomic_masses[k] * v for k, v in end_comps[forward_ion].items()])
            end_mass = sum([atomic_masses[k] * v for k, v in start_comps[backward_ion].items()])
            adjustments[forward_ion + backward_ion] = start_mass + end_mass - atomic_masses['H']

    adjustments['I'] = sum([atomic_masses[k] * v for k, v in end_comps['a'].items()])

    return adjustments


def _create_ion_type_compositions(start_comps, end_comps):
    compositions = {}
    for ion_type in start_comps:
        start_comp = start_comps[ion_type]
        end_comp = end_comps[ion_type]

        # combine the two dictionaries
        combined = deepcopy(start_comp)
        for k, v in end_comp.items():
            combined[k] = combined.get(k, 0) + v

        combined['H'] = combined.get('H', 0) - 1  # adjust for neutral mass

        # drop any zero values
        combined = {k: v for k, v in combined.items() if v != 0}

        compositions[ion_type] = combined

    for forward_ion in 'abc':
        for backward_ion in 'xyz':
            start_comp = start_comps[backward_ion]
            end_comp = end_comps[forward_ion]

            # combine the two dictionaries
            combined = deepcopy(start_comp)
            for k, v in end_comp.items():
                combined[k] = combined.get(k, 0) + v

            combined['H'] = combined.get('H', 0) - 1  # adjust for neutral mass

            # drop any zero values
            combined = {k: v for k, v in combined.items() if v != 0}

            compositions[forward_ion + backward_ion] = combined

    compositions['I'] = end_comps['a']

    return compositions


ION_TYPE_COMPOSITION_ADJUSTMENTS: Dict[str, Dict[str, int]] = \
    _create_ion_type_compositions(ION_TYPE_START_COMPOSITIONS, ION_TYPE_END_COMPOSITIONS)


# print(ION_TYPE_COMPOSITION_ADJUSTMENTS)

def _create_aa_masses(atomic_masses, aa_compositions):
    aa_masses = {}
    for aa in aa_compositions:
        aa_masses[aa] = sum([atomic_masses[k] * v for k, v in aa_compositions[aa].items()])
    return aa_masses


MONOISOTOPIC_AA_MASSES: Dict[str, float] = _create_aa_masses(ISOTOPIC_ATOMIC_MASSES, AA_COMPOSITIONS)
AVERAGE_AA_MASSES: Dict[str, float] = _create_aa_masses(AVERAGE_ATOMIC_MASSES, AA_COMPOSITIONS)
AMINO_ACIDS: Set[str] = set(AA_COMPOSITIONS.keys()) | {'B', 'Z'}
AMBIGUOUS_AMINO_ACIDS: Set[str] = {'J', 'B', 'Z', 'X'}
MASS_AMBIGUOUS_AMINO_ACIDS: Set[str] = {'B', 'Z'}

MONOISOTOPIC_ION_ADJUSTMENTS: Dict[str, float] = \
    _create_ion_adjustments(ISOTOPIC_ATOMIC_MASSES, ION_TYPE_START_COMPOSITIONS, ION_TYPE_END_COMPOSITIONS)
AVERAGE_ION_ADJUSTMENTS: Dict[str, float] = \
    _create_ion_adjustments(AVERAGE_ATOMIC_MASSES, ION_TYPE_START_COMPOSITIONS, ION_TYPE_END_COMPOSITIONS)

FORWARD_ION_TYPES: Set[str] = {'a', 'b', 'c'}
BACKWARD_ION_TYPES: Set[str] = {'x', 'y', 'z'}
INTERNAL_ION_TYPES: Set[str] = {'ax', 'ay', 'az', 'bx', 'by', 'bz', 'cx', 'cy', 'cz'}
TERMINAL_ION_TYPES: Set[str] = FORWARD_ION_TYPES | BACKWARD_ION_TYPES
IMMONIUM_ION_TYPES: Set[str] = {'I'}
VALID_ION_TYPES: Set[str] = TERMINAL_ION_TYPES | INTERNAL_ION_TYPES | IMMONIUM_ION_TYPES

AVERAGINE_RATIOS: Dict[str, float] = {'C': 4.9384, 'H': 7.7583, 'N': 1.3577, 'O': 1.4773, 'S': 0.0417}
ISOTOPIC_AVERAGINE_MASS: float = sum([v * ISOTOPIC_ATOMIC_MASSES[k] for k, v in AVERAGINE_RATIOS.items()])
AVERAGE_AVERAGINE_MASS: float = sum([v * AVERAGE_ATOMIC_MASSES[k] for k, v in AVERAGINE_RATIOS.items()])


PROTEASES: Dict[str, str] = \
    {
        'arg-c': 'R',
        'asp-n': '\\w(?=D)',
        'bnps-skatole': 'W',
        'caspase 1': '(?<=[FWYL]\\w[HAT])D(?=[^PEDQKR])',
        'caspase 2': '(?<=DVA)D(?=[^PEDQKR])',
        'caspase 3': '(?<=DMQ)D(?=[^PEDQKR])',
        'caspase 4': '(?<=LEV)D(?=[^PEDQKR])',
        'caspase 5': '(?<=[LW]EH)D',
        'caspase 6': '(?<=VE[HI])D(?=[^PEDQKR])',
        'caspase 7': '(?<=DEV)D(?=[^PEDQKR])',
        'caspase 8': '(?<=[IL]ET)D(?=[^PEDQKR])',
        'caspase 9': '(?<=LEH)D',
        'caspase 10': '(?<=IEA)D',
        'chymotrypsin high specificity': '([FY](?=[^P]))|(W(?=[^MP]))',
        'chymotrypsin low specificity': '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
        'chymotrypsin': '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
        'clostripain': 'R',
        'cnbr': 'M',
        'enterokinase': '(?<=[DE]{3})K',
        'factor xa': '(?<=[AFGILTVM][DE]G)R',
        'formic acid': 'D',
        'glutamyl endopeptidase': 'E',
        'glu-c': 'E',
        'granzyme b': '(?<=IEP)D',
        'hydroxylamine': 'N(?=G)',
        'iodosobenzoic acid': 'W',
        'lys-c': 'K',
        'lys-n': '\\w(?=K)',
        'ntcb': '\\w(?=C)',
        'pepsin ph1.3': '((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\\w[^P]))',
        'pepsin ph2.0': '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\\w[^P]))',
        'proline endopeptidase': '(?<=[HKR])P(?=[^P])',
        'proteinase k': '[AEFILTVWY]',
        'staphylococcal peptidase i': '(?<=[^E])E',
        'thermolysin': '[^DE](?=[AFILMV])',
        'thrombin': '((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
        'trypsin_full': '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
        'trypsin_exception': '((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))',
        'trypsin': '([KR](?=[^P]))',
        'trypsin/P': '([KR])',
        'non-specific': '()',
        'no-cleave': '_'
    }
