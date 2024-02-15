from typing import List, Union

# https://www.unimod.org/masses.html

isotopes = {}
isotopes["C"] = [0.9893, 0.0107]
isotopes["H"] = [0.999885, 0.000115]
isotopes["O"] = [0.99757, 0.00038, 0.00205]
isotopes["N"] = [0.99636, 0.00364]
isotopes["S"] = [0.9499, 0.0075, 0.0425, 0.0001]

aa_composition = {
    "G": {"C": 2, "H": 3, "N": 1, "O": 1, "S": 0},  # Glycine
    "A": {"C": 3, "H": 5, "N": 1, "O": 1, "S": 0},  # Alanine
    "S": {"C": 3, "H": 5, "N": 1, "O": 2, "S": 0},  # Serine
    "P": {"C": 5, "H": 7, "N": 1, "O": 1, "S": 0},  # Proline
    "V": {"C": 5, "H": 9, "N": 1, "O": 1, "S": 0},  # Valine
    "T": {"C": 4, "H": 7, "N": 1, "O": 2, "S": 0},  # Threonine
    "C": {"C": 3, "H": 5, "N": 1, "O": 1, "S": 1},  # Cysteine
    "I": {"C": 6, "H": 11, "N": 1, "O": 1, "S": 0},  # Isoleucine
    "L": {"C": 6, "H": 11, "N": 1, "O": 1, "S": 0},  # Leucine
    "N": {"C": 4, "H": 6, "N": 2, "O": 2, "S": 0},  # Asparagine
    "D": {"C": 4, "H": 5, "N": 1, "O": 3, "S": 0},  # Aspartic acid
    "Q": {"C": 5, "H": 8, "N": 2, "O": 3, "S": 0},  # Glutamine
    "K": {"C": 6, "H": 12, "N": 2, "O": 0, "S": 0},  # Lysine
    "E": {"C": 5, "H": 7, "N": 1, "O": 3, "S": 0},  # Glutamic acid
    "M": {"C": 5, "H": 9, "N": 1, "O": 1, "S": 1},  # Methionine
    "H": {"C": 6, "H": 7, "N": 3, "O": 1, "S": 0},  # Histidine
    "F": {"C": 9, "H": 9, "N": 1, "O": 1, "S": 0},  # Phenylalanine
    "R": {"C": 6, "H": 12, "N": 4, "O": 1, "S": 0},  # Arginine
    "Y": {"C": 9, "H": 9, "N": 1, "O": 2, "S": 0},  # Tyrosine
    "W": {"C": 11, "H": 10, "N": 2, "O": 0, "S": 0},  # Tryptophan
}

def get_isotopes(sequence: str, n: Union[int, float]) -> List[float]:
    pass


def get_isotope(sequence: str, n: int) -> float:
    pass