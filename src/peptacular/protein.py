import re
from typing import List


def get_peptide_indexes_in_protein(protein: str, peptide: str) -> List[int]:
    """
    Get all indexes of a peptide (substring) in a protein (string).

    Args:
        peptide (str): The peptide sequence to search for in the protein.
        protein (str): The protein sequence.

    Returns:
        list[int]: List of indexes where the peptide sequence starts in the protein sequence.
    """
    return [i.start() for i in re.finditer(peptide, protein)]


def calculate_protein_coverage(protein: str, peptides: List[int]) -> List[int]:
    """
    Calculate the coverage of a protein sequence by a list of peptides.

    The coverage is represented as a binary list where each position in the protein sequence is marked as 1 if it
    is covered by at least one peptide and 0 otherwise.

    Args:
        protein (str): The protein sequence.
        peptides (list[str]): List of peptide sequences.

    Returns:
        list[int]: A list representing the coverage of the protein sequence by the peptides. Each position in the
        list corresponds to a position in the protein sequence.
    """
    cov_arr = [0] * len(protein)
    for peptide in peptides:
        peptide_indexes = get_peptide_indexes_in_protein(protein, peptide)
        for peptide_index in peptide_indexes:
            cov_arr[peptide_index:peptide_index + len(peptide)] = [1] * len(peptide)
    return cov_arr


def calculate_protein_coverage_percent(protein: str, peptides: List[str]) -> float:
    """
    Calculates the protein coverage of a list of peptides.

    Args:
        protein: The protein sequence.
        peptides: The list of peptide sequences.

    Returns:
        The protein coverage.
    """
    cov_arr = calculate_protein_coverage(protein, peptides)

    return sum(cov_arr) / len(cov_arr)
