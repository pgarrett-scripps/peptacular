import re
from typing import List


def find_peptide_indexes(protein: str, peptide: str) -> List[int]:
    """
    Retrieves all starting indexes of a given peptide within a protein sequence.

    :param protein: The complete protein sequence in which to search.
    :type protein: str
    :param peptide: The peptide sequence to find within the protein.
    :type peptide: str
    :return: A list of starting indexes where the peptide is found in the protein sequence.
    :rtype: List[int]
    """

    if len(peptide) == 0:
        return []

    return [i.start() for i in re.finditer(peptide, protein)]


def build_coverage_array(protein: str, peptides: List[str]) -> List[int]:
    """
    Calculate the coverage of a protein sequence by a list of peptides.

    The coverage is represented as a binary list where each position in the protein sequence is marked as 1 if it
    is covered by at least one peptide and 0 otherwise.

    :param protein: The protein sequence.
    :type protein: str
    :param peptides: List of peptide sequences.
    :type peptides: List[str]

    :return: A list representing the coverage of the protein sequence by the peptides. Each position in the
             list corresponds to a position in the protein sequence.
    :rtype: List[int]
    """

    cov_arr = [0] * len(protein)
    for peptide in peptides:
        peptide_indexes = find_peptide_indexes(protein, peptide)
        for peptide_index in peptide_indexes:
            cov_arr[peptide_index:peptide_index + len(peptide)] = [1] * len(peptide)
    return cov_arr


def calculate_percent_coverage(protein: str, peptides: List[str]) -> float:
    """
    Calculates the protein coverage of a list of peptides as a percentage.

    :param protein: The protein sequence.
    :type protein: str
    :param peptides: The list of peptide sequences.
    :type peptides: List[str]

    :return: The protein coverage percentage.
    :rtype: float
    """

    cov_arr = build_coverage_array(protein, peptides)

    if len(cov_arr) == 0:
        return 0

    return sum(cov_arr) / len(cov_arr)
