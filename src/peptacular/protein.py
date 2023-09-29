import regex as re
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

    .. code-block:: python

        >>> find_peptide_indexes("PEPTIDE", "PEP")
        [0]

        >>> find_peptide_indexes("PEPTIDE", "EPT")
        [1]

        >>> find_peptide_indexes("PEPTIDE", "E")
        [1, 6]

    """

    if len(peptide) == 0:
        return []

    return [i.start() for i in re.finditer(peptide, protein, overlapped=True)]


def build_coverage_array(protein: str, peptides: List[str], accumulate: bool = False) -> List[int]:
    """
    Calculate the coverage of a protein sequence by a list of peptides.

    The coverage is represented as a binary list where each position in the protein sequence is marked as 1 if it
    is covered by at least one peptide and 0 otherwise.

    :param protein: The protein sequence.
    :type protein: str
    :param peptides: List of peptide sequences.
    :type peptides: List[str]
    :param accumulate: If True, the coverage array will be accumulated, i.e. if a position is covered by more than
                        one peptide, it will be marked as the sum of the number of peptides covering it. If False,
                        the coverage array will be binary, i.e. if a position is covered by more than one peptide,
                        it will be marked as 1.
    :type accumulate: bool

    :return: A list representing the coverage of the protein sequence by the peptides. Each position in the
             list corresponds to a position in the protein sequence.
    :rtype: List[int]

    .. code-block:: python

        >>> build_coverage_array("PEPTIDE", ["PEP"])
        [1, 1, 1, 0, 0, 0, 0]

        >>> build_coverage_array("PEPTIDE", ["PEP", "EPT"])
        [1, 1, 1, 1, 0, 0, 0]

        # If accumulate is True, overlapping indecies will be accumulated
        >>> build_coverage_array("PEPTIDE", ["PEP", "EPT"], accumulate=True)
        [1, 2, 2, 1, 0, 0, 0]

    """

    cov_arr = [0] * len(protein)
    for peptide in peptides:
        peptide_indexes = find_peptide_indexes(protein, peptide)
        for peptide_index in peptide_indexes:
            if accumulate:
                cov_arr[peptide_index:peptide_index + len(peptide)] = \
                    [x + 1 for x in cov_arr[peptide_index:peptide_index + len(peptide)]]
            else:
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

    .. code-block:: python

        >>> calculate_percent_coverage("PEPTIDE", ["PEP"])
        0.42857142857142855

        >>> calculate_percent_coverage("PEPTIDE", ["PEP", "EPT"])
        0.5714285714285714

    """

    cov_arr = build_coverage_array(protein, peptides, accumulate=False)

    if len(cov_arr) == 0:
        return 0

    return sum(cov_arr) / len(cov_arr)
