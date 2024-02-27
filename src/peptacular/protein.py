import regex as re
from typing import List

from peptacular import sequence
from peptacular.sequence import span_to_sequence_fast


def find_peptide_indexes(protein: str, peptide: str, ignore_mods: bool = False) -> List[int]:
    """
    Retrieves all starting indexes of a given peptide within a protein sequence.

    :param protein: The complete protein sequence in which to search.
    :type protein: str
    :param peptide: The peptide sequence to find within the protein.
    :type peptide: str
    :param ignore_mods: Whether to ignore modifications when searching for the peptide.
    :type ignore_mods: bool

    :return: A list of starting indexes where the peptide is found in the protein sequence.
    :rtype: List[int]

    .. code-block:: python

        # Find the starting indexes of a peptide within a protein sequence.
        >>> find_peptide_indexes("PEPTIDE", "PEP")
        [0]
        >>> find_peptide_indexes("PEPTIDE", "EPT")
        [1]
        >>> find_peptide_indexes("PEPTIDE", "E")
        [1, 6]


        # By default the function will not ignore modifications
        >>> find_peptide_indexes("[Acetyl]-PEPTIDE", "PEP")
        []
        >>> find_peptide_indexes("[Acetyl]-PEPTIDE", "[Acetyl]-PEP")
        [0]
        >>> find_peptide_indexes("PEPTIDE", "[Acetyl]-PEP")
        []
        >>> find_peptide_indexes("PEPTIDE[1.0]", "IDE")
        []
        >>> find_peptide_indexes("[Acetyl]-PEPTIDE-[Amide]", "[Acetyl]-PEPTIDE-[Amide]")
        [0]

        # If ignore_mods is set to True, the function will ignore modifications
        >>> find_peptide_indexes("[Acetyl]-PEPTIDE", "PEP", ignore_mods=True)
        [0]
        >>> find_peptide_indexes("PEPTIDE", "[Acetyl]-PEP", ignore_mods=True)
        [0]

    """

    unnmodified_peptide, peptide_mods = sequence.pop_modifications(peptide)
    unmodified_protein, protein_mods = sequence.pop_modifications(protein)

    if ignore_mods:
        return [i.start() for i in re.finditer(unnmodified_peptide, unmodified_protein, overlapped=True)]

    spans = [(i.start(), i.end(), 0) for i in re.finditer(unnmodified_peptide, unmodified_protein, overlapped=True)]
    hit_peptides = [span_to_sequence_fast(unmodified_protein, protein_mods, s) for s in spans]

    return [s[0] for s, p in zip(spans, hit_peptides) if p == peptide]


def build_coverage_array(protein: str, peptides: List[str], accumulate: bool = False, ignore_mods: bool = False) -> List[int]:
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
    :param ignore_mods: Whether to ignore modifications when searching for the peptide.
    :type ignore_mods: bool

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
        peptide_indexes = find_peptide_indexes(protein, peptide, ignore_mods=ignore_mods)
        for peptide_index in peptide_indexes:
            if accumulate:
                cov_arr[peptide_index:peptide_index + len(peptide)] = \
                    [x + 1 for x in cov_arr[peptide_index:peptide_index + len(peptide)]]
            else:
                cov_arr[peptide_index:peptide_index + len(peptide)] = [1] * len(peptide)

    return cov_arr


def calculate_percent_coverage(protein: str, peptides: List[str], ignore_mods: bool = False) -> float:
    """
    Calculates the protein coverage of a list of peptides as a percentage.

    :param protein: The protein sequence.
    :type protein: str
    :param peptides: The list of peptide sequences.
    :type peptides: List[str]
    :param ignore_mods: Whether to ignore modifications when searching for the peptide.
    :type ignore_mods: bool

    :return: The protein coverage percentage.
    :rtype: float

    .. code-block:: python

        >>> calculate_percent_coverage("PEPTIDE", ["PEP"])
        0.42857142857142855

        >>> calculate_percent_coverage("PEPTIDE", ["PEP", "EPT"])
        0.5714285714285714

    """

    cov_arr = build_coverage_array(protein, peptides, accumulate=False, ignore_mods=ignore_mods)

    if len(cov_arr) == 0:
        return 0

    return sum(cov_arr) / len(cov_arr)
