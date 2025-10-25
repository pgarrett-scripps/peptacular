from typing import Iterable

from ..proforma.annotation import (
    ProFormaAnnotation,
)
from .basic import count_residues
from .util import get_annotation_input


def is_subsequence(
    subsequence: str | ProFormaAnnotation,
    sequence: str | ProFormaAnnotation,
    order: bool = True,
) -> bool:
    """
    Checks if the input subsequence is a subsequence of the input sequence. If order is True, the subsequence must be in
    the same order as in the sequence. If order is False, the subsequence can be in any order.

    :param subsequence: The sequence or ProFormaAnnotation object, representing the subsequence.
    :type subsequence: Union[str, ProFormaAnnotation]
    :param sequence: The sequence or ProFormaAnnotation object, representing the sequence.
    :type sequence: Union[str, ProFormaAnnotation]
    :param order: If True, the subsequence must be in the same order as in the sequence.
    :type order: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: True if the subsequence is a subsequence of the sequence, False otherwise.
    :rtype: bool

    .. code-block:: python

        >>> is_subsequence('PEP', 'PEPTIDE')
        True

        >>> is_subsequence('PET', 'PEPTIDE', order=False)
        True

        >>> is_subsequence('PET', 'PEPTIDE', order=True)
        False

        >>> is_subsequence('<13C>PEP', '<13C>PEPTIDE', order=True)
        True

        >>> is_subsequence('<13C>PEP[1.0]', '<13C>PEP[1.0]TIDE', order=True)
        True

        >>> is_subsequence('<13C>PEP', '<13C>PEP[1.0]TIDE', order=True)
        False

    """

    if order is True:
        return len(find_subsequence_indices(sequence, subsequence)) != 0

    subsequence_counts = count_residues(subsequence)
    sequence_counts = count_residues(sequence)
    return all(
        subsequence_counts[aa] <= sequence_counts[aa] for aa in subsequence_counts
    )


def find_subsequence_indices(
    sequence: str | ProFormaAnnotation,
    subsequence: str | ProFormaAnnotation,
    ignore_mods: bool = False,
) -> list[int]:
    """
    Retrieves all starting indexes of a given subsequence within a sequence.

    :param sequence: The sequence or ProFormaAnnotation objectto look for the 'subsequence' in.
    :type sequence: Union[str, ProFormaAnnotation]
    :param subsequence: The sequence or ProFormaAnnotation object to be found within 'sequence'.
    :type subsequence: Union[str, ProFormaAnnotation]
    :param ignore_mods: Whether to ignore modifications.
    :type ignore_mods: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A list of starting indexes of the subsequence within the sequence.
    :rtype: List[int]

    .. code-block:: python

        # Find the starting indexes of a subsequence
        >>> find_subsequence_indices("PEPTIDE", "PEP")
        [0]
        >>> find_subsequence_indices("PEPTIDE", "EPT")
        [1]
        >>> find_subsequence_indices("PEPTIDE", "E")
        [1, 6]

        # By default the function will not ignore modifications
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE", "PEP")
        []
        >>> find_subsequence_indices("<13C>PEPTIDE", "PEP")
        []
        >>> find_subsequence_indices("<13C>PEP[1][Phospho]TIDE", "<13C>PEP[Phospho][1]")
        [0]
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE", "[Acetyl]-PEP")
        [0]
        >>> find_subsequence_indices("PEPTIDE", "[Acetyl]-PEP")
        []
        >>> find_subsequence_indices("PEPTIDE[1.0]", "IDE")
        []
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE-[Amide]", "[Acetyl]-PEPTIDE-[Amide]")
        [0]

        # If ignore_mods is set to True, the function will ignore modifications
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE", "PEP", ignore_mods=True)
        [0]
        >>> find_subsequence_indices("PEPTIDE", "[Acetyl]-PEP", ignore_mods=True)
        [0]

    """
    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequence_annot = get_annotation_input(subsequence, copy=False)
    return subsequence_annot.find_indices(other=sequence_annot, ignore_mods=ignore_mods)


def coverage(
    sequence: str | ProFormaAnnotation,
    subsequences: Iterable[str | ProFormaAnnotation],
    accumulate: bool = False,
    ignore_mods: bool = False,
    ignore_ambiguity: bool = False,
) -> list[int]:
    """
    Calculate the sequence coverage given a list of subsequecnes.

    The coverage is represented as a binary list where each position in the protein sequence is marked as 1 if it
    is covered by at least one peptide and 0 otherwise.

    :param sequence: The sequence or ProFormaAnnotation object to be covered.
    :type sequence: Union[str, ProFormaAnnotation]
    :param subsequences: The sequence's or ProFormaAnnotation object's subsequences to be used for coverage.
    :type subsequences: List[Union[str, ProFormaAnnotation]]
    :param accumulate: If True, the coverage array will be accumulated, i.e. if a position is covered by more than
                        one subsequence, it will be marked as the sum of the number of subsequences covering it.
                        If False, the coverage array will be binary, i.e. if a position is covered by more than one
                        subsequence, it will be marked as 1.
    :type accumulate: bool
    :param ignore_mods: Whether to ignore modifications when calcualting the coverage.
    :type ignore_mods: bool
    :param ignore_ambiguity: Whether to ignore ambiguity codes when calculating the coverage.
    :type ignore_ambiguity: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The sequence coverage array.
    :rtype: List[int]

    .. code-block:: python

        >>> coverage("PEPTIDE", ["PEP"])
        [1, 1, 1, 0, 0, 0, 0]

        >>> coverage("PEPTIDE", ["PEP", "EPT"])
        [1, 1, 1, 1, 0, 0, 0]

        # If accumulate is True, overlapping indecies will be accumulated
        >>> coverage("PEPTIDE", ["PEP", "EPT"], accumulate=True)
        [1, 2, 2, 1, 0, 0, 0]

        # By default ambiguity does not add to coverage
        >>> coverage("PEPTIDE", ["P(?EP)"])
        [1, 0, 0, 0, 0, 0, 0]

    """

    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequences_annot = [
        get_annotation_input(subseq, copy=False) for subseq in subsequences
    ]

    return sequence_annot.coverage(
        subsequences_annot,
        accumulate=accumulate,
        ignore_mods=ignore_mods,
        ignore_ambiguity=ignore_ambiguity,
    )


def percent_coverage(
    sequence: str | ProFormaAnnotation,
    subsequences: Iterable[str | ProFormaAnnotation],
    ignore_mods: bool = False,
    accumulate: bool = False,
    ignore_ambiguity: bool = False,
) -> float:
    """
    Calculates the coverage given a list of subsequences.

    :param sequence: The sequence or ProFormaAnnotation object to be covered.
    :type sequence: Union[str, ProFormaAnnotation]
    :param subsequences: The sequence's or ProFormaAnnotation object's subsequences to be used for coverage.
    :type subsequences: List[Union[str, ProFormaAnnotation]]
    :param ignore_mods: Whether to ignore modifications when calcualting the coverage.
    :type ignore_mods: bool
    :param accumulate: If True, the coverage array will be accumulated
    :param accumulate: bool
    :param ignore_ambiguity: Whether to ignore ambiguity codes when calculating the coverage.
    :type ignore_ambiguity: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The percent coverage.
    :rtype: float

    .. code-block:: python

        >>> round(percent_coverage("PEPTIDE", ["PEP"]), 3)
        0.429

        >>> round(percent_coverage("PEPTIDE", ["PEP", "EPT"]), 3)
        0.571

        >>> round(percent_coverage("PEPTIDE", ["PEPTIDE", "PEPTIDE"], accumulate=True), 3)
        2.0

    """

    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequences_annot = [
        get_annotation_input(subseq, copy=False) for subseq in subsequences
    ]

    return sequence_annot.percent_coverage(
        subsequences_annot,
        accumulate=accumulate,
        ignore_mods=ignore_mods,
        ignore_ambiguity=ignore_ambiguity,
    )


def modification_coverage(
    sequence: str | ProFormaAnnotation,
    subsequences: list[str | ProFormaAnnotation],
    accumulate: bool = False,
) -> dict[int, int]:
    """
    Calculate the modification coverage given a list of subsequences.

    This function identifies which modifications in the main sequence are covered by
    subsequences. It returns a dictionary where each key is a modification position/type
    (matching the format from get_mods) and each value is the number of subsequences
    that cover that modification.

    :param sequence: The sequence or ProFormaAnnotation object representing the protein.
    :type sequence: Union[str, ProFormaAnnotation]
    :param subsequences: The sequence's or ProFormaAnnotation object's subsequences (peptides) to be used for coverage.
    :type subsequences: List[Union[str, ProFormaAnnotation]]
    :param accumulate: If True, the count will accumulate for each subsequence covering the modification.
                       If False, the count will be binary (1 if covered, 0 if not).
    :type accumulate: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: Dictionary mapping modification positions/types to coverage counts.
    :rtype: Dict[Union[int, str], int]

    .. code-block:: python

        >>> modification_coverage("PEPTIDE[Phospho]", ["TIDE"])
        {6: 0}

        >>> modification_coverage("PEPTIDE[Phospho]", ["TIDE[Phospho]"])
        {6: 1}

        >>> modification_coverage("PEP[Phospho]TIDE[Methyl]", ["PEP[Phospho]", "TIDE[Methyl]"])
        {2: 1, 6: 1}

        >>> modification_coverage("PEP[Phospho]TIDE[Methyl]", ["PEP[Phospho]", "TIDE[Methyl]", "PEP[Phospho]"], accumulate=True)
        {2: 2, 6: 1}

    """
    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequence_annots = [
        get_annotation_input(subseq, copy=False) for subseq in subsequences
    ]

    return sequence_annot.modification_coverage(
        annotations=subsequence_annots, accumulate=accumulate
    )
