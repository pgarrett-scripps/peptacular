import re
from typing import Literal, Sequence, overload

from .sequence.parrallel import parallel_apply_internal


def _convert_ip2_sequence_single(sequence: str) -> str:
    """Internal function for converting a single IP2 sequence."""
    # Use regex to check if sequence starts and ends with the specified pattern
    if re.match(r"^([A-Z]|-)\..*\.([A-Z]|-)$", sequence):
        # If it matches, remove the leading and trailing characters (first and last two characters)
        sequence = sequence[2:-2]

    # Step 2: Replace () with []
    sequence = re.sub(r"\(([^)]+)\)", r"[\1]", sequence)

    # Step 3: Handle modifications at the start (can be any content, not just numbers)
    sequence = re.sub(r"^\[([^\]]+)\]", r"[\1]-", sequence)

    # Step 4: Convert consecutive modifications to use a dash
    sequence = re.sub(r"\]\[", r"]-[", sequence)

    return sequence


@overload
def convert_ip2_sequence(
    sequence: str,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> str: ...


@overload
def convert_ip2_sequence(
    sequence: Sequence[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[str]: ...


def convert_ip2_sequence(
    sequence: str | Sequence[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> str | list[str]:
    """
    Converts a IP2-Like sequence to a proforma2.0 compatible sequence.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence or list of sequences to be converted.
    :type sequence: str | Sequence[str]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: Proforma2.0 compatable sequence or list of sequences.
    :rtype: str | list[str]

    .. code-block:: python

        >>> convert_ip2_sequence('K.PEP(phospho)TIDE.K')
        'PEP[phospho]TIDE'

        >>> convert_ip2_sequence('K.(-1)PEP(phospho)TIDE.K')
        '[-1]-PEP[phospho]TIDE'

        >>> convert_ip2_sequence('K.PEPTIDE(2).K')
        'PEPTIDE[2]'

        >>> convert_ip2_sequence('K.PEPTIDE(2)(3).K')
        'PEPTIDE[2]-[3]'

        >>> convert_ip2_sequence('-.(1)PEP(phospho)TIDE(2)(3).-')
        '[1]-PEP[phospho]TIDE[2]-[3]'

        >>> convert_ip2_sequence('P')
        'P'

        >>> convert_ip2_sequence('')
        ''

        >>> convert_ip2_sequence('PEPTIDE')
        'PEPTIDE'
    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str):
        return parallel_apply_internal(
            _convert_ip2_sequence_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _convert_ip2_sequence_single(sequence)


def _convert_diann_sequence_single(sequence: str) -> str:
    """Internal function for converting a single DIANN sequence."""
    # Check if sequence starts and ends with underscores and remove them
    if sequence.startswith("_"):
        sequence = sequence[1:]

        # Check for a modification at the start of the sequence
        if re.match(r"^\[[^\]]+\]", sequence):
            sequence = re.sub(r"^\[([^\]]+)\]", r"[\1]-", sequence)

    if sequence.endswith("_"):
        sequence = sequence[:-1]

    elif re.search(r"_\[[^\]]+\]$", sequence):
        sequence = re.sub(r"_\[([^\]]+)\]$", r"-[\1]", sequence)

    return sequence


@overload
def convert_diann_sequence(
    sequence: str,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> str: ...


@overload
def convert_diann_sequence(
    sequence: Sequence[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[str]: ...


def convert_diann_sequence(
    sequence: str | Sequence[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> str | list[str]:
    """
    Converts a DIANN-Like sequence to a proforma2.0 compatible sequence.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence or list of sequences to be converted.
    :type sequence: str | Sequence[str]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: Proforma2.0 compatable sequence or list of sequences.
    :rtype: str | list[str]

    .. code-block:: python

        >>> convert_diann_sequence('_YMGTLRGC[Carbamidomethyl]LLRLYHD_')
        'YMGTLRGC[Carbamidomethyl]LLRLYHD'

        >>> convert_diann_sequence('_[Acytel]YMGTLRGC[Carbamidomethyl]LLRLYHD_')
        '[Acytel]-YMGTLRGC[Carbamidomethyl]LLRLYHD'

        >>> convert_diann_sequence('_[Acytel]YMGTLRGC[Carbamidomethyl]LLRLYHD[1.0]_[Methyl]')
        '[Acytel]-YMGTLRGC[Carbamidomethyl]LLRLYHD[1.0]-[Methyl]'

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str):
        return parallel_apply_internal(
            _convert_diann_sequence_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _convert_diann_sequence_single(sequence)


def _convert_casanovo_sequence_single(sequence: str) -> str:
    """Internal function for converting a single Casanovo sequence."""
    new_sequence_comps: list[str] = []
    in_mod = False  # Tracks if we are within a modification
    is_nterm = False  # Tracks if the current modification is at the N-terminus

    for _, char in enumerate(sequence):
        if char in {"+", "-"}:
            # Check if it's at the start (N-terminal)
            is_nterm = len(new_sequence_comps) == 0

            # Start a new modification block
            new_sequence_comps.append("[")
            new_sequence_comps.append(char)
            in_mod = True
        elif in_mod and char.isalpha():
            # End the modification block
            new_sequence_comps.append("]")

            if is_nterm:
                # Add a dash if it's an N-terminal modification
                new_sequence_comps.append("-")
                is_nterm = False

            # Add the current character and close modification
            in_mod = False
            new_sequence_comps.append(char)
        else:
            # Add regular characters
            new_sequence_comps.append(char)

    # Close any unclosed modification at the end of the sequence
    if in_mod:
        new_sequence_comps.append("]")

    return "".join(new_sequence_comps)


@overload
def convert_casanovo_sequence(
    sequence: str,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> str: ...


@overload
def convert_casanovo_sequence(
    sequence: Sequence[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[str]: ...


def convert_casanovo_sequence(
    sequence: str | Sequence[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> str | list[str]:
    """
    Converts a Casanovo sequence with modifications to a proforma2.0 compatible sequence.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence or list of sequences to be converted.
    :type sequence: str | Sequence[str]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: Proforma2.0 compatable sequence or list of sequences.
    :rtype: str | list[str]

    .. code-block:: python

        >>> convert_casanovo_sequence('+43.006P+100EPTIDE')
        '[+43.006]-P[+100]EPTIDE'

    """
    if isinstance(sequence, Sequence) and not isinstance(sequence, str):
        return parallel_apply_internal(
            _convert_casanovo_sequence_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _convert_casanovo_sequence_single(sequence)
