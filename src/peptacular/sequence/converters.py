from typing import overload
from ..constants import ParrallelMethod, ParrallelMethodLiteral
from .parrallel import parallel_apply_internal
from ..annotation import ProFormaAnnotation


def _convert_ip2_sequence_single(sequence: str) -> str:
    """Internal function for converting a single IP2 sequence."""
    # Use regex to check if sequence starts and ends with the specified pattern
    return ProFormaAnnotation.from_ip2_sequence(sequence).serialize()


@overload
def convert_ip2_sequence(
    sequence: str,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str: ...


@overload
def convert_ip2_sequence(
    sequence: list[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


def convert_ip2_sequence(
    sequence: str | list[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
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
    if isinstance(sequence, list):
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
    return ProFormaAnnotation.from_diann(sequence).serialize()


@overload
def convert_diann_sequence(
    sequence: str,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str: ...


@overload
def convert_diann_sequence(
    sequence: list[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


def convert_diann_sequence(
    sequence: str | list[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
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
    if isinstance(sequence, list):
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
    return ProFormaAnnotation.from_casanovo(sequence).serialize()


@overload
def convert_casanovo_sequence(
    sequence: str,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str: ...


@overload
def convert_casanovo_sequence(
    sequence: list[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


def convert_casanovo_sequence(
    sequence: str | list[str],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
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
    if isinstance(sequence, list):
        return parallel_apply_internal(
            _convert_casanovo_sequence_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _convert_casanovo_sequence_single(sequence)
