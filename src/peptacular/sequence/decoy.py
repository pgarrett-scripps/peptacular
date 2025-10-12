from typing import Any, Sequence, overload, Literal
from collections import defaultdict
import random

from .util import get_annotation_input
from .parrallel import parallel_apply_internal
from ..proforma.annotation import ProFormaAnnotation


# ============================================================================
# Reverse Sequence
# ============================================================================


def _reverse_sequence_single(sequence: str | ProFormaAnnotation) -> str:
    """Internal function for reversing a single sequence."""
    annot = get_annotation_input(sequence, copy=False)
    annot.sequence = annot.sequence[::-1]
    return annot.serialize()


@overload
def reverse_sequence(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> str: ...


@overload
def reverse_sequence(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[str]: ...


def reverse_sequence(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> str | list[str]:
    """
    Reverse a peptide sequence or list of sequences.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation]
    :param n_workers: Number of worker processes. If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk. If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect).
    :type method: Literal["process", "thread"] | None
    :return: The reversed sequence or list of reversed sequences.
    :rtype: str | list[str]
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _reverse_sequence_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _reverse_sequence_single(sequence)


# ============================================================================
# Shuffle Sequence
# ============================================================================


def _shuffle_sequence_single(
    sequence: str | ProFormaAnnotation,
    static_residues: str = "",
) -> str:
    """Internal function for shuffling a single sequence."""
    annot = get_annotation_input(sequence, copy=False)
    seq_str = annot.sequence
    seq_list = list(seq_str)
    indices = [i for i, aa in enumerate(seq_list) if aa not in static_residues]
    residues_to_shuffle = [seq_list[i] for i in indices]
    random.shuffle(residues_to_shuffle)
    for idx, new_aa in zip(indices, residues_to_shuffle):
        seq_list[idx] = new_aa
    annot.sequence = "".join(seq_list)
    return annot.serialize()


@overload
def shuffle_sequence(
    sequence: str | ProFormaAnnotation,
    static_residues: str = "",
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> str: ...


@overload
def shuffle_sequence(
    sequence: Sequence[str | ProFormaAnnotation],
    static_residues: str = "",
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[str]: ...


def shuffle_sequence(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    static_residues: str = "",
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> str | list[str]:
    """
    Shuffle a peptide sequence while keeping static residues in place.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation]
    :param static_residues: Amino acids to keep in their original positions.
    :type static_residues: str
    :param n_workers: Number of worker processes. If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk. If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect).
    :type method: Literal["process", "thread"] | None
    :return: The shuffled sequence or list of shuffled sequences.
    :rtype: str | list[str]
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _shuffle_sequence_single,
            sequence,
            static_residues=static_residues,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _shuffle_sequence_single(sequence, static_residues=static_residues)


# ============================================================================
# Shift Sequence
# ============================================================================


def _shift_sequence_single(
    sequence: str | ProFormaAnnotation,
    n: int,
    static_residues: str = "",
) -> str:
    """Internal function for shifting a single sequence."""
    annot = get_annotation_input(sequence, copy=False)
    seq_str = annot.sequence
    seq_list = list(seq_str)
    indices = [i for i, aa in enumerate(seq_list) if aa not in static_residues]
    residues_to_shift = [seq_list[i] for i in indices]

    try:
        n = n % len(residues_to_shift)
    except ZeroDivisionError:
        return seq_str

    shifted_residues = residues_to_shift[-n:] + residues_to_shift[:-n]
    for idx, new_aa in zip(indices, shifted_residues):
        seq_list[idx] = new_aa
    annot.sequence = "".join(seq_list)
    return annot.serialize()


@overload
def shift_sequence(
    sequence: str | ProFormaAnnotation,
    n: int,
    static_residues: str = "",
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> str: ...


@overload
def shift_sequence(
    sequence: Sequence[str | ProFormaAnnotation],
    n: int,
    static_residues: str = "",
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[str]: ...


def shift_sequence(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n: int,
    static_residues: str = "",
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> str | list[str]:
    """
    Shift a peptide sequence by n positions while keeping static residues in place.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation]
    :param n: Number of positions to shift (positive for right, negative for left).
    :type n: int
    :param static_residues: Amino acids to keep in their original positions.
    :type static_residues: str
    :param n_workers: Number of worker processes. If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk. If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect).
    :type method: Literal["process", "thread"] | None
    :return: The shifted sequence or list of shifted sequences.
    :rtype: str | list[str]
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _shift_sequence_single,
            sequence,
            n=n,
            static_residues=static_residues,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _shift_sequence_single(sequence, n=n, static_residues=static_residues)


# ============================================================================
# De Bruijn Sequence
# ============================================================================


def _build_kmer_graph(seq: str, k: int) -> dict[str, dict[str, str]]:
    """Build k-mer graph where each k-mer maps to its next k-mer and transition letter."""
    nodes: dict[str, dict[str, str]] = {}

    kmers = []
    n_kmers = len(seq) - k + 1
    for i in range(n_kmers):
        kmer = seq[i : i + k]
        kmers.append(kmer)

    for i, kmer in enumerate(kmers):
        if i == len(kmers) - 1:
            break
        next_kmer = kmers[i + 1]
        next_letter = next_kmer[-1]
        if kmer not in nodes:
            nodes[kmer] = {next_kmer: next_letter}
        else:
            nodes[kmer][next_kmer] = next_letter

    return nodes


def _randomize_graph_edges(
    graph: dict[str, dict[str, str]],
    static_residues: str,
) -> None:
    """Randomize the amino acid transitions in the graph, keeping static residues fixed."""
    all_amino_acids = set()
    for kmer in graph:
        for next_kmer in graph[kmer]:
            aa = graph[kmer][next_kmer]
            if aa != "-":
                all_amino_acids.add(aa)

    randomizable_aas = [aa for aa in all_amino_acids if aa not in static_residues]

    if not randomizable_aas:
        return

    for kmer in graph:
        for next_kmer in graph[kmer]:
            aa = graph[kmer][next_kmer]
            if aa != "-" and aa not in static_residues:
                graph[kmer][next_kmer] = random.choice(randomizable_aas)


def _reconstruct_sequence_from_graph(
    original_seq: str,
    graph: dict[str, dict[str, str]],
    k: int,
) -> str:
    """Reconstruct a sequence by traversing the randomized graph."""
    new_sequence = []

    padded_seq = "-" * k + original_seq

    kmers = []
    n_kmers = len(padded_seq) - k + 1
    for i in range(n_kmers):
        kmer = padded_seq[i : i + k]
        kmers.append(kmer)

    for i, kmer in enumerate(kmers):
        if i == len(kmers) - 1:
            break
        next_kmer = kmers[i + 1]

        if kmer in graph and next_kmer in graph[kmer]:
            letter = graph[kmer][next_kmer]
            new_sequence.append(letter)
        else:
            new_sequence.append(original_seq[i] if i < len(original_seq) else "X")

    return "".join(new_sequence)


def _debruijin_sequence_single(
    sequence: str | ProFormaAnnotation,
    k: int,
    static_residues: str = "",
) -> str:
    """Internal function for creating de Bruijn decoy from a single sequence."""
    annot = get_annotation_input(sequence, copy=False)
    seq_str = annot.sequence

    if len(seq_str) < k:
        return seq_str

    padded_seq = "-" * k + seq_str
    graph = _build_kmer_graph(padded_seq, k)
    _randomize_graph_edges(graph, static_residues)
    annot.sequence = _reconstruct_sequence_from_graph(seq_str, graph, k)
    return annot.serialize()


@overload
def debruijin_sequence(
    sequence: str | ProFormaAnnotation,
    k: int,
    static_residues: str = "",
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> str: ...


@overload
def debruijin_sequence(
    sequence: Sequence[str | ProFormaAnnotation],
    k: int,
    static_residues: str = "",
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[str]: ...


def debruijin_sequence(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    k: int,
    static_residues: str = "",
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> str | list[str]:
    """
    Generate decoy sequences using a de Bruijn graph with randomized edges.

    Builds a k-mer graph from the sequence, randomizes the amino acid transitions
    while respecting static residues, then reconstructs a sequence of the same length
    with the same k-mer patterns but different amino acids.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation]
    :param k: The k-mer size for building the graph.
    :type k: int
    :param static_residues: Amino acids to keep in their original positions.
    :type static_residues: str
    :param n_workers: Number of worker processes. If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk. If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect).
    :type method: Literal["process", "thread"] | None
    :return: The decoy sequence or list of decoy sequences.
    :rtype: str | list[str]

    .. code-block:: python

        # Single sequence
        >>> debruijin_sequence('PEPTIDE', k=3)
        'XEPXIDX'  # Randomized but maintains k-mer structure

        # With static residues
        >>> debruijin_sequence('PEPTIDE', k=3, static_residues='P')
        'PEPXIDX'  # P residues stay in place

        # Multiple sequences
        >>> debruijin_sequence(['PEPTIDE', 'SEQUENCE'], k=3)
        ['XEPXIDX', 'XEQXENX']
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _debruijin_sequence_single,
            sequence,
            k=k,
            static_residues=static_residues,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _debruijin_sequence_single(
            sequence, k=k, static_residues=static_residues
        )
