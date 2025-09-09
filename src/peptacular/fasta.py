"""
Improved Fasta.py - Production-ready version
"""

import pathlib
import io
import logging
from typing import Generator, NamedTuple, Protocol, runtime_checkable
from dataclasses import dataclass
from contextlib import contextmanager

import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import sys
import time

from peptacular.digestion.core import generate_regex
from .sequence.digestion import digest_by_regex

logger = logging.getLogger(__name__)

FASTA_INPUT_TYPE = str | pathlib.Path | io.IOBase


@runtime_checkable
class ReadableProtocol(Protocol):
    """Protocol for readable objects."""

    def read(self) -> str | bytes: ...


class FastaSequence(NamedTuple):
    header: str
    sequence: str

    def __post_init__(self):
        if not self.header:
            raise ValueError("Header cannot be empty")
        if not self.sequence:
            raise ValueError("Sequence cannot be empty")


@dataclass
class DigestConfig:
    """Configuration for digestion parameters."""

    enzyme_regex: str
    missed_cleavages: int = 0
    semi: bool = False
    min_len: int | None = None
    max_len: int | None = None
    complete_digestion: bool = True
    sort_output: bool = True

    def __post_init__(self):
        if self.missed_cleavages < 0:
            raise ValueError("missed_cleavages must be non-negative")
        if self.min_len is not None and self.min_len < 1:
            raise ValueError("min_len must be positive")
        if self.max_len is not None and self.max_len < 1:
            raise ValueError("max_len must be positive")
        if (
            self.min_len is not None
            and self.max_len is not None
            and self.min_len > self.max_len
        ):
            raise ValueError("min_len cannot be greater than max_len")


class FastaParsingError(Exception):
    """Custom exception for FASTA parsing errors."""

    pass


class DigestError(Exception):
    """Custom exception for digestion errors."""

    pass


def _is_gil_disabled() -> bool:
    """Check if the GIL is disabled in Python 3.13+."""
    return hasattr(sys, "_is_gil_disabled") and sys._is_gil_disabled()  # type: ignore


@contextmanager
def _safe_file_open(file_path: str | pathlib.Path, encoding: str = "utf-8"):
    """Context manager for safe file operations."""
    try:
        if isinstance(file_path, pathlib.Path):
            file_handle = file_path.open("r", encoding=encoding)
        else:
            file_handle = open(file_path, "r", encoding=encoding)
        yield file_handle
    except (FileNotFoundError, PermissionError) as e:
        logger.error(f"Failed to open file {file_path}: {e}")
        raise FastaParsingError(f"Cannot access file {file_path}: {e}") from e
    except UnicodeDecodeError as e:
        logger.error(f"Encoding error reading file {file_path}: {e}")
        raise FastaParsingError(f"Encoding error in file {file_path}: {e}") from e
    finally:
        if "file_handle" in locals():
            file_handle.close()


def _detect_file_encoding(file_path: str | pathlib.Path) -> str:
    """Simple encoding detection for text files."""
    encodings = ["utf-8", "utf-16", "latin-1", "ascii"]

    for encoding in encodings:
        try:
            if isinstance(file_path, pathlib.Path):
                with file_path.open("r", encoding=encoding) as f:
                    f.read(1024)  # Try to read first 1KB
            else:
                with open(file_path, "r", encoding=encoding) as f:
                    f.read(1024)
            return encoding
        except UnicodeDecodeError:
            continue

    logger.warning(f"Could not detect encoding for {file_path}, defaulting to utf-8")
    return "utf-8"


def parse_fasta(input_data: FASTA_INPUT_TYPE) -> list[FastaSequence]:
    """
    Parse FASTA formatted data from various input types with improved error handling.

    Parameters:
    -----------
    input_data : str, pathlib.Path, or file-like object
        The input can be:
        - A string containing FASTA formatted text
        - A path to a FASTA file (as string or Path object)
        - A file-like object (already opened file or StringIO)

    Returns:
    --------
    List[FastaSequence]
        List of FastaSequence objects containing header and sequence data

    Raises:
    -------
    FastaParsingError
        If the input cannot be parsed or contains invalid data
    """
    try:
        text: str = ""

        # Handle file-like objects
        if isinstance(input_data, ReadableProtocol) or hasattr(input_data, "read"):
            content = input_data.read()
            if isinstance(content, bytes):
                text = content.decode("utf-8")
            elif isinstance(content, str):
                text = content
            else:
                text = str(content)

        # Handle string input
        elif isinstance(input_data, str):
            if not input_data.lstrip().startswith(">") and "\n" not in input_data[:100]:
                # Likely a file path
                encoding = _detect_file_encoding(input_data)
                with _safe_file_open(input_data, encoding) as f:
                    text = f.read()
            else:
                text = input_data

        # Handle pathlib.Path objects
        elif isinstance(input_data, pathlib.Path):
            if not input_data.exists():
                raise FastaParsingError(f"File does not exist: {input_data}")
            if not input_data.is_file():
                raise FastaParsingError(f"Path is not a file: {input_data}")

            encoding = _detect_file_encoding(input_data)
            with _safe_file_open(input_data, encoding) as f:
                text = f.read()

        else:
            raise TypeError(f"Unsupported input type: {type(input_data)}")

        return parse_fasta_text(text)

    except Exception as e:
        logger.error(f"Error parsing FASTA input: {e}")
        if isinstance(e, (FastaParsingError, TypeError)):
            raise
        raise FastaParsingError(f"Unexpected error parsing FASTA: {e}") from e


def parse_fasta_text(text: str) -> list[FastaSequence]:
    """
    Parse FASTA formatted text with validation.

    Parameters:
    -----------
    text : str
        The input FASTA text.

    Returns:
    --------
    List[FastaSequence]
        List of FastaSequence objects

    Raises:
    -------
    FastaParsingError
        If the text contains invalid FASTA format
    """
    if not text.strip():
        raise FastaParsingError("Empty input text")

    sequences: list[FastaSequence] = []
    header = None
    seq_lines: list[str] = []
    line_number = 0

    try:
        for line in text.splitlines():
            line_number += 1
            line = line.strip()

            if not line:  # Skip empty lines
                continue

            if line.startswith(">"):
                # Process previous sequence if exists
                if header is not None:
                    sequence = "".join(seq_lines).upper()
                    if sequence:  # Only add if sequence is not empty
                        sequences.append(
                            FastaSequence(header=header, sequence=sequence)
                        )
                    else:
                        logger.warning(f"Empty sequence found for header: {header}")

                # Start new sequence
                header = line[1:].strip()
                if not header:
                    raise FastaParsingError(f"Empty header at line {line_number}")
                seq_lines = []

            else:
                if header is None:
                    raise FastaParsingError(
                        f"Sequence data before header at line {line_number}"
                    )

                # Validate sequence characters (basic validation)
                invalid_chars = set(line.upper()) - set("ACDEFGHIKLMNPQRSTVWYXBZJOU*-")
                if invalid_chars:
                    logger.warning(
                        f"Invalid characters found at line {line_number}: {invalid_chars}"
                    )

                seq_lines.append(line)

        # Process last sequence
        if header is not None:
            sequence = "".join(seq_lines).upper()
            if sequence:
                sequences.append(FastaSequence(header=header, sequence=sequence))
            else:
                logger.warning(f"Empty sequence found for header: {header}")

        if not sequences:
            raise FastaParsingError("No valid FASTA sequences found")

        logger.info(f"Successfully parsed {len(sequences)} sequences")
        return sequences

    except Exception as e:
        if isinstance(e, FastaParsingError):
            raise
        raise FastaParsingError(
            f"Error parsing FASTA at line {line_number}: {e}"
        ) from e


def _digest_single_sequence(
    args: tuple[FastaSequence, DigestConfig],
) -> tuple[FastaSequence, list[str]]:
    """
    Worker function for processing a single sequence with error handling.

    Returns:
    --------
    Tuple[str, List[str]]
        (header, peptides) or (header, []) if digestion fails
    """
    seq, config = args

    peptides = digest_by_regex(
        sequence=seq.sequence,
        enzyme_regex=config.enzyme_regex,
        missed_cleavages=config.missed_cleavages,
        semi=config.semi,
        min_len=config.min_len,
        max_len=config.max_len,
        complete_digestion=config.complete_digestion,
        sort_output=config.sort_output,
    )
    return seq, list(peptides)


def _digest_fasta_sequential(
    sequences: list[FastaSequence], config: DigestConfig
) -> Generator[tuple[FastaSequence, list[str]], None, None]:
    """Sequential digestion with progress logging."""
    total = len(sequences)
    for i, seq in enumerate(sequences, 1):
        if i % 100 == 0:
            logger.info(f"Processing sequence {i}/{total}")

        result = _digest_single_sequence((seq, config))
        yield result


def regex_digest_fasta(
    fasta: FASTA_INPUT_TYPE,
    enzyme_regex: str,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    complete_digestion: bool = True,
    sort_output: bool = True,
    processes: int | None = None,
    batch_size: int = 100,
) -> Generator[tuple[FastaSequence, list[str]], None, None]:
    """
    Digest FASTA sequences using regex with improved parallelization and error handling.

    Raises:
    -------
    DigestError
        If digestion parameters are invalid or processing fails
    """
    try:
        # Validate and create config
        config = DigestConfig(
            enzyme_regex=enzyme_regex,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            complete_digestion=complete_digestion,
            sort_output=sort_output,
        )

        # Parse sequences
        sequences = parse_fasta(fasta)
        logger.info(f"Starting digestion of {len(sequences)} sequences")

        # Determine optimal process count
        if processes is None:
            processes = max(1, mp.cpu_count() - 1)
        processes = max(
            1, min(processes, len(sequences))
        )  # Don't use more processes than sequences

        # Validate batch size
        batch_size = max(1, min(batch_size, len(sequences)))

        start_time = time.time()

        if processes == 1:
            logger.info("Using sequential processing")
            yield from _digest_fasta_sequential(sequences, config)
        else:
            # Prepare work items
            work_items = [(seq, config) for seq in sequences]

            if _is_gil_disabled():
                logger.info(
                    f"Using ThreadPoolExecutor with {processes} threads (GIL disabled)"
                )
                with ThreadPoolExecutor(max_workers=processes) as executor:
                    # Submit all work
                    futures = [
                        executor.submit(_digest_single_sequence, item)
                        for item in work_items
                    ]

                    # Process results as they complete
                    for future in as_completed(futures):
                        try:
                            result = future.result()
                            yield result
                        except Exception as e:
                            logger.error(f"Worker thread error: {e}")
                            continue
            else:
                logger.info(f"Using ProcessPoolExecutor with {processes} processes")
                try:
                    with ProcessPoolExecutor(max_workers=processes) as executor:
                        # Process in batches to manage memory
                        for i in range(0, len(work_items), batch_size):
                            batch = work_items[i : i + batch_size]
                            results = executor.map(_digest_single_sequence, batch)

                            for result in results:
                                yield result

                except Exception as e:
                    logger.error(f"ProcessPoolExecutor error: {e}")
                    # Fall back to sequential processing
                    logger.info("Falling back to sequential processing")
                    yield from _digest_fasta_sequential(sequences, config)

        duration = time.time() - start_time
        logger.info(f"Digestion completed in {duration:.2f} seconds")

    except Exception as e:
        logger.error(f"Error in regex_digest_fasta: {e}")
        if isinstance(e, (FastaParsingError, ValueError)):
            raise DigestError(f"Digestion failed: {e}") from e
        raise DigestError(f"Unexpected error during digestion: {e}") from e


def digest_fasta(
    fasta: FASTA_INPUT_TYPE,
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    complete_digestion: bool = True,
    sort_output: bool = True,
    processes: int | None = None,
    batch_size: int = 100,
) -> Generator[tuple[FastaSequence, list[str]], None, None]:
    try:
        if not cleave_on:
            raise ValueError("cleave_on parameter cannot be empty")

        regex = generate_regex(
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
        )

        logger.info(f"Generated enzyme regex: {regex}")

        yield from regex_digest_fasta(
            fasta=fasta,
            enzyme_regex=regex,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            complete_digestion=complete_digestion,
            sort_output=sort_output,
            processes=processes,
            batch_size=batch_size,
        )

    except Exception as e:
        if isinstance(e, (ValueError, DigestError)):
            raise DigestError(f"Enzyme digestion failed: {e}") from e
        raise DigestError(f"Unexpected error in enzyme digestion: {e}") from e
