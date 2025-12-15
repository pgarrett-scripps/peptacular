"""
Simplified FASTA parser - production-ready version
"""

import io
import logging
import pathlib
from typing import NamedTuple, Protocol, runtime_checkable

logger = logging.getLogger(__name__)

FASTA_INPUT_TYPE = str | pathlib.Path | io.IOBase


@runtime_checkable
class ReadableProtocol(Protocol):
    """Protocol for readable objects."""

    def read(self) -> str | bytes: ...


class FastaSequence(NamedTuple):
    header: str
    sequence: str


class FastaParsingError(Exception):
    """Custom exception for FASTA parsing errors."""

    pass


def _detect_file_encoding(file_path: str | pathlib.Path) -> str:
    """Simple encoding detection for text files."""
    encodings = ["utf-8", "utf-16", "latin-1", "ascii"]

    for encoding in encodings:
        try:
            if isinstance(file_path, pathlib.Path):
                with file_path.open("r", encoding=encoding) as f:
                    f.read(1024)
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
    Parse FASTA formatted data from various input types.

    Parameters:
    -----------
    input_data : str, pathlib.Path, or file-like object
        The input can be:
        - A string containing FASTA formatted text
        - A path to a FASTA file (as string or Path object)
        - A file-like object (already opened file or StringIO)

    Returns:
    --------
    list[FastaSequence]
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
            text = (
                content.decode("utf-8") if isinstance(content, bytes) else str(content)
            )

        # Handle string input
        elif isinstance(input_data, str):
            # Check if it's a file path
            if not input_data.lstrip().startswith(">") and "\n" not in input_data[:100]:
                encoding = _detect_file_encoding(input_data)
                with open(input_data, "r", encoding=encoding) as f:
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
            with input_data.open("r", encoding=encoding) as f:
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
    list[FastaSequence]
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

    for line in text.splitlines():
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            # Save previous sequence
            if header is not None:
                sequence = "".join(seq_lines).upper()
                if sequence:
                    sequences.append(FastaSequence(header=header, sequence=sequence))

            # Start new sequence
            header = line[1:].strip()
            if not header:
                raise FastaParsingError("Empty header found")
            seq_lines = []

        else:
            if header is None:
                raise FastaParsingError("Sequence data before header")
            seq_lines.append(line)

    # Save last sequence
    if header is not None:
        sequence = "".join(seq_lines).upper()
        if sequence:
            sequences.append(FastaSequence(header=header, sequence=sequence))

    if not sequences:
        raise FastaParsingError("No valid FASTA sequences found")

    logger.info(f"Successfully parsed {len(sequences)} sequences")
    return sequences
