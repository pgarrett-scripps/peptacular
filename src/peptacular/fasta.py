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
    def read(self) -> str | bytes: ...


class FastaSequence(NamedTuple):
    header: str
    sequence: str


def _detect_file_encoding(file_path: str | pathlib.Path) -> str:
    encodings = ["utf-8", "utf-16", "latin-1", "ascii"]

    for encoding in encodings:
        try:
            if isinstance(file_path, pathlib.Path):
                with file_path.open("r", encoding=encoding) as f:
                    f.read(1024)
            else:
                with open(file_path, encoding=encoding) as f:
                    f.read(1024)
            return encoding
        except UnicodeDecodeError:
            continue

    logger.warning(f"Could not detect encoding for {file_path}, defaulting to utf-8")
    return "utf-8"


def parse_fasta(input_data: FASTA_INPUT_TYPE) -> list[FastaSequence]:
    """
    Parse FASTA formatted data from various input types.
    """
    try:
        text: str = ""

        if isinstance(input_data, ReadableProtocol) or hasattr(input_data, "read"):
            content = input_data.read()  # type: ignore
            text = content.decode("utf-8") if isinstance(content, bytes) else str(content)

        elif isinstance(input_data, str):
            if not input_data.lstrip().startswith(">") and "\n" not in input_data[:100]:
                encoding = _detect_file_encoding(input_data)
                with open(input_data, encoding=encoding) as f:
                    text = f.read()
            else:
                text = input_data

        elif isinstance(input_data, pathlib.Path):
            if not input_data.exists():
                raise FileNotFoundError(f"File does not exist: {input_data}")
            if not input_data.is_file():
                raise FileNotFoundError(f"Path is not a file: {input_data}")

            encoding = _detect_file_encoding(input_data)
            with input_data.open("r", encoding=encoding) as f:
                text = f.read()

        else:
            raise TypeError(f"Unsupported input type: {type(input_data)}")

        return parse_fasta_text(text)

    except Exception as e:
        logger.error(f"Error parsing FASTA input: {e}")
        if isinstance(e, (FileNotFoundError, TypeError)):
            raise
        raise ValueError(f"Unexpected error parsing FASTA: {e}") from e


def parse_fasta_text(text: str) -> list[FastaSequence]:
    """
    Parse FASTA formatted text with validation.
    """
    if not text.strip():
        raise ValueError("Empty input text")

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
                raise ValueError("Empty header found")
            seq_lines = []

        else:
            if header is None:
                raise ValueError("Sequence data before header")
            seq_lines.append(line)

    # Save last sequence
    if header is not None:
        sequence = "".join(seq_lines).upper()
        if sequence:
            sequences.append(FastaSequence(header=header, sequence=sequence))

    if not sequences:
        raise ValueError("No valid FASTA sequences found")

    logger.info(f"Successfully parsed {len(sequences)} sequences")
    return sequences
