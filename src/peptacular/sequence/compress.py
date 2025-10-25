import base64
from collections import Counter
from typing import Sequence

from ..funcs import compress_with_method, decompress_with_method
from ..sequence.util import get_annotation_input


def _create_peptide_compact_format(
    peptides: Sequence[str], precision: int | None = 3
) -> str:
    """Create compact peptide representation (your original logic)"""

    peptide_counts = Counter(peptides)
    annots_to_charge: dict[str, Counter[str]] = {}

    for peptide, count in peptide_counts.items():
        annot = get_annotation_input(peptide)
        charge_suffix = annot.serialize_charge()
        annot.charge = None
        annot.charge_adducts = None
        seq = annot.serialize(precision=precision)

        if seq not in annots_to_charge:
            annots_to_charge[seq] = Counter()
        annots_to_charge[seq][charge_suffix] += count

    # Build compressed string: PEPTIDE;2:23;3:24,REPFYD;4:13;6:24,...
    compressed_parts: list[str] = []
    for seq, charge_counts in annots_to_charge.items():
        charge_strings: list[str] = []
        for charge, count in sorted(charge_counts.items()):
            charge_strings.append(f"{charge}:{count}")
        peptide_part = seq + ";" + ";".join(charge_strings)
        compressed_parts.append(peptide_part)

    return ",".join(compressed_parts)


def _parse_peptide_compact_format(compressed_str: str) -> list[str]:
    """Parse compact peptide representation (your original logic)"""
    if not compressed_str:
        return []

    peptides: list[str] = []
    for peptide_part in compressed_str.split(","):
        if not peptide_part:
            continue

        parts = peptide_part.split(";")
        seq = parts[0]

        for charge_count in parts[1:]:
            charge_str, count_str = charge_count.split(":")
            count = int(count_str)
            peptide_with_charge = f"{seq}{charge_str}"
            peptides.extend([peptide_with_charge] * count)

    return peptides


def compress_peptides(
    peptides: Sequence[str],
    precision: int | None = 3,
    compression_method: str = "gzip",
    url_safe: bool = False,
) -> str:
    """
    Compress peptides into a compact representation with additional compression

    Args:
        peptides: Sequence of peptide strings
        precision: Precision for peptide serialization (None for full precision)
        compression_method: 'gzip', 'zlib', or 'brotli'
        url_safe: If True, use URL-safe base64 encoding

    Returns:
        Compressed string with method flags
    """
    if not peptides:
        return ""

    # Validate inputs
    if compression_method not in ["gzip", "zlib", "brotli"]:
        raise ValueError("compression_method must be 'gzip', 'zlib', or 'brotli'")

    if precision is not None and (not isinstance(precision, int) or precision < 0):  # type: ignore
        raise ValueError("precision must be non-negative integer or None")

    # Create compact peptide format
    compact_str = _create_peptide_compact_format(peptides, precision)

    # Convert to bytes
    compact_bytes = compact_str.encode("utf-8")

    # Compress
    compressed_bytes = compress_with_method(compact_bytes, compression_method)

    # Create method flag
    method_flags = {"gzip": "G", "zlib": "Z", "brotli": "B"}

    encoding_flag = "U" if url_safe else "S"  # S for standard base64
    method_flag = encoding_flag + method_flags[compression_method]

    # Encode to string
    if url_safe:
        encoded = base64.urlsafe_b64encode(compressed_bytes).decode("ascii")
    else:
        encoded = base64.b64encode(compressed_bytes).decode("ascii")

    return method_flag + encoded


def decompress_peptides(compressed_str: str) -> list[str]:
    """
    Decompress peptides from compressed representation

    Args:
        compressed_str: Compressed string created by compress_peptides

    Returns:
        List of peptide strings
    """
    if not compressed_str:
        return []

    if len(compressed_str) < 3:
        raise ValueError("Invalid compressed string format")

    # Parse method flags
    encoding_flag = compressed_str[0]  # U or S
    compression_flag = compressed_str[1]  # G, Z, or B
    encoded_data = compressed_str[2:]

    # Validate flags
    if encoding_flag not in ["U", "S"]:
        raise ValueError(f"Unknown encoding flag: {encoding_flag}")

    if compression_flag not in ["G", "Z", "B"]:
        raise ValueError(f"Unknown compression flag: {compression_flag}")

    # Map flags to methods
    flag_to_method = {"G": "gzip", "Z": "zlib", "B": "brotli"}
    compression_method = flag_to_method[compression_flag]

    try:
        # Decode from base64
        if encoding_flag == "U":
            compressed_bytes = base64.urlsafe_b64decode(encoded_data)
        else:  # 'S'
            compressed_bytes = base64.b64decode(encoded_data)

        # Decompress
        compact_bytes = decompress_with_method(compressed_bytes, compression_method)

        # Convert back to string
        compact_str = compact_bytes.decode("utf-8")

        # Parse compact format
        return _parse_peptide_compact_format(compact_str)

    except Exception as e:
        raise ValueError(f"Failed to decompress peptide data: {e}")
