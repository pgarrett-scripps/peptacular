import struct
import base64
from typing import Generator, Sequence
from ..funcs import compress_with_method, decompress_with_method


def _float_to_hex(f: float) -> str:
    return format(struct.unpack("!I", struct.pack("!f", f))[0], "08x")


def _hex_to_float(s: str) -> float:
    return struct.unpack("!f", struct.pack("!I", int(s, 16)))[0]


def _encode_leading_zero(lz: int) -> str:
    if 0 <= lz < 16:
        return hex(lz)[-1]
    raise ValueError(f"Leading zero count {lz} out of range [0-15]")


def _decode_leading_zero(lz: str) -> int:
    return int(lz, 16)


def _hex_delta(a: str, b: str) -> str:
    diff = int(a, 16) - int(b, 16)
    return format(diff & 0xFFFFFFFF, "08x")


def _hex_delta_rev(a: str, b: str) -> str:
    diff = int(a, 16) + int(b, 16)
    return format(diff & 0xFFFFFFFF, "08x")


def _count_leading_zeros(s: str) -> int:
    return len(s) - len(s.lstrip("0"))


def _delta_encode_single_string(vals: Sequence[float]) -> str:
    if not vals:
        return ""

    mzs_hex = [_float_to_hex(mz) for mz in vals]
    initial_hex_value = mzs_hex[0]
    initial_hex_value_zeros = _count_leading_zeros(initial_hex_value)

    mzs_hex_deltas: list[str] = []
    leading_zeros: list[int] = []

    for i in range(1, len(mzs_hex)):
        delta = _hex_delta(mzs_hex[i], mzs_hex[i - 1])
        mzs_hex_deltas.append(delta)
        leading_zeros.append(_count_leading_zeros(delta))

    hex_delta_str = initial_hex_value.lstrip("0") + "".join(
        delta.lstrip("0") for delta in mzs_hex_deltas
    )
    leading_zero_str = _encode_leading_zero(initial_hex_value_zeros) + "".join(
        _encode_leading_zero(lz) for lz in leading_zeros
    )

    return hex_delta_str + leading_zero_str[::-1]


def _delta_decode_single_string(s: str) -> Generator[float, None, None]:
    if not s:
        return

    initial_lz = _decode_leading_zero(s[-1])
    initial_hex = "0" * initial_lz + s[: 8 - initial_lz]
    s = s[8 - initial_lz : -1]
    yield _hex_to_float(initial_hex)

    curr_value = initial_hex
    while s:
        lz = _decode_leading_zero(s[-1])
        hex_diff = "0" * lz + s[: 8 - lz]
        hex_val = _hex_delta_rev(curr_value, hex_diff)
        curr_value = hex_val
        s = s[8 - lz : -1]
        yield _hex_to_float(hex_val)


def _hex_encode(intensities: Sequence[float]) -> str:
    return "".join(_float_to_hex(intensity) for intensity in intensities)


def _hex_decode(s: str) -> Generator[float, None, None]:
    for i in range(0, len(s), 8):
        yield _hex_to_float(s[i : i + 8])


def _encode_charges(charges: Sequence[int | None]) -> str:
    """Encode charges as a compact string. None is encoded as 0, values 1-15 as hex digits."""
    if not charges:
        return ""

    result: list[str] = []
    for charge in charges:
        if charge is None:
            result.append("0")
        elif 1 <= charge <= 15:
            result.append(hex(charge)[-1])
        else:
            raise ValueError(f"Charge {charge} out of range [1-15] or None")

    return "".join(result)


def _decode_charges(s: str) -> Generator[int | None, None, None]:
    """Decode charges from compact string. '0' decodes to None, other hex digits to int."""
    if not s:
        return

    for char in s:
        val = int(char, 16)
        if val == 0:
            yield None
        else:
            yield val


def _validate_inputs(
    spectra: tuple[Sequence[float], Sequence[float]]
    | tuple[Sequence[float], Sequence[float], Sequence[int | None]],
    mz_precision: int | None = None,
    intensity_precision: int | None = None,
) -> None:
    if not isinstance(spectra, tuple) or len(spectra) not in [2, 3]:  # type: ignore
        raise ValueError(
            "spectra must be (mz_values, intensity_values) or (mz_values, intensity_values, charges) tuple"
        )

    mzs = spectra[0]
    intensities = spectra[1]

    if not isinstance(mzs, list) or not isinstance(intensities, list):
        raise ValueError("mz_values and intensity_values must be lists")

    if len(mzs) != len(intensities):
        raise ValueError(
            f"Length mismatch: {len(mzs)} m/z values vs {len(intensities)} intensities"
        )

    if len(spectra) == 3:
        charges = spectra[2]
        if not isinstance(charges, list):
            raise ValueError("charges must be a list")
        if len(charges) != len(mzs):
            raise ValueError(
                f"Length mismatch: {len(mzs)} m/z values vs {len(charges)} charges"
            )

    if mz_precision is not None and (
        not isinstance(mz_precision, int) or mz_precision < 0
    ):  # type: ignore
        raise ValueError("mz_precision must be non-negative integer or None")

    if intensity_precision is not None and (
        not isinstance(intensity_precision, int) or intensity_precision < 0
    ):  # type: ignore
        raise ValueError("intensity_precision must be non-negative integer or None")


def _encode_binary_payload(
    mz_str: str, intensity_str: str, charge_str: str = ""
) -> bytes:
    """Encode mz, intensity, and optionally charge data into binary payload."""
    mz_bytes = mz_str.encode("ascii")
    intensity_bytes = intensity_str.encode("ascii")
    charge_bytes = charge_str.encode("ascii") if charge_str else b""

    # Pack: mz_length (4 bytes) + mz_data + intensity_length (4 bytes) + intensity_data + charge_data
    return (
        struct.pack("!I", len(mz_bytes))
        + mz_bytes
        + struct.pack("!I", len(intensity_bytes))
        + intensity_bytes
        + charge_bytes
    )


def _decode_binary_payload(payload: bytes) -> tuple[str, str, str]:
    """Decode binary payload into mz, intensity, and charge strings."""
    if len(payload) < 8:
        raise ValueError("Invalid binary payload: too short")

    # Read m/z length and data
    mz_length = struct.unpack("!I", payload[:4])[0]
    if len(payload) < 4 + mz_length + 4:
        raise ValueError("Invalid binary payload: truncated m/z or intensity data")

    mz_str = payload[4 : 4 + mz_length].decode("ascii")

    # Read intensity length and data
    intensity_length = struct.unpack("!I", payload[4 + mz_length : 8 + mz_length])[0]
    if len(payload) < 8 + mz_length + intensity_length:
        raise ValueError("Invalid binary payload: truncated intensity data")

    intensity_str = payload[8 + mz_length : 8 + mz_length + intensity_length].decode(
        "ascii"
    )

    # Read charge data if present
    charge_str = ""
    if len(payload) > 8 + mz_length + intensity_length:
        charge_str = payload[8 + mz_length + intensity_length :].decode("ascii")

    return mz_str, intensity_str, charge_str


def compress_spectra(
    spectra: tuple[Sequence[float], Sequence[float]]
    | tuple[Sequence[float], Sequence[float], Sequence[int | None]],
    url_safe: bool = False,
    mz_precision: int | None = None,
    intensity_precision: int | None = None,
    compression: str = "gzip",
) -> str:
    """Compress spectra data with configurable precision and compression. Optionally include charges."""
    _validate_inputs(spectra, mz_precision, intensity_precision)

    if compression not in ["gzip", "zlib", "brotli"]:
        raise ValueError("compression must be 'gzip', 'zlib', or 'brotli'")

    mzs = spectra[0]
    intensities = spectra[1]
    charges = spectra[2] if len(spectra) == 3 else None

    if mz_precision is not None:
        mzs = [round(mz, mz_precision) for mz in mzs]

    if intensity_precision is not None:
        intensities = [
            round(intensity, intensity_precision) for intensity in intensities
        ]

    mz_str = _delta_encode_single_string(mzs) if mzs else ""
    intensity_str = _hex_encode(intensities) if intensities else ""
    charge_str = _encode_charges(charges) if charges else ""

    binary_payload = _encode_binary_payload(mz_str, intensity_str, charge_str)
    compressed_bytes = compress_with_method(binary_payload, compression)

    compression_flag = {"gzip": "G", "zlib": "Z", "brotli": "R"}[compression]

    if url_safe:
        encoded = base64.urlsafe_b64encode(compressed_bytes).decode("ascii")
        return "U" + compression_flag + encoded
    else:
        encoded = base64.b85encode(compressed_bytes).decode("ascii")
        return "B" + compression_flag + encoded


def decompress_spectra(
    compressed_str: str,
) -> (
    tuple[list[float], list[float]] | tuple[list[float], list[float], list[int | None]]
):
    """Decompress spectra data with automatic compression detection. Returns charges if present."""
    if not compressed_str:
        return [], []

    if not isinstance(compressed_str, str):  # type: ignore
        raise ValueError("compressed_str must be a string")

    if len(compressed_str) < 3:
        raise ValueError("Invalid compressed string format")

    encoding_flag = compressed_str[0]
    compression_flag = compressed_str[1]
    encoded_data = compressed_str[2:]

    if encoding_flag not in ["U", "B"]:
        raise ValueError(f"Unknown encoding method: {encoding_flag}")

    if compression_flag not in ["G", "Z", "R"]:
        raise ValueError(f"Unknown compression method: {compression_flag}")

    compression_scheme = {"G": "gzip", "Z": "zlib", "R": "brotli"}[compression_flag]

    if encoding_flag == "U":
        compressed_bytes = base64.urlsafe_b64decode(encoded_data)
    else:
        compressed_bytes = base64.b85decode(encoded_data)

    binary_payload = decompress_with_method(compressed_bytes, compression_scheme)
    mz_str, intensity_str, charge_str = _decode_binary_payload(binary_payload)

    mzs = list(_delta_decode_single_string(mz_str)) if mz_str else []
    intensities = list(_hex_decode(intensity_str)) if intensity_str else []

    if charge_str:
        charges = list(_decode_charges(charge_str))
        return mzs, intensities, charges
    else:
        return mzs, intensities
