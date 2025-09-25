import gzip
import zlib


def compress_with_method(data_bytes: bytes, method: str) -> bytes:
    """Compress data using specified method"""
    if method == "gzip":
        return gzip.compress(data_bytes, compresslevel=9)
    elif method == "zlib":
        return zlib.compress(data_bytes, level=zlib.Z_BEST_COMPRESSION)
    elif method == "brotli":
        try:
            import brotli

            return brotli.compress(data_bytes, quality=11)
        except ImportError:
            raise ImportError(
                "brotli library not available. Install with: pip install brotli"
            )
    else:
        raise ValueError(f"Unknown compression method: {method}")


def decompress_with_method(data_bytes: bytes, method: str) -> bytes:
    """Decompress data using specified method"""
    if method == "gzip":
        return gzip.decompress(data_bytes)
    elif method == "zlib":
        return zlib.decompress(data_bytes)
    elif method == "brotli":
        try:
            import brotli

            return brotli.decompress(data_bytes)
        except ImportError:
            raise ImportError("brotli library not available")
    else:
        raise ValueError(f"Unknown compression method: {method}")
