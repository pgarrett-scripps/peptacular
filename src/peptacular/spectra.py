import bisect
from collections.abc import Callable
from dataclasses import dataclass
from typing import Generator, Literal, NamedTuple, Self, Sequence, overload

from .constants import C13_NEUTRON_MASS, PROTON_MASS


class SpectrumPeak(NamedTuple):
    """Peak with m/z, intensity, and optional charge."""

    mz: float
    intensity: float
    charge: int | None = None
    ion_mobility: float | None = None

    @property
    def neutral_mass(self) -> float | None:
        """Calculate neutral mass if charge is known."""
        if self.charge is None or self.charge == 0:
            return None
        return self.mz * abs(self.charge) - self.charge * PROTON_MASS


class Spectrum:
    """Spectrum data structure with list of peaks."""

    @overload
    def __init__(self, peaks: Sequence[SpectrumPeak]) -> None: ...

    @overload
    def __init__(self, peaks: tuple[Sequence[float], Sequence[float]]) -> None: ...

    @overload
    def __init__(
        self, peaks: tuple[Sequence[float], Sequence[float], Sequence[int | None]]
    ) -> None: ...

    @overload
    def __init__(
        self,
        peaks: tuple[
            Sequence[float],
            Sequence[float],
            Sequence[int | None],
            Sequence[float | None],
        ],
    ) -> None: ...

    def __init__(self, peaks) -> None:
        """
        Initialize spectrum from:
        - List of SpectrumPeak objects
        - Tuple of (mzs, intensities), (mzs, intensities, charges), or (mzs, intensities, charges, ion_mobilities)
        - List of (mz, intensity), (mz, intensity, charge), or (mz, intensity, charge, ion_mobility) tuples
        """
        # Empty case
        if not peaks:
            self.peaks: list[SpectrumPeak] = []
            return

        # Case 1: Already SpectrumPeak objects
        if isinstance(peaks[0], SpectrumPeak):
            self.peaks = list(peaks)  # type: ignore
            return

        # Case 2: Tuple of sequences
        if isinstance(peaks, tuple) and not isinstance(peaks[0], (int, float)):
            n = len(peaks)
            if n == 2:
                mzs, intensities = peaks
                self.peaks = [
                    SpectrumPeak(float(mz), float(intensity))
                    for mz, intensity in zip(mzs, intensities, strict=True)
                ]
            elif n == 3:
                mzs, intensities, charges = peaks
                self.peaks = [
                    SpectrumPeak(
                        float(mz),
                        float(intensity),
                        int(charge) if charge is not None else None,
                    )
                    for mz, intensity, charge in zip(
                        mzs, intensities, charges, strict=True
                    )
                ]
            elif n == 4:
                mzs, intensities, charges, ion_mobilities = peaks
                self.peaks = [
                    SpectrumPeak(
                        float(mz),
                        float(intensity),
                        int(charge) if charge is not None else None,
                        float(ion_mobility) if ion_mobility is not None else None,
                    )
                    for mz, intensity, charge, ion_mobility in zip(
                        mzs, intensities, charges, ion_mobilities, strict=True
                    )
                ]
            else:
                raise ValueError("Tuple of sequences must have 2, 3, or 4 elements")
            return

        # Case 3: List of tuples
        first_item = peaks[0]
        if isinstance(first_item, tuple):
            n = len(first_item)
            if n == 2:
                self.peaks = [
                    SpectrumPeak(float(mz), float(intensity)) for mz, intensity in peaks
                ]  # type: ignore
            elif n == 3:
                self.peaks = [
                    SpectrumPeak(
                        float(mz),
                        float(intensity),
                        int(charge) if charge is not None else None,
                    )
                    for mz, intensity, charge in peaks  # type: ignore
                ]
            elif n == 4:
                self.peaks = [
                    SpectrumPeak(
                        float(mz),
                        float(intensity),
                        int(charge) if charge is not None else None,
                        float(ion_mobility) if ion_mobility is not None else None,
                    )
                    for mz, intensity, charge, ion_mobility in peaks  # type: ignore
                ]
            else:
                raise ValueError("Tuples must have 2, 3, or 4 elements")
            return

        raise ValueError(
            "peaks must be a list of SpectrumPeak objects, list of tuples, or tuple of sequences"
        )

    def __len__(self) -> int:
        return len(self.peaks)

    def __bool__(self) -> bool:
        return len(self.peaks) > 0

    @property
    def mzs(self) -> list[float]:
        return [p.mz for p in self.peaks]

    @property
    def intensities(self) -> list[float]:
        return [p.intensity for p in self.peaks]

    @property
    def charges(self) -> list[int | None]:
        return [p.charge for p in self.peaks]

    @property
    def ion_mobilities(self) -> list[float | None]:
        return [p.ion_mobility for p in self.peaks]

    def copy(self) -> Self:
        return self.__class__(
            peaks=[
                SpectrumPeak(p.mz, p.intensity, p.charge, p.ion_mobility)
                for p in self.peaks
            ]
        )

    def _filter_peaks(
        self, condition: Callable[[SpectrumPeak], bool], inplace: bool = True
    ) -> Self:
        """Filter peaks by condition function."""
        filtered = [p for p in self.peaks if condition(p)]

        if inplace:
            self.peaks = filtered
            return self
        else:
            return self.__class__(peaks=filtered[:])

    def filter_by_intensity(
        self,
        min_intensity: float | None = None,
        max_intensity: float | None = None,
        inplace: bool = True,
    ) -> Self:
        """Filter peaks by intensity range."""

        def condition(p: SpectrumPeak) -> bool:
            if min_intensity is not None and p.intensity < min_intensity:
                return False
            if max_intensity is not None and p.intensity > max_intensity:
                return False
            return True

        return self._filter_peaks(condition, inplace)

    def filter_by_mz(
        self,
        min_mz: float | None = None,
        max_mz: float | None = None,
        inplace: bool = True,
    ) -> Self:
        """Filter peaks by m/z range."""

        def condition(p: SpectrumPeak) -> bool:
            if min_mz is not None and p.mz < min_mz:
                return False
            if max_mz is not None and p.mz > max_mz:
                return False
            return True

        return self._filter_peaks(condition, inplace)

    def filter_by_charge(
        self,
        charges: list[int] | None = None,
        exclude_none: bool = False,
        inplace: bool = True,
    ) -> Self:
        """Filter peaks by charge state."""

        def condition(p: SpectrumPeak) -> bool:
            if exclude_none and p.charge is None:
                return False
            if charges is not None:
                if p.charge is None:
                    return False
                return p.charge in charges
            return True

        return self._filter_peaks(condition, inplace)

    def remove_peaks_within_mz_range(
        self, lower: float, upper: float, inplace: bool = True
    ) -> Self:
        """Remove peaks within m/z range [lower, upper]."""
        return self._filter_peaks(lambda p: not (lower <= p.mz <= upper), inplace)

    def remove_peaks_within_tolerance(
        self,
        target_mz: float,
        tolerance: float,
        tolerance_unit: Literal["ppm", "da"] = "ppm",
        inplace: bool = True,
    ) -> Self:
        """Remove peaks within tolerance of target m/z."""
        lower, upper = self._calculate_mz_window(target_mz, tolerance, tolerance_unit)
        return self.remove_peaks_within_mz_range(lower, upper, inplace)

    def _transform_intensities(
        self, transform_func: Callable[[float], float], inplace: bool = True
    ) -> Self:
        """Apply transformation function to all intensities."""
        transformed = [transform_func(p.intensity) for p in self.peaks]

        new_peaks = [
            SpectrumPeak(p.mz, trans_int, p.charge, p.ion_mobility)
            for p, trans_int in zip(self.peaks, transformed)
        ]

        if inplace:
            self.peaks = new_peaks
            return self
        else:
            return self.__class__(peaks=new_peaks)

    def normalize_intensities(
        self, method: Literal["max", "sum", "tic"] = "max", inplace: bool = True
    ) -> Self:
        """Normalize intensities by 'max', 'sum', or 'tic'."""
        if not self.peaks:
            return self if inplace else self.copy()

        intensities = self.intensities

        if method == "max":
            norm_factor = max(intensities)
        elif method in ("sum", "tic"):
            norm_factor = sum(intensities)
        else:
            raise ValueError(f"Unknown normalization method: {method}")

        if norm_factor == 0:
            return self if inplace else self.copy()

        return self._transform_intensities(lambda i: i / norm_factor, inplace)

    def _sort_peaks(
        self,
        key_func: Callable[[SpectrumPeak], float],
        reverse: bool = False,
        inplace: bool = True,
    ) -> Self:
        """Sort peaks by a key function."""
        sorted_peaks = sorted(self.peaks, key=key_func, reverse=reverse)

        if inplace:
            self.peaks = sorted_peaks
            return self
        else:
            return self.__class__(peaks=sorted_peaks[:])

    def sort_by_mz(self, inplace: bool = True) -> Self:
        """Sort peaks by m/z."""
        return self._sort_peaks(lambda p: p.mz, inplace=inplace)

    def sort_by_intensity(self, reverse: bool = True, inplace: bool = True) -> Self:
        """Sort peaks by intensity (descending by default)."""
        return self._sort_peaks(lambda p: p.intensity, reverse=reverse, inplace=inplace)

    def keep_top_n_peaks(self, n: int, inplace: bool = True) -> Self:
        """Keep only top N most intense peaks."""
        if n <= 0:
            if inplace:
                self.peaks = []
                return self
            else:
                return self.__class__(peaks=[])

        if n >= len(self.peaks):
            return self if inplace else self.copy()

        top_peaks = sorted(self.peaks, key=lambda p: p.intensity, reverse=True)[:n]

        if inplace:
            self.peaks = top_peaks
            return self
        else:
            return self.__class__(peaks=top_peaks[:])

    def has_peak_within_tolerance(
        self,
        target_mz: float,
        tolerance: float,
        tolerance_unit: Literal["ppm", "da"] = "ppm",
    ) -> bool:
        """Check if peak exists within tolerance of target m/z."""
        lower, upper = self._calculate_mz_window(target_mz, tolerance, tolerance_unit)
        return any(lower <= p.mz <= upper for p in self.peaks)

    def find_peak_within_tolerance(
        self,
        target_mz: float,
        tolerance: float,
        tolerance_unit: Literal["ppm", "da"] = "ppm",
    ) -> SpectrumPeak | None:
        """Find closest peak within tolerance of target m/z."""
        lower, upper = self._calculate_mz_window(target_mz, tolerance, tolerance_unit)

        closest_peak = None
        closest_diff = float("inf")

        for peak in self.peaks:
            if lower <= peak.mz <= upper:
                diff = abs(peak.mz - target_mz)
                if diff < closest_diff:
                    closest_diff = diff
                    closest_peak = peak

        return closest_peak

    @staticmethod
    def _calculate_mz_window(
        target_mz: float, tolerance: float, tolerance_unit: Literal["ppm", "da"]
    ) -> tuple[float, float]:
        """Calculate m/z window bounds based on tolerance."""
        if tolerance_unit == "ppm":
            delta = target_mz * tolerance / 1e6
        elif tolerance_unit == "da":
            delta = tolerance
        else:
            raise ValueError("tolerance_unit must be 'ppm' or 'da'")

        return target_mz - delta, target_mz + delta

    def sqrt_transform(self, inplace: bool = True) -> Self:
        """Apply square root transformation to intensities."""
        return self._transform_intensities(lambda i: i**0.5, inplace)

    def log_transform(self, base: float = 2.0, inplace: bool = True) -> Self:
        """Apply log transformation to intensities (adds 1 to avoid log(0))."""
        import math

        return self._transform_intensities(lambda i: math.log(i + 1, base), inplace)

    def get_total_ion_current(self) -> float:
        """Calculate total ion current (sum of all intensities)."""
        return sum(p.intensity for p in self.peaks)

    @property
    def has_charges(self) -> bool:
        """Check if any peaks have assigned charges."""
        return any(p.charge is not None for p in self.peaks)

    def __str__(self) -> str:
        return f"Spectrum with {len(self.peaks)} peaks"

    def __repr__(self) -> str:
        return f"Spectrum(peaks=[{', '.join(repr(p) for p in self.peaks)}])"

    def __iter__(self) -> Generator[SpectrumPeak, None, None]:
        for peak in self.peaks:
            yield peak

    # more dunder
    def __getitem__(self, index: int) -> SpectrumPeak:
        return self.peaks[index]

    def __eq__(self, other: object) -> bool:
        """Check equality based on peaks."""
        if not isinstance(other, Spectrum):
            return NotImplemented
        return self.peaks == other.peaks

    def __contains__(self, item: SpectrumPeak | float) -> bool:
        """Check if peak or m/z exists in spectrum."""
        if isinstance(item, SpectrumPeak):
            return item in self.peaks
        elif isinstance(item, (int, float)):
            return any(abs(p.mz - item) < 1e-6 for p in self.peaks)
        return False

    def __add__(self, other: "Spectrum") -> "Spectrum":
        """Combine two spectra (merge peaks)."""
        if not isinstance(other, Spectrum):
            return NotImplemented
        return self.__class__(peaks=self.peaks + other.peaks)

    def __mul__(self, scalar: float) -> "Spectrum":
        """Scale all intensities by a scalar."""
        if not isinstance(scalar, (int, float)):
            return NotImplemented
        return self._transform_intensities(lambda i: i * scalar, inplace=False)

    def __rmul__(self, scalar: float) -> "Spectrum":
        """Right multiplication (scalar * spectrum)."""
        return self.__mul__(scalar)

    def __reversed__(self) -> Generator[SpectrumPeak, None, None]:
        """Iterate peaks in reverse order."""
        for peak in reversed(self.peaks):
            yield peak

    def merge_peaks_within_tolerance(
        self,
        tolerance: float,
        tolerance_unit: Literal["ppm", "da"] = "ppm",
        inplace: bool = True,
    ) -> Self:
        """
        Merge peaks within tolerance using bisect for binary search.
        """
        if not self.peaks:
            return self if inplace else self.copy()

        spectrum = self if inplace else self.copy()

        # Sort peaks by m/z
        mz_sorted_peaks = sorted(spectrum.peaks, key=lambda p: p.mz)
        mz_values = [p.mz for p in mz_sorted_peaks]  # For binary search

        # Create intensity-sorted indices
        intensity_idx = sorted(
            range(len(mz_sorted_peaks)),
            key=lambda i: mz_sorted_peaks[i].intensity,
            reverse=True,
        )

        merged = [False] * len(mz_sorted_peaks)
        merged_peaks: list[SpectrumPeak] = []

        for base_idx in intensity_idx:
            if merged[base_idx]:
                continue

            base_peak = mz_sorted_peaks[base_idx]
            base_mz = base_peak.mz

            # Calculate tolerance window
            if tolerance_unit == "ppm":
                tol = base_mz * tolerance * 1e-6
            else:
                tol = tolerance

            lower = base_mz - tol
            upper = base_mz + tol

            # Binary search for range
            left_idx = bisect.bisect_left(mz_values, lower)
            right_idx = bisect.bisect_right(mz_values, upper)

            # Merge peaks in range
            merged_intensity = 0.0
            merged_charge = base_peak.charge
            has_inconsistent_charge = False

            for i in range(left_idx, right_idx):
                if not merged[i]:
                    peak = mz_sorted_peaks[i]
                    merged_intensity += peak.intensity
                    merged[i] = True

                    # Check charge consistency
                    if merged_charge is not None and peak.charge is not None:
                        if merged_charge != peak.charge:
                            has_inconsistent_charge = True

            # Create merged peak
            final_charge = None if has_inconsistent_charge else merged_charge
            merged_peaks.append(SpectrumPeak(base_mz, merged_intensity, final_charge))

        # Sort by m/z
        merged_peaks.sort(key=lambda p: p.mz)

        spectrum.peaks = merged_peaks
        return spectrum


@dataclass
class Ms1Spectrum(Spectrum):
    retention_time: float | None = None
    ms1_scan_number: int | None = None

    def __init__(
        self,
        peaks,  # type: ignore
        retention_time: float | None = None,
        ms1_scan_number: int | None = None,
    ):
        super().__init__(peaks)
        self.retention_time = retention_time
        self.ms1_scan_number = ms1_scan_number

    def copy(self) -> Self:
        return self.__class__(
            peaks=[
                SpectrumPeak(p.mz, p.intensity, p.charge, p.ion_mobility)
                for p in self.peaks
            ],
            retention_time=self.retention_time,
            ms1_scan_number=self.ms1_scan_number,
        )


@dataclass
class Ms2Spectrum(Ms1Spectrum):
    """MS2 spectrum with precursor information."""

    precursor_mz: float | None = None
    precursor_charge: int | None = None
    ms2_scan_number: int | None = None
    activation_method: str | None = None
    isolation_window: tuple[float, float] | None = None

    def __init__(
        self,
        peaks,  # type: ignore
        retention_time: float | None = None,
        ms1_scan_number: int | None = None,
        precursor_mz: float | None = None,
        precursor_charge: int | None = None,
        ms2_scan_number: int | None = None,
        activation_method: str | None = None,
        isolation_window: tuple[float, float] | None = None,
    ):
        super().__init__(
            peaks=peaks,
            retention_time=retention_time,
            ms1_scan_number=ms1_scan_number,
        )
        self.precursor_mz = precursor_mz
        self.precursor_charge = precursor_charge
        self.ms2_scan_number = ms2_scan_number
        self.activation_method = activation_method
        self.isolation_window = isolation_window

    def copy(self) -> Self:
        return self.__class__(
            peaks=[
                SpectrumPeak(p.mz, p.intensity, p.charge, p.ion_mobility)
                for p in self.peaks
            ],
            retention_time=self.retention_time,
            ms1_scan_number=self.ms1_scan_number,
            precursor_mz=self.precursor_mz,
            precursor_charge=self.precursor_charge,
            ms2_scan_number=self.ms2_scan_number,
            activation_method=self.activation_method,
            isolation_window=self.isolation_window,
        )

    def remove_precursor_peak(
        self,
        tolerance: float,
        tolerance_unit: Literal["ppm", "da"] = "ppm",
        remove_isotopes: bool = True,
        max_isotopes: int = 10,
        inplace: bool = True,
    ) -> Self:
        """Remove precursor peak and optionally its isotope envelope."""
        if self.precursor_mz is None:
            return self if inplace else self.copy()

        spectrum = self if inplace else self.copy()

        # Remove main precursor peak
        spectrum.remove_peaks_within_tolerance(
            self.precursor_mz, tolerance, tolerance_unit, inplace=True
        )

        # Remove isotopes if requested and charge is known
        if (
            remove_isotopes
            and self.precursor_charge is not None
            and self.precursor_charge != 0
        ):
            isotopic_spacing = 1.003355 / abs(self.precursor_charge)

            # Remove forward isotopes
            for i in range(1, max_isotopes + 1):
                isotope_mz = self.precursor_mz + (i * isotopic_spacing)
                if not spectrum.has_peak_within_tolerance(
                    isotope_mz, tolerance, tolerance_unit
                ):
                    break
                spectrum.remove_peaks_within_tolerance(
                    isotope_mz, tolerance, tolerance_unit, inplace=True
                )

            # Remove backward isotopes
            for i in range(1, max_isotopes + 1):
                isotope_mz = self.precursor_mz - (i * isotopic_spacing)
                if isotope_mz < 0:
                    break
                if not spectrum.has_peak_within_tolerance(
                    isotope_mz, tolerance, tolerance_unit
                ):
                    break
                spectrum.remove_peaks_within_tolerance(
                    isotope_mz, tolerance, tolerance_unit, inplace=True
                )

        return spectrum
