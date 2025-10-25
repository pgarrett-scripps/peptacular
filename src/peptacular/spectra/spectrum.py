from collections.abc import Callable
from dataclasses import dataclass
from typing import Generator, Literal, Self, Sequence

from ..constants import C13_NEUTRON_MASS
from ..isotope import IsotopeLookup
from .compress import compress_spectra, decompress_spectra
from .decon import SpectrumPeak, deconvolute


class Spectrum:
    """Spectrum data structure with list of peaks."""

    def __init__(
        self,
        peaks: Sequence[SpectrumPeak]
        | tuple[Sequence[float], Sequence[float]]
        | tuple[Sequence[float], Sequence[float], Sequence[int | None]]
        | Sequence[tuple[float, float]]
        | Sequence[tuple[float, float, int | None]],
    ) -> None:
        """
        Initialize spectrum from:
        - List of Peak objects
        - Tuple of (mzs, intensities) or (mzs, intensities, charges)
        - List of (mz, intensity) or (mz, intensity, charge) tuples
        """
        # Case 1: Tuple of sequences - (mzs, intensities) or (mzs, intensities, charges)
        if isinstance(peaks, tuple) and not isinstance(peaks[0], (int, float)):
            if len(peaks) == 2:
                mzs, intensities = peaks
                if len(mzs) != len(intensities):
                    raise ValueError("mzs and intensities must be the same length")
                self.peaks: list[SpectrumPeak] = [
                    SpectrumPeak(mz, intensity)
                    for mz, intensity in zip(mzs, intensities)
                ]
            elif len(peaks) == 3:
                mzs, intensities, charges = peaks
                if not (len(mzs) == len(intensities) == len(charges)):
                    raise ValueError(
                        "mzs, intensities, and charges must be the same length"
                    )
                self.peaks: list[SpectrumPeak] = [
                    SpectrumPeak(mz, intensity, charge)
                    for mz, intensity, charge in zip(mzs, intensities, charges)
                ]
            else:
                raise ValueError(
                    "If peaks is a tuple of sequences, it must be (mzs, intensities) or (mzs, intensities, charges)"
                )

        # Case 2: List of Peak objects or list of tuples
        else:
            first_item = peaks[0] if peaks else None

            # Check if it's already Peak objects
            if isinstance(first_item, SpectrumPeak):
                self.peaks = list(peaks)

            # Check if it's tuples
            elif isinstance(first_item, tuple):
                if len(first_item) == 2:
                    self.peaks = [
                        SpectrumPeak(mz, intensity) for mz, intensity in peaks
                    ]
                elif len(first_item) == 3:
                    self.peaks = [
                        SpectrumPeak(mz, intensity, charge)
                        for mz, intensity, charge in peaks
                    ]
                else:
                    raise ValueError(
                        "Tuples must be (mz, intensity) or (mz, intensity, charge)"
                    )

            # Empty list
            elif first_item is None:
                self.peaks = []

            else:
                raise ValueError(
                    "peaks must be a list of Peak objects, list of tuples, or tuple of sequences"
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

    def copy(self) -> Self:
        return self.__class__(
            peaks=[SpectrumPeak(p.mz, p.intensity, p.charge) for p in self.peaks]
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

    def normalize_intensities(
        self, method: Literal["max", "sum", "tic"] = "max", inplace: bool = True
    ) -> Self:
        """Normalize intensities by 'max', 'sum', or 'tic'."""
        if not self.peaks:
            return self if inplace else self.copy()

        intensities = [p.intensity for p in self.peaks]

        if method == "max":
            max_intensity = max(intensities)
            normalized = (
                [i / max_intensity for i in intensities]
                if max_intensity > 0
                else intensities
            )
        elif method in ("sum", "tic"):
            total_intensity = sum(intensities)
            normalized = (
                [i / total_intensity for i in intensities]
                if total_intensity > 0
                else intensities
            )
        else:
            raise ValueError(f"Unknown normalization method: {method}")

        if inplace:
            self.peaks = [
                SpectrumPeak(p.mz, norm_int, p.charge)
                for p, norm_int in zip(self.peaks, normalized)
            ]
            return self
        else:
            new_peaks = [
                SpectrumPeak(p.mz, norm_int, p.charge)
                for p, norm_int in zip(self.peaks, normalized)
            ]
            return self.__class__(peaks=new_peaks)

    def sort_by_mz(self, inplace: bool = True) -> Self:
        """Sort peaks by m/z."""
        sorted_peaks = sorted(self.peaks, key=lambda p: p.mz)

        if inplace:
            self.peaks = sorted_peaks
            return self
        else:
            return self.__class__(peaks=sorted_peaks[:])

    def sort_by_intensity(self, reverse: bool = True, inplace: bool = True) -> Self:
        """Sort peaks by intensity (descending by default)."""
        sorted_peaks = sorted(self.peaks, key=lambda p: p.intensity, reverse=reverse)

        if inplace:
            self.peaks = sorted_peaks
            return self
        else:
            return self.__class__(peaks=sorted_peaks[:])

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
        transformed = [p.intensity**0.5 for p in self.peaks]

        if inplace:
            self.peaks = [
                SpectrumPeak(p.mz, trans_int, p.charge)
                for p, trans_int in zip(self.peaks, transformed)
            ]
            return self
        else:
            new_peaks = [
                SpectrumPeak(p.mz, trans_int, p.charge)
                for p, trans_int in zip(self.peaks, transformed)
            ]
            return self.__class__(peaks=new_peaks)

    def log_transform(self, base: float = 2.0, inplace: bool = True) -> Self:
        """Apply log transformation to intensities (adds 1 to avoid log(0))."""
        import math

        transformed = [math.log(p.intensity + 1, base) for p in self.peaks]

        if inplace:
            self.peaks = [
                SpectrumPeak(p.mz, trans_int, p.charge)
                for p, trans_int in zip(self.peaks, transformed)
            ]
            return self
        else:
            new_peaks = [
                SpectrumPeak(p.mz, trans_int, p.charge)
                for p, trans_int in zip(self.peaks, transformed)
            ]
            return self.__class__(peaks=new_peaks)

    def get_total_ion_current(self) -> float:
        """Calculate total ion current (sum of all intensities)."""
        return sum(p.intensity for p in self.peaks)

    @property
    def has_charges(self) -> bool:
        """Check if any peaks have assigned charges."""
        return any(p.charge is not None for p in self.peaks)

    def deconvolute(
        self,
        tolerance: float = 50,
        tolerance_type: Literal["ppm", "da"] = "ppm",
        charge_range: tuple[int, int] = (1, 3),
        max_left_decrease: float = 0.6,
        max_right_decrease: float = 0.9,
        isotope_mass: float = C13_NEUTRON_MASS,
        isotope_lookup: IsotopeLookup | None = None,
        inplace: bool = False,
    ) -> Self:
        """Deconvolute spectrum to identify isotopic envelopes."""
        if self.has_charges:
            raise ValueError(
                "Spectrum already has assigned charges; cannot deconvolute."
            )

        if not inplace:
            return self.copy().deconvolute(
                tolerance=tolerance,
                tolerance_type=tolerance_type,
                charge_range=charge_range,
                max_left_decrease=max_left_decrease,
                max_right_decrease=max_right_decrease,
                isotope_mass=isotope_mass,
                isotope_lookup=isotope_lookup,
                inplace=True,
            )

        peak_tupls = [(p.mz, p.intensity) for p in self.peaks]

        dpeaks = deconvolute(
            peaks=peak_tupls,
            tolerance=tolerance,
            tolerance_type=tolerance_type,
            charge_range=charge_range,
            max_left_decrease=max_left_decrease,
            max_right_decrease=max_right_decrease,
            isotope_mass=isotope_mass,
            isotope_lookup=isotope_lookup,
        )

        # update charges in original peaks
        mzs = [p.base_peak.mz for p in dpeaks]
        ints = [p.base_peak.intensity for p in dpeaks]
        charges = [p.charge for p in dpeaks]

        # make new peaks list friom dpeaks
        new_peaks: list[SpectrumPeak] = []
        for mz, intensity, charge in zip(mzs, ints, charges):
            new_peaks.append(SpectrumPeak(mz, intensity, charge))

        # sort
        new_peaks.sort(key=lambda p: p.mz)

        self.peaks = new_peaks

        return self

    def __str__(self) -> str:
        return f"Spectrum with {len(self.peaks)} peaks"

    def __repr__(self) -> str:
        return f"Spectrum(peaks=[{', '.join(repr(p) for p in self.peaks)}])"

    def compress(
        self,
        compression: str = "gzip",
        mz_precision: int | None = None,
        intensity_precision: int | None = None,
        url_safe: bool = False,
    ) -> str:
        """Compress spectrum to binary payload."""
        if self.has_charges:
            return compress_spectra(
                spectra=(self.mzs, self.intensities, self.charges),
                compression=compression,
                mz_precision=mz_precision,
                intensity_precision=intensity_precision,
                url_safe=url_safe,
            )
        else:
            return compress_spectra(
                spectra=(self.mzs, self.intensities),
                compression=compression,
                mz_precision=mz_precision,
                intensity_precision=intensity_precision,
                url_safe=url_safe,
            )

    @staticmethod
    def decompress(compressed_str: str) -> "Spectrum":
        """Decompress binary payload to Spectrum."""

        elems = decompress_spectra(compressed_str)

        if len(elems) == 2:
            mzs, intensities = elems
            return Spectrum(peaks=(mzs, intensities))
        elif len(elems) == 3:
            mzs, intensities, charges = elems
            return Spectrum(peaks=(mzs, intensities, charges))
        else:
            raise ValueError(
                "Decompressed data must be (mzs, intensities) or (mzs, intensities, charges)"
            )

        # dunder iter method  so i can use for peak in spectrum:

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
        new_peaks = [
            SpectrumPeak(p.mz, p.intensity * scalar, p.charge) for p in self.peaks
        ]
        return self.__class__(peaks=new_peaks)

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
        Merge peaks within tolerance by combining them into the most abundant peak.

        Iterates through peaks in order of decreasing intensity, merging nearby peaks
        into the most abundant one and removing the merged peaks.

        Args:
            tolerance: Tolerance value for merging
            tolerance_unit: Unit of tolerance ('ppm' or 'da')
            inplace: Whether to modify spectrum in place

        Returns:
            Modified spectrum
        """
        if not self.peaks:
            return self if inplace else self.copy()

        spectrum = self if inplace else self.copy()

        # Sort peaks by intensity (descending)
        sorted_peaks = sorted(spectrum.peaks, key=lambda p: p.intensity, reverse=True)

        merged_peaks: list[SpectrumPeak] = []
        merged_indices: set[int] = set()

        # Create index mapping for original peaks
        peak_to_idx = {id(peak): i for i, peak in enumerate(sorted_peaks)}

        for i, base_peak in enumerate(sorted_peaks):
            if i in merged_indices:
                continue

            # Calculate tolerance window for this peak
            lower, upper = self._calculate_mz_window(
                base_peak.mz, tolerance, tolerance_unit
            )

            # Find all peaks within tolerance
            merged_intensity = base_peak.intensity
            merged_charge = base_peak.charge

            for j, other_peak in enumerate(sorted_peaks[i + 1 :], start=i + 1):
                if j in merged_indices:
                    continue

                if lower <= other_peak.mz <= upper:
                    # Merge this peak
                    merged_intensity += other_peak.intensity
                    merged_indices.add(j)

                    # Keep charge if consistent, otherwise set to None
                    if merged_charge is not None and other_peak.charge is not None:
                        if merged_charge != other_peak.charge:
                            merged_charge = None

            # Create merged peak with combined intensity
            merged_peak = SpectrumPeak(base_peak.mz, merged_intensity, merged_charge)
            merged_peaks.append(merged_peak)

        # Sort merged peaks by m/z
        merged_peaks.sort(key=lambda p: p.mz)

        spectrum.peaks = merged_peaks
        return spectrum


@dataclass
class Ms1Spectrum(Spectrum):
    retention_time: float | None = None
    ms1_scan_number: int | None = None
    ion_mobility: float | None = None

    def __init__(
        self,
        peaks,
        retention_time: float | None = None,
        ms1_scan_number: int | None = None,
        ion_mobility: float | None = None,
    ):
        super().__init__(peaks)
        self.retention_time = retention_time
        self.ms1_scan_number = ms1_scan_number
        self.ion_mobility = ion_mobility

    def copy(self) -> Self:
        return self.__class__(
            peaks=[SpectrumPeak(p.mz, p.intensity, p.charge) for p in self.peaks],
            retention_time=self.retention_time,
            ms1_scan_number=self.ms1_scan_number,
            ion_mobility=self.ion_mobility,
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
        peaks,
        retention_time: float | None = None,
        ms1_scan_number: int | None = None,
        ion_mobility: float | None = None,
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
            ion_mobility=ion_mobility,
        )
        self.precursor_mz = precursor_mz
        self.precursor_charge = precursor_charge
        self.ms2_scan_number = ms2_scan_number
        self.activation_method = activation_method
        self.isolation_window = isolation_window

    def copy(self) -> Self:
        return self.__class__(
            peaks=[SpectrumPeak(p.mz, p.intensity, p.charge) for p in self.peaks],
            retention_time=self.retention_time,
            ms1_scan_number=self.ms1_scan_number,
            ion_mobility=self.ion_mobility,
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


def parse_pymzml_spectrum(spectrum) -> Spectrum:
    """Parse pymzML spectrum object into Spectrum."""

    mzs = []
    intensities = []
    for mz, intensity in spectrum.peaks("centroided"):
        mzs.append(float(mz))
        intensities.append(float(intensity))

    if spectrum.ms_level == 1:
        return Ms1Spectrum(
            retention_time=float(spectrum.scan_time_in_minutes() or 0),
            ms1_scan_number=spectrum.ID,
            ion_mobility=float(spectrum.get("ion mobility drift time (ms)", 0) or 0)
            if spectrum.get("ion mobility drift time (ms)")
            else None,
            peaks=[
                SpectrumPeak(mz, intensity) for mz, intensity in zip(mzs, intensities)
            ],
        )
    elif spectrum.ms_level == 2:
        return Ms2Spectrum(
            retention_time=float(spectrum.scan_time_in_minutes() or 0),
            ms1_scan_number=int(spectrum.get("ms1 scan number", 0) or 0)
            if spectrum.get("ms1 scan number")
            else None,
            ion_mobility=float(spectrum.get("ion mobility drift time (ms)", 0) or 0)
            if spectrum.get("ion mobility drift time (ms)")
            else None,
            precursor_mz=float(spectrum.selected_precursors[0]["mz"])
            if spectrum.selected_precursors and "mz" in spectrum.selected_precursors[0]
            else None,
            precursor_charge=int(spectrum.selected_precursors[0]["charge"])
            if spectrum.selected_precursors
            and "charge" in spectrum.selected_precursors[0]
            and isinstance(spectrum.selected_precursors[0]["charge"], int)
            else None,
            ms2_scan_number=spectrum.ID,
            activation_method=spectrum.get("activation method"),
            isolation_window=(
                float(spectrum.get("isolation window target m/z", 0) or 0)
                - float(spectrum.get("isolation window lower offset", 0) or 0),
                float(spectrum.get("isolation window target m/z", 0) or 0)
                + float(spectrum.get("isolation window upper offset", 0) or 0),
            )
            if all(
                k in spectrum
                for k in [
                    "isolation window target m/z",
                    "isolation window lower offset",
                    "isolation window upper offset",
                ]
            )
            else None,
            peaks=[
                SpectrumPeak(mz, intensity) for mz, intensity in zip(mzs, intensities)
            ],
        )
    else:
        return Spectrum((mzs, intensities))


def parse_mzml(
    mzml_file: str, ms_level: int | None = None
) -> Generator[Spectrum, None, None]:
    """Parse mzML file and yield Spectrum objects."""
    import pymzml

    for spectrum in pymzml.run.Reader(mzml_file):
        if spectrum.ms_level != ms_level and ms_level is not None:
            continue

        spec = parse_pymzml_spectrum(spectrum)
        yield spec
