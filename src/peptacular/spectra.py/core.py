from dataclasses import dataclass, field
from typing import Literal, Self
from collections.abc import Callable


@dataclass
class Spectrum:
    """A simple data class to hold spectrum information."""

    mzs: list[float] = field(default_factory=list)
    intensities: list[float] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Validate that mzs and intensities have the same length."""
        if len(self.mzs) != len(self.intensities):
            raise ValueError(
                f"mzs and intensities must have the same length. "
                f"Got {len(self.mzs)} mzs and {len(self.intensities)} intensities."
            )

    def __len__(self) -> int:
        """Return the number of peaks in the spectrum."""
        return len(self.mzs)

    def __bool__(self) -> bool:
        """Return True if spectrum has peaks, False otherwise."""
        return len(self.mzs) > 0

    def copy(self) -> Self:
        """Return a deep copy of the spectrum."""
        return self.__class__(
            mzs=self.mzs.copy(),
            intensities=self.intensities.copy(),
        )

    def _filter_peaks(
        self, condition: Callable[[float, float], bool], inplace: bool = True
    ) -> Self:
        """
        Internal method to filter peaks based on a condition function.

        :param condition: Function that takes (mz, intensity) and returns True to keep the peak
        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The filtered spectrum (self if inplace=True, new spectrum if inplace=False)
        """
        filtered = [
            (mz, intensity)
            for mz, intensity in zip(self.mzs, self.intensities)
            if condition(mz, intensity)
        ]

        if inplace:
            if filtered:
                self.mzs, self.intensities = zip(*filtered)  # type: ignore
                self.mzs = list(self.mzs)
                self.intensities = list(self.intensities)
            else:
                self.mzs = []
                self.intensities = []
            return self
        else:
            if filtered:
                mzs, intensities = zip(*filtered)  # type: ignore
                return self.__class__(mzs=list(mzs), intensities=list(intensities))
            else:
                return self.__class__(mzs=[], intensities=[])

    def filter_by_intensity(
        self,
        min_intensity: float | None = None,
        max_intensity: float | None = None,
        inplace: bool = True,
    ) -> Self:
        """
        Filter peaks by intensity range.

        :param min_intensity: Minimum intensity (inclusive). None means no lower bound.
        :param max_intensity: Maximum intensity (inclusive). None means no upper bound.
        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The filtered spectrum

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 50, 100])
            >>> spec.filter_by_intensity(min_intensity=20)
            >>> spec.intensities
            [50, 100]

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 50, 100])
            >>> filtered = spec.filter_by_intensity(min_intensity=20, max_intensity=75, inplace=False)
            >>> filtered.intensities
            [50]
            >>> spec.intensities  # Original unchanged
            [10, 50, 100]
        """

        def condition(mz: float, intensity: float) -> bool:
            if min_intensity is not None and intensity < min_intensity:
                return False
            if max_intensity is not None and intensity > max_intensity:
                return False
            return True

        return self._filter_peaks(condition, inplace)

    def filter_by_mz(
        self,
        min_mz: float | None = None,
        max_mz: float | None = None,
        inplace: bool = True,
    ) -> Self:
        """
        Filter peaks by m/z range.

        :param min_mz: Minimum m/z (inclusive). None means no lower bound.
        :param max_mz: Maximum m/z (inclusive). None means no upper bound.
        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The filtered spectrum

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 50, 100])
            >>> spec.filter_by_mz(min_mz=150.0, max_mz=250.0)
            >>> spec.mzs
            [200.0]
        """

        def condition(mz: float, intensity: float) -> bool:
            if min_mz is not None and mz < min_mz:
                return False
            if max_mz is not None and mz > max_mz:
                return False
            return True

        return self._filter_peaks(condition, inplace)

    def remove_peaks_within_mz_range(
        self, lower: float, upper: float, inplace: bool = True
    ) -> Self:
        """
        Remove peaks within the specified m/z range [lower, upper].

        :param lower: Lower bound of m/z range (inclusive)
        :param upper: Upper bound of m/z range (inclusive)
        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The filtered spectrum

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 50, 100])
            >>> spec.remove_peaks_within_mz_range(150.0, 250.0)
            >>> spec.mzs
            [100.0, 300.0]
        """
        return self._filter_peaks(lambda mz, _: not (lower <= mz <= upper), inplace)

    def remove_peaks_within_tolerance(
        self,
        target_mz: float,
        tolerance: float,
        tolerance_unit: Literal["ppm", "da"] = "ppm",
        inplace: bool = True,
    ) -> Self:
        """
        Remove peaks within the specified m/z tolerance of target_mz.

        :param target_mz: Target m/z value
        :param tolerance: Tolerance value
        :param tolerance_unit: Unit of tolerance ('ppm' or 'da')
        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The filtered spectrum

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 50, 100])
            >>> spec.remove_peaks_within_tolerance(200.0, tolerance=50, tolerance_unit='da')
            >>> spec.mzs
            [100.0, 300.0]
        """
        lower, upper = self._calculate_mz_window(target_mz, tolerance, tolerance_unit)
        return self.remove_peaks_within_mz_range(lower, upper, inplace)

    def normalize_intensities(
        self, method: Literal["max", "sum", "tic"] = "max", inplace: bool = True
    ) -> Self:
        """
        Normalize intensities.

        :param method: Normalization method:
            - 'max': Scale to maximum intensity of 1.0
            - 'sum'/'tic': Scale so intensities sum to 1.0
        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The normalized spectrum

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 50, 100])
            >>> spec.normalize_intensities(method='max')
            >>> spec.intensities
            [0.1, 0.5, 1.0]

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 20, 70])
            >>> spec.normalize_intensities(method='sum')
            >>> spec.intensities
            [0.1, 0.2, 0.7]
        """
        if not self.intensities:
            return self if inplace else self.copy()

        if method == "max":
            max_intensity = max(self.intensities)
            if max_intensity > 0:
                normalized = [i / max_intensity for i in self.intensities]
            else:
                normalized = self.intensities.copy()
        elif method in ("sum", "tic"):
            total_intensity = sum(self.intensities)
            if total_intensity > 0:
                normalized = [i / total_intensity for i in self.intensities]
            else:
                normalized = self.intensities.copy()
        else:
            raise ValueError(f"Unknown normalization method: {method}")

        if inplace:
            self.intensities = normalized
            return self
        else:
            return self.__class__(mzs=self.mzs.copy(), intensities=normalized)

    def sort_by_mz(self, inplace: bool = True) -> Self:
        """
        Sort peaks by m/z.

        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The sorted spectrum

        .. code-block:: python

            >>> spec = Spectrum([300.0, 100.0, 200.0], [10, 50, 100])
            >>> spec.sort_by_mz()
            >>> spec.mzs
            [100.0, 200.0, 300.0]
        """
        sorted_peaks = sorted(zip(self.mzs, self.intensities))

        if inplace:
            if sorted_peaks:
                self.mzs, self.intensities = zip(*sorted_peaks)  # type: ignore
                self.mzs = list(self.mzs)
                self.intensities = list(self.intensities)
            return self
        else:
            if sorted_peaks:
                mzs, intensities = zip(*sorted_peaks)  # type: ignore
                return self.__class__(mzs=list(mzs), intensities=list(intensities))
            else:
                return self.__class__(mzs=[], intensities=[])

    def sort_by_intensity(self, reverse: bool = True, inplace: bool = True) -> Self:
        """
        Sort peaks by intensity.

        :param reverse: If True, sort in descending order (highest intensity first)
        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The sorted spectrum

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 100, 50])
            >>> spec.sort_by_intensity()
            >>> spec.intensities
            [100, 50, 10]
        """
        sorted_peaks = sorted(
            zip(self.mzs, self.intensities), key=lambda x: x[1], reverse=reverse
        )

        if inplace:
            if sorted_peaks:
                self.mzs, self.intensities = zip(*sorted_peaks)  # type: ignore
                self.mzs = list(self.mzs)
                self.intensities = list(self.intensities)
            return self
        else:
            if sorted_peaks:
                mzs, intensities = zip(*sorted_peaks)  # type: ignore
                return self.__class__(mzs=list(mzs), intensities=list(intensities))
            else:
                return self.__class__(mzs=[], intensities=[])

    def keep_top_n_peaks(self, n: int, inplace: bool = True) -> Self:
        """
        Keep only the top N most intense peaks.

        :param n: Number of peaks to keep
        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The filtered spectrum

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0, 400.0], [10, 100, 50, 75])
            >>> spec.keep_top_n_peaks(2)
            >>> spec.intensities
            [100, 75]
        """
        if n <= 0:
            if inplace:
                self.mzs = []
                self.intensities = []
                return self
            else:
                return self.__class__(mzs=[], intensities=[])

        if n >= len(self.mzs):
            return self if inplace else self.copy()

        top_peaks = sorted(
            zip(self.mzs, self.intensities), key=lambda x: x[1], reverse=True
        )[:n]

        if inplace:
            self.mzs, self.intensities = zip(*top_peaks)  # type: ignore
            self.mzs = list(self.mzs)
            self.intensities = list(self.intensities)
            return self
        else:
            mzs, intensities = zip(*top_peaks)  # type: ignore
            return self.__class__(mzs=list(mzs), intensities=list(intensities))

    def has_peak_within_tolerance(
        self,
        target_mz: float,
        tolerance: float,
        tolerance_unit: Literal["ppm", "da"] = "ppm",
    ) -> bool:
        """
        Check if there is any peak within the specified m/z tolerance of target_mz.

        :param target_mz: Target m/z value
        :param tolerance: Tolerance value
        :param tolerance_unit: Unit of tolerance ('ppm' or 'da')
        :return: True if a peak exists within tolerance, False otherwise

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 50, 100])
            >>> spec.has_peak_within_tolerance(200.5, tolerance=1.0, tolerance_unit='da')
            True
            >>> spec.has_peak_within_tolerance(250.0, tolerance=1.0, tolerance_unit='da')
            False
        """
        lower, upper = self._calculate_mz_window(target_mz, tolerance, tolerance_unit)
        return any(lower <= mz <= upper for mz in self.mzs)

    def find_peak_within_tolerance(
        self,
        target_mz: float,
        tolerance: float,
        tolerance_unit: Literal["ppm", "da"] = "ppm",
    ) -> tuple[float, float] | None:
        """
        Find the closest peak within tolerance to target_mz.

        :param target_mz: Target m/z value
        :param tolerance: Tolerance value
        :param tolerance_unit: Unit of tolerance ('ppm' or 'da')
        :return: (mz, intensity) tuple of closest peak, or None if no peak found

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 50, 100])
            >>> spec.find_peak_within_tolerance(200.5, tolerance=1.0, tolerance_unit='da')
            (200.0, 50)
        """
        lower, upper = self._calculate_mz_window(target_mz, tolerance, tolerance_unit)

        closest_peak = None
        closest_diff = float("inf")

        for mz, intensity in zip(self.mzs, self.intensities):
            if lower <= mz <= upper:
                diff = abs(mz - target_mz)
                if diff < closest_diff:
                    closest_diff = diff
                    closest_peak = (mz, intensity)

        return closest_peak

    @staticmethod
    def _calculate_mz_window(
        target_mz: float, tolerance: float, tolerance_unit: Literal["ppm", "da"]
    ) -> tuple[float, float]:
        """
        Calculate m/z window bounds based on tolerance.

        :param target_mz: Target m/z value
        :param tolerance: Tolerance value
        :param tolerance_unit: Unit of tolerance ('ppm' or 'da')
        :return: (lower_bound, upper_bound) tuple
        """
        if tolerance_unit == "ppm":
            delta = target_mz * tolerance / 1e6
        elif tolerance_unit == "da":
            delta = tolerance
        else:
            raise ValueError("tolerance_unit must be 'ppm' or 'da'")

        return target_mz - delta, target_mz + delta

    def sqrt_transform(self, inplace: bool = True) -> Self:
        """
        Apply square root transformation to intensities.
        Useful for reducing the influence of very intense peaks.

        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The transformed spectrum

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0], [4, 16, 100])
            >>> spec.sqrt_transform()
            >>> spec.intensities
            [2.0, 4.0, 10.0]
        """
        transformed = [i**0.5 for i in self.intensities]

        if inplace:
            self.intensities = transformed
            return self
        else:
            return self.__class__(mzs=self.mzs.copy(), intensities=transformed)

    def log_transform(self, base: float = 2.0, inplace: bool = True) -> Self:
        """
        Apply log transformation to intensities.
        Adds 1 to avoid log(0).

        :param base: Base of the logarithm
        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The transformed spectrum

        .. code-block:: python

            >>> import math
            >>> spec = Spectrum([100.0, 200.0], [1, 7])
            >>> spec.log_transform(base=2)
            >>> [round(i, 2) for i in spec.intensities]
            [1.0, 3.0]
        """
        import math

        transformed = [math.log(i + 1, base) for i in self.intensities]

        if inplace:
            self.intensities = transformed
            return self
        else:
            return self.__class__(mzs=self.mzs.copy(), intensities=transformed)

    def get_total_ion_current(self) -> float:
        """
        Calculate the total ion current (sum of all intensities).

        :return: Sum of all intensities

        .. code-block:: python

            >>> spec = Spectrum([100.0, 200.0, 300.0], [10, 50, 40])
            >>> spec.get_total_ion_current()
            100
        """
        return sum(self.intensities)


@dataclass
class Ms1Spectrum(Spectrum):
    """A data class to hold MS1 spectrum information."""

    retention_time: float | None = None
    ms1_scan_number: int | None = None
    ion_mobility: float | None = None

    def copy(self) -> Self:
        """Return a deep copy of the spectrum."""
        return self.__class__(
            mzs=self.mzs.copy(),
            intensities=self.intensities.copy(),
            retention_time=self.retention_time,
            ms1_scan_number=self.ms1_scan_number,
            ion_mobility=self.ion_mobility,
        )


@dataclass
class Ms2Spectrum(Ms1Spectrum):
    """A data class to hold MS2 spectrum information."""

    precursor_mz: float | None = None
    precursor_charge: int | None = None
    ms2_scan_number: int | None = None
    activation_method: str | None = None
    isolation_window: tuple[float, float] | None = None  # (lower, upper)

    def copy(self) -> Self:
        """Return a deep copy of the spectrum."""
        return self.__class__(
            mzs=self.mzs.copy(),
            intensities=self.intensities.copy(),
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
        """
        Remove precursor peak and optionally its isotope envelope.

        :param tolerance: Tolerance for peak removal
        :param tolerance_unit: Unit of tolerance ('ppm' or 'da')
        :param remove_isotopes: If True, also remove isotopic peaks
        :param max_isotopes: Maximum number of isotopes to check (forward and backward)
        :param inplace: If True, modify this spectrum. If False, return a new spectrum.
        :return: The filtered spectrum

        .. code-block:: python

            >>> spec = Ms2Spectrum([100.0, 500.0, 500.5, 600.0], [10, 50, 30, 100],
            ...                     precursor_mz=500.0, precursor_charge=2)
            >>> spec.remove_precursor_peak(tolerance=20, tolerance_unit='ppm')
            >>> spec.mzs
            [100.0, 600.0]
        """
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
