import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Literal, Sequence
import multiprocessing as mp
from .spectrum import Ms1Spectrum, Ms2Spectrum, Spectrum


@dataclass
class Hill:
    mz_profile: list[float | None] = field(default_factory=lambda: [])
    intensity_profile: list[float] = field(default_factory=lambda: [])
    rt_profile: list[float] = field(default_factory=lambda: [])
    _scan_start: int = 0

    @property
    def mz(self) -> float:
        """Return the m/z value at the maximum intensity (apex)."""
        if not self.intensity_profile or not self.mz_profile:
            return 0.0
        apex_idx = self.intensity_profile.index(max(self.intensity_profile))
        mz_at_apex = self.mz_profile[apex_idx]
        # If apex m/z is None, find nearest non-None value
        if mz_at_apex is None:
            # Search outward from apex
            for offset in range(1, len(self.mz_profile)):
                if apex_idx - offset >= 0:
                    left_mz = self.mz_profile[apex_idx - offset]
                    if left_mz is not None:
                        return left_mz
                if apex_idx + offset < len(self.mz_profile):
                    right_mz = self.mz_profile[apex_idx + offset]
                    if right_mz is not None:
                        return right_mz
            return 0.0
        return mz_at_apex

    @property
    def scan_start(self) -> int:
        """Return the starting scan index."""
        return self._scan_start

    @property
    def scan_end(self) -> int:
        """Return the ending scan index."""
        return self._scan_start + len(self.intensity_profile) - 1

    @property
    def scan_apex(self) -> int:
        """Return the scan index of the maximum intensity."""
        apex_idx = self.intensity_profile.index(max(self.intensity_profile))
        return self._scan_start + apex_idx

    @property
    def rt_start(self) -> float:
        """Return the starting retention time."""
        return self.rt_profile[0] if self.rt_profile else 0.0

    @property
    def rt_end(self) -> float:
        """Return the ending retention time."""
        return self.rt_profile[-1] if self.rt_profile else 0.0

    @property
    def rt_apex(self) -> float:
        """Return the retention time at maximum intensity."""
        apex_idx = self.intensity_profile.index(max(self.intensity_profile))
        return self.rt_profile[apex_idx]

    @property
    def n_scans(self) -> int:
        """Return the number of scans in the hill."""
        return len(self.intensity_profile)

    @property
    def skipped_scans(self) -> int:
        """Return the number of skipped scans in the hill. (where intensity is zero)"""
        return sum(1 for intensity in self.intensity_profile if intensity == 0.0)

    @property
    def intensity_sum(self) -> float:
        """Return the sum of intensities in the hill."""
        return sum(self.intensity_profile)

    @property
    def intensity_max(self) -> float:
        """Return the maximum intensity in the hill."""
        return max(self.intensity_profile) if self.intensity_profile else 0.0

    @property
    def intensity(self) -> float:
        """Return the sum of intensities in the hill. (alias for intensity_sum)"""
        return self.intensity_max

    @property
    def scan_list(self) -> list[int]:
        """Return the list of scan indices for the hill."""
        return [self._scan_start + i for i in range(len(self.intensity_profile))]

    @property
    def intensity_list(self) -> list[float]:
        """Return the intensity profile list for the hill."""
        return self.intensity_profile

    @property
    def rt_list(self) -> list[float]:
        """Return the retention time profile list for the hill."""
        return self.rt_profile

    def draw_profile(
        self, height: int = 10, width: int | None = 30, log_scale: bool = False
    ) -> str:
        """
        Draw a text-based histogram of the intensity profile.

        Args:
            height: Number of rows in the output (default 10)
            width: Number of columns (default None = use actual profile length)
            log_scale: If True, use log scale for y-axis (default False)

        Returns:
            String representation of the histogram
        """
        if not self.intensity_profile:
            return ""

        profile = self.intensity_profile
        if width and width < len(profile):
            # Downsample if width specified
            step = len(profile) / width
            profile = [
                max(self.intensity_profile[int(i * step) : int((i + 1) * step)])
                for i in range(width)
            ]

        max_intensity = max(profile)
        if max_intensity == 0:
            return ""

        # Apply log scaling if requested
        if log_scale:
            # Add small epsilon to avoid log(0)
            epsilon = max_intensity * 1e-10
            profile = [math.log10(v + epsilon) for v in profile]
            max_val = math.log10(max_intensity + epsilon)
        else:
            max_val = max_intensity

        # Normalize to height
        normalized = [int(v / max_val * height) for v in profile]

        # Build histogram from top to bottom
        lines: list[str] = []
        for row in range(height, 0, -1):
            line = ""
            for val in normalized:
                line += "█" if val >= row else " "
            lines.append(line)

        return "\n".join(lines)


@dataclass
class _ActiveHill:
    """Internal class to track hills being built."""

    mz: float
    rt_start: float
    scan_start: int
    mz_profile: list[float | None] = field(default_factory=lambda: [])
    intensity_profile: list[float] = field(default_factory=lambda: [])
    rt_profile: list[float] = field(default_factory=lambda: [])
    last_scan_seen: int = 0
    mz_sum: float = 0.0
    mz_count: int = 0
    mean_mz: float = 0.0
    current_bin: int = 0


def _finalize_hill(hill: _ActiveHill, min_scans: int) -> Hill | None:
    """Convert an active hill to a finalized hill if it meets criteria."""
    if not hill.intensity_profile:
        return None

    # Trim zeros from start (shouldn't happen, but safety check)
    start_idx = 0
    while (
        start_idx < len(hill.intensity_profile)
        and hill.intensity_profile[start_idx] == 0.0
    ):
        start_idx += 1

    # Trim zeros from end
    end_idx = len(hill.intensity_profile) - 1
    while end_idx >= 0 and hill.intensity_profile[end_idx] == 0.0:
        end_idx -= 1

    # Check if we have any valid data left
    if start_idx > end_idx:
        return None

    # Trim the profiles
    mz_profile = hill.mz_profile[start_idx : end_idx + 1]
    intensity_profile = hill.intensity_profile[start_idx : end_idx + 1]
    rt_profile = hill.rt_profile[start_idx : end_idx + 1]

    profile_len = len(intensity_profile)
    if profile_len < min_scans:
        return None

    return Hill(
        mz_profile=mz_profile,
        intensity_profile=intensity_profile,
        rt_profile=rt_profile,
        _scan_start=hill.scan_start + start_idx,
    )


def _remove_stale_hills(
    active_hills: list[_ActiveHill],
    mz_bins: dict[int, list[int]],
    scan_idx: int,
    max_gap: int,
    min_scans: int,
) -> list[Hill]:
    """Remove hills that haven't been seen in max_gap scans and finalize them."""
    if not active_hills:
        return []

    if max_gap > 0 and scan_idx <= max_gap:
        return []

    hills_to_remove: list[int] = []
    finalized_hills: list[Hill] = []

    cutoff_scan = scan_idx - max_gap
    for hill_idx, hill in enumerate(active_hills):
        if hill.last_scan_seen < cutoff_scan:
            finalized = _finalize_hill(hill, min_scans)
            if finalized:
                finalized_hills.append(finalized)
            hills_to_remove.append(hill_idx)

    if not hills_to_remove:
        return finalized_hills

    # Remove hills from bins
    for hill_idx in hills_to_remove:
        hill = active_hills[hill_idx]
        if hill_idx in mz_bins[hill.current_bin]:
            mz_bins[hill.current_bin].remove(hill_idx)

    # Build new index mapping
    hills_to_remove_set = set(hills_to_remove)
    new_idx_map = {}
    new_idx = 0
    for old_idx in range(len(active_hills)):
        if old_idx not in hills_to_remove_set:
            new_idx_map[old_idx] = new_idx
            new_idx += 1

    # Remove hills
    for hill_idx in reversed(hills_to_remove):
        active_hills.pop(hill_idx)

    # Update bin indices
    for bin_idx in list(mz_bins.keys()):
        old_indices = mz_bins[bin_idx]
        new_indices: list[int] = [
            new_idx_map[idx] for idx in old_indices if idx in new_idx_map
        ]
        if new_indices:
            mz_bins[bin_idx] = new_indices
        else:
            del mz_bins[bin_idx]

    return finalized_hills


def _process_spectrum(
    spectrum: Spectrum,
    scan_idx: int,
    active_hills: list[_ActiveHill],
    mz_bins: dict[int, list[int]],
    tolerance: float,
    tolerance_type: str,
    bin_size_inv: float,
    min_mz: float,
    max_mz: float,
) -> None:
    """Process a single spectrum and update active hills."""
    # Get retention time
    if isinstance(spectrum, (Ms1Spectrum, Ms2Spectrum)):
        rt = (
            spectrum.retention_time
            if spectrum.retention_time is not None
            else float(scan_idx)
        )
    else:
        rt = float(scan_idx)


    # Sort peaks once
    sorted_peaks = sorted(spectrum.peaks, key=lambda p: p.mz)

    use_ppm = tolerance_type == "ppm"

    matched_hills: set[int] = set()
    tolerance_mult = tolerance / 1e6 if use_ppm else tolerance

    # Process each peak
    for peak in sorted_peaks:

        if peak.mz < min_mz:
            continue

        if peak.mz > max_mz:
            break

        peak_mz = peak.mz
        peak_intensity = peak.intensity

        # Calculate search radius
        if use_ppm:
            search_radius = peak_mz * tolerance_mult
        else:
            search_radius = tolerance_mult
        delta_for_bins = search_radius * 2

        # Get search bins inline
        min_bin = int((peak_mz - delta_for_bins) * bin_size_inv)
        max_bin = int((peak_mz + delta_for_bins) * bin_size_inv)

        # Find best matching hill
        best_hill_idx = None
        best_distance = float("inf")

        for bin_idx in range(min_bin, max_bin + 1):
            if bin_idx not in mz_bins:
                continue

            for hill_idx in mz_bins[bin_idx]:
                if hill_idx in matched_hills:
                    continue

                hill = active_hills[hill_idx]
                mean_mz = hill.mean_mz

                # Calculate tolerance
                if use_ppm:
                    delta = mean_mz * tolerance_mult
                else:
                    delta = tolerance_mult

                # Check distance
                mz_distance = abs(peak_mz - mean_mz)
                if mz_distance <= delta and mz_distance < best_distance:
                    best_distance = mz_distance
                    best_hill_idx = hill_idx

        # Update or create hill
        if best_hill_idx is not None:
            hill = active_hills[best_hill_idx]
            old_bin = hill.current_bin

            # Fill in any gaps with interpolated values
            gap = scan_idx - hill.last_scan_seen - 1
            if gap > 0:
                # Get last intensity and interpolate to current
                last_intensity = (
                    hill.intensity_profile[-1] if hill.intensity_profile else 0.0
                )
                last_rt = hill.rt_profile[-1] if hill.rt_profile else rt

                # Linear interpolation for both intensity and RT
                intensity_step = (peak_intensity - last_intensity) / (gap + 1)
                rt_step = (rt - last_rt) / (gap + 1)

                for i in range(gap):
                    interpolated_intensity = last_intensity + intensity_step * (i + 1)
                    interpolated_rt = last_rt + rt_step * (i + 1)
                    hill.mz_profile.append(None)  # No m/z for missing peaks
                    hill.intensity_profile.append(interpolated_intensity)
                    hill.rt_profile.append(interpolated_rt)

            # Add current peak
            hill.mz_profile.append(peak_mz)
            hill.intensity_profile.append(peak_intensity)
            hill.rt_profile.append(rt)
            hill.last_scan_seen = scan_idx

            # Update running average for m/z
            hill.mz_sum += peak_mz
            hill.mz_count += 1
            hill.mean_mz = hill.mz_sum / hill.mz_count

            matched_hills.add(best_hill_idx)

            # Update bin if needed
            new_bin = int(hill.mean_mz * bin_size_inv)
            if old_bin != new_bin:
                # Safety check before removal
                if old_bin in mz_bins and best_hill_idx in mz_bins[old_bin]:
                    mz_bins[old_bin].remove(best_hill_idx)
                if new_bin not in mz_bins:
                    mz_bins[new_bin] = []
                mz_bins[new_bin].append(best_hill_idx)
                hill.current_bin = new_bin
        else:
            new_bin = int(peak_mz * bin_size_inv)
            new_hill = _ActiveHill(
                mz=peak_mz,
                rt_start=rt,
                scan_start=scan_idx,
                mz_profile=[peak_mz],
                intensity_profile=[peak_intensity],
                rt_profile=[rt],
                last_scan_seen=scan_idx,
                mz_sum=peak_mz,
                mz_count=1,
                mean_mz=peak_mz,
                current_bin=new_bin,
            )
            hill_idx = len(active_hills)
            active_hills.append(new_hill)
            if new_bin not in mz_bins:
                mz_bins[new_bin] = []
            mz_bins[new_bin].append(hill_idx)


def detect_hills(
    spectra: Sequence[Spectrum],
    tolerance: float = 8,
    tolerance_type: Literal["ppm", "da"] = "ppm",
    min_scans: int = 3,
    max_gap: int = 1,
    split_hills: bool = True,
    min_mz: float = 0.0,
    max_mz: float = float("inf"),
) -> list[Hill]:
    """
    Detect hills in a sequence of spectra.

    Args:
        spectra: Sequence of Spectrum objects
        tolerance: Tolerance for m/z matching
        tolerance_type: 'ppm' or 'da'
        min_scans: Minimum scans for a valid hill
        max_gap: Maximum gap before finalizing a hill
        mz_window_size: Window size for calculating mean m/z

    Returns:
        List of detected hills
    """
    if not spectra:
        return []

    # Check for mixed spectrum types
    spectrum_types = set(type(s) for s in spectra)
    if len(spectrum_types) > 1:
        raise ValueError(
            f"Mixed spectrum types detected: {spectrum_types}. "
            "All spectra must be of the same type."
        )

    # Initialize
    mz_bins: defaultdict[int, list[int]] = defaultdict(list)
    active_hills: list[_ActiveHill] = []
    finalized_hills: list[Hill] = []

    # Determine check frequency
    if max_gap == 0:
        check_frequency = 1
    elif max_gap < 10:
        check_frequency = max_gap + 1
    else:
        check_frequency = max(10, max_gap // 2)

    bin_size = 0.1  # m/z bin size for hill tracking
    bin_size_inv = 1.0 / bin_size if bin_size > 0 else 1.0

    # Process spectra
    for scan_idx, spectrum in enumerate(spectra):
        _process_spectrum(
            spectrum,
            scan_idx,
            active_hills,
            mz_bins,
            tolerance,
            tolerance_type,
            bin_size_inv,
            min_mz,
            max_mz
        )

        # Periodic cleanup
        should_check = (scan_idx % check_frequency == 0) or (
            scan_idx == len(spectra) - 1
        )
        if should_check:
            newly_finalized = _remove_stale_hills(
                active_hills, mz_bins, scan_idx, max_gap, min_scans
            )
            finalized_hills.extend(newly_finalized)

    # Finalize remaining
    for hill in active_hills:
        finalized = _finalize_hill(hill, min_scans)
        if finalized:
            finalized_hills.append(finalized)

    # split hills using watershed
    if not split_hills:
        return finalized_hills

    return finalized_hills


def _worker_detect_hills(args: tuple) -> list[Hill]:
    """Worker function for multiprocessing hill detection."""
    (
        spectra,
        min_mz,
        max_mz,
        tolerance,
        tolerance_type,
        min_scans,
        max_gap,
        split_hills,
        chunk_min_mz,
        chunk_max_mz,
    ) = args
    
    hills = detect_hills(
        spectra=spectra,
        tolerance=tolerance,
        tolerance_type=tolerance_type,
        min_scans=min_scans,
        max_gap=max_gap,
        split_hills=split_hills,
        min_mz=min_mz,
        max_mz=max_mz,
    )
    
    # Remove hills within 2*tolerance of either boundary
    filtered_hills = []
    for hill in hills:
        hill_mz = hill.mz
        
        if tolerance_type == "ppm":
            lower_boundary = chunk_min_mz + 4 * tolerance * chunk_min_mz / 1e6
            upper_boundary = chunk_max_mz - 4 * tolerance * chunk_max_mz / 1e6
        else:  # da
            lower_boundary = chunk_min_mz + 4 * tolerance
            upper_boundary = chunk_max_mz - 4 * tolerance
        
        # Keep hill only if it's not near boundaries
        if lower_boundary <= hill_mz <= upper_boundary:
            filtered_hills.append(hill)
    
    return filtered_hills


def _deduplicate_hills(
    hills: list[Hill],
) -> list[Hill]:
    """
    Deduplicate hills that represent the same feature across multiple chunks.
    
    Since each peak can only be added to a single hill, if two hills have the
    same first (m/z, intensity) tuple, they are duplicates from overlapping chunks.
    
    Args:
        hills: List of hills (possibly with duplicates)
        tolerance: Unused (kept for API compatibility)
        tolerance_type: Unused (kept for API compatibility)
        
    Returns:
        Deduplicated list of hills
    """
    if not hills:
        return []
    
    # Use a dict to track first peak -> best hill
    # Key: (first_mz, first_intensity, first_scan)
    # Value: index of the hill with highest intensity_sum
    first_peak_to_hill: dict[tuple[float | None, float, int], int] = {}
    
    for hill_idx, hill in enumerate(hills):
        # Get first peak (m/z might be None for interpolated values, intensity always present)
        first_mz = hill.mz_profile[0]
        first_intensity = hill.intensity_profile[0]
        first_scan = hill.scan_start
        
        # Create unique key for this first peak
        key = (first_mz, first_intensity, first_scan)
        
        # If we've seen this peak before, keep the hill with higher total intensity
        if key in first_peak_to_hill:
            existing_idx = first_peak_to_hill[key]
            if hill.intensity_sum > hills[existing_idx].intensity_sum:
                first_peak_to_hill[key] = hill_idx
        else:
            first_peak_to_hill[key] = hill_idx
    
    # Return unique hills (in original order)
    unique_indices = set(first_peak_to_hill.values())
    return [hill for idx, hill in enumerate(hills) if idx in unique_indices]


def detect_hills_parallel(
    spectra: Sequence[Spectrum],
    tolerance: float = 8,
    tolerance_type: Literal["ppm", "da"] = "ppm",
    min_scans: int = 3,
    max_gap: int = 1,
    split_hills: bool = True,
    n_processes: int | None = None,
    global_min_mz: float = 0.0,
    global_max_mz: float = float("inf"),
    method: Literal["process", "thread"] = "process",
) -> list[Hill]:
    """
    Detect hills in parallel by splitting the m/z range into overlapping chunks.
    
    Args:
        spectra: Sequence of Spectrum objects
        tolerance: Tolerance for m/z matching
        tolerance_type: 'ppm' or 'da'
        min_scans: Minimum scans for a valid hill
        max_gap: Maximum gap before finalizing a hill
        split_hills: Whether to apply watershed splitting
        n_processes: Number of processes (default: CPU count)
        global_min_mz: Global minimum m/z to consider
        global_max_mz: Global maximum m/z to consider
        
    Returns:
        Deduplicated list of detected hills
    """
    if not spectra:
        return []
    
    if n_processes is None:
        n_processes = mp.cpu_count()
    
    # Find actual m/z range from data
    actual_min_mz = float("inf")
    actual_max_mz = 0.0
    
    for spectrum in spectra:
        if spectrum.peaks:
            spectrum_min = min(p.mz for p in spectrum.peaks)
            spectrum_max = max(p.mz for p in spectrum.peaks)
            actual_min_mz = min(actual_min_mz, spectrum_min)
            actual_max_mz = max(actual_max_mz, spectrum_max)
    
    # Apply global bounds
    actual_min_mz = max(actual_min_mz, global_min_mz)
    actual_max_mz = min(actual_max_mz, global_max_mz)
    actual_min_mz -= (tolerance * actual_min_mz / 1e6 if tolerance_type == "ppm" else tolerance)
    actual_max_mz += (tolerance * actual_max_mz / 1e6 if tolerance_type == "ppm" else tolerance)
    

    if actual_min_mz >= actual_max_mz:
        return []
    
    # Calculate overlap based on tolerance
    # 10x tolerance as specified
    if tolerance_type == "ppm":
        # Use average m/z for ppm calculation
        avg_mz = (actual_min_mz + actual_max_mz) / 2
        overlap = 20 * tolerance * avg_mz / 1e6
    else:  # da
        overlap = 20 * tolerance
    
    # Create m/z boundaries with overlap
    mz_range = actual_max_mz - actual_min_mz
    chunk_size = mz_range / n_processes
    
    # Build overlapping chunks
    chunks: list[tuple[float, float]] = []
    for i in range(n_processes):
        chunk_min = actual_min_mz + i * chunk_size - (overlap if i > 0 else 0)
        chunk_max = actual_min_mz + (i + 1) * chunk_size + (overlap if i < n_processes - 1 else 0)
        
        # Clamp to actual bounds
        chunk_min = max(chunk_min, actual_min_mz)
        chunk_max = min(chunk_max, actual_max_mz)
        
        chunks.append((chunk_min, chunk_max))
    
    # Prepare worker arguments (now includes chunk boundaries for filtering)
    worker_args = [
        (
            spectra,
            chunk_min,
            chunk_max,
            tolerance,
            tolerance_type,
            min_scans,
            max_gap,
            split_hills,
            chunk_min,  # Original chunk boundary for filtering
            chunk_max,  # Original chunk boundary for filtering
        )
        for chunk_min, chunk_max in chunks
    ]
    
    # Run in parallel
    if method == "process":
        with mp.Pool(processes=n_processes) as pool:
            results = pool.map(_worker_detect_hills, worker_args)
    elif method == "thread":
        from concurrent.futures import ThreadPoolExecutor
        
        with ThreadPoolExecutor(max_workers=n_processes) as executor:
            results = list(executor.map(_worker_detect_hills, worker_args))
    else:
        raise ValueError(f"Invalid method: {method}. Use 'process' or 'thread'.")
    
    # Combine results
    all_hills: list[Hill] = []
    for hills in results:
        all_hills.extend(hills)
    
    # Deduplicate 
    print(f"Combining results from {n_processes} processes, total hills before deduplication: {len(all_hills)}")
    deduplicated_hills = _deduplicate_hills(all_hills)
    
    return deduplicated_hills