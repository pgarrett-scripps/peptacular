from collections.abc import Callable
from dataclasses import dataclass, field
import bisect

from .decon.score import IsotopicPatternScore, score_isotopic_pattern
from .hill import Hill

from .hill import Hill
from .decon.deconvolution import construct_graph, navigate_left, navigate_right
from .decon.dclass import Graph
from ..isotope import IsotopeLookup, estimate_isotopic_distribution
from ..constants import C13_NEUTRON_MASS, PROTON_MASS
from ..isotope import IsotopeLookup


@dataclass
class IsotopeFeature:
    """Represents a feature with isotopic peaks."""

    monoisotopic_hill: Hill
    charge: int
    isotope_hills: list[Hill] = field(default_factory=lambda: [])  # [M+1, M+2, ...]
    isotope_positions: list[int] = field(
        default_factory=lambda: []
    )  # [1, 2, ...] corresponding to each hill

    @property
    def hills(self) -> list[Hill]:
        """All hills including monoisotopic and isotopes."""
        return [self.monoisotopic_hill] + self.isotope_hills

    @property
    def observed_mz(self) -> float:
        """M/z of monoisotopic peak."""
        return self.monoisotopic_hill.mz

    @property
    def observed_neutral_mass(self) -> float | None:
        """Calculate neutral mass if charge is known."""
        if self.charge == 0:
            return None
        return self.observed_mz * abs(self.charge) - self.charge * PROTON_MASS

    @property
    def rt_apex(self) -> float:
        """RT at apex of monoisotopic peak."""
        return self.monoisotopic_hill.rt_apex

    @property
    def intensity(self) -> float:
        """Max intensity of monoisotopic peak."""
        return self.monoisotopic_hill.intensity_max

    @property
    def total_intensity(self) -> float:
        """Total intensity across all isotope peaks."""
        total = self.monoisotopic_hill.intensity_sum
        for hill in self.isotope_hills:
            total += hill.intensity_sum
        return total

    @property
    def n_isotopes(self) -> int:
        """Number of isotope peaks detected (including monoisotopic)."""
        return 1 + len(self.isotope_hills)

    def overlaps_with_hill(self, hill: Hill) -> bool:
        """Check if this feature uses the given hill."""
        if hill is self.monoisotopic_hill:
            return True
        return any(h is hill for h in self.isotope_hills)

    @property
    def hill_similarity(self) -> float | None:
        """Calculate similarity score based on RT apex alignment."""
        if not self.isotope_hills:
            return None

        rt_apexes = [self.monoisotopic_hill.rt_apex] + [
            h.rt_apex for h in self.isotope_hills
        ]

        # Calculate standard deviation
        mean_rt = sum(rt_apexes) / len(rt_apexes)
        variance = sum((rt - mean_rt) ** 2 for rt in rt_apexes) / len(rt_apexes)
        std_dev = variance**0.5

        # Handle perfect alignment (std_dev = 0)
        if std_dev == 0:
            return 1.0

        # Convert to score: lower std_dev = higher score
        scale = 0.05
        return min(1.0, max(0.0, (1.0 / (1.0 + std_dev / scale))))

    @property
    def isotope_profile(self) -> list[tuple[int, float]]:
        """
        Get the isotope profile as a list of (position, intensity) tuples.

        Returns:
            List of (isotope_position, apex_intensity) tuples, sorted by position.
        """
        all_hills = [self.monoisotopic_hill] + self.isotope_hills
        all_positions = [0] + self.isotope_positions
        # Use intensity_max for relative peak height comparison
        return [
            (pos, hill.intensity_max) for pos, hill in zip(all_positions, all_hills)
        ]

    @property
    def elution_profile(self) -> tuple[int, int, list[float]]:
        """
        Get the elution profile across all isotope hills.

        Returns:
            Tuple of (min_scan, max_scan, intensity_per_scan) where intensity_per_scan
            is a list of summed intensities for each scan in the range.
        """
        all_hills = [self.monoisotopic_hill] + self.isotope_hills
        min_scan = min(hill.scan_start for hill in all_hills)
        max_scan = max(hill.scan_end for hill in all_hills)
        scan_range = max_scan - min_scan + 1

        # Accumulate intensity for each scan
        scan_intensities = [0.0] * scan_range

        for hill in all_hills:
            for i, intensity in enumerate(hill.intensity_profile):
                scan_idx = hill.scan_start + i - min_scan
                if 0 <= scan_idx < scan_range:
                    scan_intensities[scan_idx] += intensity

        return (min_scan, max_scan, scan_intensities)

    def draw_isotope_profile(
        self, height: int = 10, width: int = 20, offset: int = 0
    ) -> str:
        """
        Draw a text-based bar chart of the summed intensity for each isotope peak.

        Args:
            height: Number of rows in the output (default 10)
            width: Width of each bar in characters (default 20)

        Returns:
            String representation of the isotope profile
        """
        profile = self.isotope_profile

        if not profile:
            return ""

        positions = [pos for pos, _ in profile]
        intensities = [intensity for _, intensity in profile]

        max_intensity = max(intensities)
        if max_intensity == 0:
            return ""

        # Build the chart
        lines: list[str] = []

        # Draw bars from top to bottom
        for row in range(height, 0, -1):
            line = ""
            for intensity in intensities:
                normalized = int(intensity / max_intensity * height)
                bar_width = width // len(intensities)

                if normalized >= row:
                    bar = "█" * bar_width
                else:
                    bar = " " * bar_width

                line += bar
            lines.append(line)

        # Add labels at bottom with offset correction
        label_line = ""
        for pos in positions:
            bar_width = width // len(intensities)
            # Adjust position based on offset
            corrected_pos = pos - offset

            if corrected_pos == 0:
                label = "M"
            elif corrected_pos > 0:
                label = f"M+{corrected_pos}"
            else:
                label = f"M{corrected_pos}"  # Negative, so it'll show as M-1, M-2, etc.

            padding = (bar_width - len(label)) // 2
            label_line += (
                " " * padding + label + " " * (bar_width - padding - len(label))
            )
        lines.append(label_line)

        return "\n".join(lines)

    def draw_elution_profile(
        self, height: int = 10, width: int | None = 40, show_individual: bool = False
    ) -> str:
        """
        Draw a text-based elution profile showing intensity over scans.

        Args:
            height: Number of rows in the output (default 10)
            width: Number of columns (default 40, None = use actual scan range)
            show_individual: If True, show individual hill profiles stacked (default False)

        Returns:
            String representation of the elution profile
        """
        min_scan, max_scan, scan_intensities = self.elution_profile

        if not scan_intensities:
            return ""

        # Downsample if width is specified and smaller than scan range
        if width and width < len(scan_intensities):
            step = len(scan_intensities) / width
            downsampled = []
            for i in range(width):
                start_idx = int(i * step)
                end_idx = int((i + 1) * step)
                downsampled.append(
                    max(scan_intensities[start_idx:end_idx])
                    if start_idx < end_idx
                    else 0.0
                )
            scan_intensities = downsampled

        max_intensity = max(scan_intensities) if scan_intensities else 0.0
        if max_intensity == 0:
            return ""

        # Normalize to height
        normalized = [int(v / max_intensity * height) for v in scan_intensities]

        # Build histogram from top to bottom
        lines: list[str] = []
        for row in range(height, 0, -1):
            line = ""
            for val in normalized:
                line += "█" if val >= row else " "
            lines.append(line)

        # Add scan range info at bottom (removed RT range since it's scan-based)
        scan_range = max_scan - min_scan + 1
        info_line = f"Scans: {min_scan}-{max_scan} ({scan_range} scans, RT: {self.rt_start:.2f}-{self.rt_end:.2f}min)"
        lines.append(info_line)

        return "\n".join(lines)

    def draw_hills(self, height: int = 10, width: int | None = None) -> str:
        """
        Draw individual hill profiles stacked vertically, aligned on the same scan axis.

        Args:
            height: Number of rows per hill (default 10)
            width: Number of columns (default None = use actual scan range)

        Returns:
            String representation of the individual hills aligned by scan
        """
        all_hills = [self.monoisotopic_hill] + self.isotope_hills
        all_labels = ["M"] + [f"M+{pos}" for pos in self.isotope_positions]

        # Get overall scan range
        min_scan, max_scan, _ = self.elution_profile
        scan_range = max_scan - min_scan + 1

        # Use width or actual scan range
        display_width = width if width is not None else scan_range

        lines: list[str] = []

        for label, hill in zip(all_labels, all_hills):
            # Create intensity array aligned to overall scan range
            aligned_intensities = [0.0] * scan_range

            for i, intensity in enumerate(hill.intensity_profile):
                scan_idx = hill.scan_start + i - min_scan
                if 0 <= scan_idx < scan_range:
                    aligned_intensities[scan_idx] = intensity

            # Downsample if width is specified and smaller than scan range
            if width and width < scan_range:
                step = scan_range / width
                downsampled: list[float] = []
                for i in range(width):
                    start_idx = int(i * step)
                    end_idx = int((i + 1) * step)
                    downsampled.append(
                        max(aligned_intensities[start_idx:end_idx])
                        if start_idx < end_idx
                        else 0.0
                    )
                aligned_intensities = downsampled

            # Normalize and draw
            max_intensity = max(aligned_intensities) if aligned_intensities else 0.0

            if max_intensity > 0:
                normalized = [
                    int(v / max_intensity * height) for v in aligned_intensities
                ]

                # Build histogram
                hill_lines: list[str] = []
                for row in range(height, 0, -1):
                    line = ""
                    for val in normalized:
                        line += "█" if val >= row else " "
                    hill_lines.append(line)

                # Add label and chart
                lines.append(
                    f"{label} (scans {hill.scan_start}-{hill.scan_end}. apex {hill.scan_apex}), mz: {hill.mz:.4f}:"
                )
                lines.extend(hill_lines)
                lines.append("")  # Empty line between hills
            else:
                lines.append(f"{label}: (no intensity)")
                lines.append("")

        # Add axis labels at bottom
        axis_label = f"Scans: {min_scan} {'─' * (display_width - 20)} {max_scan}"
        lines.append(axis_label)

        return "\n".join(lines)

    @property
    def rt_start(self) -> float:
        """Earliest RT across all isotope hills."""
        all_hills = [self.monoisotopic_hill] + self.isotope_hills
        return min(hill.rt_start for hill in all_hills)

    @property
    def rt_end(self) -> float:
        """Latest RT across all isotope hills."""
        all_hills = [self.monoisotopic_hill] + self.isotope_hills
        return max(hill.rt_end for hill in all_hills)


"""
@dataclass
class IsotopicPatternScore:
    combined_score: float
    offset: int
    bhattacharyya: float
    cosine_similarity: float
    ratio_score: float
    coverage: float
    missed_penalty: float
    normalized_theo: list[float]
"""


@dataclass
class ScoredIsotopeFeature:
    """An isotope feature with theoretical isotopic pattern scoring."""

    feature: IsotopeFeature
    score: IsotopicPatternScore

    @property
    def monoisotopic_mz(self) -> float:
        """Monoisotopic m/z corrected for neutron offset."""
        return self.feature.observed_mz - (
            self.score.offset * C13_NEUTRON_MASS / self.feature.charge
        )

    @property
    def neutral_mass(self) -> float | None:
        """Neutral mass from corrected monoisotopic m/z."""
        if self.feature.charge == 0:
            return None
        return (
            self.monoisotopic_mz * abs(self.feature.charge)
            - self.feature.charge * PROTON_MASS
        )

    def draw_isotope_profile(self, height: int = 10, width: int = 20) -> str:
        """
        Draw a text-based bar chart of the summed intensity for each isotope peak,
        adjusted for the scoring offset.

        Args:
            height: Number of rows in the output (default 10)
            width: Width of each bar in characters (default 20)
        Returns:
            String representation of the isotope profile
        """
        return self.feature.draw_isotope_profile(
            height=height, width=width, offset=self.score.offset
        )

    def draw_theoretical_profile(self, height: int = 10, width: int = 20) -> str:
        """
        Draw a text-based bar chart of the theoretical isotopic pattern.

        Args:
            height: Number of rows in the output (default 10)
            width: Width of each bar in characters (default 20)
        Returns:
            String representation of the theoretical isotope profile
        """
        intensities = self.score.normalized_theo
        if not intensities:
            return ""

        max_intensity = max(intensities)
        if max_intensity == 0:
            return ""

        # Build the chart
        lines: list[str] = []

        # prepend 0 for M-1 if offset > 0
        if self.score.offset > 0:
            intensities = [0.0] * self.score.offset + intensities
        # append 0 for M+1 if offset < 0
        elif self.score.offset < 0:
            intensities = intensities + [0.0] * abs(self.score.offset)

        # Draw bars from top to bottom
        for row in range(height, 0, -1):
            line = ""
            for intensity in intensities:
                normalized = int(intensity / max_intensity * height)
                bar_width = width // len(intensities)

                if normalized >= row:
                    bar = "█" * bar_width
                else:
                    bar = " " * bar_width

                line += bar
            lines.append(line)

        # Add labels at bottom with offset correction
        label_line = ""
        for pos in range(len(intensities)):
            bar_width = width // len(intensities)
            # Adjust position based on offset
            corrected_pos = pos - self.score.offset

            if corrected_pos == 0:
                label = "M"
            elif corrected_pos > 0:
                label = f"M+{corrected_pos}"
            else:
                label = f"M{corrected_pos}"  # Negative, so it'll show as M-1, M-2, etc.

            padding = (bar_width - len(label)) // 2
            label_line += (
                " " * padding + label + " " * (bar_width - padding - len(label))
            )
        lines.append(label_line)

        return "\n".join(lines)

    # Delegate other properties to feature
    def __getattr__(self, name: str):
        """Delegate attribute access to the wrapped feature."""
        try:
            return getattr(self.score, name)
        except AttributeError:
            return getattr(self.feature, name)


def score_isotope_feature(
    feature: IsotopeFeature,
    isotope_lookup: IsotopeLookup | None = None,
    min_intensity: float = 0.0,
    offset_range: int = 2,
) -> ScoredIsotopeFeature | None:
    """
    Score an isotope feature against theoretical distribution.

    Returns:
        ScoredIsotopeFeature if scoring succeeds, None otherwise
    """
    if feature.charge == 0 or feature.observed_neutral_mass is None:
        return None

    if len(feature.isotope_hills) < 1:  # Need at least M+1
        return None

    theoretical_dist = (
        isotope_lookup.get_isotope_pattern(feature.observed_neutral_mass)
        if isotope_lookup is not None
        else estimate_isotopic_distribution(
            feature.observed_neutral_mass,
            max_isotopes=20,
            output_masses_for_neutron_offset=False,
            min_abundance_threshold=0.005,
            use_neutron_count=True,
            is_abundance_sum=True,
        )
    )

    observed: list[float] = [hill.intensity_max for hill in feature.hills]
    theo: list[float] = [abundance for _, abundance in theoretical_dist]

    score = score_isotopic_pattern(
        observed, theo, min_intensity=min_intensity, offset_range=offset_range
    )

    return ScoredIsotopeFeature(feature=feature, score=score)


def create_hill_filter(
    max_apex_scan_diff: int,
) -> Callable[[Hill, Hill, int, float], bool]:
    """
    Create a hill filter function with specific parameters.

    Args:
        max_apex_scan_diff: Maximum allowed difference between apex scans (None = no limit)
        require_overlap: If True, require scan ranges to overlap

    Returns:
        Filter function compatible with construct_graph
    """

    def filter_func(h1: Hill, h2: Hill, charge: int, mz_diff: float) -> bool:
        """Filter function to require scan overlap and/or apex proximity between hills."""

        # Check apex scan proximity if specified
        apex_diff = abs(h1.scan_apex - h2.scan_apex)
        if apex_diff > max_apex_scan_diff:
            return False

        return True

    return filter_func


def detect_isotope_features(
    hills: list[Hill],
    neutron_mass: float = C13_NEUTRON_MASS,
    charge_carrier: float = PROTON_MASS,
    min_charge: int = 1,
    max_charge: int = 7,
    ppm_tolerance: float = 5.0,
    max_left_decrease: float = 0.75,
    max_right_decrease: float = 0.9,
    isotope_lookup: IsotopeLookup | None = None,
    offset_range: int = 1,
    debug: bool = False,
) -> list[ScoredIsotopeFeature]:
    """
    Detect isotope features using graph-based deconvolution approach.

    This uses the same algorithm as spectrum deconvolution but operates on Hill objects,
    requiring scan overlap for isotope connections.

    Args:
        hills: List of Hill objects
        neutron_mass: Mass of neutron for isotope spacing
        charge_carrier: Mass of charge carrier (proton)
        min_charge: Minimum charge state to consider
        max_charge: Maximum charge state to consider
        ppm_tolerance: m/z tolerance in ppm
        max_left_decrease: Max intensity drop going left in graph
        max_right_decrease: Max intensity drop going right in graph
        isotope_lookup: Optional isotope pattern lookup for scoring
        debug: Print debug information

    Returns:
        List of IsotopeFeature objects
    """
    if isotope_lookup is None:
        isotope_lookup = IsotopeLookup()

    if debug:
        print(f"Building graph from {len(hills)} hills...")

    hill_filter = create_hill_filter(
        max_apex_scan_diff=5,
    )

    # Construct graph with scan overlap requirement
    graph: Graph[Hill] = construct_graph(
        peaks=hills,
        tolerance=ppm_tolerance,
        tolerance_type="ppm",
        charge_range=(min_charge, max_charge),
        isotope_mass=neutron_mass,
        edge_filter=hill_filter,  # Require scan overlap
    )

    if debug:
        n_edges = len(graph.edge_index_map())
        print(f"Graph built with {len(hills)} nodes and {n_edges} edges")

    # Sort hills by intensity descending
    sorted_indices = sorted(
        range(len(hills)), key=lambda i: hills[i].intensity_max, reverse=True
    )

    features: list[IsotopeFeature] = []

    for idx, hill_idx in enumerate(sorted_indices):
        # Skip if already used
        if graph[hill_idx].seen:
            continue

        if debug and idx < 3:
            hill = hills[hill_idx]
            print(
                f"\nProcessing hill {idx}: m/z={hill.mz:.4f}, RT={hill.rt_apex:.2f}, scans={hill.scan_start}-{hill.scan_end}"
            )

        # Try each charge state
        results: dict[int, tuple[list[int], list[Hill]]] = {}

        for charge in range(max_charge, min_charge - 1, -1):
            # Navigate left and right to find isotope envelope
            left_indices = navigate_left(
                graph,
                hill_idx,
                charge,
                max_left_decrease,
            )
            right_indices = navigate_right(
                graph,
                hill_idx,
                charge,
                max_right_decrease,
            )

            # Combine and sort
            envelope_indices = sorted(set(left_indices + right_indices))
            envelope_hills = [graph[i].value for i in envelope_indices]

            results[charge] = (envelope_indices, envelope_hills)

            if debug and idx < 3:
                print(f"  Charge {charge}: Found {len(envelope_hills)} isotope peaks")

        # Select best charge state
        # Priority: most isotopes, then highest total intensity, then lowest charge
        def score_result(charge: int) -> tuple[int, float, int]:
            _, hills_list = results[charge]
            total_intensity = sum(h.intensity_sum for h in hills_list)
            return (
                len(hills_list),
                total_intensity,
                -charge,
            )  # Negative charge for ascending sort

        best_charge = max(results, key=score_result)
        best_indices, best_hills = results[best_charge]

        if len(best_hills) == 1:
            # Single peak - charge is uncertain, but we keep it
            single_peak = True
        else:
            single_peak = False

        # Mark as seen
        for i in best_indices:
            graph[i].seen = True

        # Create feature
        if len(best_hills) > 0:
            # Sort hills by m/z to ensure proper order
            sorted_hills = sorted(best_hills, key=lambda h: h.mz)

            feature = IsotopeFeature(
                monoisotopic_hill=sorted_hills[0],
                charge=best_charge if not single_peak else 0,  # 0 = unknown charge
                isotope_hills=sorted_hills[1:],
                isotope_positions=list(range(1, len(sorted_hills))),
            )
            features.append(feature)

    if debug:
        print(f"\nTotal features detected: {len(features)}")

    min_hill_intensity = min(h.intensity_max for h in hills) if hills else 0.0

    scored_features: list[ScoredIsotopeFeature] = []
    for feature in features:
        scored = score_isotope_feature(
            feature,
            isotope_lookup=isotope_lookup,
            min_intensity=min_hill_intensity / 2,
            offset_range=offset_range,
        )
        if scored is not None:
            scored_features.append(scored)

    return scored_features
