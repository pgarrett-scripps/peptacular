from collections.abc import Iterable, Sequence
import math
from enum import StrEnum
from typing import Any, NamedTuple

from .fragment import Fragment
from .sequence.sequence import strip_mods


class ToleranceType(StrEnum):
    PPM = "ppm"
    TH = "th"


class MatchMode(StrEnum):
    CLOSEST = "closest"
    LARGEST = "largest"
    ALL = "all"


def get_matched_indices(
    mz_spectrum1: Sequence[float],
    mz_spectrum2: Sequence[float],
    tolerance_value: float = 0.1,
    tolerance_type: ToleranceType = ToleranceType.PPM,
) -> list[tuple[int, int] | None]:
    """
    Matches two m/z spectra based on a specified tolerance value and type.

    :param mz_spectrum1: list of m/z values for the first spectrum.
    :type mz_spectrum1: list[float]
    :param mz_spectrum2: list of m/z values for the second spectrum.
    :type mz_spectrum2: list[float]
    :param tolerance_value: Tolerance value for matching. Default is 0.1.
    :type tolerance_value: float
    :param tolerance_type: Type of tolerance ('ppm' or 'th'). Default is 'ppm'.
    :type tolerance_type: str

    :raises ValueError: If the provided tolerance type is not 'ppm' or 'th'.

    :return: list of index pairs representing matched peaks between the two spectra.
    :rtype: list[tuple[int, int]]

    .. code-block:: python

        >>> get_matched_indices([100, 200, 300], [100, 200, 300])
        [(0, 1), (1, 2), (2, 3)]

        >>> get_matched_indices([100.1, 250, 300, 400], [100, 200, 300], 1, 'th')
        [(0, 1), None, (2, 3), None]

        >>> get_matched_indices([100, 100, 100], [100, 100, 100])
        [(0, 3), (0, 3), (0, 3)]

    """

    if tolerance_value <= 0:
        raise ValueError("Tolerance value must be positive")

    if not mz_spectrum1 or not mz_spectrum2:
        return [None] * len(mz_spectrum1)

    # Convert to list for faster indexing if needed
    mz2_list = (
        list(mz_spectrum2) if not isinstance(mz_spectrum2, list) else mz_spectrum2
    )

    def _binary_search_left(arr: list[float], target: float) -> int:
        """Find leftmost position where target could be inserted to keep array sorted."""
        left, right = 0, len(arr)
        while left < right:
            mid = (left + right) // 2
            if arr[mid] < target:
                left = mid + 1
            else:
                right = mid
        return left

    def _binary_search_right(arr: list[float], target: float) -> int:
        """Find rightmost position where target could be inserted to keep array sorted."""
        left, right = 0, len(arr)
        while left < right:
            mid = (left + right) // 2
            if arr[mid] <= target:
                left = mid + 1
            else:
                right = mid
        return left

    def _calculate_tolerance_offset(mz_value: float) -> float:
        """Calculate tolerance offset based on type."""
        match tolerance_type:
            case ToleranceType.TH:
                return tolerance_value
            case ToleranceType.PPM:
                return mz_value * tolerance_value / 1e6

    indices: list[tuple[int, int] | None] = []

    for mz1 in mz_spectrum1:
        tolerance_offset = _calculate_tolerance_offset(mz1)
        mz1_start = mz1 - tolerance_offset
        mz1_end = mz1 + tolerance_offset

        # Use binary search for efficient range finding
        start_idx = _binary_search_left(mz2_list, mz1_start)
        end_idx = _binary_search_right(mz2_list, mz1_end)

        # Check if we found any matches
        if start_idx >= len(mz2_list) or start_idx >= end_idx:
            indices.append(None)
        else:
            indices.append((start_idx, end_idx))

    return indices


def match_spectra(
    fragments: Sequence[float],
    mz_spectra: Sequence[float],
    tolerance_value: float,
    tolerance_type: ToleranceType = ToleranceType.PPM,
    mode: MatchMode = MatchMode.CLOSEST,
    intensity_spectra: Sequence[float] | None = None,
) -> list[list[int]]:
    """
    Matches two m/z spectra based on a specified tolerance value and type.

    :param fragments: list of m/z values for the first spectrum.
    :type fragments: list[float]
    :param mz_spectra: list of m/z values for the second spectrum.
    :type mz_spectra: list[float]
    :param intensity_spectra: list of intensity values for the second spectrum. Required if mode is 'largest'.
    :type intensity_spectra: Union[list[float], None]
    :param tolerance_value: Tolerance value for matching. Default is 0.1.
    :type tolerance_value: float
    :param tolerance_type: Type of tolerance ('ppm' or 'th'). Default is 'pmp'.
    :type tolerance_type: str
    :param mode: Mode of matching ('closest' or 'largest' or 'all'). Default is 'closest'.
    :type mode: str

    :raises ValueError: If the provided tolerance type is not 'ppm' or 'th', or if the provided mode is
    not 'closest' or 'largest'.

    :return: list of lists containing matched indices. Empty list for no matches.
    :rtype: list[list[int]]

    .. code-block:: python

        >>> match_spectra([100, 200, 300], [100, 200, 300], 1, 'th')
        [[0], [1], [2]]

        >>> match_spectra([100.1, 250, 300, 400], [100, 200, 300], 1, 'th')
        [[0], [], [2], []]

        >>> f = [100.1, 250, 300, 400]
        >>> s = [100, 100.1, 200, 300]
        >>> match_spectra(f, s, 1, 'th')
        [[1], [], [3], []]

        >>> match_spectra(f, s, 1, 'th', mode='all')
        [[0, 1], [], [3], []]

        >>> match_spectra(f, s, 1, 'th', mode='largest', intensity_spectra=[10, 2, 3, 4])
        [[0], [], [3], []]

    """

    if tolerance_type not in [ToleranceType.PPM, ToleranceType.TH]:
        raise ValueError('Invalid tolerance type. Must be "ppm" or "th"')

    if mode not in [MatchMode.CLOSEST, MatchMode.LARGEST, MatchMode.ALL]:
        raise ValueError('Invalid mode. Must be "closest" or "largest" or "all"')

    results: list[list[int]] = []
    for i, indexes in enumerate(
        get_matched_indices(fragments, mz_spectra, tolerance_value, tolerance_type)
    ):
        if indexes is None:
            results.append([])
        else:
            match mode:
                case MatchMode.ALL:
                    results.append(list(range(indexes[0], indexes[1])))

                case MatchMode.CLOSEST:
                    mz_diffs = [
                        abs(fragments[i] - mz_spectra[idx])
                        for idx in range(indexes[0], indexes[1])
                    ]
                    min_index = mz_diffs.index(min(mz_diffs))
                    results.append([indexes[0] + min_index])

                case MatchMode.LARGEST:
                    if intensity_spectra is None:
                        raise ValueError(
                            'Intensity spectra must be provided when mode is "largest"'
                        )
                    intensities = intensity_spectra[indexes[0] : indexes[1]]
                    largest_index = intensities.index(max(intensities))
                    results.append([indexes[0] + largest_index])

    return results


class FragmentMatch(NamedTuple):
    """
    Represents a match between a theoretical fragment and an experimental spectrum.

    :ivar fragment: Theoretical fragment object.
    :ivar mz: Experimental m/z value.
    :ivar intensity: Intensity of the experimental m/z value.
    """

    fragment: Fragment | None
    mz: float
    intensity: float

    @property
    def error(self) -> float:
        """
        The error between the theoretical and experimental m/z values.

        :return: error between the theoretical and experimental m/z values.
        :rtype: float
        """
        return self.theo_mz - self.mz

    @property
    def error_ppm(self) -> float:
        """
        The error between the theoretical and experimental m/z values in parts-per-million (ppm).

        :return: error between the theoretical and experimental m/z values in parts-per-million (ppm).
        :rtype: float
        """
        if self.fragment is None:
            return 0.0

        return self.error / self.fragment.mz * 1e6

    @property
    def charge(self) -> int:
        if self.fragment is None:
            return 0

        return abs(self.fragment.charge)

    @property
    def ion_type(self) -> str:
        if self.fragment is None:
            return ""

        return self.fragment.ion_type

    @property
    def start(self) -> int:
        if self.fragment is None:
            return 0

        return self.fragment.start

    @property
    def end(self) -> int:
        if self.fragment is None:
            return 0

        return self.fragment.end

    @property
    def monoisotopic(self) -> bool:
        if self.fragment is None:
            return False

        return self.fragment.monoisotopic

    @property
    def isotope(self) -> int:
        if self.fragment is None:
            return 0
        return self.fragment.isotope

    @property
    def loss(self) -> float:
        if self.fragment is None:
            return 0.0
        return self.fragment.loss

    @property
    def sequence(self) -> str:
        if self.fragment is None:
            return ""
        return self.fragment.sequence

    @property
    def theo_mz(self) -> float:
        if self.fragment is None:
            return 0.0
        return self.fragment.mz

    @property
    def internal(self) -> bool:
        if self.fragment is None:
            return False
        return self.fragment.internal

    @property
    def label(self) -> str:
        if self.fragment is None:
            return ""
        return self.fragment.label

    @property
    def parent_sequence(self) -> str:
        if self.fragment is None:
            return ""
        return self.fragment.parent_sequence

    @property
    def number(self) -> str:
        if self.fragment is None:
            return "0"
        return self.fragment.number

    def to_dict(self) -> dict[str, Any]:
        """
        Converts the FragmentMatch object to a dictionary.

        :return: Dictionary representation of the FragmentMatch object.
        :rtype: Dict[str, Any]
        """

        return {
            "mz": self.mz,
            "intensity": self.intensity,
            "error": self.error,
            "error_ppm": self.error_ppm,
            "charge": self.charge,
            "ion_type": self.ion_type,
            "start": self.start,
            "end": self.end,
            "monoisotopic": self.monoisotopic,
            "isotope": self.isotope,
            "loss": self.loss,
            "sequence": self.sequence,
            "theo_mz": self.theo_mz,
            "internal": self.internal,
            "label": self.label,
            "number": self.number,
        }


def get_fragment_matches(
    fragments: Sequence[Fragment],
    mz_spectra: Sequence[float],
    intensity_spectra: Sequence[float],
    tolerance_value: float,
    tolerance_type: ToleranceType = ToleranceType.PPM,
    mode: MatchMode = MatchMode.ALL,
) -> list[FragmentMatch]:
    """
    Computes the fragment matches for a given set of fragments and an experimental spectrum.

    :param fragments: list of theoretical fragments.
    :type fragments: list[Fragment]
    :param mz_spectra: list of m/z values for the experimental spectrum.
    :type mz_spectra: list[float]
    :param intensity_spectra: list of intensity values for the experimental spectrum.
    :type intensity_spectra: list[float]
    :param tolerance_value: Tolerance value for matching. Default is 0.1.
    :type tolerance_value: float
    :param tolerance_type: Type of tolerance ('ppm' or 'th'). Default is 'ppm'.
    :type tolerance_type: str

    :param mode: Mode of matching ('all', 'closest' or 'largest'). Default is 'all'.
    :type mode: str

    :return: list of fragment matches.
    :rtype: list[FragmentMatch]
    """

    if mode not in [MatchMode.ALL, MatchMode.CLOSEST, MatchMode.LARGEST]:
        raise ValueError('Invalid mode. Must be "all", "closest" or "largest"')

    # sort fragments by mass
    fragments = sorted(fragments, key=lambda x: x.mz)
    fragment_spectrum = [f.mz for f in fragments]

    mz_spectra, intensity_spectra = zip(  # type: ignore
        *sorted(zip(mz_spectra, intensity_spectra), key=lambda x: x[0])
    )

    indices = match_spectra(
        fragment_spectrum,
        mz_spectra,
        tolerance_value,
        tolerance_type,
        mode,
        intensity_spectra,
    )

    fragment_matches: list[FragmentMatch] = []
    for i, index_list in enumerate(indices):
        if not index_list:  # Empty list means no matches
            continue

        for j in index_list:
            fragment_matches.append(
                FragmentMatch(fragments[i], mz_spectra[j], intensity_spectra[j])
            )

    return fragment_matches


def get_match_coverage(
    fragment_matches: Sequence[FragmentMatch],
) -> dict[str, list[int]]:
    """
    Returns a dictionary of fragment coverage for each fragment type / charge.

    :param fragment_matches: list of fragments.
    :type fragment_matches: list[Fragment]

    :return: Dictionary of fragment coverage.
    :rtype: dict[str, list[int]]
    """
    cov: dict[str, list[int]] = {}

    if not fragment_matches:
        return cov

    unmod_sequence = strip_mods(fragment_matches[0].parent_sequence)

    for frag in fragment_matches:
        label = f"{'+' * frag.charge}{frag.ion_type}"
        if label not in cov:
            cov[label] = [0] * len(unmod_sequence)
        for i in range(frag.start, frag.end):
            cov[label][i] += 1

    return cov


def _get_monoisotopic_label(label: str) -> str:
    return label.replace("*", "")


def filter_missing_mono_isotope(
    fragment_matches: list[FragmentMatch],
) -> list[FragmentMatch]:
    """
    Filters out fragment matches that do not have a monoisotopic peak.

    :param fragment_matches: list of fragment matches.
    :type fragment_matches: list[FragmentMatch]

    :return: list of fragment matches with monoisotopic peaks.
    :rtype: list[FragmentMatch]
    """

    mono_labels = set(f.label for f in fragment_matches if f.isotope == 0)
    return [
        f for f in fragment_matches if _get_monoisotopic_label(f.label) in mono_labels
    ]


def filter_skipped_isotopes(
    fragment_matches: list[FragmentMatch],
) -> list[FragmentMatch]:
    """
    Filters out fragment matches that have skipped isotopes - removes isotopes after a gap.

    :param fragment_matches: list of fragment matches.
    :type fragment_matches: list[FragmentMatch]

    :return: list of fragment matches without isotopes after gaps.
    :rtype: list[FragmentMatch]
    """

    # Group fragments by their base identity (same fragment, different isotopes)
    fragment_groups: dict[str, list[FragmentMatch]] = {}

    for match in fragment_matches:
        if match.fragment is None:
            continue

        # Create a key that identifies the same fragment regardless of isotope
        # Remove isotope markers from label to group isotopes of same fragment
        base_label = match.fragment.label.replace("*", "")  # Remove isotope markers
        key = f"{base_label}_{match.fragment.charge}_{match.fragment.start}_{match.fragment.end}_{match.fragment.loss}"

        if key not in fragment_groups:
            fragment_groups[key] = []
        fragment_groups[key].append(match)

    filtered_matches: list[FragmentMatch] = []

    for group in fragment_groups.values():
        if len(group) == 1:
            # Only one isotope, keep it
            filtered_matches.extend(group)
            continue

        # Sort by isotope number
        group.sort(key=lambda x: x.fragment.isotope if x.fragment else 0)

        # Find the first gap and remove everything after it
        keep_matches: list[FragmentMatch] = []
        expected_isotope = 0

        for match in group:
            if match.fragment is None:
                continue

            current_isotope = match.fragment.isotope

            if current_isotope == expected_isotope:
                # No gap, keep this match
                keep_matches.append(match)
                expected_isotope += 1
            else:
                # Gap found - don't keep this or any subsequent isotopes
                break

        filtered_matches.extend(keep_matches)

    return filtered_matches


def filter_losses(
    fragment_matches: list[FragmentMatch],
) -> list[FragmentMatch]:
    """
    Filters out fragment matches that do not have a base peak.

    :param fragment_matches: list of fragment matches.
    :type fragment_matches: list[FragmentMatch]

    :return: list of fragment matches with base peaks.
    :rtype: list[FragmentMatch]
    """

    new_matches: list[FragmentMatch] = []
    for fragment_match in fragment_matches:
        if fragment_match.loss != 0.0:
            for other in fragment_matches:
                if (
                    other.ion_type == fragment_match.ion_type
                    and other.start == fragment_match.start
                    and other.end == fragment_match.end
                    and other.charge == fragment_match.charge
                    and other.loss == 0.0
                ):
                    new_matches.append(fragment_match)
                    break
        else:
            new_matches.append(fragment_match)

    return new_matches


def _binomial_probability(n: int, k: int, p: float) -> float:
    """
    Computes the binomial probability P(X=k) for given parameters.

    :param n: Number of trials.
    :type n: int
    :param k: Number of successes.
    :type k: int
    :param p: Probability of success in a single trial.
    :type p: int

    :return: Binomial probability P(X=k).
    :rtype: float
    """

    return math.comb(n, k) * (p**k) * ((1 - p) ** (n - k))


def _estimate_probability_of_random_match(
    error_tolerance: float,
    mz_spectrum: Sequence[float],
    tolerance_type: ToleranceType,
    min_mz: float | None = None,
    max_mz: float | None = None,
) -> float:
    """
    Estimate the probability of a random match between two peaks based on error tolerance and the experimental spectrum.

    :param error_tolerance: Tolerance value for matching peaks.
    :type error_tolerance: float
    :param mz_spectrum: list of m/z values from the experimental spectrum.
    :type mz_spectrum: list[float]
    :param tolerance_type: Type of tolerance ('ppm' or 'th').
    :type tolerance_type: str

    :return: Estimated probability of a random match.
    :rtype: float
    """

    if min_mz is None:
        min_mz = min(mz_spectrum)

    if max_mz is None:
        max_mz = max(mz_spectrum)

    # Calculate the spectrum range
    spectrum_range = max_mz - min_mz

    # Calculate the error tolerance range
    if tolerance_type == ToleranceType.PPM:
        average_mz = sum(mz_spectrum) / len(mz_spectrum)
        error_tolerance_range = average_mz * error_tolerance / 1e6
    elif tolerance_type == ToleranceType.TH:
        error_tolerance_range = error_tolerance
    else:
        raise ValueError("Invalid tolerance type")

    # Estimate the probability
    num_bins = spectrum_range / error_tolerance_range

    return len(mz_spectrum) / num_bins


def binomial_score(
    fragments: Sequence[Fragment | float],
    mz_spectra: Sequence[float],
    tolerance_value: float,
    tolerance_type: ToleranceType,
    min_mz: float | None = None,
    max_mz: float | None = None,
) -> float:
    """
    Computes a score based on binomial probability for a given set of fragments and an experimental spectrum.

    :param min_mz:
    :param fragments: list of theoretical fragments.
    :type fragments: Union[list[Fragment], list[float]]
    :param mz_spectra: list of m/z values from the experimental spectrum.
    :type mz_spectra: list[float]
    :param tolerance_value: Tolerance value for matching fragments to the spectrum.
    :type tolerance_value: float
    :param tolerance_type: Type of tolerance ('ppm' or 'th').
    :type tolerance_type: str
    :param min_mz: Minimum m/z value for the spectrum.
    :type min_mz: Union[float, None]
    :param max_mz: Maximum m/z value for the spectrum.
    :type max_mz: Union[float, None]

    :return: Score based on binomial probability.
    :rtype: float

    .. code-block:: python

            >>> binomial_score([100, 200, 300], [100, 300, 400], 1, 'th')
            0.000297

            >>> binomial_score([103, 203, 303], [100, 300, 400], 1, 'th')
            0.970299

    """

    fragments_mzs: list[float] = []
    for fragment in fragments:
        if isinstance(fragment, Fragment):
            fragments_mzs.append(fragment.mz)
        else:
            fragments_mzs.append(fragment)

    indexes = get_matched_indices(
        fragments_mzs, mz_spectra, tolerance_value, tolerance_type
    )
    fragment_matches = [i is not None for i in indexes]
    matches = sum(fragment_matches)

    # Number of matched fragments (successes)
    k = matches

    # Total number of trials is the number of theoretical fragments
    n = len(fragments)

    # Estimate the probability of a random match
    p_success = _estimate_probability_of_random_match(
        tolerance_value, mz_spectra, tolerance_type, min_mz, max_mz
    )

    # Compute the score based on binomial probability
    score = _binomial_probability(n, k, p_success)

    return score


def bimodal_score_from_matches(
    fragment_matches: Sequence[FragmentMatch],
    fragments: Sequence[Fragment | float],
    mz_spectra: Sequence[float],
    tolerance_value: float,
    tolerance_type: ToleranceType,
    min_mz: float | None = None,
    max_mz: float | None = None,
) -> float:
    """
    Compute the binomial (bimodal) score using an existing set of FragmentMatch objects.

    Parameters:
    - fragment_matches: sequence of FragmentMatch (matches already found).
    - fragments: full list of theoretical fragments used as the trials (order/identity expected to be the same objects).
    - mz_spectra: experimental m/z values (used to estimate random-match probability).
    - tolerance_value, tolerance_type: used for estimating random-match probability.
    - min_mz, max_mz: optional m/z range for probability estimation.

    Returns:
    - binomial probability P(X = k) where k is the number of unique matched theoretical fragments
      and n is the total number of theoretical fragments.
    """

    n = len(fragments)
    if n == 0:
        return 0.0

    # Map theoretical fragments by identity for fast lookup
    frag_id_map = {id(f): idx for idx, f in enumerate(fragments)}

    matched_indices: set[int] = set()
    for fm in fragment_matches:
        if fm.fragment is None:
            continue
        idx = frag_id_map.get(id(fm.fragment))
        if idx is not None:
            matched_indices.add(idx)

    k = len(matched_indices)

    # Estimate probability of a random match given the experimental spectrum and tolerance
    p_success = _estimate_probability_of_random_match(
        tolerance_value, mz_spectra, tolerance_type, min_mz, max_mz
    )

    # Compute and return binomial probability P(X = k)
    return _binomial_probability(n, k, p_success)


def get_matched_intensity_percentage(
    fragment_matches: Sequence[FragmentMatch], intensities: Sequence[float]
) -> float:
    """
    Calculates the proportion of matched intensity to total intensity.

    :param fragment_matches: list of fragment matches.
    :type fragment_matches: list[FragmentMatch]
    :param intensities: list of intensities from the experimental spectrum.
    :type intensities: list[float]

    :return: Proportion of matched intensity to total intensity.
    :rtype: float
    """

    # group matches by mz
    matches: dict[float, FragmentMatch] = {f.mz: f for f in fragment_matches}
    matched_intensity = sum(f.intensity for f in matches.values())
    total_intensity = sum(intensities)

    if total_intensity == 0:
        return 0

    return matched_intensity / total_intensity


class Scorer:
    def __init__(
        self,
        experimental_spectra: tuple[Sequence[float], Sequence[float]],
        fragments: Sequence[Fragment],
        tolerance_type: ToleranceType,
        tolerance: float,
        match_mode: MatchMode,
        filter_fragments_with_iso_gap: bool = False,
        filter_fragments_without_mono: bool = False,
        remove_duplicate_matches: bool = True,
        filter_losses_without_base: bool = False,
    ) -> None:
        self.mz_spectra = experimental_spectra[0]
        self.intensities = experimental_spectra[1]
        self.fragments = fragments
        self.tolerance_type = tolerance_type
        self.tolerance = tolerance
        self.match_mode = match_mode
        self.filter_fragments_with_iso_gap = filter_fragments_with_iso_gap
        self.filter_fragments_without_mono = filter_fragments_without_mono
        self.remove_duplicate_matches = remove_duplicate_matches
        self.filter_losses_without_base = filter_losses_without_base

        self.fragment_matches = get_fragment_matches(
            self.fragments,
            self.mz_spectra,
            self.intensities,
            self.tolerance,
            self.tolerance_type,
            self.match_mode,
        )

        if self.remove_duplicate_matches:
            seen_mzs: set[float] = set()
            new_fragment_matches: list[FragmentMatch] = []
            # sort fragments by priority
            self.fragment_matches.sort(key=lambda x: x.fragment.priority if x.fragment else 0, reverse=True)
            for fm in self.fragment_matches:
                if fm.mz in seen_mzs:
                    continue
                seen_mzs.add(fm.mz)
                new_fragment_matches.append(fm)

            self.fragment_matches = new_fragment_matches

        if self.filter_fragments_with_iso_gap:
            self.fragment_matches = filter_skipped_isotopes(self.fragment_matches)

        if self.filter_fragments_without_mono:
            self.fragment_matches = filter_missing_mono_isotope(self.fragment_matches)

        if self.filter_losses_without_base:
            self.fragment_matches = filter_losses(self.fragment_matches)

    def _get_fragment_matches_by_direction(
        self, forward: bool, include_losses: bool, include_isotopes: bool
    ) -> list[FragmentMatch]:
        return self._get_fragment_matches_by_ion_type(
            ["a", "b", "c"] if forward else ["x", "y", "z"],
            include_losses,
            include_isotopes,
        )

    def _get_fragment_matches_by_ion_type(
        self,
        ion_types: str | Iterable[str] | None,
        include_losses: bool,
        include_isotopes: bool,
    ) -> list[FragmentMatch]:
        if ion_types is None:
            return [
                f
                for f in self.fragment_matches
                if (include_losses or f.loss == 0.0)
                and (include_isotopes or not f.isotope == 0)
            ]
        if isinstance(ion_types, str):
            ion_types = [ion_types]
        return [
            f
            for f in self.fragment_matches
            if f.ion_type in ion_types
            and (include_losses or f.loss == 0.0)
            and (include_isotopes or not f.isotope == 0)
        ]

    def _get_fragments_by_direction(
        self,
        forward: bool = True,
        include_losses: bool = True,
        include_isotopes: bool = True,
    ) -> Sequence[Fragment]:
        return self._get_fragments_by_ion_type(
            ["a", "b", "c"] if forward else ["x", "y", "z"],
            include_losses,
            include_isotopes,
        )

    def _get_fragments_by_ion_type(
        self,
        ion_types: str | Iterable[str] | None,
        include_losses: bool = True,
        include_isotopes: bool = True,
    ) -> Sequence[Fragment]:
        if ion_types is None:
            return [
                f
                for f in self.fragments
                if (include_losses or f.loss == 0.0)
                and (include_isotopes or not f.isotope == 0)
            ]
        if isinstance(ion_types, str):
            ion_types = [ion_types]
        return [
            f
            for f in self.fragments
            if f.ion_type in ion_types
            and (include_losses or f.loss == 0.0)
            and (include_isotopes or not f.isotope == 0)
        ]

    def get_match_coverage(
        self,
        ion_types: str | Iterable[str] | None = None,
        include_losses: bool = True,
        include_isotopes: bool = True,
    ) -> dict[str, list[int]]:
        return get_match_coverage(
            self._get_fragment_matches_by_ion_type(
                ion_types, include_losses, include_isotopes
            )
        )

    def matched_intensity(
        self,
        ion_types: str | Iterable[str] | None = None,
        include_losses: bool = True,
        include_isotopes: bool = True,
    ) -> float:
        return get_matched_intensity_percentage(
            self._get_fragment_matches_by_ion_type(
                ion_types, include_losses, include_isotopes
            ),
            self.intensities,
        )

    def bimodal_score(
        self,
        ion_types: str | Iterable[str] | None,
        include_losses: bool = True,
        include_isotopes: bool = True,
    ) -> float:
        return bimodal_score_from_matches(
            self._get_fragment_matches_by_ion_type(
                ion_types, include_losses, include_isotopes
            ),
            self._get_fragments_by_ion_type(
                ion_types, include_losses, include_isotopes
            ),
            self.mz_spectra,
            self.tolerance,
            self.tolerance_type,
        )
