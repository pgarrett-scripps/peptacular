import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum

from .annotation import Fragment


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

    # Convert to list for faster indexing if needed (ensure elements are float)
    mz2_list: list[float] = [float(x) for x in mz_spectrum2]

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
    fragments: Sequence[float | Fragment],
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

    fragments_mzs: list[float] = []
    for fragment in fragments:
        if isinstance(fragment, Fragment):
            fragments_mzs.append(fragment.mz)
        else:
            fragments_mzs.append(fragment)

    results: list[list[int]] = []
    for i, indexes in enumerate(
        get_matched_indices(fragments_mzs, mz_spectra, tolerance_value, tolerance_type)
    ):
        if indexes is None:
            results.append([])
        else:
            match mode:
                case MatchMode.ALL:
                    results.append(list(range(indexes[0], indexes[1])))

                case MatchMode.CLOSEST:
                    mz_diffs = [
                        abs(fragments_mzs[i] - mz_spectra[idx])
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


@dataclass(frozen=True)
class FragmentMatch:
    """
    Represents a match between a theoretical fragment and an experimental spectrum.
    """

    fragment: Fragment
    obs_mz: float
    intensity: float

    @property
    def error(self) -> float:
        """Error between theoretical and experimental m/z values (Da)."""
        return self.fragment.mz - self.obs_mz

    @property
    def error_ppm(self) -> float:
        """Error between theoretical and experimental m/z values (ppm)."""
        return (self.error / self.fragment.mz) * 1e6

    @property
    def mz(self) -> float:
        """Theoretical m/z from the fragment."""
        return self.fragment.mz

    def __str__(self) -> str:
        return f"{self.fragment} | Obs m/z: {self.obs_mz:.4f} | Δ: {self.error_ppm:.2f} ppm | I: {self.intensity:.2e}"


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

    mz_spectra, intensity_spectra = zip(
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
