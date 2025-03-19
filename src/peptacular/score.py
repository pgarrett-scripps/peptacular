import math
from dataclasses import dataclass
from typing import List, Tuple, Union, Dict, Optional
from peptacular.fragmentation import Fragment
from peptacular.sequence.sequence_funcs import strip_mods


def get_matched_indices(mz_spectrum1: List[float], mz_spectrum2: List[float], tolerance_value: float = 0.1,
                        tolerance_type: str = 'ppm') -> List[Tuple[int, int]]:
    """
    Matches two m/z spectra based on a specified tolerance value and type.

    :param mz_spectrum1: List of m/z values for the first spectrum.
    :type mz_spectrum1: List[float]
    :param mz_spectrum2: List of m/z values for the second spectrum.
    :type mz_spectrum2: List[float]
    :param tolerance_value: Tolerance value for matching. Default is 0.1.
    :type tolerance_value: float
    :param tolerance_type: Type of tolerance ('ppm' or 'th'). Default is 'ppm'.
    :type tolerance_type: str

    :raises ValueError: If the provided tolerance type is not 'ppm' or 'th'.

    :return: List of index pairs representing matched peaks between the two spectra.
    :rtype: List[Tuple[int, int]]

    .. code-block:: python

        >>> get_matched_indices([100, 200, 300], [100, 200, 300])
        [(0, 1), (1, 2), (2, 3)]

        >>> get_matched_indices([100.1, 250, 300, 400], [100, 200, 300], 1, 'th')
        [(0, 1), None, (2, 3), None]

        >>> get_matched_indices([100, 100, 100], [100, 100, 100])
        [(0, 3), (0, 3), (0, 3)]

    """

    if tolerance_type not in ['ppm', 'th']:
        raise ValueError('Invalid tolerance type. Must be "ppm" or "th"')

    indices = []
    mz2_start_index = 0
    for mz1 in mz_spectrum1:

        if mz2_start_index >= len(mz_spectrum2):  # if we've exhausted mz_spectrum2
            indices.append(None)
            continue

        tolerance_offset = tolerance_value if tolerance_type == 'th' else mz1 * tolerance_value / 1e6
        mz1_start = mz1 - tolerance_offset
        mz1_end = mz1 + tolerance_offset

        # Find the first mz2 that is greater than or equal to mz1_start
        while mz2_start_index < len(mz_spectrum2) and mz_spectrum2[mz2_start_index] < mz1_start:
            mz2_start_index += 1

        # If we've reached the end of mz_spectrum2, break
        if mz2_start_index >= len(mz_spectrum2):
            indices.append(None)
            continue

        # Find the last mz2 that is less than or equal to mz1_end
        mz2_end_index = mz2_start_index
        while mz2_end_index < len(mz_spectrum2) and mz_spectrum2[mz2_end_index] <= mz1_end:
            mz2_end_index += 1

        # Decrement to get the actual last index (because while loop exits when condition is no longer met)
        if mz2_end_index - 1 < mz2_start_index:
            indices.append(None)
        else:
            indices.append((mz2_start_index, mz2_end_index))

    return indices


def match_spectra(fragments: List[float], mz_spectra: List[float], tolerance_value: float,
                  tolerance_type: str = 'ppm', mode: str = 'closest',
                  intensity_spectra: Optional[Union[List[float]]] = None) -> List[Union[int, None]]:
    """
    Matches two m/z spectra based on a specified tolerance value and type.

    :param fragments: List of m/z values for the first spectrum.
    :type fragments: List[float]
    :param mz_spectra: List of m/z values for the second spectrum.
    :type mz_spectra: List[float]
    :param intensity_spectra: List of intensity values for the second spectrum. Required if mode is 'largest'.
    :type intensity_spectra: Union[List[float], None]
    :param tolerance_value: Tolerance value for matching. Default is 0.1.
    :type tolerance_value: float
    :param tolerance_type: Type of tolerance ('ppm' or 'th'). Default is 'ppm'.
    :type tolerance_type: str
    :param mode: Mode of matching ('closest' or 'largest' or 'all'). Default is 'closest'.
    :type mode: str

    :raises ValueError: If the provided tolerance type is not 'ppm' or 'th', or if the provided mode is
    not 'closest' or 'largest'.
    
    :return: List of index pairs representing matched peaks between the two spectra.
    :rtype: List[int]
    
    .. code-block:: python

            >>> match_spectra([100, 200, 300], [100, 200, 300], 1, 'th')
            [0, 1, 2]

            >>> match_spectra([100.1, 250, 300, 400], [100, 200, 300], 1, 'th')
            [0, None, 2, None]

            >>> f = [100.1, 250, 300, 400]
            >>> s = [100, 100.1, 200, 300]
            >>> match_spectra(f, s, 1, 'th')
            [1, None, 3, None]

            >>> match_spectra(f, s, 1, 'th', mode='all')
            [[0, 1], None, [3], None]

            >>> match_spectra(f, s, 1, 'th', mode='largest', intensity_spectra=[10, 2, 3, 4])
            [0, None, 3, None]

    """

    if tolerance_type not in ['ppm', 'th']:
        raise ValueError('Invalid tolerance type. Must be "ppm" or "th"')

    if mode not in ['closest', 'largest', 'all']:
        raise ValueError('Invalid mode. Must be "closest" or "largest" or "all"')

    results = []
    for i, indexes in enumerate(get_matched_indices(fragments, mz_spectra, tolerance_value, tolerance_type)):
        if indexes is None:
            results.append(None)
        else:

            if mode == 'all':
                results.append(list(range(indexes[0], indexes[1])))

            if mode == 'closest':
                mz_diffs = [abs(fragments[i] - mz_spectra[idx]) for idx in range(indexes[0], indexes[1])]
                min_index = mz_diffs.index(min(mz_diffs))

                results.append(indexes[0] + min_index)

            if mode == 'largest':
                intensities = intensity_spectra[indexes[0]:indexes[1]]
                largest_index = intensities.index(max(intensities))

                results.append(indexes[0] + largest_index)

    return results


@dataclass(frozen=True)
class FragmentMatch:
    """
    Represents a match between a theoretical fragment and an experimental spectrum.

    :ivar fragment: Theoretical fragment object.
    :ivar mz: Experimental m/z value.
    :ivar intensity: Intensity of the experimental m/z value.
    """
    fragment: Optional[Fragment]
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
        return self.error / self.fragment.mz * 1e6

    @property
    def charge(self) -> int:
        return abs(self.fragment.charge)

    @property
    def ion_type(self) -> str:
        return self.fragment.ion_type

    @property
    def start(self) -> int:
        return self.fragment.start

    @property
    def end(self) -> int:
        return self.fragment.end

    @property
    def monoisotopic(self) -> bool:
        return self.fragment.monoisotopic

    @property
    def isotope(self) -> int:
        return self.fragment.isotope

    @property
    def loss(self) -> float:
        return self.fragment.loss

    @property
    def sequence(self) -> str:
        return self.fragment.sequence

    @property
    def theo_mz(self) -> float:
        return self.fragment.mz

    @property
    def internal(self) -> bool:
        return self.fragment.internal

    @property
    def label(self) -> str:
        return self.fragment.label

    @property
    def parent_sequence(self) -> str:
        return self.fragment.parent_sequence

    @property
    def number(self) -> int:
        return self.fragment.number


    def to_dict(self) -> Dict:
        """
        Converts the FragmentMatch object to a dictionary.

        :return: Dictionary representation of the FragmentMatch object.
        :rtype: Dict[str, Any]
        """

        if self.fragment is not None:

            return {
                'mz': self.mz,
                'intensity': self.intensity,
                'error': self.error,
                'error_ppm': self.error_ppm,
                'charge': self.charge,
                'ion_type': self.ion_type,
                'start': self.start,
                'end': self.end,
                'monoisotopic': self.monoisotopic,
                'isotope': self.isotope,
                'loss': self.loss,
                'sequence': self.sequence,
                'theo_mz': self.theo_mz,
                'internal': self.internal,
                'label': self.fragment.label,
                'number': self.number
            }


        return {
            'mz': self.mz,
            'intensity': self.intensity,
            'error': 0,
            'error_ppm': 0,
            'charge': 0,
            'ion_type': '',
            'start': 0,
            'end': 0,
            'monoisotopic': True,
            'isotope': 0,
            'loss': 0,
            'sequence': '',
            'theo_mz': 0,
            'internal': False,
            'label': '',
            'number': 0
        }


def get_fragment_matches(fragments: List[Fragment], mz_spectra: List[float],
                         intensity_spectra: List[float], tolerance_value: float,
                         tolerance_type: str = 'ppm', mode: str = 'all') -> List[FragmentMatch]:
    """
    Computes the fragment matches for a given set of fragments and an experimental spectrum.

    :param fragments: List of theoretical fragments.
    :type fragments: List[Fragment]
    :param mz_spectra: List of m/z values for the experimental spectrum.
    :type mz_spectra: List[float]
    :param intensity_spectra: List of intensity values for the experimental spectrum.
    :type intensity_spectra: List[float]
    :param tolerance_value: Tolerance value for matching. Default is 0.1.
    :type tolerance_value: float
    :param tolerance_type: Type of tolerance ('ppm' or 'th'). Default is 'ppm'.
    :type tolerance_type: str

    :param mode: Mode of matching ('all', 'closest' or 'largest'). Default is 'all'.
    :type mode: str

    :return: List of fragment matches.
    :rtype: List[FragmentMatch]
    """

    if mode not in ['all', 'closest', 'largest']:
        raise ValueError('Invalid mode. Must be "all", "closest" or "largest"')

    # sort fragments by mass
    fragments.sort(key=lambda x: x.mz)
    fragment_spectrum = [f.mz for f in fragments]

    # sort peaks
    mz_spectra, intensity_spectra = zip(*sorted(zip(mz_spectra, intensity_spectra), key=lambda x: x[0]))

    indices = match_spectra(fragment_spectrum, mz_spectra, tolerance_value, tolerance_type, mode, intensity_spectra)

    fragment_matches = []
    for i, index in enumerate(indices):
        if index is None:
            continue
        if isinstance(index, int):
            fragment_matches.append(
                FragmentMatch(fragments[i], mz_spectra[index], intensity_spectra[index]))
        elif isinstance(index, List):
            for j in index:
                fragment_matches.append(
                    FragmentMatch(fragments[i], mz_spectra[j], intensity_spectra[j]))

    return fragment_matches


def get_match_coverage(fragment_matches: List[FragmentMatch]) -> Dict[str, List[int]]:
    """
    Returns a dictionary of fragment coverage for each fragment type / charge.

    :param fragment_matches: List of fragments.
    :type fragment_matches: List[Fragment]

    :return: Dictionary of fragment coverage.
    :rtype: Dict[str, List[int]]
    """
    cov = {}

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
    return label.replace('*', '')


def _remove_isotope_from_label(label: str) -> str:
    if label.endswith('*'):
        return label[:-1]
    return ''


def _add_isotope_to_label(label: str) -> str:
    return label + '*'


def filter_missing_mono_isotope(fragment_matches: List[FragmentMatch]) -> List[FragmentMatch]:
    """
    Filters out fragment matches that do not have a monoisotopic peak.

    :param fragment_matches: List of fragment matches.
    :type fragment_matches: List[FragmentMatch]

    :return: List of fragment matches with monoisotopic peaks.
    :rtype: List[FragmentMatch]
    """

    mono_labels = set(f.label for f in fragment_matches if f.isotope == 0)
    return [f for f in fragment_matches if _get_monoisotopic_label(f.label) in mono_labels]


def filter_skipped_isotopes(fragment_matches: List[FragmentMatch]) -> List[FragmentMatch]:
    """
    Filters out fragment matches that have skipped isotopes.

    :param fragment_matches: List of fragment matches.
    :type fragment_matches: List[FragmentMatch]

    :return: List of fragment matches with monoisotopic peaks.
    :rtype: List[FragmentMatch]
    """

    # remove peaks any peaks after a missing isotope
    labels = set([f.label for f in fragment_matches])
    return [f for f in fragment_matches if
            _remove_isotope_from_label(f.label) in labels or _add_isotope_to_label(f.label) in labels]


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

    return math.comb(n, k) * (p ** k) * ((1 - p) ** (n - k))


def _estimate_probability_of_random_match(error_tolerance: float, mz_spectrum: List[float],
                                          tolerance_type: str = 'ppm', min_mz: float = None,
                                          max_mz: float = None) -> float:
    """
    Estimate the probability of a random match between two peaks based on error tolerance and the experimental spectrum.

    :param error_tolerance: Tolerance value for matching peaks.
    :type error_tolerance: float
    :param mz_spectrum: List of m/z values from the experimental spectrum.
    :type mz_spectrum: List[float]
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
    if tolerance_type == 'ppm':
        average_mz = sum(mz_spectrum) / len(mz_spectrum)
        error_tolerance_range = average_mz * error_tolerance / 1e6
    else:  # 'th'
        error_tolerance_range = error_tolerance

    # Estimate the probability
    num_bins = spectrum_range / error_tolerance_range

    return len(mz_spectrum) / num_bins


def binomial_score(fragments: Union[List[Fragment], List[float]], mz_spectra: List[float],
                   tolerance_value: float, tolerance_type='ppm', min_mz: float = None,
                   max_mz: float = None) -> float:
    """
    Computes a score based on binomial probability for a given set of fragments and an experimental spectrum.

    :param min_mz:
    :param fragments: List of theoretical fragments.
    :type fragments: Union[List[Fragment], List[float]]
    :param mz_spectra: List of m/z values from the experimental spectrum.
    :type mz_spectra: List[float]
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

    if fragments and isinstance(fragments[0], Fragment):
        fragments = [fragment.mz for fragment in fragments]

    indexes = get_matched_indices(fragments, mz_spectra, tolerance_value, tolerance_type)
    fragment_matches = [i is not None for i in indexes]
    matches = sum(fragment_matches)

    # Number of matched fragments (successes)
    k = matches

    # Total number of trials is the number of theoretical fragments
    n = len(fragments)

    # Estimate the probability of a random match
    p_success = _estimate_probability_of_random_match(tolerance_value, mz_spectra, tolerance_type, min_mz, max_mz)

    # Compute the score based on binomial probability
    score = _binomial_probability(n, k, p_success)

    return score


def get_matched_intensity_percentage(fragment_matches: List[FragmentMatch], intensities: List[float]) -> float:
    """
    Calculates the proportion of matched intensity to total intensity.

    :param fragment_matches: List of fragment matches.
    :type fragment_matches: List[FragmentMatch]
    :param intensities: List of intensities from the experimental spectrum.
    :type intensities: List[float]

    :return: Proportion of matched intensity to total intensity.
    :rtype: float
    """

    # group matches by mz
    matches = {f.mz: f for f in fragment_matches}
    matched_intensity = sum(f.intenisty for f in matches.values())
    total_intensity = sum(intensities)

    if total_intensity == 0:
        return 0

    return matched_intensity / total_intensity
