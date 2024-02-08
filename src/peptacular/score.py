import math
from dataclasses import dataclass
from typing import List, Tuple, Union, Dict, Any
from peptacular.fragment import Fragment


def match_spectra_range(mz_spectrum1: List[float], mz_spectrum2: List[float], tolerance_value: float = 0.1,
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

        >>> match_spectra_range([100, 200, 300], [100, 200, 300])
        [(0, 1), (1, 2), (2, 3)]

        >>> match_spectra_range([100.1, 250, 300, 400], [100, 200, 300], 1, 'th')
        [(0, 1), None, (2, 3), None]

        >>> match_spectra_range([100, 100, 100], [100, 100, 100])
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


def match_spectra_best(mz_spectrum1: List[float], mz_spectrum2: List[float], tolerance_value: float = 0.1,
                       tolerance_type: str = 'ppm') -> List[Union[int, None]]:
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
    :rtype: List[int]
    
    .. code-block:: python

            >>> match_spectra_best([100, 200, 300], [100, 200, 300])
            [0, 1, 2]

            >>> match_spectra_best([100.1, 250, 300, 400], [100, 200, 300], 1, 'th')
            [0, None, 2, None]

            >>> match_spectra_best([100.1, 250, 300, 400], [100, 100.1, 200, 300], 1, 'th')
            [1, None, 3, None]

    """

    if tolerance_type not in ['ppm', 'th']:
        raise ValueError('Invalid tolerance type. Must be "ppm" or "th"')

    indices = [None] * len(mz_spectrum1)
    j = 0
    for i, mz1 in enumerate(mz_spectrum1):
        tolerance_offset = tolerance_value if tolerance_type == 'th' else mz1 * tolerance_value / 1e6
        mz1_min = mz1 - tolerance_offset
        mz1_max = mz1 + tolerance_offset

        best_match_index = None
        best_match_distance = float('inf')

        while j < len(mz_spectrum2) and mz_spectrum2[j] <= mz1_max:
            if mz1_min <= mz_spectrum2[j]:
                distance = abs(mz1 - mz_spectrum2[j])
                if distance < best_match_distance:
                    best_match_distance = distance
                    best_match_index = j
            j += 1

        indices[i] = best_match_index

    return indices


@dataclass(frozen=True)
class FragmentMatch:
    """
    Represents a match between a theoretical fragment and an experimental spectrum.

    :ivar fragment: Theoretical fragment object.
    :ivar mz: Experimental m/z value.
    :ivar intensity: Intensity of the experimental m/z value.
    """
    fragment: Fragment
    mz: float
    intensity: float

    @property
    def error(self):
        """
        The error between the theoretical and experimental m/z values.

        :return: error between the theoretical and experimental m/z values.
        :rtype: float
        """
        return self.fragment.mz - self.mz

    @property
    def error_ppm(self):
        """
        The error between the theoretical and experimental m/z values in parts-per-million (ppm).

        :return: error between the theoretical and experimental m/z values in parts-per-million (ppm).
        :rtype: float
        """
        return self.error / self.fragment.mz * 1e6

    def to_dict(self) -> Dict[str, Any]:
        """
        Converts the FragmentMatch object to a dictionary.

        :return: Dictionary representation of the FragmentMatch object.
        :rtype: Dict[str, Any]
        """
        return {
            'fragment': self.fragment.to_dict(),
            'mz': self.mz,
            'intensity': self.intensity,
            'error': self.error,
            'error_ppm': self.error_ppm
        }


def compute_fragment_matches(fragments: List[Fragment], mz_spectrum: List[float],
                             intensity_spectrum: List[float], tolerance_value: float = 0.1,
                             tolerance_type: str = 'ppm') -> List[FragmentMatch]:
    """
    Computes the fragment matches for a given set of fragments and an experimental spectrum.

    :param fragments: List of theoretical fragments.
    :type fragments: List[Fragment]
    :param mz_spectrum: List of m/z values for the experimental spectrum.
    :type mz_spectrum: List[float]
    :param intensity_spectrum: List of intensity values for the experimental spectrum.
    :type intensity_spectrum: List[float]
    :param tolerance_value: Tolerance value for matching. Default is 0.1.
    :type tolerance_value: float
    :param tolerance_type: Type of tolerance ('ppm' or 'th'). Default is 'ppm'.
    :type tolerance_type: str

    :return: List of fragment matches.
    :rtype: List[FragmentMatch]
    """

    # sort fragments by mass
    fragments.sort(key=lambda x: x.mz)
    fragment_spectrum = [f.mz for f in fragments]

    # sort peaks
    mz_spectrum, intensity_spectrum = zip(*sorted(zip(mz_spectrum, intensity_spectrum), key=lambda x: x[0]))

    indices = match_spectra_range(fragment_spectrum, mz_spectrum, tolerance_value, tolerance_type)

    fragment_matches = []
    for i, index in enumerate(indices):
        if index is None:
            continue
        for j in range(*index):
            fragment_matches.append(
                FragmentMatch(fragments[i], mz_spectrum[j], intensity_spectrum[j]))

    return fragment_matches


# TODO: Doesnt work with multiple charges or internal fragments
def hyper_score(fragments: Union[List[Fragment], Dict[str, List[float]]], mz_spectrum: List[float],
                intensity_spectrum: List[float], tolerance_value=0.1, tolerance_type='ppm',
                filter_by='intensity') -> float:
    """
    Computes the hyperscore for a given set of fragments and an experimental spectrum.

    :param fragments: List of theoretical fragments
    :type fragments: List[Fragment] or Dict[str, List[float]]
    :param mz_spectrum: List of m/z values from the experimental spectrum.
    :type mz_spectrum: List[float]
    :param intensity_spectrum: List of intensity values corresponding to the m/z values in mz_spectrum.
    :type intensity_spectrum: List[float]
    :param tolerance_value: Tolerance value for matching fragments to the spectrum.
    :type tolerance_value: float
    :param tolerance_type: Type of tolerance ('ppm' or 'th').
    :type tolerance_type: str
    :param filter_by: How to filter the matched fragments ('intensity' or 'error').
    :type filter_by: str

    :return: Computed hyper score.
    :rtype: float

    .. code-block:: python

            >>> hyper_score({'a': [100, 200, 300]}, [100, 300, 400], [1000, 2000, 1500])
            3.0

    """

    @dataclass(frozen=True)
    class MockFragment:
        mz: float
        ion_type: str

    if isinstance(fragments, Dict):
        fragments = [MockFragment(mz, ion_type) for ion_type in fragments for mz in fragments[ion_type]]

    max_intensity = max(intensity_spectrum)
    intensity_spectrum = [intensity / max_intensity for intensity in intensity_spectrum]
    # Compute the fragment matches using the provided function
    fragment_matches = compute_fragment_matches(fragments, mz_spectrum, intensity_spectrum, tolerance_value,
                                                tolerance_type)

    if len(fragment_matches) == 0:
        return 0

    if filter_by == 'intensity':
        fragment_matches.sort(key=lambda x: x.intensity, reverse=False)
    elif filter_by == 'error':
        fragment_matches.sort(key=lambda x: abs(x.error), reverse=True)

    fragment_matches = {match.fragment: match for match in fragment_matches}

    frag_nums = set()
    for frag in fragment_matches:
        key = (frag.ion_type, frag.parent_number)
        if key not in frag_nums:
            frag_nums.add(key)

    fragment_matches = fragment_matches.values()

    # Calculate the dot product using matched fragments
    dot_product = sum(match.intensity for match in fragment_matches)

    # Count the number of b and y ions
    Na = sum(1 for ion_type, _ in frag_nums if ion_type == 'a')
    Nb = sum(1 for ion_type, _ in frag_nums if ion_type == 'b')
    Nc = sum(1 for ion_type, _ in frag_nums if ion_type == 'c')
    Nx = sum(1 for ion_type, _ in frag_nums if ion_type == 'x')
    Ny = sum(1 for ion_type, _ in frag_nums if ion_type == 'y')
    Nz = sum(1 for ion_type, _ in frag_nums if ion_type == 'z')

    # Compute the hyper score
    score = dot_product * math.factorial(Na) * math.factorial(Nb) * math.factorial(Nc) * math.factorial(Nx) * \
            math.factorial(Ny) * math.factorial(Nz)

    return score


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
                                          tolerance_type: str = 'ppm') -> float:
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

    # Calculate the spectrum range
    spectrum_range = max(mz_spectrum) - min(mz_spectrum)

    # Calculate the error tolerance range
    if tolerance_type == 'ppm':
        average_mz = sum(mz_spectrum) / len(mz_spectrum)
        error_tolerance_range = average_mz * error_tolerance / 1e6
    else:  # 'th'
        error_tolerance_range = error_tolerance

    # Estimate the probability
    num_bins = spectrum_range / error_tolerance_range

    return len(mz_spectrum) / num_bins


def binomial_score(fragments: Union[List[Fragment], List[float]], mz_spectrum: List[float],
                   intensity_spectrum: List[float],
                   tolerance_value=0.1, tolerance_type='ppm') -> float:
    """
    Computes a score based on binomial probability for a given set of fragments and an experimental spectrum.

    :param fragments: List of theoretical fragments.
    :type fragments: List[Fragment]
    :param mz_spectrum: List of m/z values from the experimental spectrum.
    :type mz_spectrum: List[float]
    :param intensity_spectrum: List of intensity values corresponding to the m/z values in mz_spectrum.
    :type intensity_spectrum: List[float]
    :param tolerance_value: Tolerance value for matching fragments to the spectrum.
    :type tolerance_value: float
    :param tolerance_type: Type of tolerance ('ppm' or 'th').
    :type tolerance_type: str

    :return: Score based on binomial probability.
    :rtype: float
    """

    @dataclass(frozen=True)
    class MockFragment:
        mz: float

    if isinstance(fragments, Dict):
        fragments = [MockFragment(mz) for mz in fragments]

    # Compute the fragment matches using the provided function
    fragment_matches = compute_fragment_matches(fragments, mz_spectrum, intensity_spectrum, tolerance_value,
                                                tolerance_type)

    if len(fragment_matches) == 0:
        return 1.0

    fragment_matches = {match.fragment: match for match in fragment_matches}.values()

    # Number of matched fragments (successes)
    k = len(fragment_matches)

    # Total number of trials is the number of theoretical fragments
    n = len(fragments)

    # Estimate the probability of a random match
    p_success = _estimate_probability_of_random_match(tolerance_value, mz_spectrum, tolerance_type)

    # Compute the score based on binomial probability
    score = _binomial_probability(n, k, p_success)

    return score
