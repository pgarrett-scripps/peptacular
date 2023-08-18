import math
from dataclasses import dataclass
from typing import List
from .fragment import Fragment


# TODO: Validation


def match_spectra(mz_spectrum1, mz_spectrum2, tolerance_value=0.1, tolerance_type='ppm'):
    """
    Matches two m/z spectra based on a specified tolerance value and type.

    :param mz_spectrum1: List of m/z values for the first spectrum.
    :param mz_spectrum2: List of m/z values for the second spectrum.
    :param tolerance_value: Tolerance value for matching. Default is 0.1.
    :param tolerance_type: Type of tolerance ('ppm' or 'th'). Default is 'ppm'.

    :return: List of index pairs representing matched peaks between the two spectra.

    :raises ValueError: If the provided tolerance type is not 'ppm' or 'th'.
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


def compute_fragment_matches(fragments: List[Fragment], mz_spectrum, intensity_spectrum, tolerance_value=0.1,
                             tolerance_type='ppm'):
    """
    Computes the fragment matches for a given set of fragments and an experimental spectrum.
    :param fragments:  A list of Fragment objects.
    :param mz_spectrum:  A list of m/z values.
    :param intensity_spectrum:  A list of intensity values corresponding to the m/z values in mz_spectrum.
    :param tolerance_value:  The tolerance value for matching fragments to the spectrum.
    :param tolerance_type:  The type of tolerance ('ppm' or 'th').
    :return:  A list of FragmentMatch objects.
    """

    # sort fragments by mass
    fragments.sort(key=lambda x: x.mz)
    fragment_spectrum = [f.mz for f in fragments]
    indices = match_spectra(fragment_spectrum, mz_spectrum, tolerance_value, tolerance_type)

    fragment_matches = []
    for i, index in enumerate(indices):
        if index is None:
            continue
        for j in range(*index):
            fragment_matches.append(
                FragmentMatch(fragments[i], mz_spectrum[j], intensity_spectrum[j]))

    return fragment_matches


def hyper_score(fragments: List[Fragment], mz_spectrum: List[float], intensity_spectrum: List[float],
                tolerance_value=0.1, tolerance_type='ppm', filter_by='intensity') -> float:
    """
    Computes the hyperscore for a given set of fragments and an experimental spectrum.

    :param fragments: List of theoretical fragments.
    :param mz_spectrum: List of m/z values from the experimental spectrum.
    :param intensity_spectrum: List of intensity values corresponding to the m/z values in mz_spectrum.
    :param tolerance_value: Tolerance value for matching fragments to the spectrum.
    :param tolerance_type: Type of tolerance ('ppm' or 'th').
    :param filter_by: How to filter the matched fragments ('intensity' or 'error').

    :return: Computed hyper score.
    """

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

    fragment_matches = {match.fragment: match for match in fragment_matches}.values()

    # Calculate the dot product using matched fragments
    dot_product = sum(match.intensity for match in fragment_matches)

    # Count the number of b and y ions
    Na = sum(1 for match in fragment_matches if match.fragment.ion_type == 'a')
    Nb = sum(1 for match in fragment_matches if match.fragment.ion_type == 'b')
    Nc = sum(1 for match in fragment_matches if match.fragment.ion_type == 'c')
    Nx = sum(1 for match in fragment_matches if match.fragment.ion_type == 'x')
    Ny = sum(1 for match in fragment_matches if match.fragment.ion_type == 'y')
    Nz = sum(1 for match in fragment_matches if match.fragment.ion_type == 'z')

    # Compute the hyper score
    score = dot_product * math.factorial(Na) * math.factorial(Nb) * math.factorial(Nc) * \
            math.factorial(Nx) * math.factorial(Ny) * math.factorial(Nz)

    return score


def binomial_probability(n: int, k: int, p: float) -> float:
    """
    Computes the binomial probability P(X=k) for given parameters.

    :param n: Number of trials.
    :param k: Number of successes.
    :param p: Probability of success in a single trial.

    :return: Binomial probability P(X=k).
    """

    return math.comb(n, k) * (p ** k) * ((1 - p) ** (n - k))


def estimate_probability_of_random_match(error_tolerance: float, mz_spectrum: List[float],
                                         tolerance_type: str = 'ppm') -> float:
    """
    Estimate the probability of a random match between two peaks based on error tolerance and the experimental spectrum.

    :param error_tolerance: Tolerance value for matching peaks.
    :param mz_spectrum: List of m/z values from the experimental spectrum.
    :param tolerance_type: Type of tolerance ('ppm' or 'th').

    :return: Estimated probability of a random match.
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


def binomial_score(fragments: List[Fragment], mz_spectrum: List[float],
                   intensity_spectrum: List[float],
                   tolerance_value=0.1, tolerance_type='ppm') -> float:
    """
    Computes a score based on binomial probability for a given set of fragments and an experimental spectrum.

    :param fragments: List of theoretical fragments.
    :param mz_spectrum: List of m/z values from the experimental spectrum.
    :param intensity_spectrum: List of intensity values corresponding to the m/z values in mz_spectrum.
    :param tolerance_value: Tolerance value for matching fragments to the spectrum.
    :param tolerance_type: Type of tolerance ('ppm' or 'th').
    :param p_success: Probability of a successful match in a single trial (default is 0.5).

    :return: Score based on binomial probability.
    """

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
    p_success = estimate_probability_of_random_match(tolerance_value, mz_spectrum, tolerance_type)

    # Compute the score based on binomial probability
    score = binomial_probability(n, k, p_success)

    return score
