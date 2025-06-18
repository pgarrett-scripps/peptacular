"""
Isotope.py - A module for calculating isotopic distributions of molecules.
"""

import sys
from typing import List, Dict, Tuple, Union, Optional
import warnings

from . import constants
from .types import ChemComposition
from .chem.chem_calc import estimate_comp
from .mass_calc import chem_mass


def isotopic_distribution(
    chemical_formula: ChemComposition,
    max_isotopes: Optional[int] = None,
    min_abundance_threshold: Optional[float] = None,
    distribution_resolution: Union[int] = 5,
    use_neutron_count: bool = False,
    conv_min_abundance_threshold: Optional[float] = None,
    distribution_abundance: float = 1.0,
    is_abundance_sum: bool = False,
    output_masses_for_neutron_offset: bool = False,
    neutron_mass: float = constants.NEUTRON_MASS,
    precision: Optional[int] = None,
) -> List[Tuple[float, float]]:
    """
    Calculate the isotopic distribution for a given formula.

    :param chemical_formula: The chemical formula. Can contain non integer element counts (will be rounded).
    :type chemical_formula: Dict[str, int]
    :param max_isotopes: The maximum number of isotopes to track during convolution. For very complex compositions,
                         it can be useful to bound this value to avoid memory/speed issues.
                         If None, no limit is applied. It is not recommended to use this to limit the number of
                         returned isotopes since it could affect the acuracy of the abundance values. Instead,
                         set `min_abundance_threshold` to limit the number of isotopes returned.
    :type max_isotopes: Optional[float]
    :param min_abundance_threshold: The minimum abundance of an isotope to keep. If None, no threshold is applied.
    :type min_abundance_threshold: Optional[float]
    :param distribution_resolution: The resolution of the distribution. If None, no rounding is applied. Applied
                                    during convolution to limit the number of isotopes tracked. For example,
                                    setting this value to the mass-precision of a mass spectrometer will result in
                                    similar distributions to those seen in the mass spectrometer.
    :type distribution_resolution: Optional[float]
    :param use_neutron_count: Whether to use neutron offsets instead of masses.
    :type use_neutron_count: bool
    :param conv_min_abundance_threshold: The minimum abundance of an isotope to keep during convolution.
                                         If None, no threshold is applied. Works in tandem with `max_isotopes` to limit
                                         the total number of isotopes tracked during convolution. Use
                                         `min_abundance_threshold` to limit the number of isotopes returned.
    :type conv_min_abundance_threshold: float
    :param distribution_abundance: The abundance of the distribution. If None, this will be set to 1.
    :type distribution_abundance: float
    :param is_abundance_sum: if False (Default), the largest peak will be equal to distribution_abundance and all
                                other peaks will be scaled by this value. If True, the total abundances will sum to
                                distribution_abundance.
    :type is_abundance_sum: bool
    :param output_masses_for_neutron_offset: Whether to output masses for when using neutron offsets.
    :type output_masses_for_neutron_offset: bool
    :param neutron_mass: The mass of a neutron. Default is the value is pt.NEUTRON_MASS. Might want to set to C13
                         neutron mass (1.003350), the Averagine peptide neutron mass (1.002856).
    :type neutron_mass: float
    :param precision: The precision of the output. Applied at the end of the function (not during convolution).
    :type precision: Optional[int]

    :return: The isotopic distribution.
    :rtype: List[Tuple[float, float]]

    .. code-block:: python

        # Example usage
        >>> formula = {'C': 12, 'H': 6, 'N': 3}
        >>> isotopic_distribution(formula, 3, 0.0, 5)
        [(192.05617, 1.0), (193.05321, 0.010959894014211729), (193.05952, 0.1297887395127868)]

        >>> formula = {'C': 12, 'H': 6, 'N': 3, 'Li': 0}
        >>> isotopic_distribution(formula, 3, 0.0, 5)
        [(192.05617, 1.0), (193.05321, 0.010959894014211729), (193.05952, 0.1297887395127868)]

        >>> formula = {'C': 100, 'H': 100, 'N': 100, 'O': 100}
        >>> isotopic_distribution(formula, 3, 0.0, 5)
        [(4300.58136, 0.9245794392523361), (4301.58471, 1.0), (4302.58807, 0.5353785504902455)]

        # Example usage
        >>> formula = {'C': 100, 'H': 150, 'N': 15, 'O': 20}
        >>> isotopic_distribution(formula, 3, 0.0, 5)
        [(1881.11815, 0.9245794392523363), (1882.1215, 1.0), (1883.12486, 0.5353785504902455)]

        # Example usage
        >>> formula = {'C': 100, 'H': 150, 'N': 15, 'O': 20}
        >>> isotopic_distribution(formula, 3, 0.0, 5, True)
        [(0, 0.8611463538706062), (1, 1.0), (2, 0.6108892554305453)]

        # Example usage with negative values
        >>> formula = {'C': 100, 'Li': 10, 'N': 15, 'O': 20}
        >>> isotopic_distribution(formula, 3, 0.0, 5, True)
        [(-1, 0.5556664041575184), (0, 1.0), (1, 0.8123504851488393)]

        # example with float values
        >>> formula = {'C': 100.2, 'Li': 10, 'N': 15, 'O': 20}
        >>> isotopic_distribution(formula, 3, 0.0, 5, True,precision=3)
        [(-1, 0.556), (0, 1.0), (1, 0.812)]

        # example with abundance sum: peaks will be scalled by distribution_abundance.
        #With the largest peak beign equal to distribution_abundance
        >>> formula = {'C': 12, 'H': 6, 'N': 3}
        >>> isotopic_distribution(formula, 3, 0.0, 5, distribution_abundance=100, precision=3)
        [(192.056, 100.0), (193.053, 1.096), (193.06, 12.979)]

        # example with abundance sum (summed intensity), total abundances will sum to distribution_abundance
        >>> formula = {'C': 12, 'H': 6, 'N': 3}
        >>> isotopic_distribution(formula, 3, 0.0, 5, distribution_abundance=100, precision=3, is_abundance_sum=True)
        [(192.056, 87.662), (193.053, 0.961), (193.06, 11.378)]

        >>> formula = {'C': 20, 'H': 20.5, 'N': 20, 'O': 20}
        >>> isotopic_distribution(formula, 5, 0.0, 0, True, output_masses_for_neutron_offset=True, precision=3)
        [(860.625, 1.0), (861.633, 0.299), (862.642, 0.085), (863.651, 0.016), (864.659, 0.003)]

        >>> formula = {'C': 20, 'H': 20, 'N': 20, 'O': 20}
        >>> isotopic_distribution(formula, 5, 0.0, 0, False, output_masses_for_neutron_offset=True, precision=3)
        [(860.0, 1.0), (861.0, 0.299), (862.0, 0.085), (863.0, 0.016), (864.0, 0.003)]

        >>> formula = {'C': 20, 'H': 20, 'N': 20, 'O': 20}
        >>> isotopic_distribution(formula, 5, 0.0, 0, True, output_masses_for_neutron_offset=True, precision=3)
        [(860.116, 1.0), (861.125, 0.299), (862.134, 0.085), (863.142, 0.016), (864.151, 0.003)]

        # Example usage with floating values
        >>> isotopic_distribution({'C': 100, 'H': 150.5, 'N': 15, 'O': 20}, 3, 0.0, 5, False)
        [(1881.622062516115, 0.9245794392523363), (1882.625412516115, 1.0), (1883.6287725161149, 0.5353785504902455)]

        # Example usage with floating values
        >>> isotopic_distribution({'C': -100, 'H': 150.5, 'N': 15, 'O': 20}, 3, 0.0, 5, False)
        Traceback (most recent call last):
        ValueError: Negative values are not allowed in the chemical formula.

    """

    if min_abundance_threshold is None:
        min_abundance_threshold = 0.0

    if max_isotopes is None:
        max_isotopes = sys.maxsize

    # Just use these to correct mass
    electron_count = chemical_formula.pop("e", 0)
    proton_count = chemical_formula.pop("p", 0)
    neutron_count = chemical_formula.pop("n", 0)
    particle_mass_offset = (
        (proton_count * constants.PROTON_MASS)
        + (neutron_count * constants.NEUTRON_MASS)
        + (electron_count * constants.ELECTRON_MASS)
    )

    # check if any values are negative:
    if any(v < 0 for v in chemical_formula.values()):
        raise ValueError("Negative values are not allowed in the chemical formula.")

    # Remove 0 values from the chemical formula
    for elem in list(chemical_formula.keys()):
        if chemical_formula[elem] == 0:
            del chemical_formula[elem]

    delta_mass = 0.0
    if not all(isinstance(v, int) for v in chemical_formula.values()):
        prior_mass = chem_mass(chemical_formula)
        chemical_formula = _fix_chemical_formula(chemical_formula, False)
        post_mass = chem_mass(chemical_formula)
        delta_mass = prior_mass - post_mass

    formula_mass = chem_mass(chemical_formula)

    # raise warning
    if delta_mass != 0.0:
        warnings.warn(
            f"The chemical formula has a mass difference of {delta_mass} Da. This is likely due to floating point errors. The mass will be corrected for this."
        )

    total_distribution = {0: 1.0}  # Start with a base distribution
    for element, count in chemical_formula.items():
        elemental_distribution = _calculate_elemental_distribution(
            element, count, use_neutron_count
        )
        total_distribution = _convolve_distributions(
            total_distribution,
            elemental_distribution,
            max_isotopes,
            conv_min_abundance_threshold,
            distribution_resolution,
        )

    # Normalize abundances before filtering to ensure consistent abundance calculation
    max_abundance = max(total_distribution.values())
    normalized_distribution = [
        (mass, abundance / max_abundance)
        for mass, abundance in sorted(total_distribution.items())
        if abundance / max_abundance >= min_abundance_threshold
    ]

    if delta_mass != 0.0:
        if not use_neutron_count:
            normalized_distribution = [
                (mass + delta_mass + particle_mass_offset, abundance)
                for mass, abundance in normalized_distribution
            ]
        elif output_masses_for_neutron_offset:
            normalized_distribution = [
                (mass + delta_mass + particle_mass_offset, abundance)
                for mass, abundance in normalized_distribution
            ]

    if output_masses_for_neutron_offset and use_neutron_count:
        normalized_distribution = [
            (formula_mass + offset * neutron_mass, abundance)
            for offset, abundance in normalized_distribution
        ]

    return _scale_isotope_abundances(
        normalized_distribution, distribution_abundance, is_abundance_sum, precision
    )


def merge_isotopic_distributions(
    *distributions: List[Tuple[float, float]], precision: Optional[float] = None
) -> List[Tuple[float, float]]:
    """
    Merge multiple isotopic distributions into one.

    :param distributions: The isotopic distributions to merge.
    :type distributions: List[Tuple[float, float]]
    :param precision: The precision of the output.
    :type precision: Optional[float]

    :return: The merged isotopic distribution.
    :rtype: List[Tuple[float, float]]

    .. code-block:: python

        # Example usage
        >>> d1 = [(1.0, 0.5), (2.0, 0.5)]
        >>> d2 = [(1.0, 0.5), (2.0, 0.5)]
        >>> merge_isotopic_distributions(d1, d2)
        [(1.0, 1.0), (2.0, 1.0)]

    """
    merged_distribution = {}
    for distribution in distributions:
        for mass, abundance in distribution:
            if precision is not None:
                mass = round(mass, precision)
            # Merge the distributions
            if mass in merged_distribution:
                merged_distribution[mass] += abundance
            else:
                merged_distribution[mass] = abundance

    return sorted(merged_distribution.items(), key=lambda x: x[0])


def estimate_isotopic_distribution(
    neutral_mass: float,
    max_isotopes: Union[int, None] = None,
    min_abundance_threshold: Union[float, None] = None,
    distribution_resolution: Union[int, None] = 5,
    use_neutron_count: bool = False,
    conv_min_abundance_threshold: Union[float, None] = None,
    distribution_abundance: float = 1.0,
    is_abundance_sum: bool = False,
    output_masses_for_neutron_offset: bool = False,
    neutron_mass: float = constants.NEUTRON_MASS,
    precision: Optional[int] = None,
) -> List[Tuple[float, float]]:
    """
    Predict the isotopic distribution of a molecule using the averagine model.

    :param neutral_mass: The total neutral mass of the molecule.
    :type neutral_mass: float
    :param max_isotopes: The maximum number of isotopes to keep. If None, no limit is applied.
    :type max_isotopes: int
    :param min_abundance_threshold: The minimum abundance of an isotope to keep. If None, no threshold is applied.
    :type min_abundance_threshold: float
    :param distribution_resolution: The resolution of the distribution. If None, no rounding is applied.
    :type distribution_resolution: int
    :param use_neutron_count: Whether to use neutron offsets instead of masses.
    :type use_neutron_count: bool
    :param conv_min_abundance_threshold: The minimum abundance of an isotope to keep during convolution. If None,
    no threshold is applied.
    :type conv_min_abundance_threshold: float
    :param distribution_abundance: The abundance of the distribution. If None, no abundance is applied.
    :type distribution_abundance: float
    :param is_abundance_sum: Whether to sum the abundances to 1.
    :type is_abundance_sum: bool
    :param output_masses_for_neutron_offset: Whether to output masses for neutron offsets.
    :type output_masses_for_neutron_offset: bool
    :param neutron_mass: The mass of a neutron.
    :type neutron_mass: float
    :param precision: The precision of the output.
    :type precision: Optional[int]

    :return: A list of tuples with (mass, abundance) representing the predicted isotopic distribution.
    :rtype: List[Tuple[float, float]]

    .. code-block:: python

        # Example usage
        >>> estimate_isotopic_distribution(800, 3, 0.0, 5, precision=3)
        [(800.0, 1.0), (801.003, 0.389), (802.007, 0.074)]

        # Example usage
        >>> estimate_isotopic_distribution(800, 3, 0.0, 5, True, precision=3)
        [(0, 1.0), (1, 0.437), (2, 0.116)]

    """
    # Calculate the total number of each atom in the molecule based on its molecular mass
    total_atoms = estimate_comp(neutral_mass)

    distributions = isotopic_distribution(
        total_atoms,
        max_isotopes,
        min_abundance_threshold,
        distribution_resolution,
        use_neutron_count,
        conv_min_abundance_threshold,
        distribution_abundance,
        is_abundance_sum,
        output_masses_for_neutron_offset,
        neutron_mass,
        precision,
    )

    return distributions


def _convolve_distributions(
    dist1: Dict[float, float],
    dist2: Dict[float, float],
    max_isotopes: Union[int, None],
    min_abundance_threshold: Union[float, None],
    distribution_resolution: Union[int, None],
) -> Dict[float, float]:
    """
    Convolve two distributions to calculate the distribution of their sum.

    :param dist1: The first distribution.
    :type dist2: Dict[float, float]
    :param dist2: The second distribution.
    :type dist2: Dict[float, float]
    :param max_isotopes: The maximum number of isotopes to keep. If None, no limit is applied.
    :type max_isotopes: Union[int, None]
    :param min_abundance_threshold: The minimum abundance of an isotope to keep. If None, no threshold is applied.
    :type min_abundance_threshold: Union[float, None]
    :param distribution_resolution: The resolution of the distribution. If None, no rounding is applied.
    :type distribution_resolution: Union[int, None]

    :return: A convolved isotopic distribution as a mass-to-abundance mapping.
    :rtype: Dict[float, float]

    .. code-block:: python

        # Example usage
        >>> d1 = {1.0: 0.5, 2.0: 0.5}
        >>> d2 = {1.0: 0.5, 2.0: 0.5}
        >>> _convolve_distributions(d1, d2, 3, 0.0, 5)
        {3.0: 0.5, 2.0: 0.25, 4.0: 0.25}

    """

    if min_abundance_threshold is None:
        min_abundance_threshold = 0.0

    if max_isotopes is None:
        max_isotopes = sys.maxsize

    result = {}
    for mass1, abundance1 in dist1.items():
        for mass2, abundance2 in dist2.items():
            new_mass = mass1 + mass2
            if distribution_resolution is not None:
                new_mass = round(new_mass, distribution_resolution)
            new_abundance = abundance1 * abundance2
            if new_abundance >= min_abundance_threshold:
                if new_mass in result:
                    result[new_mass] += new_abundance
                else:
                    result[new_mass] = new_abundance

    # Apply max isotopes limit and sort by abundance
    if max_isotopes != sys.maxsize:
        sorted_result = sorted(result.items(), key=lambda x: x[1], reverse=True)
        # Retain only the top `max_isotopes` isotopes based on abundance
        return dict(sorted_result[:max_isotopes])

    return result


def _calculate_elemental_distribution(
    element: str,
    count: int,
    use_neutron_count: bool,
    min_abundance_threshold: float = 10e-9,
) -> Dict[float, float]:
    """
    Calculate the isotopic distribution for an element.

    :param element: The element.
    :type element: str
    :param count: The number of atoms.
    :type count: int
    :param use_neutron_count: Whether to use neutron offsets instead of masses.
    :type use_neutron_count: bool

    :return: The isotopic distribution.
    :rtype: Counter[float, float]

    .. code-block:: python

        # Example usage
        >>> _calculate_elemental_distribution('C', 2, False)
        {24.0: 0.9787144899999999, 25.00335483507: 0.02117102, 26.00670967014: 0.00011448999999999998}

        # Example using neutron count
        >>> _calculate_elemental_distribution('C', 2, True)
        {0: 0.9787144899999999, 1: 0.02117102, 2: 0.00011448999999999998}

    """

    if use_neutron_count is True:
        isotopes = constants.ATOMIC_SYMBOL_TO_ISOTOPE_NEUTRON_OFFSETS_AND_ABUNDANCES[
            element
        ]
    else:
        isotopes = constants.ATOMIC_SYMBOL_TO_ISOTOPE_MASSES_AND_ABUNDANCES[element]

    # Start with a distribution for an element not present (mass=0, abundance=1)
    distribution = {0: 1.0}
    for _ in range(count):
        # Update the distribution by convolving it with the isotopes' distribution each time
        isotope_distribution = dict(isotopes)
        distribution = _convolve_distributions(
            distribution, isotope_distribution, None, min_abundance_threshold, None
        )
    return distribution


def _fix_chemical_formula(
    chemical_formula: Dict[str, float], add_hydrogens: bool = True
) -> Dict[str, int]:
    """
    Fix a chemical formula, by rounding the atom counts to the nearest integer and adding hydrogen atoms to reach the
    correct molecular mass.

    :param chemical_formula: The chemical formula.
    :type chemical_formula: Dict[str, float]
    :param add_hydrogens: Whether to add hydrogen atoms to reach the correct molecular mass.
    :type add_hydrogens: bool
    :return: The fixed chemical formula.
    :rtype: Dict[str, int]

    .. code-block:: python

        # Example usage
        >>> _fix_chemical_formula({'C': 12.0, 'H': 6.0, 'N': 3.0})
        {'C': 12, 'H': 6, 'N': 3}

        # Example usage
        >>> _fix_chemical_formula({'C': 12.1, 'H': 6.0, 'N': 3.0})
        {'C': 12, 'H': 7, 'N': 3}

        # Example usage
        >>> _fix_chemical_formula({'C': 12.1, 'H': 6.0, 'N': 3.9})
        {'C': 12, 'H': 6, 'N': 4}

    """

    starting_mass = chem_mass(chemical_formula)

    # get the floor of the total atoms
    total_atoms = {k: round(v) for k, v in chemical_formula.items()}

    if add_hydrogens:
        if "H" not in total_atoms:
            total_atoms["H"] = 0

        # add hydrogen's till the molecular mass is reached
        total_atoms["H"] += int(
            (starting_mass - chem_mass(total_atoms))
            / constants.ISOTOPIC_ATOMIC_MASSES["H"]
        )

    return total_atoms


def _scale_isotope_abundances(
    isotopes: List[Tuple[float, float]],
    distribution_abundance: float,
    is_abundance_sum: bool,
    precision: Optional[int] = None,
) -> List[Tuple[float, float]]:
    """
    Scale the abundances of the isotopes to the distribution abundance.
    """

    if is_abundance_sum is True:
        # isotope abudnances should sum to 1
        total_abundance = sum(abundance for _, abundance in isotopes)
        isotopes = [(mass, abundance / total_abundance) for mass, abundance in isotopes]

    # multiply the abundances by the abundance factor
    isotopes = [
        (mass, abundance * distribution_abundance) for mass, abundance in isotopes
    ]

    if precision is not None:
        isotopes = [
            (round(mass, precision), round(abundance, precision))
            for mass, abundance in isotopes
        ]

    return isotopes
