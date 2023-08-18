from typing import List

from .constants import MONO_ISOTOPIC_ATOMIC_MASSES, AVERAGE_ATOMIC_MASSES, AVERAGE_AA_MASSES, MONO_ISOTOPIC_AA_MASSES, \
    ION_ADJUSTMENTS, UWPR_MONO_ISOTOPIC_ATOMIC_MASSES, UWPR_AVERAGE_AA_MASSES, UWPR_AVERAGE_ATOMIC_MASSES, \
    UWPR_MONO_ISOTOPIC_AA_MASSES
from .sequence import parse_modifications, strip_modifications
from .util import validate_ion_type


def filter_by_mass(sequences: List[str], min_mass: float = 0, max_mass: float = float('inf'),
                   charge: int = 0, ion_type: str = 'y', monoisotopic: bool = True) -> List[str]:
    """
    Filters a list of sequences by mass.

    :param sequences: A list of sequences to be filtered.
    :type sequences: list[str]
    :param min_mass: Minimum mass of the returned sequences. Default is 0.
    :type min_mass: float
    :param max_mass: Maximum mass of the returned sequences. Default is infinity.
    :type max_mass: float
    :param charge: Peptide charge. Default is 0.
    :type charge: int
    :param ion_type: Ion type ('a', 'b', 'c', 'x', 'y', 'z'). Default is 'y'.
    :type ion_type: str
    :param monoisotopic: If true, uses monoisotopic masses. Defaults to True.
    :type monoisotopic: bool
    :return: A list of sequences that fall within the specified mass range.
    :rtype: list[str]
    """

    return [sequence for sequence in sequences if min_mass <=
            calculate_mass(sequence=sequence, charge=charge, ion_type=ion_type, monoisotopic=monoisotopic) <= max_mass]


def calculate_mass(sequence: str, charge: int = 0, ion_type: str = 'y', monoisotopic: bool = True,
                   uwpr_mass: bool = False) -> float:
    """
    Calculate the mass of a peptide sequence.

    :param sequence: Peptide sequence, possibly including modifications.
    :type sequence: str
    :param charge: Peptide charge. Default is 0.
    :type charge: int
    :param ion_type: Ion type ('a', 'b', 'c', 'x', 'y', 'z'). Default is 'y'.
    :type ion_type: str
    :param monoisotopic: If true, uses monoisotopic masses. Defaults to True.
    :type monoisotopic: bool
    :param uwpr_mass: If true, uses uwpr masses. If false, uses Pyteomics masses. Defaults to False.
    :type uwpr_mass: bool
    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.
    :return: The calculated mass of the peptide sequence.
    :rtype: float
    """

    validate_ion_type(ion_type=ion_type)

    # Select mass set based on monoisotopic flag
    atomic_masses, aa_masses = (MONO_ISOTOPIC_ATOMIC_MASSES, MONO_ISOTOPIC_AA_MASSES) \
        if monoisotopic is True else (AVERAGE_ATOMIC_MASSES, AVERAGE_AA_MASSES)

    if uwpr_mass is True:
        atomic_masses, aa_masses = (UWPR_MONO_ISOTOPIC_ATOMIC_MASSES, UWPR_MONO_ISOTOPIC_AA_MASSES) \
            if monoisotopic is True else (UWPR_AVERAGE_ATOMIC_MASSES, UWPR_AVERAGE_AA_MASSES)

    # Parse modifications and strip them from sequence
    mods = parse_modifications(sequence=sequence)
    stripped_sequence = strip_modifications(sequence=sequence)

    # Calculate mass
    mass = sum(aa_masses[aa] for aa in stripped_sequence)
    mass += sum(float(value) for value in mods.values())
    mass += (charge * atomic_masses['PROTON'])

    return mass + ION_ADJUSTMENTS[ion_type]


def calculate_mz(sequence: str, charge: int = 0, ion_type: str = 'y',
                 monoisotopic: bool = True, uwpr_mass: bool = False) -> float:
    """
    Calculate the m/z (mass-to-charge ratio) of a peptide sequence.

    :param sequence: Peptide sequence which may include modifications.
    :type sequence: str
    :param charge: Charge of the peptide. Default is 0.
    :type charge: int
    :param ion_type: Ion type, can be one of 'a', 'b', 'c', 'x', 'y', 'z'. Default is 'y'.
    :type ion_type: str
    :param monoisotopic: If true, uses monoisotopic masses. Defaults to True.
    :type monoisotopic: bool
    :param uwpr_mass: If true, uses uwpr masses. If false, uses Pyteomics masses. Defaults to False.
    :type uwpr_mass: bool
    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.
    :return: The m/z of the peptide sequence. If charge is 0, returns the mass.
    :rtype: float
    """

    mass = calculate_mass(sequence=sequence, charge=charge, ion_type=ion_type,
                          monoisotopic=monoisotopic, uwpr_mass=uwpr_mass)
    return mass if charge == 0 else mass / charge
