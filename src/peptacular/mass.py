from typing import List

from .constants import MONO_ISOTOPIC_ATOMIC_MASSES, AVERAGE_ATOMIC_MASSES, AVERAGE_AA_MASSES, MONO_ISOTOPIC_AA_MASSES
from .sequence import parse_modified_sequence, strip_modifications, sequence_generator


def filter_by_mass(sequences: List[str], min_mass: float = None, max_mass: float = None) -> List[str]:
    """
    Filters a list of sequences by mass, returning only those that fall within the specified range.

    Args:
        sequences (list[str]): A list of sequences to be filtered.
        min_mass (float): Minimum mass of the returned sequences.
        max_mass (float): Maximum mass of the returned sequences.

    Returns:
        list[str]: A list of sequences that fall within the specified mass range.
    """

    min_mass = min_mass or 0
    max_mass = max_mass or float('inf')

    return [sequence for sequence in sequences if min_mass <= calculate_mass(sequence) <= max_mass]


def calculate_mass(sequence: str, charge=0, ion_type: str = 'y', monoisotopic=True) -> float:
    """
    Calculate the mass of a peptide sequence considering possible modifications, charge and ion type.

    Args:
        sequence (str): Peptide sequence, possibly including modifications.
        charge (int, optional): Peptide charge. Default is 0.
        ion_type (str, optional): Ion type ('a', 'b', 'c', 'x', 'y', 'z'). Default is 'y'.
        monoisotopic (bool, optional): If true, uses monoisotopic masses. Defaults to True.

    Returns:
        float: The calculated mass of the peptide sequence.
    """
    # Select mass set based on monoisotopic flag
    atomic_masses, aa_masses = (MONO_ISOTOPIC_ATOMIC_MASSES, MONO_ISOTOPIC_AA_MASSES) if monoisotopic else \
        (AVERAGE_ATOMIC_MASSES, AVERAGE_AA_MASSES)

    # Parse modifications and strip them from sequence
    mods = parse_modified_sequence(sequence)
    stripped_sequence = strip_modifications(sequence)

    # Calculate mass
    mass = sum(aa_masses[aa] for aa in stripped_sequence) + sum(float(value) for value in mods.values()) + \
           (charge * atomic_masses['PROTON'])

    # Adjust mass based on ion type
    ion_adjustments = {
        'a': -(atomic_masses['CARBON'] + atomic_masses['OXYGEN']),
        'c': atomic_masses['HYDROGEN'] * 3 + atomic_masses['NITROGEN'],
        'x': atomic_masses['CARBON'] + atomic_masses['OXYGEN'] * 2,
        'y': atomic_masses['HYDROGEN'] * 2 + atomic_masses['OXYGEN'],
        'z': atomic_masses['OXYGEN'] - atomic_masses['NITROGEN'] - atomic_masses['HYDROGEN']
    }

    return mass + ion_adjustments.get(ion_type, 0)


def calculate_mz(sequence: str, charge=0, ion_type: str = 'y', monoisotopic=True) -> float:
    """
    Calculate the m/z (mass-to-charge ratio) of a peptide sequence.

    Args:
        sequence (str): Peptide sequence which may include modifications.
        charge (int, optional): Charge of the peptide. Default is 0.
        ion_type (str, optional): Ion type, can be one of 'a', 'b', 'c', 'x', 'y', 'z'. Default is 'y'.
        monoisotopic (bool, optional): If true, uses monoisotopic masses. Defaults to True.

    Returns:
        float: The m/z of the peptide sequence. If charge is 0, returns the mass.
    """
    # Calculate mass
    mass = calculate_mass(sequence, charge, ion_type, monoisotopic)

    # Return mass if charge is zero; otherwise, return m/z
    return mass if charge == 0 else mass / charge


def fragment_sequence(sequence: str, types=('b', 'y'), max_charge=1, monoisotopic=True):
    """
    Generates fragments of a given amino acid sequence based on the specified ion types and maximum charge.
    This function is a generator that yields the m/z (mass-to-charge ratio) of each fragment.

    Args:
        sequence (str): The amino acid sequence to be fragmented.
        types (tuple, optional): A tuple containing the types of ions to be considered in the fragmentation.
            Each ion type is represented by a single character. Defaults to ('b', 'y').
        max_charge (int, optional): The maximum charge to consider for the fragments. Defaults to 1.
        monoisotopic (bool, optional): Indicates whether to calculate monoisotopic mass. Defaults to True.

    Yields:
        float: The calculated m/z for the fragment.
    """
    # Loop through the ion types and charges
    for ion_type in types:
        for charge in range(1, max_charge + 1):
            # If the ion type is in 'xyz' or 'abc', yield the series
            if ion_type in 'xyzabc':
                yield from fragment_series(sequence, ion_type=ion_type, charge=charge, monoisotopic=monoisotopic)


def fragment_series(sequence: str, ion_type='y', charge=1, monoisotopic=True):
    """
    Generates fragment series based on ion type and charge.

    Args:
        sequence (str): The amino acid sequence to be fragmented.
        ion_type (str, optional): The type of ions to be considered in the fragmentation.
            Each ion type is represented by a single character. Defaults to 'y'.
        charge (int, optional): The charge to consider for the fragments. Defaults to 1.
        monoisotopic (bool, optional): Indicates whether to calculate monoisotopic mass. Defaults to True.

    Yields:
        float: The calculated m/z for the fragment.
    """
    # Check if ion type is forward or reverse
    forward = ion_type in 'xyz'

    # Loop through the sequence to generate fragments and yield m/z
    for pep in sequence_generator(sequence, forward=forward):
        yield calculate_mz(pep, charge=charge, ion_type=ion_type, monoisotopic=monoisotopic)
