from typing import List

import numpy as np

from .constants import *
from peptacular.sequence import parse_modified_sequence, strip_modifications, sequence_generator, \
    identify_cleavage_sites


# TODO: Remove numpy dependency and separate index functions into separate project


def filter_by_mass(sequences: list[str], min_mass: float = None, max_mass: float = None) -> List[str]:
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


def calculate_mass_array(sequence: str, monoisotopic: bool = True):
    if monoisotopic is True:
        aa_masses = MONO_ISOTOPIC_AA_MASSES
    else:
        aa_masses = AVERAGE_AA_MASSES

    sequence_mods = parse_modified_sequence(sequence)
    stripped_sequence = strip_modifications(sequence)

    mass_arr = np.zeros(len(stripped_sequence), dtype=np.float32)
    for i, aa in enumerate(stripped_sequence):

        mod_mass = 0
        if i in sequence_mods:
            mod_mass += float(sequence_mods[i])

        if i == 0 and -1 in sequence_mods:
            mod_mass += float(sequence_mods[-1])

        mass_arr[i] = np.float32(aa_masses[aa] + mod_mass)

    mass_arr = np.array(mass_arr, dtype=np.float32)

    return mass_arr


def create_ion_table(sequence: str, max_len: int = 50, ion_type: str = 'y', charge: int = 0,
                     enzyme: str = 'non-specific'):
    if ion_type == 'y':
        sequence = sequence[::-1]

    sequence_cumulative_sum = np.cumsum(calculate_mass_array(sequence))
    len_sequence = len(sequence)

    # Add padding zeros to make our life easier when computing window sums
    sequence_cumulative_sum = np.pad(sequence_cumulative_sum, (1, 0), 'constant')

    # Window start and end indices
    window_start = np.arange(len_sequence)[:, None]
    window_end = window_start + np.arange(1, max_len + 1)

    # Limit window_end indices to array bounds
    window_end = np.minimum(window_end, len_sequence)
    # Compute window sums
    ion_table = sequence_cumulative_sum[window_end] - sequence_cumulative_sum[window_start]

    # Add (HYDROGEN * 2 + OXYGEN) to the computed mass, only where window_end > window_start
    if ion_type == 'y':
        ion_table[window_end > window_start] += (MONO_ISOTOPIC_ATOMIC_MASSES['HYDROGEN'] * 2 +
                                                 MONO_ISOTOPIC_ATOMIC_MASSES['OXYGEN']) + \
                                                (MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * charge)
    elif ion_type == 'b':
        ion_table[window_end > window_start] += (MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * charge)
    else:
        raise ValueError("Invalid ion type. Must be either 'y' or 'b'")

    if charge != 0:
        ion_table = ion_table / charge

    # set invalid indexes to 0
    for i, j in enumerate(range(len(ion_table) - 1, len(ion_table) - min(max_len, len_sequence) - 1, -1)):
        ion_table[j, i + 1:] = 0

    if len_sequence < max_len:
        ion_table[:, len_sequence:] = 0

    if enzyme != 'non-specific':
        filter_ion_table(ion_table, sequence, enzyme)

    return ion_table


def filter_ion_table(ion_table: np.ndarray, sequence: str, enzyme: str):
    enzyme_sites = set(identify_cleavage_sites(sequence, enzyme))

    for i in range(len(ion_table)):
        if i not in enzyme_sites:
            ion_table[i, :] = 0


def get_sorted_ion_table_indexes(matrix: np.ndarray) -> np.ndarray:
    # Get 2D indexes without flattening
    # Apply mass filters if provided
    indexes_2d = np.argwhere(matrix)

    # Sort by values at these indexes
    sorted_2d_indexes = indexes_2d[np.argsort(matrix[indexes_2d[:, 0], indexes_2d[:, 1]], kind='mergesort')]
    return sorted_2d_indexes
