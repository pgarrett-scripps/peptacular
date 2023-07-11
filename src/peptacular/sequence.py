import re
from copy import deepcopy
from typing import Dict, List, Any, Tuple, Union

import numpy as np
import regex as reg

from .constants import PROTEASES, MONO_ISOTOPIC_ATOMIC_MASSES, \
    MONO_ISOTOPIC_AA_MASSES, AVERAGE_AA_MASSES, AVERAGE_ATOMIC_MASSES

from .util import check_parentheses


def convert_ip2_mod_to_uniprot_mod(sequence: str, mod_dict: Dict[str, Any]):
    """
    Convert a peptide sequence with IP2-style modifications to UniProt-style modifications.

    This function takes a peptide sequence in string format and a dictionary of IP2-style
    modifications and their corresponding UniProt-style modifications, and returns a
    modified peptide sequence with the IP2-style modifications replaced with their
    corresponding UniProt-style modifications.
    """

    for ip2_mod, uniprot_mod in mod_dict.items():
        ip2_mod = f'({ip2_mod})'
        uniprot_mod = f'({uniprot_mod})'
        sequence = sequence.replace(ip2_mod, uniprot_mod)

    return sequence


def convert_to_ms2_pip_style(sequence: str) -> str:
    """
    Convert a peptide sequence to MS2PIP-style format.

    This function takes a peptide sequence in string format and converts it to a modified
    format that is used by the MS2PIP prediction software. The modified format consists of
    a series of location and modification pairs, separated by a pipe symbol ('|').

    Example:
    >>> convert_to_ms2_pip_style('PEPTIDEK')
    '2|Phospho|6|Lysine'
    """

    mod_dict = parse_modified_sequence(sequence)
    locs, mods = [], []
    for loc in mod_dict:
        mod = mod_dict[loc]
        ms2_pip_style_loc = loc + 1
        if ms2_pip_style_loc == len(strip_modifications(sequence)):
            ms2_pip_style_loc = -1
        locs.append(ms2_pip_style_loc)
        mods.append(mod)

    return '|'.join([f'{loc}|{mod}' for loc, mod in zip(locs, mods)])


def parse_modified_sequence(sequence: str) -> Dict[int, Any]:
    """
    This function reads a peptide sequence with modifications and returns a dictionary
    with the modification values indexed by the position of the modified amino acid.
    """

    if check_parentheses(sequence) is False:
        raise ValueError(f'Incorrect modification notation in peptide sequence : {sequence}!')

    matches = re.finditer(r'\(([^)]*)\)', sequence)

    modifications = {}
    mod_offset = 0

    # Loop through the substrings
    for match in matches:
        modifications[match.start() - mod_offset - 1] = match.group()[1:-1]
        mod_offset += match.end() - match.start()
    return modifications


def create_modified_sequence(sequence: str, modifications: Dict[int, Any]) -> str:
    """
    Creates a modified amino acid sequence from an unmodified amino acid sequence and a dictionary of modifications.

    The modifications are specified as a dictionary where the keys are the indices of the modified amino acids
    in the peptide sequence, and the values are the modifications to apply at those indices. The modifications are
    added to the peptide sequence in descending order of index, so that the indices remain valid after each
    modification is applied.
    """

    modified_sequence = []
    prev_index = 0

    # Sort the modifications by index in descending order
    for i, mod in sorted(modifications.items()):

        if i < -1 or i >= len(sequence):
            raise ValueError(f'Index of modification: {i} is invalid for peptide sequence: {sequence}')

        # Insert the modification into the modified peptide sequence
        modified_sequence.append(sequence[prev_index: i + 1])
        modified_sequence.append(f"({mod})")
        prev_index = i + 1

    modified_sequence.append(sequence[prev_index:])
    return ''.join(modified_sequence)


def strip_modifications(sequence: str) -> str:
    """
    Removes any non-amino-acid characters from the given peptide sequence.
    """

    if check_parentheses(sequence) is False:
        raise ValueError(f'Incorrect modification notation in peptide sequence : {sequence}!')

    unmodified_sequence = re.sub(r'\([^)]*\)', '', sequence)
    return unmodified_sequence


def get_left_semi_sequences(sequence: str, min_len: int = None, max_len: int = None, max_semi_offset: int = None):
    """
    Returns a set of left substrings of string `sequence` that have lengths between `min_len` and `max_len`.

    Example:
    >>> get_left_semi_sequences("abc", 1, 2)
    {'a', 'ab'}
    """

    if min_len is None:
        min_len = 1
    if max_len is None:
        max_len = len(sequence) - 1

    cnt = 0
    sequences = set()
    for i in range(min_len, min(max_len, len(sequence) - 1) + 1)[::-1]:
        sequences.add(sequence[0:i])
        cnt += 1

        if max_semi_offset is not None and cnt >= max_semi_offset:
            break

    return sequences


def get_right_semi_sequences(sequence: str, min_len: int = None, max_len: int = None, max_semi_offset: int = None):
    """
    Returns a set of right substrings of string `sequence` that have lengths between `min_len` and `max_len`.

    Example:
    >>> get_right_semi_sequences("abc", 1, 2)
    {'c', 'bc'}
    """
    if min_len is None:
        min_len = 1
    if max_len is None:
        max_len = len(sequence) - 1
    end = min(max_len, len(sequence) - 1)

    cnt = 0
    sequences = set()
    for i in range(min_len, end + 1)[::-1]:
        sequences.add(sequence[-i:])
        cnt += 1

        if max_semi_offset is not None and cnt >= max_semi_offset:
            break

    return sequences


def get_semi_sequences(sequence: str, min_len: int = None, max_len: int = None, max_semi_offset: int = None):
    """
    Returns a set of all semi enzymatic amino acid sequences of string `sequence` that have lengths
    between `min_len` and `max_len`.

    Example:
    >>> get_semi_sequences("abc", 1, 2)
    {'a', 'ab', 'c', 'bc'}
    """

    return get_left_semi_sequences(sequence, min_len, max_len, max_semi_offset).union(
        get_right_semi_sequences(sequence, min_len, max_len, max_semi_offset))


def get_non_enzymatic_sequences(sequence: str, min_len: int = None, max_len: int = None):
    """
    Returns a set of all non-enzymatic amino acid sequences of string `sequence` that have lengths
    between `min_len` and `max_len`.

    Example:
    >>> get_non_enzymatic_sequences("abc", 1, 2)
    {'a', 'b', 'c', 'ab', 'bc'}
    """

    if min_len is None:
        min_len = 1

    if max_len is None:
        max_len = len(sequence)

    substrings = set()
    for i in range(len(sequence)):
        for j in range(i + min_len, min(i + max_len + 1, len(sequence) + 1)):
            substrings.add(sequence[i:j])
    return substrings


def add_static_mods(sequence: str, mod_map: Dict[str, Any]):
    """
    Add static modifications to an amino acid sequence.

    Args:
        sequence (str): Original amino acid sequence.
        mod_map (Dict[str, float]): Dictionary mapping amino acids to the mass of their modifications.

    Returns:
        str: Modified amino acid sequence.
    """

    stripped_sequence = strip_modifications(sequence)
    sequence_mod_map = parse_modified_sequence(sequence)
    for aa, mass in mod_map.items():
        indexes = [i for i, x in enumerate(stripped_sequence) if x == aa]
        for index in indexes:
            sequence_mod_map[index] = mass

    return create_modified_sequence(stripped_sequence, sequence_mod_map)


def rec_mod_builder(mods: Dict[str, float], sequence: str, i: int, mod_dict: Dict[int, float], max_mod_count: int):
    """
    Recursive function to generate all possible combinations of modifications on an amino acid sequence.

    Args:
        mods (Dict[str, float]): Dictionary mapping amino acids to the mass of their modifications.
        sequence (str): amino acid sequence to modify.
        i (int): Current index on the amino acid sequence.
        mod_dict (Dict[int, float]): Dictionary mapping the indices of the amino acids to the mass of their
        modifications.
        max_mod_count (int): Maximum number of modifications allowed on the amino acid sequence.

    Yields:
        Dict[int, float]: Modified dictionary mapping the indices of the amino acid sequence to the mass of their
        modifications.
    """

    if i == len(sequence):
        yield mod_dict
        return
    original_mod_dict = deepcopy(mod_dict)
    for aa, mass in mods.items():

        if len(mod_dict) >= max_mod_count:
            break

        if sequence[i] == aa and i not in mod_dict:
            mod_dict[i] = mass
            yield from rec_mod_builder(mods, sequence, i + 1, mod_dict, max_mod_count)

    yield from rec_mod_builder(mods, sequence, i + 1, original_mod_dict, max_mod_count)


def add_variable_mods(sequence: str, mod_map: Dict[str, float], max_mods: int) -> List[str]:
    """
    Add variable modifications to an amino acid sequence.

    Args:
        sequence (str): Original amino acid sequence.
        mod_map (Dict[str, float]): Dictionary mapping amino acids to the mass of their modifications.
        max_mods (int): Maximum number of modifications allowed on the peptide.

    Returns:
        List[str]: List of all possible modified peptide sequences, each represented as a string.
    """

    original_mods = parse_modified_sequence(sequence)
    stripped_sequence = strip_modifications(sequence)
    mods = rec_mod_builder(mod_map, stripped_sequence, 0, original_mods, max_mods + len(original_mods))
    return [create_modified_sequence(stripped_sequence, mod) for mod in mods]


def calculate_mass(sequence: str, charge=0, ion_type: str = 'y', monoisotopic=True) -> float:
    """
    Calculate the mass of a peptide sequence.

    Args:
        sequence (str): Peptide sequence. The sequence may include modifications.
        charge (int, optional): Charge of the peptide. Default is 0.
        ion_type (str, optional): Ion type. Can be either 'b' or 'y'. Default is 'y'.

    Returns:
        float: The mass of the peptide sequence.
    """

    if monoisotopic is True:
        atomic_masses = MONO_ISOTOPIC_ATOMIC_MASSES
        aa_masses = MONO_ISOTOPIC_AA_MASSES
    else:
        atomic_masses = AVERAGE_ATOMIC_MASSES
        aa_masses = AVERAGE_AA_MASSES

    mods = parse_modified_sequence(sequence)
    stripped_sequence = strip_modifications(sequence)

    mass = sum(aa_masses[aa] for aa in stripped_sequence)
    mass += sum(float(value) for value in mods.values())
    mass += (charge * atomic_masses['PROTON'])

    if ion_type == 'a':
        mass -= (atomic_masses['CARBON'] + atomic_masses['OXYGEN'])
    elif ion_type == 'b':
        pass
    elif ion_type == 'c':
        mass += (atomic_masses['HYDROGEN'] * 3 + atomic_masses['NITROGEN'])
    elif ion_type == 'x':
        mass += (atomic_masses['CARBON'] + atomic_masses['OXYGEN'] * 2)
    elif ion_type == 'y':
        mass += (atomic_masses['HYDROGEN'] * 2 + atomic_masses['OXYGEN'])
    elif ion_type == 'z':
        mass += (atomic_masses['OXYGEN'] - atomic_masses['NITROGEN'] - atomic_masses['HYDROGEN'])

    return mass


def calculate_mz(sequence: str, charge=0, ion_type: str = 'y', monoisotopic=True) -> float:
    """
    Calculate the m/z (mass-to-charge ratio) of a peptide sequence.

    Args:
        sequence (str): Peptide sequence. The sequence may include modifications.
        charge (int, optional): Charge of the peptide. Default is 0.
        ion_type (str, optional): Ion type. Can be either 'b' or 'y'. Default is 'y'.

    Returns:
        float: The m/z of the peptide sequence.
    """

    mass = calculate_mass(sequence, charge, ion_type, monoisotopic)
    if charge == 0:
        return mass
    return mass / charge


def fragment_sequence(sequence: str, types=('b', 'y'), max_charge=1, monoisotopic=True):
    """
    Generates fragments of a given amino acid sequence based on the specified ion types and maximum charge.
    This function is a generator that yields the m/z (mass-to-charge ratio) of each fragment.

    Args:
        sequence (str): The amino acid sequence to be fragmented.
        types (tuple, optional): A tuple containing the types of ions to be considered in the fragmentation.
            Each ion type is represented by a single character. Defaults to ('b', 'y').
        max_charge (int, optional): The maximum charge to consider for the fragments. Defaults to 1.

    Yields:
        float: The calculated m/z for the fragment.
    """
    if any([ion in 'xyz' for ion in types]):
        for pep in sequence_generator(sequence, forward=True):
            for ion_type in types:
                if ion_type not in 'xyz':
                    continue
                for charge in range(1, max_charge + 1):
                    yield calculate_mz(pep, charge=charge, ion_type=ion_type, monoisotopic=monoisotopic)

    if any([ion in 'abc' for ion in types]):
        for pep in sequence_generator(sequence, forward=False):
            for ion_type in types:
                if ion_type not in 'abc':
                    continue
                for charge in range(1, max_charge + 1):
                    yield calculate_mz(pep, charge=charge, ion_type=ion_type)


def fragment_series(sequence: str, ion_type='y', charge=1, monoisotopic=True):
    if ion_type in 'xyz':
        for pep in sequence_generator(sequence, forward=True):
            yield calculate_mz(pep, charge=charge, ion_type=ion_type, monoisotopic=monoisotopic)

    if ion_type in 'abc':
        for pep in sequence_generator(sequence, forward=False):
            yield calculate_mz(pep, charge=charge, ion_type=ion_type, monoisotopic=monoisotopic)


def sequence_generator(sequence: str, forward: bool):
    """
    Generates amino acid sequences either from the front or the back.

    Args:
        sequence (str): The initial amino acid sequence.
        forward (bool): If True, start generating from the front; else, start from the back.

    Returns:
        generator: A generator yielding amino acid sequences.
    """
    if forward:
        return _sequence_generator_front(sequence)
    else:
        return _sequence_generator_back(sequence)


def _sequence_generator_back(sequence: str):
    """
    Generates amino acid sequences starting from the back of the sequence.

    Args:
        sequence (str): The initial amino acid sequence.

    Yields:
        str: The sequence with the back amino acid (and modification) removed.

    Returns:
        None: If the modification starts at the beginning of the sequence.
    """
    while sequence:
        index = 0
        if sequence[-1] == ')':
            start_of_modification = sequence.rfind('(')

            if start_of_modification == 0:
                return None

            index = start_of_modification

        yield sequence
        sequence = sequence[:index - 1]


def _sequence_generator_front(sequence: str):
    """
    Generates sub-sequences starting from the front of the sequence.

    Args:
        sequence (str): The initial sequence.

    Yields:
        str: The sequence with the front amino acid (and modification) removed.

    Returns:
        None: If the modification ends at the end of the sequence.
    """
    if sequence[0] == '(':
        yield sequence

        end_of_modification = sequence.find(')')
        sequence = sequence[end_of_modification + 2:]

        if sequence[0] == '(':
            end_of_modification = sequence.find(')')
            sequence = sequence[end_of_modification + 1:]

    while sequence:
        index = 0
        if len(sequence) > 1 and sequence[1] == '(':
            end_of_modification = sequence.find(')')

            if end_of_modification == len(sequence):
                return None

            index = end_of_modification

        yield sequence
        sequence = sequence[index + 1:]


def reverse_sequence(sequence: str):
    """
    Reverses the sequence, while preserving the position of any modifications.

    Args:
        sequence (str): The amino acid sequence to be reversed.

    Returns:
        str: The reversed sequence with modifications preserved.
    """
    mods = parse_modified_sequence(sequence)
    stripped_sequence = strip_modifications(sequence)[::-1]
    mod_reverse = {len(stripped_sequence) - 1 - k if k != -1 else -1: v for k, v in mods.items()}
    return create_modified_sequence(stripped_sequence, mod_reverse)


def shift_sequence_left(sequence: str, shift: int):
    """
    Shifts the sequence to the left by a given number of positions, while preserving the position of any modifications.

    Args:
        sequence (str): The sequence to be shifted.
        shift (int): The number of positions to shift the sequence to the left.

    Returns:
        str: The shifted sequence with modifications preserved.
    """
    mods = parse_modified_sequence(sequence)
    stripped_sequence = strip_modifications(sequence)
    seq_len = len(stripped_sequence)

    # Take modulus to ensure proper wrapping around
    effective_shift = shift % seq_len

    shifted_stripped_sequence = stripped_sequence[effective_shift:] + stripped_sequence[:effective_shift]

    # Shift mod positions, taking care to avoid negative keys
    shifted_sequence_mods = {(k - effective_shift) % seq_len if k != -1 else -1: v for k, v in mods.items()}

    return create_modified_sequence(shifted_stripped_sequence, shifted_sequence_mods)


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


def create_ion_table(sequence: str, max_len: int = 50, ion_type: str = 'y', charge: int = 0, enzyme: str = 'non-specific'):

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


def identify_cleavage_sites(protein_sequence: str, enzyme_regex: str):
    """
    Identifies cleavage sites in a protein sequence, based on enzyme_regex string. Since cleavages sites can only occur
    before or after a residues, cleavage_offset provides a way for the user to specify where the cleavage should occur
    in the matching regex.
    """

    if enzyme_regex in PROTEASES:
        enzyme_regex = PROTEASES[enzyme_regex]

    enzyme_sites = []
    for site in reg.finditer(enzyme_regex, protein_sequence, overlapped=True):
        enzyme_sites.append(site.span(0))
    return [site[0] + 1 for site in enzyme_sites]


def _get_spans(start_site: int, future_sites: List[int], missed_cleavages: int, last_sequence_index: int) -> \
        List[Tuple[int, int, int]]:
    """
    The function computes the spans by taking the start_site and the first missed_cleavages number of elements
    of future_sites and creating a tuple of start and end indices. If there are not enough elements in next_sites to
    fill the missed_cleavages, the end index of the tuple is set to the last index of the sequence, last_sequence_index.
    """

    spans = []
    for i, next_site in enumerate(future_sites[:missed_cleavages + 1]):
        spans.append((start_site, next_site, i))

    if len(spans) != missed_cleavages + 1:
        spans.append((start_site, last_sequence_index, len(spans)))

    return spans


def _digest_sequence(protein_sequence: str, enzyme_sites: List[int], missed_cleaves: int, min_len: int,
                     max_len: int) -> List[str]:
    """
    Digests a protein sequence according to the enzyme regex string and number of missed cleavages.
    """

    if len(enzyme_sites) == 0:
        if min_len <= len(protein_sequence) <= max_len:
            return [protein_sequence]
        else:
            return []

    peptide_spans = []

    cur_site = 0
    curr_enzyme_index = -1
    while True:
        future_start_sites = enzyme_sites[curr_enzyme_index + 1:]
        spans = _get_spans(cur_site, future_start_sites, missed_cleaves, len(protein_sequence))
        peptide_spans.extend(spans)
        curr_enzyme_index += 1

        if curr_enzyme_index >= len(enzyme_sites):
            break

        cur_site = enzyme_sites[curr_enzyme_index]

    peptides = []
    for span in peptide_spans:
        span_len = span[1] - span[0]
        if min_len <= span_len <= max_len:
            peptides.append(protein_sequence[span[0]:span[1]])
    return peptides


def digest_sequence(sequence: str, enzyme_regexes: Union[List[str], str], missed_cleavages: int, min_len: int,
                    max_len: int, semi_enzymatic: bool) -> List[str]:
    """
    A function that digests a given amino acid sequence based on the provided positive and negative regular expressions and
    the number of missed cleavages. It returns a list of tuples containing the digested_sequences and the number of missed
    cleavages. The returned digested_sequences will be filtered based on the provided minimum and maximum length.
    """

    if isinstance(enzyme_regexes, str):
        enzyme_regexes = [enzyme_regexes]
    digested_sequences = []
    for enzyme_regex in enzyme_regexes:

        if enzyme_regex == PROTEASES['non-specific']:
            digested_sequences.extend(
                [peptide for peptide in get_non_enzymatic_sequences(sequence, min_len, max_len)])

        else:
            cleavage_sites = identify_cleavage_sites(sequence, enzyme_regex)
            sequences = _digest_sequence(sequence, cleavage_sites, missed_cleavages, min_len, max_len)
            digested_sequences.extend(sequences)

            if semi_enzymatic is True:
                for sequence in sequences:
                    digested_sequences.extend(get_semi_sequences(sequence, min_len, max_len))

    return digested_sequences
