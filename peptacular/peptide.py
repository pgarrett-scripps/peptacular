import regex as reg

from . import constants


def parse_modified_peptide(peptide_sequence: str) -> dict[int, int]:
    """
    This function reads a peptide sequence with modifications and returns a dictionary
    with the modification values indexed by the position of the modified amino acid.

    :param peptide_sequence: The peptide sequence to read, with modifications indicated
                             by opening parenthesis followed by the modification value,
                             followed by a closing parenthesis, before the modified amino acid.
    :type peptide_sequence: str
    :return: A dictionary with the modification values indexed by the position of the modified
             amino acid.
    :rtype: dict[int, int]
    :raises ValueError: If the peptide sequence contains incorrect modification notation.
    """
    starting_par_count, ending_par_count = peptide_sequence.count('('), peptide_sequence.count(')')
    if starting_par_count != ending_par_count:
        raise ValueError(f'Incorrect modification notation in peptide sequence : {peptide_sequence}!')
    mod_count = ending_par_count

    modifications = {}

    # Initialize the index counter to 0
    index = -1

    # Loop through the substrings
    for s in peptide_sequence.split('('):  # Split the sequence into a list of substrings at each opening parenthesis
        # Check if the substring contains a closing parenthesis
        if ')' in s:
            # Split the substring at the closing parenthesis
            mod, peptide = s.split(')')

            # Convert the modification value to an integer
            mod = int(mod)

            # Add the modification to the dictionary, using the index as the key
            modifications[index] = mod

            # Increment the index by the length of the peptide substring
            index += len(peptide)
        else:
            # Increment the index by the length of the substring
            index += len(s)

    if len(modifications) != mod_count:
        raise ValueError(f'Incorrect modification notation in peptide sequence : {peptide_sequence}!')

    return modifications


def create_modified_peptide(unmodified_sequence: str, modifications: dict[int, int]) -> str:
    """
    Creates a modified peptide sequence from an unmodified peptide sequence and a dictionary of modifications.

    The modifications are specified as a dictionary where the keys are the indices of the modified amino acids
    in the peptide sequence, and the values are the modifications to apply at those indices. The modifications are
    added to the peptide sequence in descending order of index, so that the indices remain valid after each
    modification is applied.

    :param unmodified_sequence: The unmodified peptide sequence
    :param modifications: The modifications to apply to the peptide sequence
    :return: The modified peptide sequence
    :raises ValueError: If the index of a modification is invalid for the given peptide sequence
    """

    modified_sequence = []
    prev_index = 0

    # Sort the modifications by index in descending order
    for i, mod in sorted(modifications.items()):

        if i < 0 or i >= len(unmodified_sequence):
            raise ValueError(f'Index of modification: {i} is invalid for peptide sequence: {unmodified_sequence}')

        # Insert the modification into the modified peptide sequence
        modified_sequence.append(unmodified_sequence[prev_index: i + 1])
        modified_sequence.append(f"({mod})")
        prev_index = i + 1

    modified_sequence.append(unmodified_sequence[prev_index:])
    return ''.join(modified_sequence)


def strip_modifications(peptide_sequence: str) -> str:
    """
    Removes any non-amino-acid characters from the given peptide sequence.

    Args:
        peptide_sequence: The peptide sequence to be stripped of modifications.

    Returns:
        The peptide sequence with all non-amino-acid characters removed.
    """
    return ''.join([c for c in peptide_sequence if c in constants.amino_acids])


def get_modification(peptide_sequence: str, modification_regex: str, modification: float, regex_offset=0):
    unmodified_peptide_sequence = strip_modifications(peptide_sequence)
    sequon_sites = {}
    for site in reg.finditer(modification_regex, unmodified_peptide_sequence, overlapped=True):
        sequon_sites[site.span(0)[0] + regex_offset] = modification
    return sequon_sites
