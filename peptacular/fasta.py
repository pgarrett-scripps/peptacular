from io import StringIO, TextIOWrapper
from typing import IO


def fasta_to_dict(fasta_io: IO[str]):
    """
    A function that reads a FASTA file from a file-like object, and returns a dictionary where the keys are the locus
    (protein id) and the values are a dictionary containing the sequence and the description of the protein.
    Empty lines will be ignored and the locus will be cleaned of leading ">", and any whitespace will be removed.
    The sequence of each protein will be concatenated.

    Parameters:
        fasta_io (StringIO | TextIOWrapper): A file-like object containing FASTA data

    Returns:
        dict: A dictionary where the keys are the locus (protein id) and the values are a dictionary containing
        the sequence and the description of the protein
    """

    locus_to_sequence_map = {}
    locus = None
    for line in fasta_io:
        if line == "":
            continue
        elif line[0] == ">":  # new protein
            line_elems = line.rstrip().split(" ")
            locus = line_elems[0].replace(">", "")
            description = " ".join(line_elems[1:])
            locus_to_sequence_map[locus] = {'sequence': "", 'description': description}
        else:  # protein sequence
            locus_to_sequence_map[locus]['sequence'] += line.rstrip()
    return locus_to_sequence_map
