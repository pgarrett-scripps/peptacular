from typing import Tuple, List, Any, Set

import regex as reg

from .peptide import get_non_enzymatic_sequences, get_semi_sequences


def identify_cleavage_sites(protein_sequence: str, enzyme_regex: str, cleavage_offset=1):
    """
    Identifies cleavage sites in a protein sequence, based on enzyme_regex string. Since cleavages sites can only occur
    before or after a residues, cleavage_offset provides a way for the user to specify where the cleavage should occur
    in the matching regex.

    For Trypsin:
    identify_cleavage_sites('PEPKTIDEPERPTIDE', '([KR])([^P])', 1):
        (4, 16)

            site        site
            |           |
        PEPKTIDEPERPTIDRES
    """

    enzyme_sites = []
    for site in reg.finditer(enzyme_regex, protein_sequence, overlapped=True):
        enzyme_sites.append(site.span(0))
    return [site[0] + cleavage_offset for site in enzyme_sites]


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


def digest_protein_sequence(protein_sequence: str, enzyme_sites: List[int], missed_cleaves: int, min_len: int,
                            max_len: int) -> \
        List[Tuple[str, int]]:
    """
    Digests a protein sequence according to the enzyme regex string and number of missed cleavages.
    """

    if len(enzyme_sites) == 0:
        if min_len <= len(protein_sequence) <= max_len:
            return [(protein_sequence, 0)]
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
            peptides.append((protein_sequence[span[0]:span[1]], span[2]))
    return peptides


def flatten(nested_list: List[List[Any]]) -> List[Any]:
    """
    returns a flattened versions of the nested list
    """

    return [e for sublist in nested_list for e in sublist]


def combine_site_regexes(protein_sequence: str, positive_site_info: List[Tuple[str, int]],
                         negative_site_info: List[Tuple[str, int]]) -> List[int]:
    """
    The function uses the identify_cleavage_sites() function to find the positive and negative sites in the protein
    sequence using the regex strings and offsets provided in positive_site_info and negative_site_info.
    It then takes the set difference of the two sets of sites, resulting in a set of sites that are only
    found in the positive set and not the negative set. Finally, it returns the list of sites sorted in
    ascending order.
    """

    positive_sites = [identify_cleavage_sites(protein_sequence, site_info[0], site_info[1]) for site_info in
                      positive_site_info]
    negative_sites = [identify_cleavage_sites(protein_sequence, site_info[0], site_info[1]) for site_info in
                      negative_site_info]
    positive_sites = set(flatten(positive_sites))
    negative_sites = set(flatten(negative_sites))
    return sorted(list(positive_sites - negative_sites))


def filter_peptides(peptides: List[str]) -> Set[Tuple[str, int]]:
    peptide_dict = {}
    for (peptide, peptide_type) in peptides:
        found_peptide_type = peptide_dict.setdefault(peptide, peptide_type)

        if found_peptide_type == -1:
            peptide_dict[peptide] = peptide_type
        elif peptide_type != -1:
            continue
        else:
            if peptide_type < found_peptide_type:
                peptide_dict[peptide] = peptide_type

    return {(peptide, peptide_dict[peptide]) for peptide in peptide_dict}


def digest_protein(protein_sequence: str, enzyme_regexes: Tuple[List[Tuple[str, int]], List[Tuple[str, int]]],
                   missed_cleavages: int, min_len: int, max_len: int, non_enzymatic: bool, semi_enzymatic: bool) -> \
        List[Tuple[str, int]]:
    """
    A function that digests a given protein sequence based on the provided positive and negative regular expressions and
    the number of missed cleavages. It returns a list of tuples containing the peptides and the number of missed
    cleavages. The returned peptides will be filtered based on the provided minimum and maximum length.
    """

    if non_enzymatic is True:
        return [(peptide, -1) for peptide in get_non_enzymatic_sequences(protein_sequence, min_len, max_len)]

    sites = combine_site_regexes(protein_sequence, enzyme_regexes[0], enzyme_regexes[1])
    peptides = digest_protein_sequence(protein_sequence, sites, missed_cleavages, min_len, max_len)

    all_semi_peptides = []
    if semi_enzymatic is True:
        for (peptide, peptide_type) in peptides:
            semi_peptides = get_semi_sequences(peptide, min_len, max_len)
            semi_peptides.remove(peptide)

            semi_peptides = [(semi_peptide, -1) for semi_peptide in semi_peptides]
            all_semi_peptides.extend(semi_peptides)
    peptides.extend(all_semi_peptides)
    return peptides
