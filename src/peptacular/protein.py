import copy
import re
from collections import defaultdict
from dataclasses import field, dataclass
from typing import Tuple, List, Any, Set, Dict, Generator, Union

import numpy as np
import pandas as pd
import regex as reg

from .peptide import get_non_enzymatic_sequences, get_semi_sequences, strip_modifications


def get_peptide_indexes_in_protein(substring, string):
    return [i.start() for i in re.finditer(substring, string)]


def calculate_protein_coverage(protein_sequence: str, peptides: list[str]):
    cov_arr = [0] * len(protein_sequence)
    for peptide in peptides:
        peptide_indexes = get_peptide_indexes_in_protein(peptide, protein_sequence)
        for peptide_index in peptide_indexes:
            cov_arr[peptide_index:peptide_index + len(peptide)] = [1] * len(peptide)
    return cov_arr


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


@dataclass
class PeptideProteinRelationship:
    protein_sequence_to_protein_id: Dict[str, int] = field(default_factory=lambda: dict())
    protein_locus_to_protein_id: Dict[str, int] = field(default_factory=lambda: dict())
    peptide_sequence_to_peptide_id: Dict[str, int] = field(default_factory=lambda: dict())

    protein_id_to_protein_sequence: Dict[int, str] = field(default_factory=lambda: dict())
    protein_id_to_protein_locus: Dict[int, str] = field(default_factory=lambda: dict())
    peptide_id_to_peptide_sequence: Dict[int, str] = field(default_factory=lambda: dict())

    protein_id_to_peptide_ids: Dict[int, Set[int]] = field(default_factory=lambda: defaultdict(set))
    peptide_id_to_protein_ids: Dict[int, Set[int]] = field(default_factory=lambda: defaultdict(set))

    def _get_id(self, item: str, seq_to_id: Dict[str, int], id_to_seq: Dict[int, str]) -> int:
        if item not in seq_to_id:
            id = len(seq_to_id)
            seq_to_id[item] = id
            id_to_seq[id] = item
        return seq_to_id[item]

    def add(self, protein_sequence: str, protein_locus: str, peptide_sequences: List[str]) -> None:
        protein_id = self._get_id(protein_locus, self.protein_locus_to_protein_id, self.protein_id_to_protein_locus)
        self.protein_sequence_to_protein_id[protein_sequence] = protein_id
        self.protein_id_to_protein_sequence[protein_id] = protein_sequence

        for peptide_sequence in peptide_sequences:
            peptide_id = self._get_id(peptide_sequence, self.peptide_sequence_to_peptide_id,
                                      self.peptide_id_to_peptide_sequence)

            self.protein_id_to_peptide_ids[protein_id].add(peptide_id)
            self.peptide_id_to_protein_ids[peptide_id].add(protein_id)

    def get_protein_ids(self, peptide_input: str) -> List[int]:
        if peptide_input in self.peptide_sequence_to_peptide_id:
            peptide_id = self.peptide_sequence_to_peptide_id.get(peptide_input, None)
        elif peptide_input in self.peptide_id_to_protein_ids:
            peptide_id = peptide_input
        else:
            peptide_id = None

        return [protein_id for protein_id in
                self.peptide_id_to_protein_ids[peptide_id]] if peptide_id is not None else []

    def get_proteins(self, peptide_input: str, output: str = 'p') -> List[Union[str, int]]:
        if output == 'p':
            return [self.protein_id_to_protein_sequence[protein_id] for protein_id in
                    self.get_protein_ids(peptide_input)]
        elif output == 'l':
            return [self.protein_id_to_protein_locus[protein_id] for protein_id in self.get_protein_ids(peptide_input)]
        elif output == 'i':
            return [protein_id for protein_id in self.get_protein_ids(peptide_input)]
        else:
            raise TypeError(f'Unsupported output: {output}')

    def get_peptide_ids(self, protein_input: str) -> List[int]:
        if protein_input in self.protein_sequence_to_protein_id:
            protein_id = self.protein_sequence_to_protein_id.get(protein_input, None)
        elif protein_input in self.protein_locus_to_protein_id:
            protein_id = self.protein_locus_to_protein_id.get(protein_input, None)
        else:
            protein_id = None
        return [peptide_id for peptide_id in
                self.protein_id_to_peptide_ids[protein_id]] if protein_id is not None else []

    def get_peptides(self, protein_input: str, output: str = 'p') -> List[Union[str, int]]:
        if output == 'p':
            return [self.peptide_id_to_peptide_sequence[peptide_id] for peptide_id in
                    self.get_peptide_ids(protein_input)]
        elif output == 'i':
            return [peptide_id for peptide_id in self.get_peptide_ids(protein_input)]
        else:
            raise TypeError(f'Unsupported output: {output}')

    def sort_protein_ids_by_peptide_count(self):
        return sorted(list(self.protein_ids), key=lambda i: len(self.get_peptides(i)), reverse=True)

    @property
    def peptides(self) -> Generator[str, None, None]:
        return (peptide for peptide in self.peptide_sequence_to_peptide_id.keys())

    def peptide(self, peptide_id: int) -> str:
        return self.peptide_id_to_peptide_sequence.get(peptide_id)

    @property
    def peptide_ids(self) -> Generator[int, Any, None]:
        return (peptide_id for peptide_id in self.peptide_id_to_peptide_sequence.keys())

    def peptide_id(self, peptide: str) -> int:
        return self.peptide_sequence_to_peptide_id.get(peptide)

    @property
    def protein_ids(self) -> Generator[int, Any, None]:
        return (protein_id for protein_id in self.protein_id_to_protein_sequence.keys())

    def protein_id(self, protein_input: str) -> int:
        if protein_input in self.protein_locus_to_protein_id:
            return self.protein_locus_to_protein_id.get(protein_input)
        elif protein_input in self.protein_sequence_to_protein_id:
            return self.protein_sequence_to_protein_id.get(protein_input)
        else:
            return None

    @property
    def proteins(self) -> Generator[str, None, None]:
        return (protein for protein in self.protein_sequence_to_protein_id.keys())

    def protein(self, protein_id: int) -> str:
        return self.protein_id_to_protein_sequence.get(protein_id)

    @property
    def loci(self) -> Generator[str, None, None]:
        return (locus for locus in self.protein_locus_to_protein_id.keys())

    def locus(self, protein_id: int) -> str:
        return self.protein_id_to_protein_locus.get(protein_id)

    @property
    def num_proteins(self) -> int:
        return len(self.protein_sequence_to_protein_id)

    @property
    def num_peptides(self) -> int:
        return len(self.peptide_sequence_to_peptide_id)

    def group_proteins(self):
        grouped = {}
        peptides = copy.deepcopy(self.peptide_id_to_protein_ids)
        for prot, peps in sorted(self.protein_id_to_peptide_ids.items(), key=lambda x: -len(x[1])):
            if not grouped:
                grouped[prot] = peps
                continue

            matches = set.intersection(*[self.peptide_id_to_protein_ids[p] for p in peps])
            matches = [m for m in matches if m in grouped.keys()]

            # If the entry is unique:
            if not matches:
                grouped[prot] = peps
                continue

            # Create new entries from subsets:
            for match in matches:
                new_prot = set()

                if isinstance(match, Set):
                    new_prot.update(match)
                elif isinstance(match, List):
                    new_prot.update(set(match))
                elif isinstance(match, int):
                    new_prot.add(match)

                new_prot.add(prot)
                new_prot = frozenset(new_prot)

                # Update grouped proteins:
                grouped[new_prot] = grouped.pop(match)

                # Update peptides:
                for pep in grouped[new_prot]:
                    peptides[pep].remove(match)
                    if prot in peptides[pep]:
                        peptides[pep].remove(prot)

                    peptides[pep].add(new_prot)

        return grouped, peptides


def map_peptide_protein_relationship(protein_sequences: List[str],
                                     protein_locuses: List[str],
                                     enzyme_regexes: Tuple[List[Tuple[str, int]], List[Tuple[str, int]]],
                                     missed_cleavages: int, min_len: int, max_len: int, non_enzymatic: bool,
                                     semi_enzymatic: bool) -> PeptideProteinRelationship:
    r = PeptideProteinRelationship()
    for protein, locus in zip(protein_sequences, protein_locuses):
        peptide_type_pairs = digest_protein(protein, enzyme_regexes, missed_cleavages, min_len, max_len, non_enzymatic,
                                            semi_enzymatic)
        r.add(protein, locus, [pp[0] for pp in peptide_type_pairs])
    return r


def pd_digest(protein_sequences: List[str],
              protein_locuses: List[str],
              enzyme_regexes: Tuple[List[Tuple[str, int]], List[Tuple[str, int]]],
              missed_cleavages: int, min_len: int, max_len: int, non_enzymatic: bool,
              semi_enzymatic: bool) -> [pd.DataFrame, pd.DataFrame]:
    protein_data = {'protein_sequence': [], 'protein_locus': [], 'protein_id': [], 'peptide_ids': []}

    peptide_to_id = {}
    peptide_id_to_protein_ids = {}
    protein_id = 0
    for protein, locus in zip(protein_sequences, protein_locuses):
        peptide_type_pairs = digest_protein(protein, enzyme_regexes, missed_cleavages, min_len, max_len, non_enzymatic,
                                            semi_enzymatic)

        peptides = set([pp[0] for pp in peptide_type_pairs])
        peptide_ids = set()
        for peptide in peptides:
            if peptide not in peptide_to_id:
                peptide_id = len(peptides)
                peptide_to_id[peptide] = peptide_id
            peptide_id = peptide_to_id[peptide]
            peptide_ids.add(peptide_id)

            peptide_id_to_protein_ids.setdefault(peptide_id, set()).add(protein_id)

        protein_data['protein_sequence'].append(protein)
        protein_data['protein_locus'].append(locus)
        protein_data['protein_id'].append(protein_id)
        protein_data['peptide_ids'].append(peptide_ids)

        protein_id += 1

    peptide_data = {'peptide_sequence': [], 'peptide_id': [], 'protein_ids': []}
    peptide_id_to_peptide = {v: k for k, v in peptide_to_id.items()}
    for peptide_id, protein_ids in peptide_id_to_protein_ids.items():
        peptide_data['peptide_sequence'].append(peptide_id_to_peptide[peptide_id])
        peptide_data['peptide_id'].append(peptide_id)
        peptide_data['protein_ids'].append(protein_ids)

    return pd.DataFrame(protein_data), pd.DataFrame(peptide_data)


def get_enzyme_specificity(sequence: str, enzyme_regexes: Tuple[List[Tuple[str, int]], List[Tuple[str, int]]]) -> int:
    unmod_sequence = strip_modifications(sequence).replace(".", "")
    sites = combine_site_regexes(protein_sequence=unmod_sequence, positive_site_info=enzyme_regexes[0],
                                 negative_site_info=enzyme_regexes[1])

    specificity = 0
    if 1 in sites:
        specificity += 1
    elif len(unmod_sequence) - 1 in sites:
        specificity += 1

    return specificity
