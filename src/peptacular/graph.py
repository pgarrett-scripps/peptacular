import random
from collections import Counter
from typing import List, Dict, Iterable

from peptacular.sequence import build_kmers, split_sequence, calculate_sequence_length
from peptacular.term.modification import get_n_term_modification, get_c_term_modification


class DebruijnGraph:

    def __init__(self, k: int, static_residues: List[str] = None, sequences: List[str] = None, seed: int = None):
        self.graph = {}
        self.k = k
        self.static_residues = [] if static_residues is None else static_residues
        self.seed = seed

        self.aa_counter = Counter()

        if sequences:
            for sequence in sequences:
                self.add_sequence(sequence)

    def add_sequence(self, sequence: str) -> None:
        update_debruijn_graph(self.graph, sequence, self.k)
        components = split_sequence(sequence)
        self.aa_counter.update(components)

    def reset(self) -> None:
        self.graph = {}

    def randomize_graph(self) -> None:
        aa_choices = list(self.aa_counter.keys())
        aa_weights = list(self.aa_counter.values())
        randomize_debruijn_graph(self.graph, aa_choices, aa_weights, self.static_residues)

    def construct_sequence(self, sequence: str) -> str:
        return construct_sequence_with_debruijn_graph(sequence, self.graph)


class MarkovChain:

    def __init__(self, k: int, sequences: List[str] = None, seed: int = None):
        self.chain = {}
        self.k = k
        self.seed = seed

        self.aa_counter = Counter()

        if sequences:
            for sequence in sequences:
                self.add_sequence(sequence)

    def add_sequence(self, sequence: str) -> None:
        update_markov_chain(self.chain, sequence, self.k)
        components = split_sequence(sequence)
        self.aa_counter.update(components)

    def reset(self) -> None:
        self.chain = {}

    def construct_sequence(self, max_length: int = None, start_key: str = None) -> str:
        return construct_sequence_with_markov_chain(self.chain, self.k, max_length, start_key, self.seed)


def update_markov_chain(chain: Dict[str, Counter[str]], sequence: str, k: int) -> None:
    """
    Adds a sequence to a Markov chain.

    :param chain: The Markov chain to add the sequence to.
    :type chain: dict
    :param sequence: The sequence to add to the chain.
    :type sequence: str
    :param k: The length of the kmers.
    :type k: int

    :return: None

    .. code-block:: python

        >>> chain = {'--': Counter({'P': 1}), '-P': Counter({'E': 1}), 'PE': Counter({'P': 1})}
        >>> update_markov_chain(chain, "PEP", 2)
        >>> chain
        {'--': Counter({'P': 2}), '-P': Counter({'E': 2}), 'PE': Counter({'P': 2}), 'EP': Counter({'#': 1})}

    """

    if get_n_term_modification(sequence) is not None or get_c_term_modification(sequence) is not None:
        raise ValueError("N or C term modification not supported")

    sequence = '-' * k + sequence + '#'
    for kmer in build_kmers(sequence, k + 1):
        components = split_sequence(kmer)
        kmer, aa = ''.join(components[:-1]), components[-1]
        if kmer not in chain:
            chain[kmer] = Counter()
        chain[kmer][aa] += 1


def construct_markov_chain(sequences: List[str], k: int) -> Dict[str, Counter[str]]:
    """
    Constructs a Markov chain from a list of sequences.

    :param sequences: A list of sequences to construct the chain from.
    :type sequences: list
    :param k: The length of the kmers.
    :type k: int

    :return: A dictionary representing the Markov chain with each k-1-mer mapped to a set of succeeding characters.

    .. code-block:: python

        >>> construct_markov_chain(["PEP"], 2)
        {'--': Counter({'P': 1}), '-P': Counter({'E': 1}), 'PE': Counter({'P': 1}), 'EP': Counter({'#': 1})}

    """
    chain = {}
    for sequence in sequences:
        update_markov_chain(chain, sequence, k)
    return chain


def update_debruijn_graph(graph: Dict[str, Dict[str, str]], sequence: str, k: int) -> None:
    """
    Adds a sequence to a De Bruijn graph.

    :param graph: The De Bruijn graph to add the sequence to.
    :type graph: dict
    :param sequence: The sequence to add to the graph.
    :type sequence: str
    :param k: The length of the kmers.
    :type k: int

    :return: None

    .. code-block:: python

        >>> graph = {'--': {'P': 'P'}, '-P': {'E': 'E'}, 'PE': {'P': 'P'}}
        >>> update_debruijn_graph(graph, "PEP", 2)
        >>> graph
        {'--': {'P': 'P'}, '-P': {'E': 'E'}, 'PE': {'P': 'P'}}

    """

    if get_n_term_modification(sequence) is not None or get_c_term_modification(sequence) is not None:
        raise ValueError("N or C term modification not supported")

    sequence = '-' * k + sequence
    for kmer in build_kmers(sequence, k + 1):
        comps = split_sequence(kmer)
        kmer, aa = ''.join(comps[:-1]), comps[-1]
        if kmer not in graph:
            graph[kmer] = {}
        graph[kmer][aa] = aa  # Add the last character of the suffix


def construct_bruijn_graph(sequences: List[str], k: int) -> Dict[str, Dict[str, str]]:
    """
    Constructs a traditional De Bruijn graph from a list of sequences.

    :param sequences: A list of sequences to construct the graph from.
    :type sequences: list
    :param k: The length of the kmers.
    :type k: int

    :raises ValueError: If the input sequences contain N or C term modifications.

    :return: A dictionary representing the De Bruijn graph with each k-1-mer mapped to a set of succeeding characters.


    .. code-block:: python

        >>> construct_bruijn_graph(["PEP"], 2)
        {'--': {'P': 'P'}, '-P': {'E': 'E'}, 'PE': {'P': 'P'}}

        >>> construct_bruijn_graph(["PE[2]P"], 2)
        {'--': {'P': 'P'}, '-P': {'E[2]': 'E[2]'}, 'PE[2]': {'P': 'P'}}

    """

    graph = {}
    for sequence in sequences:
        update_debruijn_graph(graph, sequence, k)
    return graph


def randomize_debruijn_graph(graph: Dict[str, Dict[str, str]], aa_choices: List[str], aa_weights: List[float] = None,
                             static_aa: Iterable = None, seed: int = None) -> None:
    """
    Randomizes a De Bruijn graph.

    :param graph: The De Bruijn graph to randomize.
    :type graph: dict
    :param aa_choices: A list of valid amino acids.
    :type aa_choices: list
    :param aa_weights: A list of amino acid frequencies.
    :type aa_weights: list
    :param static_aa: A list of amino acids that should not be randomized.
    :type static_aa: list
    :param seed: The seed to use for the random number generator.
    :type seed: int

    :return: None


    .. code-block:: python

        >>> graph = construct_bruijn_graph(["PEP"], 2)
        >>> randomize_debruijn_graph(graph, ['T', 'R'], seed=1)
        >>> graph
        {'--': {'P': 'T'}, '-P': {'E': 'R'}, 'PE': {'P': 'R'}}

    """
    if seed is not None:
        random.seed(seed)

    static_aa = set(static_aa) if static_aa is not None else set()

    for kmer in graph:
        for aa in graph[kmer]:
            if graph[kmer][aa] in static_aa:
                continue
            graph[kmer][aa] = random.choices(aa_choices, weights=aa_weights, k=1)[0]


def construct_sequence_with_debruijn_graph(sequence: str, graph: Dict[str, Dict[str, str]]) -> str:
    """
    Constructs a sequence from a De Bruijn graph.

    :param sequence: The input sequence.
    :type sequence: str
    :param graph: The De Bruijn graph to construct the sequence from.
    :type graph: dict

    :return: The constructed sequence.
    :rtype: str

    .. code-block:: python

        >>> graph = construct_bruijn_graph(["PEP"], 2)
        >>> construct_sequence_with_debruijn_graph("PEP", graph)
        'PEP'

        >>> graph = construct_bruijn_graph(["PEP"], 2)
        >>> randomize_debruijn_graph(graph, ['T', 'R'], seed=1)
        >>> construct_sequence_with_debruijn_graph("PEP", graph)
        'TRR'

    """

    k = calculate_sequence_length(list(graph.keys())[0])

    sequence = '-' * k + sequence
    kmers = list(build_kmers(sequence, k))

    new_sequence = ''
    for i, kmer in enumerate(kmers):
        if i == len(kmers) - 1:
            break
        next_aa = sequence[i + 2]
        new_sequence += graph[kmer][next_aa]

    return new_sequence


def construct_sequence_with_markov_chain(graph: Dict[str, Counter], k: int, max_length: int = None,
                                         start_key: str = None, seed: int = None) -> str:
    """
    Constructs a random sequence from a Markov Chain.

    :param graph: The De Bruijn graph to construct the sequence from.
    :type graph: dict
    :param k: The length of the kmers.
    :type k: int
    :param max_length: The length of the sequence to construct. If None, the sequence will be constructed until no more
    transitions are possible.
    :type max_length: int
    :param start_key: The starting key for the sequence. If None, the sequence will start with a string of '-'*k.
    :type start_key: str
    :param seed: The seed to use for the random number generator.
    :type seed: int

    :return: The constructed sequence.
    :rtype: str

    .. code-block:: python

        >>> graph = construct_markov_chain(["PEPTIDE"], 2)
        >>> construct_sequence_with_markov_chain(graph, 2)
        'PEPTIDE'

    """

    if seed is not None:
        random.seed(seed)

    if start_key is None:
        start_key = '-' * k

    if max_length is None:
        max_length = float('inf')

    if start_key not in graph:
        raise ValueError("Start key not in graph")

    new_sequence = ''

    while True:
        choices = list(graph[start_key].keys())
        weights = list(graph[start_key].values())
        next_aa = random.choices(choices, weights=weights, k=1)[0]

        if next_aa == '#' or calculate_sequence_length(new_sequence) >= max_length:
            break

        new_sequence += next_aa
        start_key = start_key[1:] + next_aa

    return new_sequence
