import pickle
import random
from collections import Counter, defaultdict
from typing import Dict, Sequence, Set, Tuple


class DeBruijnGraph:
    """
    A De Bruijn graph for analyzing and generating sequences based on k-mer statistics.

    Nodes represent (k-1)-mers, edges represent k-mers.
    """

    def __init__(self, k: int):
        self.k = k
        self.nodes: Set[str] = set()
        self.edges: Dict[Tuple[str, str], int] = defaultdict(int)
        self.kmers: Set[str] = set()
        self.amino_acid_counts: Counter = Counter()
        self.total_amino_acids: int = 0

        # Add these for performance
        self._outgoing_cache: Dict[
            str, list
        ] = {}  # node -> [(next_node, cumulative_prob), ...]
        self._cache_built: bool = False

    def _build_transition_cache(self):
        """Build a cache of cumulative probabilities for fast random selection."""
        self._outgoing_cache.clear()

        # Group edges by source node
        node_edges = defaultdict(list)
        for (from_node, to_node), count in self.edges.items():
            node_edges[from_node].append((to_node, count))

        # Convert to cumulative probabilities for binary search
        for node, edges in node_edges.items():
            total = sum(count for _, count in edges)
            cumulative = 0.0
            cumulative_probs = []
            for next_node, count in edges:
                cumulative += count / total
                cumulative_probs.append((next_node, cumulative))
            self._outgoing_cache[node] = cumulative_probs

        self._cache_built = True

    def _get_next_node_fast(self, current_node: str) -> str | None:
        """Fast weighted random selection using cached cumulative probabilities."""
        if not self._cache_built:
            self._build_transition_cache()

        transitions = self._outgoing_cache.get(current_node)
        if not transitions:
            return None

        # Random selection using cumulative probabilities
        r = random.random()
        for next_node, cumulative_prob in transitions:
            if r <= cumulative_prob:
                return next_node

        # Fallback (shouldn't happen due to floating point)
        return transitions[-1][0]

    def get_amino_acid_frequencies(self) -> Dict[str, float]:
        """Return the percentage of each amino acid observed."""
        if self.total_amino_acids == 0:
            return {}
        return {
            aa: count / self.total_amino_acids
            for aa, count in self.amino_acid_counts.items()
        }

    def get_edge_probability(self, from_node: str, to_node: str) -> float:
        """
        Get the probability of transitioning from from_node to to_node.
        Returns 0.0 if the edge doesn't exist.
        """
        total_from_node = sum(
            count for (f, t), count in self.edges.items() if f == from_node
        )
        if total_from_node == 0:
            return 0.0
        return self.edges.get((from_node, to_node), 0) / total_from_node

    def get_outgoing_edges(self, node: str) -> Dict[str, float]:
        """
        Get all outgoing edges from a node with their probabilities.
        Returns dict of {to_node: probability}
        """
        outgoing = {
            to_node: count
            for (from_node, to_node), count in self.edges.items()
            if from_node == node
        }
        total = sum(outgoing.values())
        if total == 0:
            return {}
        return {to_node: count / total for to_node, count in outgoing.items()}

    def get_kmers(self) -> Set[str]:
        """Return the set of all k-mers observed."""
        return self.kmers.copy()

    def get_random_kmer(self) -> str:
        """Get a random k-mer to start sequence generation."""
        if not self.kmers:
            raise ValueError("No k-mers available")
        return random.choice(list(self.kmers))

    def generate_sequence(
        self,
        length: int,
        start_kmer: str | None = None,
        ending_chars: str | None = None,
    ) -> str:
        """
        Generate a random sequence of the specified length using the graph's transition probabilities.

        Args:
            length: Desired length of the generated sequence
            start_kmer: Optional k-mer to start with. If None, a random k-mer is chosen.
            ending_chars: Optional characters that sequences should end with

        Returns:
            Generated sequence string
        """
        if length < self.k:
            raise ValueError(f"Length must be at least k={self.k}")

        # Build cache on first use
        if not self._cache_built:
            self._build_transition_cache()

        # Start with a k-mer
        if start_kmer is None:
            current_kmer = self.get_random_kmer()
        else:
            if start_kmer not in self.kmers:
                raise ValueError(f"Start k-mer '{start_kmer}' not found in graph")
            current_kmer = start_kmer

        # Use list for faster concatenation, convert to string at end
        sequence_parts = [current_kmer]
        current_length = len(current_kmer)
        current_node = current_kmer[1:]  # (k-1)-mer suffix

        # Extend the sequence
        while current_length < length:
            next_node = self._get_next_node_fast(current_node)

            if next_node is None:
                # Dead end - restart with a random k-mer
                current_kmer = self.get_random_kmer()
                sequence_parts.append(current_kmer)
                current_length += len(current_kmer)
                current_node = current_kmer[1:]
            else:
                # Add the last character of the next node
                new_char = next_node[-1]
                sequence_parts.append(new_char)
                current_length += 1
                current_node = next_node

        seq = "".join(sequence_parts)[:length]

        # if ending_chars is specified, ensure the sequence ends with one of them
        if ending_chars is not None and seq[-1] not in ending_chars:
            seq = seq[:-1] + random.choice(ending_chars)

        return seq

    def generate_sequence_with_ending(
        self,
        length: int,
        ending_chars: str,
        start_kmer: str | None = None,
        max_search_size: int | None = None,
    ) -> str:
        """
        Generate a sequence that ends with one of the specified amino acid characters.

        Args:
            length: Desired length of the final sequence
            ending_chars: String of amino acid characters that the sequence should end with (any one of them)
            start_kmer: Optional k-mer to start with. If None, a random k-mer is chosen.
            max_search_size: Maximum length to search for an ending character.
                            Defaults to length + 10.

        Returns:
            Generated sequence string of the specified length ending with one of the ending_chars
        """
        if length < self.k:
            raise ValueError(f"Length must be at least k={self.k}")

        if max_search_size is None:
            max_search_size = length * 5

        # Generate a longer sequence to search for ending character
        extended_sequence = self.generate_sequence(max_search_size, start_kmer)

        # Find the last occurrence of any ending character
        last_position = -1
        for i in range(len(extended_sequence) - 1, -1, -1):
            if extended_sequence[i] in ending_chars:
                last_position = i
                break

        if last_position == -1:
            # No ending character found, just use regular generation
            return self.generate_sequence(length, start_kmer)

        # Take sequence up to and including the ending character
        result = extended_sequence[: last_position + 1]

        # If too short, pad with regular generation
        if len(result) < length:
            additional_needed = length - len(result)
            additional_seq = self.generate_sequence(additional_needed)
            result = additional_seq + result

        # If too long, trim from the beginning
        if len(result) > length:
            result = result[-length:]

        return result

    def save(self, filepath: str):
        """Save the De Bruijn graph to a file."""
        with open(filepath, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def load(filepath: str) -> "DeBruijnGraph":
        """Load a De Bruijn graph from a file."""
        with open(filepath, "rb") as f:
            graph = pickle.load(f)
        return graph


def generate_debruijn_graph(sequences: Sequence[str], k: int = 3) -> DeBruijnGraph:
    """
    Generate a De Bruijn graph from the given sequences and kmer size.

    The graph stores:
    - Amino acid frequencies across all sequences
    - Transition probabilities for each edge
    - Set of all k-mers for sequence generation

    Args:
        k: Size of k-mers
        sequences: Input sequences to build the graph from

    Returns:
        DeBruijnGraph object with statistical properties of input sequences
    """
    if k < 2:
        raise ValueError("k must be at least 2")

    graph = DeBruijnGraph(k)

    for sequence in sequences:
        if len(sequence) < k:
            continue  # Skip sequences shorter than k

        # Count amino acids
        for aa in sequence:
            graph.amino_acid_counts[aa] += 1
            graph.total_amino_acids += 1

        # Extract k-mers and build graph
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i : i + k]
            graph.kmers.add(kmer)

            # Node is (k-1)-mer prefix, edge goes to (k-1)-mer suffix
            from_node = kmer[:-1]  # First k-1 characters
            to_node = kmer[1:]  # Last k-1 characters

            graph.nodes.add(from_node)
            graph.nodes.add(to_node)
            graph.edges[(from_node, to_node)] += 1

    return graph


from multiprocessing import Pool, cpu_count


def _generate_single_sequence(args):
    """Helper function for multiprocessing - generates a single sequence."""
    graph, length, start_kmer, ending_chars, seed = args
    # Set random seed for this process to ensure different sequences
    if seed is not None:
        random.seed(seed)
    return graph.generate_sequence(length, start_kmer, ending_chars)


def generate_sequences_parallel(
    graph: DeBruijnGraph,
    n_sequences: int,
    length_range: Tuple[int, int] = (6, 40),
    start_kmer: str | None = None,
    ending_chars: str | None = None,
    n_processes: int | None = None,
) -> list[str]:
    """
    Generate multiple sequences in parallel using multiprocessing.

    Args:
        graph: DeBruijnGraph object to use for generation
        n_sequences: Number of sequences to generate
        length_range: Tuple of (min_length, max_length) for random sequence lengths
        start_kmer: Optional k-mer to start with. If None, random k-mers are chosen.
        ending_chars: Optional characters that sequences should end with
        n_processes: Number of processes to use. If None, uses cpu_count().

    Returns:
        List of generated sequences
    """
    if n_processes is None:
        n_processes = cpu_count()

    # Create arguments for each sequence generation with random lengths
    # Use different seeds for each sequence to ensure randomness
    args_list = [
        (
            graph,
            random.randint(length_range[0], length_range[1]),
            start_kmer,
            ending_chars,
            random.randint(0, 2**32 - 1),
        )
        for _ in range(n_sequences)
    ]

    # Generate sequences in parallel
    with Pool(processes=n_processes) as pool:
        sequences = pool.map(_generate_single_sequence, args_list)

    return sequences
