from typing import List


def check_parentheses(text):
    stack = []

    for i, char in enumerate(text):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if not stack:
                return False
            else:
                stack.pop()

    return len(stack) == 0

def calculate_protein_coverage(protein_sequence: str, peptide_sequences: List[str]) -> float:
    """
    Calculates the protein coverage of a list of peptides.

    Args:
        protein_sequence: The protein sequence.
        peptide_sequences: The list of peptide sequences.

    Returns:
        The protein coverage.
    """
    covered = set()
    for peptide_sequence in peptide_sequences:
        for i in range(len(protein_sequence) - len(peptide_sequence) + 1):
            if protein_sequence[i:i + len(peptide_sequence)] == peptide_sequence:
                covered.update(range(i, i + len(peptide_sequence)))

    return len(covered) / len(protein_sequence)