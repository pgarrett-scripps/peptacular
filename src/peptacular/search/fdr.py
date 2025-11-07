from collections.abc import Sequence


def calculate_fdr(
    scores: Sequence[float], target_flags: Sequence[bool], reverse_score: bool = False
) -> list[float]:
    """
    Calculate FDR values given target and decoy scores.

    Args:
        scores: Sequence of scores for each entry
        target_flags: Boolean flags indicating if each entry is a target (True) or decoy (False)
        reverse_score: If True, lower scores are better; if False, higher scores are better

    Returns:
        List of FDR values corresponding to each entry in the original order
    """
    # Create list of (score, is_target, original_index) tuples
    indexed_data = list(zip(scores, target_flags, range(len(scores))))

    # Sort by score (reverse if higher scores are better)
    indexed_data.sort(key=lambda x: x[0], reverse=not reverse_score)

    # Calculate cumulative FDR
    cum_targets = 0
    cum_decoys = 0
    fdrs: list[float] = []

    for _, is_target, _ in indexed_data:
        if is_target:
            cum_targets += 1
        else:
            cum_decoys += 1

        fdr = (cum_decoys + 1) / cum_targets if cum_targets > 0 else 1.0
        fdrs.append(fdr)

    # Convert to q-values (minimum FDR from this point forward)
    qvalues_sorted: list[float] = []
    min_fdr = 1.0
    for fdr in reversed(fdrs):
        min_fdr = min(min_fdr, fdr)
        qvalues_sorted.append(min_fdr)
    qvalues_sorted.reverse()

    # Map q-values back to original order
    qvalues_with_index = [
        (orig_idx, qval) for (_, _, orig_idx), qval in zip(indexed_data, qvalues_sorted)
    ]
    qvalues_with_index.sort(key=lambda x: x[0])  # Sort by original index

    return [qval for _, qval in qvalues_with_index]
