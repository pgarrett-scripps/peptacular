from dataclasses import dataclass
import math


@dataclass
class IsotopicPatternScore:
    combined_score: float
    offset: int
    bhattacharyya: float
    cosine_similarity: float
    ratio_score: float
    coverage: float
    missed_penalty: float
    normalized_theo: list[float]

    def __str__(self) -> str:
        # simple
        return self.__repr__()

    def __repr__(self) -> str:
        return (
            f"IsoScore(comb={self.combined_score:.4f}, offset={self.offset}, "
            f"bhatt={self.bhattacharyya:.4f}, cosine={self.cosine_similarity:.4f}, "
            f"ratio={self.ratio_score:.4f}, cov={self.coverage:.4f}, "
            f"missed={self.missed_penalty:.4f})"
        )


def score_isotopic_pattern(
    observed: list[float],
    theoretical: list[float],
    min_intensity: float = 0.0,
    offset_range: int = 1,
    offset_zero_bonus: float = 0.02,  # Prefer offset=0 when scores are close
    min_score_threshold: float = 0.1,  # Don't accept terrible matches
) -> IsotopicPatternScore:
    """
    Score isotopic pattern match between observed and theoretical patterns.

    Args:
        observed: List of observed intensities
        theoretical: List of theoretical abundances
        min_intensity: Minimum intensity threshold for detection
        offset_range: Check offsets from -offset_range to +offset_range
        offset_zero_bonus: Score bonus for offset=0 to prefer perfect alignment
        min_score_threshold: Return offset=0 if all scores below this threshold

    Returns:
        IsotopicPatternScore containing all scoring metrics

    Notes:
        - offset = 0: obs[0] aligns with theo[0] (perfect alignment)
        - offset < 0: pad observed with |offset| zeros at start (we missed first |offset| isotopes)
        - offset > 0: skip first offset observed peaks as noise
    """
    if not observed or not theoretical:
        return IsotopicPatternScore(0.0, 0, 0.0, 0.0, 0.0, 0.0, 1.0, [])

    best_score = -1.0
    best_offset = 0
    best_bhattacharyya = 0.0
    best_cosine = 0.0
    best_ratio = 0.0
    best_coverage = 0.0
    best_missed_penalty = 1.0
    best_normalized_theo = []

    # Store all results for comparison
    all_results: dict[
        int, tuple[float, float, float, float, float, float, list[float]]
    ] = {}

    # Try different offsets
    for offset in range(-offset_range, offset_range + 1):
        # Keep theoretical pattern intact, adjust observed based on offset
        theo_aligned = theoretical[:]

        if offset == 0:
            # Perfect alignment
            obs_aligned = observed[:]
        elif offset < 0:
            # Pad observed with zeros at the start (we missed first |offset| isotopes)
            obs_aligned = [0.0] * abs(offset) + observed
        else:  # offset > 0
            # Skip first offset observed peaks (treat as noise)
            obs_aligned = observed[offset:]

        # Make arrays same length by padding with zeros at the end
        max_len = max(len(obs_aligned), len(theo_aligned))
        obs_aligned = obs_aligned + [0.0] * (max_len - len(obs_aligned))
        theo_aligned = theo_aligned + [0.0] * (max_len - len(theo_aligned))

        # Scale theoretical to observed max
        max_obs = max(obs_aligned) if obs_aligned else 0.0
        max_theo = max(theo_aligned) if theo_aligned else 0.0
        scaled_theo = [
            t * (max_obs / max_theo) if max_theo > 0 else t for t in theo_aligned
        ]

        # Create masks
        theo_detectable_mask = [t >= min_intensity for t in scaled_theo]
        obs_detectable_mask = [o > 0.0 for o in obs_aligned]
        detectable_mask = [
            t or o for t, o in zip(theo_detectable_mask, obs_detectable_mask)
        ]

        # Apply mask
        scaled_theo = [t if m else 0.0 for t, m in zip(scaled_theo, detectable_mask)]

        # Normalize both to sum to 1
        sum_obs = sum(obs_aligned)
        sum_theo = sum(scaled_theo)
        obs_normalized = [o / sum_obs if sum_obs > 0 else o for o in obs_aligned]
        theo_normalized = [t / sum_theo if sum_theo > 0 else t for t in scaled_theo]

        # Calculate metrics

        # 1. Bhattacharyya coefficient
        bhattacharyya = sum(
            math.sqrt(o * t) for o, t in zip(obs_normalized, theo_normalized)
        )

        # 2. Cosine similarity
        obs_norm = math.sqrt(sum(o**2 for o in obs_normalized))
        theo_norm = math.sqrt(sum(t**2 for t in theo_normalized))
        dot_product = sum(o * t for o, t in zip(obs_normalized, theo_normalized))
        cosine_sim = (
            dot_product / (obs_norm * theo_norm)
            if obs_norm > 0 and theo_norm > 0
            else 0.0
        )

        # 3. Ratio score
        min_sum = sum(min(o, t) for o, t in zip(obs_normalized, theo_normalized))
        max_sum = sum(max(o, t) for o, t in zip(obs_normalized, theo_normalized))
        ratio_score = min_sum / max_sum if max_sum > 0 else 0.0

        # 4. Coverage
        num_expected = sum(theo_detectable_mask)
        num_detected = sum(
            1 for o, t in zip(obs_detectable_mask, theo_detectable_mask) if o and t
        )
        coverage = num_detected / num_expected if num_expected > 0 else 0.0

        # 5. Missed penalty
        missed_intensity = sum(
            t
            for t, td, od in zip(scaled_theo, theo_detectable_mask, obs_detectable_mask)
            if td and not od
        )
        total_theo_intensity = sum(
            t for t, td in zip(scaled_theo, theo_detectable_mask) if td
        )
        missed_penalty = (
            missed_intensity / total_theo_intensity if total_theo_intensity > 0 else 1.0
        )

        # Combined score
        combined = (
            (0.4 * bhattacharyya)
            + (0.3 * cosine_sim)
            + (0.2 * ratio_score)
            + (0.1 * coverage)
        )
        combined *= 1.0 - missed_penalty

        # Apply offset=0 bonus to prefer perfect alignment when scores are close
        if offset == 0:
            combined += offset_zero_bonus

        # Store results
        all_results[offset] = (
            combined,
            bhattacharyya,
            cosine_sim,
            ratio_score,
            coverage,
            missed_penalty,
            theo_normalized,
        )

        # Track best result
        if combined > best_score:
            best_score = combined
            best_offset = offset
            best_bhattacharyya = bhattacharyya
            best_cosine = cosine_sim
            best_ratio = ratio_score
            best_coverage = coverage
            best_missed_penalty = missed_penalty
            best_normalized_theo = theo_normalized

    # If all scores are terrible, default to offset=0
    if best_score < min_score_threshold:
        best_offset = 0
        if 0 in all_results:
            (
                best_score,
                best_bhattacharyya,
                best_cosine,
                best_ratio,
                best_coverage,
                best_missed_penalty,
                best_normalized_theo,
            ) = all_results[0]
            # Remove the bonus for reporting
            best_score -= offset_zero_bonus

    # Remove offset bonus from final score for reporting
    if best_offset == 0:
        best_score -= offset_zero_bonus
        best_score = max(0.0, best_score)  # Don't go negative

    return IsotopicPatternScore(
        combined_score=float(best_score),
        offset=int(best_offset),
        bhattacharyya=float(best_bhattacharyya),
        cosine_similarity=float(best_cosine),
        ratio_score=float(best_ratio),
        coverage=float(best_coverage),
        missed_penalty=float(best_missed_penalty),
        normalized_theo=best_normalized_theo,
    )
