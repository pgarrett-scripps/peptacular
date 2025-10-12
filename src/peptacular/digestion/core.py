from __future__ import annotations

from typing import Generator, Iterable, Sequence

import regex as re

from ..util import get_regex_match_indices
from ..spans import (
    build_left_semi_spans,
    build_non_enzymatic_spans,
    build_right_semi_spans,
    build_spans,
)
from ..proforma.dclasses import SPAN_TYPE

from .constants import PROTEASES_COMPILED
from .types import (
    DigestReturnType,
    EnzymeConfig,
    DigestGenerator,
    DigestProtocol,
    DigestReturnTypeLiterals,
)


def _return_digested_sequences(
    annotation: DigestProtocol,
    spans: Iterable[SPAN_TYPE],
    return_type: DigestReturnTypeLiterals | DigestReturnType,
) -> DigestGenerator:
    """Helper function to return digested sequences in various formats."""
    return_type_enum = DigestReturnType(return_type)

    # Check if we can use fast path (no modifications)
    has_mods = annotation.has_mods() if hasattr(annotation, "has_mods") else True

    match return_type_enum:
        case DigestReturnType.SPAN:
            return (span for span in spans)
        case DigestReturnType.ANNOTATION:
            return (annotation.slice(span[0], span[1], inplace=False) for span in spans)
        case DigestReturnType.STR:
            if not has_mods:
                # Fast path: no mods, just slice the string directly
                seq = annotation.stripped_sequence
                return (seq[span[0] : span[1]] for span in spans)
            return (
                annotation.slice(span[0], span[1], inplace=False).serialize()
                for span in spans
            )
        case DigestReturnType.STR_SPAN:
            if not has_mods:
                # Fast path: no mods, just slice the string directly
                seq = annotation.stripped_sequence
                return ((seq[span[0] : span[1]], span) for span in spans)
            return (
                (annotation.slice(span[0], span[1], inplace=False).serialize(), span)
                for span in spans
            )
        case DigestReturnType.ANNOTATION_SPAN:
            return (
                (annotation.slice(span[0], span[1], inplace=False), span)
                for span in spans
            )
        case _:
            raise ValueError(f"Unsupported return type: {return_type}")


def get_left_semi_enzymatic_sequences(
    annotation: DigestProtocol,
    min_len: int | None = None,
    max_len: int | None = None,
    return_type: (
        DigestReturnTypeLiterals | DigestReturnType
    ) = DigestReturnType.ANNOTATION,
) -> DigestGenerator:
    """Get left semi-enzymatic sequences from a ProForma annotation."""
    span = (0, len(annotation), 0)
    spans = build_left_semi_spans(span=span, min_len=min_len, max_len=max_len)
    return _return_digested_sequences(annotation, spans, return_type)


def get_right_semi_enzymatic_sequences(
    annotation: DigestProtocol,
    min_len: int | None = None,
    max_len: int | None = None,
    return_type: (
        DigestReturnTypeLiterals | DigestReturnType
    ) = DigestReturnType.ANNOTATION,
) -> DigestGenerator:
    """Get right semi-enzymatic sequences from a ProForma annotation."""
    s = (0, len(annotation), 0)
    spans = build_right_semi_spans(span=s, min_len=min_len, max_len=max_len)
    return _return_digested_sequences(annotation, spans, return_type)


def get_semi_enzymatic_sequences(
    annotation: DigestProtocol,
    min_len: int | None = None,
    max_len: int | None = None,
    return_type: (
        DigestReturnTypeLiterals | DigestReturnType
    ) = DigestReturnType.ANNOTATION,
) -> DigestGenerator:
    """Get both left and right semi-enzymatic sequences from a ProForma annotation."""
    yield from get_left_semi_enzymatic_sequences(
        annotation, min_len=min_len, max_len=max_len, return_type=return_type
    )
    yield from get_right_semi_enzymatic_sequences(
        annotation, min_len=min_len, max_len=max_len, return_type=return_type
    )


def get_non_enzymatic_sequences(
    annotation: DigestProtocol,
    min_len: int | None = None,
    max_len: int | None = None,
    return_type: (
        DigestReturnTypeLiterals | DigestReturnType
    ) = DigestReturnType.ANNOTATION,
) -> DigestGenerator:
    """Get non-enzymatic sequences from a ProForma annotation."""
    s = (0, len(annotation), 0)
    spans = build_non_enzymatic_spans(span=s, min_len=min_len, max_len=max_len)
    return _return_digested_sequences(annotation, spans, return_type)


def get_cleavage_sites(
    annotation: DigestProtocol, enzyme_regex: str
) -> Generator[int, None, None]:
    """Get cleavage sites for a given enzyme regex."""
    enzyme_compiled_regex = PROTEASES_COMPILED.get(
        enzyme_regex, re.compile(enzyme_regex)
    )
    return get_regex_match_indices(
        input_str=annotation.stripped_sequence, regex_str=enzyme_compiled_regex
    )


def _convert_to_aa_set(aa_keys: str | None) -> set[str]:
    X = set("GASPVTCILJNDQKEMHFRYWUO")
    B = set("DN")
    J = set("IL")
    Z = set("EQ")

    def make_aa_set(aa_keys: str) -> set[str]:
        aas: set[str] = set()
        for char in aa_keys:
            if char == "X":
                aas.update(X)
            elif char == "B":
                aas.update(B)
            elif char == "J":
                aas.update(J)
            elif char == "Z":
                aas.update(Z)
            else:
                aas.add(char)

        return aas

    if aa_keys is None:
        return set()

    if "-" in aa_keys:
        if aa_keys.count("-") != 1:
            raise ValueError(
                "Amino acid keys with '-' must contain exactly one '-' character."
            )
        first_part, second_part = aa_keys.split("-")
        return make_aa_set(first_part) - make_aa_set(second_part)

    return make_aa_set(aa_keys)


def generate_regex(
    cleave_on: str | None = None,
    restrict_before: str | None = None,
    restrict_after: str | None = None,
    cterminal: bool = True,
) -> str:
    """Generate regex for proteolytic cleavage sites."""

    if cleave_on is None:
        return ""

    if not cleave_on:
        return ""

    # Escape special regex characters
    escaped_cleave = "".join(re.escape(char) for char in _convert_to_aa_set(cleave_on))

    parts: list[str] = []

    if cterminal:
        # Cleave after residue: lookbehind for cleavage residue
        parts.append(f"(?<=[{escaped_cleave}])")

        # Restrict after: what follows cleavage point
        if restrict_after:
            escaped_restrict = "".join(
                re.escape(char) for char in _convert_to_aa_set(restrict_after)
            )
            parts.append(f"(?=[^{escaped_restrict}])")

        # Restrict before: what precedes cleavage residue
        if restrict_before:
            escaped_restrict = "".join(
                re.escape(char) for char in _convert_to_aa_set(restrict_before)
            )
            # Need to check what comes before the cleavage residue
            parts.insert(0, f"(?<=[^{escaped_restrict}][{escaped_cleave}])")
            # Remove the simple lookbehind since we now have the combined one
            parts = [parts[0]] + parts[2:]
    else:
        # Cleave before residue: lookahead for cleavage residue
        parts.append(f"(?=[{escaped_cleave}])")

        # Restrict before: what precedes cleavage point
        if restrict_before:
            escaped_restrict = "".join(
                re.escape(char) for char in _convert_to_aa_set(restrict_before)
            )
            parts.insert(0, f"(?<=[^{escaped_restrict}])")

        # Restrict after: what follows cleavage residue
        if restrict_after:
            escaped_restrict = "".join(
                re.escape(char) for char in _convert_to_aa_set(restrict_after)
            )
            # Replace simple lookahead with one that checks after the residue
            parts[-1] = f"(?=[{escaped_cleave}][^{escaped_restrict}])"

    return "".join(parts)


def digest_annotation_by_aa(
    annotation: DigestProtocol,
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    complete_digestion: bool = True,
    sort_output: bool = True,
    return_type: (
        DigestReturnTypeLiterals | DigestReturnType
    ) = DigestReturnType.ANNOTATION,
) -> DigestGenerator:
    return digest_annotation_by_regex(
        annotation=annotation,
        enzyme_regex=generate_regex(
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
        ),
        missed_cleavages=missed_cleavages,
        semi=semi,
        min_len=min_len,
        max_len=max_len,
        complete_digestion=complete_digestion,
        sort_output=sort_output,
        return_type=return_type,
    )


def digest_annotation_by_regex(
    annotation: DigestProtocol,
    enzyme_regex: str,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    complete_digestion: bool = True,
    sort_output: bool = True,
    return_type: (
        DigestReturnTypeLiterals | DigestReturnType
    ) = DigestReturnType.ANNOTATION,
) -> DigestGenerator:
    all_spans: set[SPAN_TYPE] = set()
    if not complete_digestion:
        all_spans.add((0, len(annotation), 0))

    cleavage_sites: list[int] = []
    cleavage_sites.extend(
        list(get_cleavage_sites(annotation, enzyme_regex=enzyme_regex))
    )

    spans = build_spans(
        max_index=len(annotation),
        enzyme_sites=cleavage_sites,
        missed_cleavages=missed_cleavages,
        min_len=min_len,
        max_len=max_len,
        semi=semi,
    )
    all_spans.update(spans)

    spans = all_spans
    if sort_output:
        spans = sorted(spans, key=lambda x: (x[0], x[1], x[2]))

    return _return_digested_sequences(annotation, spans, return_type)


def sequential_digest_annotation(
    annotation: DigestProtocol,
    enzyme_configs: Sequence[EnzymeConfig],
    min_len: int | None = None,
    max_len: int | None = None,
    return_type: (
        DigestReturnTypeLiterals | DigestReturnType
    ) = DigestReturnType.ANNOTATION_SPAN,
) -> DigestGenerator:
    """Perform sequential digestion with multiple enzymes."""
    digested_anot_spans: list[tuple[DigestProtocol, SPAN_TYPE]] = []

    for i, enzyme_config in enumerate(enzyme_configs):
        if i == 0:
            digested_anot_spans = list(
                digest_annotation_by_regex(
                    annotation=annotation,
                    enzyme_regex=enzyme_config.enzyme_regex,
                    missed_cleavages=enzyme_config.missed_cleavages,
                    semi=enzyme_config.semi_enzymatic,
                    min_len=min_len,
                    max_len=None,
                    sort_output=True,
                    complete_digestion=enzyme_config.complete_digestion,
                    return_type=DigestReturnType.ANNOTATION_SPAN,
                )  # type: ignore (could add overloads but would really just add bloat)
            )
        else:
            if len(digested_anot_spans) == 0:
                break

            sequential_digested_anot_spans: list[tuple[DigestProtocol, SPAN_TYPE]] = []

            for anot, span in digested_anot_spans:
                _digested_anot_spans = list(
                    digest_annotation_by_regex(
                        annotation=anot,
                        enzyme_regex=enzyme_config.enzyme_regex,
                        missed_cleavages=enzyme_config.missed_cleavages,
                        semi=enzyme_config.semi_enzymatic,
                        min_len=min_len,
                        max_len=None,
                        sort_output=True,
                        complete_digestion=enzyme_config.complete_digestion,
                        return_type=DigestReturnType.ANNOTATION_SPAN,
                    )
                )

                fixed_digested_anot_spans: list[tuple[DigestProtocol, SPAN_TYPE]] = []
                for digested_anot, digested_span in _digested_anot_spans:  # type: ignore
                    fixed_digested_span = (  # type: ignore
                        span[0] + digested_span[0],  # type: ignore
                        span[0] + digested_span[1],  # type: ignore
                        span[2],  # type: ignore
                    )
                    fixed_digested_anot_spans.append(
                        (digested_anot, fixed_digested_span)
                    )  # type: ignore

                sequential_digested_anot_spans.extend(fixed_digested_anot_spans)

            digested_anot_spans = sequential_digested_anot_spans

    if max_len is not None:
        digested_anot_spans = [
            (anot, span)
            for anot, span in digested_anot_spans
            if span[1] - span[0] <= max_len
        ]

    return _return_digested_sequences(
        annotation, [span for _, span in digested_anot_spans], return_type
    )
