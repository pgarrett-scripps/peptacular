



from typing import Generator, Iterable, List, Literal, Optional, Tuple, Union

from peptacular.constants import PROTEASES_COMPILED
from peptacular.util import get_regex_match_indices

from ..spans import build_left_semi_spans, build_non_enzymatic_spans, build_right_semi_spans, build_spans

from ..proforma_dataclasses import Span
from .annot_fragmentation import ProFormaAnnotationFragmentation


DigestReturnType = Literal["str", "annotation", "span", "str-span", "annotation-span"]

DIGEST_RETURN_TYPING = Union[
    Generator[str, None, None],
    Generator['ProFormaAnnotationDigestion', None, None],
    Generator[Span, None, None],
    Generator[Tuple[str, Span], None, None],
    Generator[Tuple['ProFormaAnnotationDigestion', Span], None, None],
]

class ProFormaAnnotationDigestion(ProFormaAnnotationFragmentation):
    """
    Fragmentation Methods
    """

    def _return_digested_sequences(
        self, spans: Iterable[Span], return_type: DigestReturnType
    ) -> DIGEST_RETURN_TYPING:
        """
        Helper function to return the digested sequences as strings or ProFormaAnnotations.
        """

        if return_type == "span":
            return (span for span in spans)

        if return_type == "annotation":
            return (self.slice(span[0], span[1], inplace=False) for span in spans)

        if return_type == "str":
            return (
                self.slice(span[0], span[1], inplace=False).serialize()
                for span in spans
            )

        if return_type == "str-span":
            return (
                (self.slice(span[0], span[1], inplace=False).serialize(), span)
                for span in spans
            )

        if return_type == "annotation-span":
            return (
                (self.slice(span[0], span[1], inplace=False), span) for span in spans
            )

        raise ValueError(f"Unsupported return type: {return_type}")

    def get_left_semi_enzymatic_sequences(
        self,
        min_len: Optional[int] = None,
        max_len: Optional[int] = None,
        return_type: DigestReturnType = "str",
    ) -> DIGEST_RETURN_TYPING:
        span = (0, len(self), 0)
        spans = build_left_semi_spans(span=span, min_len=min_len, max_len=max_len)
        return self._return_digested_sequences(spans, return_type)


    def get_right_semi_enzymatic_sequences(
        self,
        min_len: Optional[int] = None,
        max_len: Optional[int] = None,
        return_type: DigestReturnType = "str",
    ) -> DIGEST_RETURN_TYPING:
        s = (0, len(self), 0)
        spans = build_right_semi_spans(span=s, min_len=min_len, max_len=max_len)

        return self._return_digested_sequences(spans, return_type)


    def get_semi_enzymatic_sequences(
        self,
        min_len: Optional[int] = None,
        max_len: Optional[int] = None,
        return_type: DigestReturnType = "str",
    ) -> DIGEST_RETURN_TYPING:
        yield from self.get_left_semi_enzymatic_sequences(
            min_len=min_len, max_len=max_len, return_type=return_type
        )
        yield from self.get_right_semi_enzymatic_sequences(
            min_len=min_len, max_len=max_len, return_type=return_type
        )


    def get_non_enzymatic_sequences(
        self,
        min_len: Optional[int] = None,
        max_len: Optional[int] = None,
        return_type: DigestReturnType = "str",
    ) -> DIGEST_RETURN_TYPING:
        s = (0, len(self), 0)
        spans = build_non_enzymatic_spans(span=s, min_len=min_len, max_len=max_len)
        return self._return_digested_sequences(spans, return_type)


    def get_cleavage_sites(
        self, enzyme_regex: str
    ) -> Generator[int, None, None]:
        enzyme_regex = PROTEASES_COMPILED.get(enzyme_regex, enzyme_regex)
        return get_regex_match_indices(
            input_str=self.stripped_sequence, regex_str=enzyme_regex
        )


    def digest(
        self,
        enzyme_regex: Union[List[str], str],
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: Optional[int] = None,
        max_len: Optional[int] = None,
        complete_digestion: bool = True,
        return_type: DigestReturnType = "str",
        sort_output: bool = True,
    ) -> DIGEST_RETURN_TYPING:

        if isinstance(enzyme_regex, str):
            enzyme_regex = [enzyme_regex]

        all_spans = set()
        if not complete_digestion:
            all_spans.add((0, len(self), 0))

        cleavage_sites = []
        for regex in enzyme_regex:
            cleavage_sites.extend(
                list(self.get_cleavage_sites(enzyme_regex=regex))
            )

        spans = build_spans(
            max_index=len(self),
            enzyme_sites=cleavage_sites,
            missed_cleavages=missed_cleavages,
            min_len=min_len,
            max_len=max_len,
            semi=semi,
        )
        all_spans.update(spans)

        if sort_output:
            all_spans = sorted(all_spans, key=lambda x: (x[0], x[1], x[2]))

        return self._return_digested_sequences(all_spans, return_type)
    

    def sequential_digest(
    self,
    enzyme_configs: List[EnzymeConfig],
    min_len: Optional[int] = None,
    max_len: Optional[int] = None,
    return_type: DigestReturnType = "str",
) -> DIGEST_RETURN_TYPING:

    # Start with the whole sequence
    digested_anot_spans = []

    # Apply each enzyme sequentially
    for i, enzyme_config in enumerate(enzyme_configs):
        if i == 0:
            # First enzyme digests the original sequence
            digested_anot_spans = list(
                digest(
                    sequence=annotation,
                    enzyme_regex=enzyme_config.regex,
                    missed_cleavages=enzyme_config.missed_cleavages,
                    semi=enzyme_config.semi_enzymatic,
                    min_len=min_len,
                    max_len=None,
                    return_type="annotation-span",
                    complete_digestion=enzyme_config.complete_digestion,
                )
            )
        else:
            # Break early if previous enzyme digestion yielded nothing
            if len(digested_anot_spans) == 0:
                break

            sequential_digested_anot_spans = []

            # Apply this enzyme to each fragment from previous digestion
            for anot, span in digested_anot_spans:
                _digested_anot_spans = list(
                    digest(
                        anot,
                        enzyme_regex=enzyme_config.regex,
                        missed_cleavages=enzyme_config.missed_cleavages,
                        semi=enzyme_config.semi_enzymatic,
                        min_len=min_len,
                        max_len=None,
                        return_type="annotation-span",
                        complete_digestion=enzyme_config.complete_digestion,
                    )
                )

                # Fix span to be in reference to original sequence
                for j, digested_anot_span in enumerate(_digested_anot_spans):
                    digested_span = digested_anot_span[1]
                    fixed_digested_span = (
                        span[0] + digested_span[0],
                        span[0] + digested_span[1],
                        span[2],
                    )
                    _digested_anot_spans[j] = (
                        digested_anot_span[0],
                        fixed_digested_span,
                    )

                sequential_digested_anot_spans.extend(_digested_anot_spans)

            digested_anot_spans = sequential_digested_anot_spans

    # Apply max_len filter if specified
    if max_len is not None:
        digested_anot_spans = [
            (anot, span)
            for anot, span in digested_anot_spans
            if span[1] - span[0] <= max_len
        ]

    # Format and return the results
    return _return_digested_sequences(
        annotation, [span for anot, span in digested_anot_spans], return_type
    )