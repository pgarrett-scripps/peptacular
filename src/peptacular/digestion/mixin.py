from __future__ import annotations

import re
from typing import Generator, overload

from ..spans import Span
from .core import (
    digest_annotation_by_aa,
    digest_annotation_by_regex,
    get_cleavage_sites,
    generate_regex,
    left_semi_spans,
    nonspecific_spans,
    right_semi_spans,
    semi_spans,
    sequential_digest_annotation,
)
from .types import (
    DigestProtocol,
    EnzymeConfig,
)


class DigestionMixin(DigestProtocol):
    """
    Mixin to add digestion methods to ProFormaAnnotation.
    All methods return Span objects only.
    """

    __slots__ = ()

    def left_semi_spans(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Get left semi-enzymatic sequences (N-terminus fixed)."""
        return left_semi_spans(self, min_len, max_len)

    def right_semi_spans(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Get right semi-enzymatic sequences (C-terminus fixed)."""
        return right_semi_spans(self, min_len, max_len)

    def semi_spans(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Get all semi-enzymatic sequences."""
        return semi_spans(self, min_len, max_len)

    def nonspecific_spans(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Get all non-enzymatic sequences (all possible subsequences)."""
        return nonspecific_spans(self, min_len, max_len)

    def cleavage_sites(
        self,
        enzyme: str | re.Pattern[str],
    ) -> Generator[int, None, None]:
        # Call the underlying function
        return get_cleavage_sites(self, enzyme)

    def simple_cleavage_sites(
        self,
        cleave_on: str,
        restrict_before: str = "",
        restrict_after: str = "",
        cterminal: bool = True,
    ) -> Generator[int, None, None]:
        """Get cleavage sites using simple amino acid rules."""
        enzyme_regex = generate_regex(
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
        )
        return self.cleavage_sites(enzyme_regex)

    def digest(
        self,
        enzyme: str,
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Digest this annotation using a regex pattern."""
        return digest_annotation_by_regex(
            annotation=self,
            enzyme_regex=enzyme,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )

    def simple_digest(
        self,
        cleave_on: str,
        restrict_before: str = "",
        restrict_after: str = "",
        cterminal: bool = True,
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Digest this annotation with specified enzyme parameters."""
        return digest_annotation_by_aa(
            annotation=self,
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
        )

    def sequential_digest(
        self,
        enzyme_configs: list[EnzymeConfig],
        min_len: int | None = None,
        max_len: int | None = None,
    ) -> Generator[Span, None, None]:
        """Perform sequential digestion with multiple enzymes."""
        return sequential_digest_annotation(self, enzyme_configs, min_len, max_len)
