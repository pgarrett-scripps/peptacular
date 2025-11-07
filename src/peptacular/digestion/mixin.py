from __future__ import annotations

from typing import Generator, Literal, Self, Sequence, overload

from ..proforma.dclasses import Span
from .core import (
    digest_annotation_by_aa,
    digest_annotation_by_regex,
    get_cleavage_sites,
    get_left_semi_enzymatic_sequences,
    get_non_enzymatic_sequences,
    get_right_semi_enzymatic_sequences,
    get_semi_enzymatic_sequences,
    sequential_digest_annotation,
)
from .types import (
    DigestGenerator,
    DigestProtocol,
    DigestReturnType,
    DigestReturnTypeAnnotationLiteral,
    DigestReturnTypeAnnotationSpanLiteral,
    DigestReturnTypeLiterals,
    DigestReturnTypeSpanLiteral,
    DigestReturnTypeStrLiteral,
    DigestReturnTypeStrSpanLiteral,
    EnzymeConfig,
)


class DigestionMixin(DigestProtocol):
    """
    Mixin to add digestion methods to ProFormaAnnotation.
    Preserves all type overloads from the standalone functions.
    """
    __slots__ = ()

    # Overloads for get_left_semi_enzymatic_sequences
    @overload
    def get_left_semi_enzymatic_sequences(  # type: ignore[reportOverlappingOverloads]
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: DigestReturnTypeSpanLiteral | Literal[DigestReturnType.SPAN],
    ) -> Generator[Span, None, None]: ...

    @overload
    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeAnnotationLiteral | Literal[DigestReturnType.ANNOTATION]
        ),
    ) -> Generator[Self, None, None]: ...

    @overload
    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: DigestReturnTypeStrLiteral | Literal[DigestReturnType.STR],
    ) -> Generator[str, None, None]: ...

    @overload
    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeStrSpanLiteral | Literal[DigestReturnType.STR_SPAN]
        ),
    ) -> Generator[tuple[str, Span], None, None]: ...

    @overload
    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeAnnotationSpanLiteral
            | Literal[DigestReturnType.ANNOTATION_SPAN]
        ),
    ) -> Generator[tuple[Self, Span], None, None]: ...

    @overload
    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeLiterals | DigestReturnType
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeLiterals | DigestReturnType
        ) = DigestReturnType.ANNOTATION,
    ) -> DigestGenerator:
        return get_left_semi_enzymatic_sequences(self, min_len, max_len, return_type)

    # Overloads for get_right_semi_enzymatic_sequences
    @overload
    def get_right_semi_enzymatic_sequences(  # type: ignore[reportOverlappingOverloads]
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeSpanLiteral | Literal[DigestReturnType.SPAN]
        ) = DigestReturnType.SPAN,
    ) -> Generator[Span, None, None]: ...

    @overload
    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeAnnotationLiteral | Literal[DigestReturnType.ANNOTATION]
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    @overload
    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeStrLiteral | Literal[DigestReturnType.STR]
        ) = DigestReturnType.STR,
    ) -> Generator[str, None, None]: ...

    @overload
    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeStrSpanLiteral | Literal[DigestReturnType.STR_SPAN]
        ) = DigestReturnType.STR_SPAN,
    ) -> Generator[tuple[str, Span], None, None]: ...

    @overload
    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeAnnotationSpanLiteral
            | Literal[DigestReturnType.ANNOTATION_SPAN]
        ) = DigestReturnType.ANNOTATION_SPAN,
    ) -> Generator[tuple[Self, Span], None, None]: ...

    @overload
    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeLiterals | DigestReturnType
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeLiterals | DigestReturnType
        ) = DigestReturnType.ANNOTATION,
    ) -> DigestGenerator:
        return get_right_semi_enzymatic_sequences(self, min_len, max_len, return_type)

    # Overloads for get_semi_enzymatic_sequences
    @overload
    def get_semi_enzymatic_sequences(  # type: ignore[reportOverlappingOverloads]
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeSpanLiteral | Literal[DigestReturnType.SPAN]
        ) = DigestReturnType.SPAN,
    ) -> Generator[Span, None, None]: ...

    @overload
    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeAnnotationLiteral | Literal[DigestReturnType.ANNOTATION]
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    @overload
    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeStrLiteral | Literal[DigestReturnType.STR]
        ) = DigestReturnType.STR,
    ) -> Generator[str, None, None]: ...

    @overload
    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeStrSpanLiteral | Literal[DigestReturnType.STR_SPAN]
        ) = DigestReturnType.STR_SPAN,
    ) -> Generator[tuple[str, Span], None, None]: ...

    @overload
    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeAnnotationSpanLiteral
            | Literal[DigestReturnType.ANNOTATION_SPAN]
        ) = DigestReturnType.ANNOTATION_SPAN,
    ) -> Generator[tuple[Self, Span], None, None]: ...

    @overload
    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        return_type: (
            DigestReturnTypeLiterals | DigestReturnType
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        return_type: (
            DigestReturnTypeLiterals | DigestReturnType
        ) = DigestReturnType.ANNOTATION,
    ) -> DigestGenerator:
        return get_semi_enzymatic_sequences(self, min_len, max_len, return_type)

    # Overloads for get_non_enzymatic_sequences
    @overload
    def get_non_enzymatic_sequences(  # type: ignore[reportOverlappingOverloads]
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeSpanLiteral | Literal[DigestReturnType.SPAN]
        ) = DigestReturnType.SPAN,
    ) -> Generator[Span, None, None]: ...

    @overload
    def get_non_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeAnnotationLiteral | Literal[DigestReturnType.ANNOTATION]
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    @overload
    def get_non_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeStrLiteral | Literal[DigestReturnType.STR]
        ) = DigestReturnType.STR,
    ) -> Generator[str, None, None]: ...

    @overload
    def get_non_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeStrSpanLiteral | Literal[DigestReturnType.STR_SPAN]
        ) = DigestReturnType.STR_SPAN,
    ) -> Generator[tuple[str, Span], None, None]: ...

    @overload
    def get_non_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeAnnotationSpanLiteral
            | Literal[DigestReturnType.ANNOTATION_SPAN]
        ) = DigestReturnType.ANNOTATION_SPAN,
    ) -> Generator[tuple[Self, Span], None, None]: ...

    @overload
    def get_non_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        return_type: (
            DigestReturnTypeLiterals | DigestReturnType
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    def get_non_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        return_type: (
            DigestReturnTypeLiterals | DigestReturnType
        ) = DigestReturnType.ANNOTATION,
    ) -> DigestGenerator:
        return get_non_enzymatic_sequences(self, min_len, max_len, return_type)

    def get_cleavage_sites(self, enzyme_regex: str) -> Generator[int, None, None]:
        """Get cleavage sites for a given enzyme_regex regex."""
        return get_cleavage_sites(self, enzyme_regex)

    @overload
    def regex_digest(  # type: ignore[reportOverlappingOverloads]
        self,
        enzyme_regex: str,
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        complete_digestion: bool = True,
        sort_output: bool = True,
        return_type: (
            DigestReturnTypeSpanLiteral | Literal[DigestReturnType.SPAN]
        ) = DigestReturnType.SPAN,
    ) -> Generator[Span, None, None]: ...

    @overload
    def regex_digest(
        self,
        enzyme_regex: str,
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        complete_digestion: bool = True,
        sort_output: bool = True,
        return_type: (
            DigestReturnTypeAnnotationLiteral | Literal[DigestReturnType.ANNOTATION]
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    @overload
    def regex_digest(
        self,
        enzyme_regex: str,
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        complete_digestion: bool = True,
        sort_output: bool = True,
        return_type: (
            DigestReturnTypeStrLiteral | Literal[DigestReturnType.STR]
        ) = DigestReturnType.STR,
    ) -> Generator[str, None, None]: ...

    @overload
    def regex_digest(
        self,
        enzyme_regex: str,
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        complete_digestion: bool = True,
        sort_output: bool = True,
        return_type: (
            DigestReturnTypeStrSpanLiteral | Literal[DigestReturnType.STR_SPAN]
        ) = DigestReturnType.STR_SPAN,
    ) -> Generator[tuple[str, Span], None, None]: ...

    @overload
    def regex_digest(
        self,
        enzyme_regex: str,
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        complete_digestion: bool = True,
        sort_output: bool = True,
        return_type: (
            DigestReturnTypeAnnotationSpanLiteral
            | Literal[DigestReturnType.ANNOTATION_SPAN]
        ) = DigestReturnType.ANNOTATION_SPAN,
    ) -> Generator[tuple[Self, Span], None, None]: ...

    @overload
    def regex_digest(
        self,
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
    ) -> Generator[Self, None, None]: ...

    def regex_digest(
        self,
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
        """Digest this annotation with specified enzymes and parameters."""
        return digest_annotation_by_regex(
            annotation=self,
            enzyme_regex=enzyme_regex,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            complete_digestion=complete_digestion,
            sort_output=sort_output,
            return_type=return_type,
        )

    @overload
    def digest(  # type: ignore[reportOverlappingOverloads]
        self,
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
            DigestReturnTypeSpanLiteral | Literal[DigestReturnType.SPAN]
        ) = DigestReturnType.SPAN,
    ) -> Generator[Span, None, None]: ...

    @overload
    def digest(
        self,
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
            DigestReturnTypeAnnotationLiteral | Literal[DigestReturnType.ANNOTATION]
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    @overload
    def digest(
        self,
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
            DigestReturnTypeStrLiteral | Literal[DigestReturnType.STR]
        ) = DigestReturnType.STR,
    ) -> Generator[str, None, None]: ...

    @overload
    def digest(
        self,
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
            DigestReturnTypeStrSpanLiteral | Literal[DigestReturnType.STR_SPAN]
        ) = DigestReturnType.STR_SPAN,
    ) -> Generator[tuple[str, Span], None, None]: ...

    @overload
    def digest(
        self,
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
            DigestReturnTypeAnnotationSpanLiteral
            | Literal[DigestReturnType.ANNOTATION_SPAN]
        ) = DigestReturnType.ANNOTATION_SPAN,
    ) -> Generator[tuple[Self, Span], None, None]: ...

    @overload
    def digest(
        self,
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
    ) -> Generator[Self, None, None]: ...

    def digest(
        self,
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
        """Digest this annotation with specified enzymes and parameters."""
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
            complete_digestion=complete_digestion,
            sort_output=sort_output,
            return_type=return_type,
        )

    # Overloads for sequential_digest_annotation
    @overload
    def sequential_digest(  # type: ignore[reportOverlappingOverloads]
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeSpanLiteral | Literal[DigestReturnType.SPAN]
        ) = DigestReturnType.SPAN,
    ) -> Generator[Span, None, None]: ...

    @overload
    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeAnnotationLiteral | Literal[DigestReturnType.ANNOTATION]
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    @overload
    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeStrLiteral | Literal[DigestReturnType.STR]
        ) = DigestReturnType.STR,
    ) -> Generator[str, None, None]: ...

    @overload
    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeStrSpanLiteral | Literal[DigestReturnType.STR_SPAN]
        ) = DigestReturnType.STR_SPAN,
    ) -> Generator[tuple[str, Span], None, None]: ...

    @overload
    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeAnnotationSpanLiteral
            | Literal[DigestReturnType.ANNOTATION_SPAN]
        ) = DigestReturnType.ANNOTATION_SPAN,
    ) -> Generator[tuple[Self, Span], None, None]: ...

    @overload
    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeLiterals | DigestReturnType
        ) = DigestReturnType.ANNOTATION,
    ) -> Generator[Self, None, None]: ...

    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None = None,
        max_len: int | None = None,
        *,
        return_type: (
            DigestReturnTypeLiterals | DigestReturnType
        ) = DigestReturnType.ANNOTATION,
    ) -> DigestGenerator:
        """Perform sequential digestion with multiple enzymes."""
        return sequential_digest_annotation(
            self, enzyme_configs, min_len, max_len, return_type
        )
