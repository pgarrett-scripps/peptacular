from __future__ import annotations

from typing import Generator, Iterable, Self, Sequence, Union, overload, Literal

from ..proforma.dclasses import SPAN_TYPE
from .types import EnzymeConfig, DigestGenerator, DigestProtocol, ReturnTypeLiteral
from .core import (
    digest_annotation,
    sequential_digest_annotation,
    get_left_semi_enzymatic_sequences,
    get_right_semi_enzymatic_sequences,
    get_semi_enzymatic_sequences,
    get_non_enzymatic_sequences,
    get_cleavage_sites,
)


class DigestionMixin(DigestProtocol):
    """
    Mixin to add digestion methods to ProFormaAnnotation.
    Preserves all type overloads from the standalone functions.
    """

    # Overloads for get_left_semi_enzymatic_sequences
    @overload
    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['span'],
    ) -> Generator[SPAN_TYPE, None, None]: ...

    @overload
    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['annotation'],
    ) -> Generator[Self, None, None]: ...

    @overload
    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['str'],
    ) -> Generator[str, None, None]: ...

    @overload
    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['str-span'],
    ) -> Generator[tuple[str, SPAN_TYPE], None, None]: ...

    @overload
    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['annotation-span'],
    ) -> Generator[tuple[Self, SPAN_TYPE], None, None]: ...

    def get_left_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        return_type: ReturnTypeLiteral = 'annotation',
    ) -> DigestGenerator:
        return get_left_semi_enzymatic_sequences(self, min_len, max_len, return_type)

    # Overloads for get_right_semi_enzymatic_sequences
    @overload
    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['span'],
    ) -> Generator[SPAN_TYPE, None, None]: ...

    @overload
    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['annotation'],
    ) -> Generator[Self, None, None]: ...

    @overload
    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['str'],
    ) -> Generator[str, None, None]: ...

    @overload
    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['str-span'],
    ) -> Generator[tuple[str, SPAN_TYPE], None, None]: ...

    @overload
    def get_right_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['annotation-span'],
    ) -> Generator[tuple[Self, SPAN_TYPE], None, None]: ...

    def get_right_semi_enzymatic_sequences(
        self: DigestProtocol,
        min_len: int | None = None,
        max_len: int | None = None,
        return_type: ReturnTypeLiteral = 'annotation',
    ) -> DigestGenerator:
        return get_right_semi_enzymatic_sequences(self, min_len, max_len, return_type)

    # Overloads for get_semi_enzymatic_sequences
    @overload
    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['span'],
    ) -> Generator[SPAN_TYPE, None, None]: ...

    @overload
    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['annotation'],
    ) -> Generator[Self, None, None]: ...

    @overload
    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['str'],
    ) -> Generator[str, None, None]: ...

    @overload
    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['str-span'],
    ) -> Generator[tuple[str, SPAN_TYPE], None, None]: ...

    @overload
    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['annotation-span'],
    ) -> Generator[tuple[Self, SPAN_TYPE], None, None]: ...

    def get_semi_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        return_type: ReturnTypeLiteral = 'annotation',
    ) -> DigestGenerator:
        return get_semi_enzymatic_sequences(self, min_len, max_len, return_type)

    # Overloads for get_non_enzymatic_sequences
    @overload
    def get_non_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['span'],
    ) -> Generator[SPAN_TYPE, None, None]: ...

    @overload
    def get_non_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['annotation'],
    ) -> Generator[Self, None, None]: ...

    @overload
    def get_non_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['str'],
    ) -> Generator[str, None, None]: ...

    @overload
    def get_non_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['str-span'],
    ) -> Generator[tuple[str, SPAN_TYPE], None, None]: ...

    @overload
    def get_non_enzymatic_sequences(
        self,
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['annotation-span'],
    ) -> Generator[tuple[Self, SPAN_TYPE], None, None]: ...

    def get_non_enzymatic_sequences(
        self,
        min_len: int | None = None,
        max_len: int | None = None,
        return_type: ReturnTypeLiteral = 'annotation',
    ) -> DigestGenerator:
        return get_non_enzymatic_sequences(self, min_len, max_len, return_type)

    def get_cleavage_sites(
        self,
        enzyme_regex: str
    ) -> Generator[int, None, None]:
        """Get cleavage sites for a given enzyme regex."""
        return get_cleavage_sites(self, enzyme_regex)

    @overload
    def digest(
        self,  # Remove type annotation - it's inferred
        enzyme_regex: Union[Iterable[str], str],
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: Union[int, None] = None,
        max_len: Union[int, None] = None,
        complete_digestion: bool = True,
        sort_output: bool = True,
        *,
        return_type: Literal["span"],
    ) -> Generator[SPAN_TYPE, None, None]: ...

    @overload
    def digest(
        self,
        enzyme_regex: Union[Iterable[str], str],
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: Union[int, None] = None,
        max_len: Union[int, None] = None,
        complete_digestion: bool = True,
        sort_output: bool = True,
        *,
        return_type: Literal["annotation"],
    ) -> Generator[Self, None, None]: ...

    @overload
    def digest(
        self,
        enzyme_regex: Union[Iterable[str], str],
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: Union[int, None] = None,
        max_len: Union[int, None] = None,
        complete_digestion: bool = True,
        sort_output: bool = True,
        *,
        return_type: Literal["str"],
    ) -> Generator[str, None, None]: ...

    @overload
    def digest(
        self,
        enzyme_regex: Union[Iterable[str], str],
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: Union[int, None] = None,
        max_len: Union[int, None] = None,
        complete_digestion: bool = True,
        sort_output: bool = True,
        *,
        return_type: Literal["str-span"],
    ) -> Generator[tuple[str, SPAN_TYPE], None, None]: ...

    @overload
    def digest(
        self,
        enzyme_regex: Union[Iterable[str], str],
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: Union[int, None] = None,
        max_len: Union[int, None] = None,
        complete_digestion: bool = True,
        sort_output: bool = True,
        return_type: Literal['annotation-span'] = 'annotation-span',
    ) -> Generator[tuple[Self, SPAN_TYPE], None, None]: ...

    def digest(
        self,
        enzyme_regex: Union[Iterable[str], str],
        missed_cleavages: int = 0,
        semi: bool = False,
        min_len: Union[int, None] = None,
        max_len: Union[int, None] = None,
        complete_digestion: bool = True,
        sort_output: bool = True,
        return_type: ReturnTypeLiteral = "annotation",
    ) -> DigestGenerator:
        """Digest this annotation with specified enzymes and parameters."""
        return digest_annotation(
            annotation=self,
            enzyme_regex=enzyme_regex,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            complete_digestion=complete_digestion,
            sort_output=sort_output,
            return_type=return_type
        )

    # Overloads for sequential_digest_annotation
    @overload
    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['span'],
    ) -> Generator[SPAN_TYPE, None, None]: ...

    @overload
    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['annotation'],
    ) -> Generator[Self, None, None]: ...

    @overload
    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['str'],
    ) -> Generator[str, None, None]: ...

    @overload
    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['str-span'],
    ) -> Generator[tuple[str, SPAN_TYPE], None, None]: ...

    @overload
    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None,
        max_len: int | None,
        return_type: Literal['annotation-span'],
    ) -> Generator[tuple[Self, SPAN_TYPE], None, None]: ...

    def sequential_digest(
        self,
        enzyme_configs: Sequence[EnzymeConfig],
        min_len: int | None = None,
        max_len: int | None = None,
        return_type: ReturnTypeLiteral = 'annotation',
    ) -> DigestGenerator:
        """Perform sequential digestion with multiple enzymes."""
        return sequential_digest_annotation(self, enzyme_configs, min_len, max_len, return_type)