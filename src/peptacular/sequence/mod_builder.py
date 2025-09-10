from typing import Generator, Iterable, Mapping
from ..proforma.annotation import ProFormaAnnotation
from .util import get_annotation_input
from ..mod import Mod

MOD_BUILDER_INPUT_TYPE = Mapping[str, Iterable[str | float | int | Mod]]


def build_mods(
    sequence: ProFormaAnnotation | str,
    nterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_static: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_static: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_static: MOD_BUILDER_INPUT_TYPE | None = None,
    nterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    cterm_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    internal_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    labile_variable: MOD_BUILDER_INPUT_TYPE | None = None,
    max_variable_mods: int = 2,
    use_regex: bool = False,
    precision: int | None = None,
    include_plus: bool = False,
) -> Generator[str, None, None]:
    annotation = get_annotation_input(sequence, copy=False)
    for annot in annotation.build_mods(
        nterm_static=nterm_static,
        cterm_static=cterm_static,
        internal_static=internal_static,
        labile_static=labile_static,
        nterm_variable=nterm_variable,
        cterm_variable=cterm_variable,
        internal_variable=internal_variable,
        labile_variable=labile_variable,
        max_variable_mods=max_variable_mods,
        use_regex=use_regex,
    ):
        yield annot.serialize(precision=precision, include_plus=include_plus)
