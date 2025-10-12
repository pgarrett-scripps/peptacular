from ..proforma.annotation import ProFormaAnnotation


def predict_rt(sequence: str | ProFormaAnnotation) -> float:
    pass


def predict_ccs(sequence: str | ProFormaAnnotation, charge: int | None = None) -> float:
    pass


def predict_msms_spectra(
    peptides: list[str | ProFormaAnnotation],
    charge: int | None = None,
) -> tuple[list[float], list[float]]:
    pass
