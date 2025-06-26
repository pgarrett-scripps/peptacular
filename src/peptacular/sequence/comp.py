from typing import *

from .util import get_annotation_input, override_annotation_properties

from ..proforma.annotation import ProFormaAnnotation
from ..proforma_dataclasses import ChemComposition, Mod, ModValue


def comp(
    sequence: Union[str, ProFormaAnnotation],
    ion_type: str = "p",
    estimate_delta: bool = False,
    charge: Optional[int] = None,
    isotope: int = 0,
    charge_adducts: Optional[List[Mod]] = None,
    isotope_mods: Optional[List[ModValue]] = None,
    use_isotope_on_mods: bool = False,
) -> ChemComposition:
    """
    Calculates the elemental composition of a peptide sequence, including modifications,
    and optionally estimates the composition based on the delta mass from modifications.

    :param sequence: A sequence or ProFormaAnnotation.
    :type sequence: str | ProFormaAnnotation
    :param ion_type: The type of ion. Default is 'p'.
    :type ion_type: str
    :param estimate_delta: If True, estimate the composition based on the delta mass from modifications.
    Default is False.
    :type estimate_delta: bool
    :param charge: The charge state of the ion. Default is None.
    :type charge: int | None
    :param isotope: The number of Neutrons to add/subtract from the final mass. Default is 0.
    :type isotope: int
    :param charge_adducts: The charge adducts. Default is None.
    :type charge_adducts: int | None
    :param isotope_mods: The isotope modifications. Default is None.
    :type isotope_mods: List[Mod] | None
    :param use_isotope_on_mods: If True, use the isotope on the modifications. Default is False.
    :type use_isotope_on_mods: bool

    :raises ValueError: If delta_mass is nonzero and estimate_delta is False, indicating an unaccounted modification.

    :return: The elemental composition of the peptide sequence.
    :rtype: Dict[str, int | float]

    .. code-block:: python

        # Calculate the elemental composition of a peptide sequence.
        >>> comp('PEPTIDE')
        {'C': 34, 'H': 53, 'N': 7, 'O': 15}

        >>> comp('PEPTIDE')
        {'C': 34, 'H': 53, 'N': 7, 'O': 15}

        >>> comp('PEPTIDE[1.0]', estimate_delta=True)['C']
        34.04446833455479

        >>> comp('PEPTIDE[1.0]', estimate_delta=False)  # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        ValueError: Non-zero delta mass (1.0) encountered without estimation enabled for sequence 'PEPTIDE[1.00]'.

    """

    annotation = get_annotation_input(sequence=sequence, copy=True)

    override_annotation_properties(
        annotation=annotation,
        charge=charge,
        charge_adducts=charge_adducts,
        isotope_mods=isotope_mods,
    )

    return annotation.comp(
        ion_type=ion_type,
        estimate_delta=estimate_delta,
        isotope=isotope,
        use_isotope_on_mods=use_isotope_on_mods,
        inplace=True,
    )


def comp_mass(
    sequence: Union[str, ProFormaAnnotation],
    ion_type: str = "p",
    charge: Optional[int] = None,
    isotope: int = 0,
    charge_adducts: Optional[List[Mod]] = None,
    isotope_mods: Optional[List[ModValue]] = None,
    use_isotope_on_mods: bool = False,
) -> Tuple[ChemComposition, float]:

    annotation = get_annotation_input(sequence, copy=True)

    override_annotation_properties(
        annotation,
        charge=charge,
        charge_adducts=charge_adducts,
        isotope_mods=isotope_mods,
    )

    return annotation.comp_mass(
        ion_type=ion_type,
        isotope=isotope,
        use_isotope_on_mods=use_isotope_on_mods,
        precision=None,
    )


def condense_to_mass_mods(
    sequence: Union[str, ProFormaAnnotation],
    include_plus: bool = False,
    precision: Optional[int] = None,
) -> str:
    """
    Converts all modifications in a sequence to their mass equivalents by calculating
    the mass difference between modified and unmodified segments.

    :param sequence: The sequence or ProFormaAnnotation object to convert.
    :type sequence: Union[str, ProFormaAnnotation]
    :param include_plus: Whether to include the plus sign for positive mass modifications.
    :type include_plus: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The sequence with all modifications converted to mass modifications.
    :rtype: str

    .. code-block:: python

        >>> condense_to_mass_mods('PEP[Phospho]TIDE', include_plus=False, precision=3)
        'PEP[79.966]TIDE'

        >>> condense_to_mass_mods('PEP[Phospho]TIDE', include_plus=True, precision=3)
        'PEP[+79.966]TIDE'

        >>> condense_to_mass_mods('[Acetyl]-PEPTIDE', precision=3)
        '[42.011]-PEPTIDE'

        >>> condense_to_mass_mods('PEPTIDE-[Amidated]', precision=3)
        'PEPTIDE-[-0.984]'

        >>> condense_to_mass_mods('<13C>PEP[Phospho]TIDE', precision=3)
        'P[5.017]E[5.017]P[84.983]T[4.013]I[6.020]D[4.013]E[5.017]'

    """
    annotation = get_annotation_input(sequence=sequence, copy=True)

    return annotation.condense_to_mass_mods(
        include_plus=include_plus,
        inplace=True,
    ).serialize(include_plus=include_plus, precision=precision)
