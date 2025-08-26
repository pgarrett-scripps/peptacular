from .util import get_annotation_input, override_annotation_properties

from ..proforma import ProFormaAnnotation
from ..constants import IonType, IonTypeLiteral


def mass(
    sequence: str | ProFormaAnnotation,
    charge: int | None = None,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
) -> float:
    """
    Calculate the mass of an amino acid 'sequence'.

    :param sequence: A sequence or ProFormaAnnotation.
    :type sequence: str | ProFormaAnnotation
    :param charge: The charge state, default is None.
    :type charge: int | None
    :param ion_type: The ion type. Default is IonType.PRECURSOR.
    :type ion_type: IonTypeLiteral | IonType
    :param monoisotopic: If True, use monoisotopic mass else use average mass. Default is True.
    :type monoisotopic: bool
    :param isotope: The number of Neutrons to add/subtract from the final mass. Default is 0.
    :type isotope: int
    :param loss: The loss to add/subtract to the final mass. Default is 0.0.
    :type loss: float
    :param charge_adducts: The charge adducts. Default is None.
    :type charge_adducts: int | None
    :param isotope_mods: The isotope modifications. Default is None.
    :type isotope_mods: List[Mod] | None
    :param use_isotope_on_mods:
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None

    :raise ValueError: If the ion type is not supported.
    :raise UnknownAminoAcidError: If an unknown amino acid is encountered.
    :raise AmbiguousAminoAcidError: If an ambiguous amino acid is encountered.

    :return: The mass of the sequence.
    :rtype: float

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> mass('PEPTIDE', precision=3)
        799.36

        >>> mass('<12C>PEPTIDE', precision=3)
        799.36

        >>> mass('<13C>PEPTIDE', precision=3)
        833.474

        >>> mass('<13C>PEPT[10]IDE', precision=3)
        843.474

        >>> mass('<13C><[Formula:[13C6]H20]@T>PEPTIDE', precision=3)
        931.651

        >>> mass('<[10]@T>PEPTIDE', precision=3)
        809.36

        >>> mass('<[10]@N-Term>PEPTIDE', precision=3)
        809.36

        >>> mass('<13C><[10]@T>PEPTIDE', precision=3)
        843.474

        >>> mass('<[Formula:[13C6]H20]@T>PEPTIDE', precision=3)
        897.537

        # Calculate the b-ion mass of a peptide sequence.
        >>> mass('PEPTIDE', ion_type='b', precision=3)
        781.349

        # Calulate the average mass of a peptide sequence.
        >>> mass('PEPTIDE', monoisotopic=False, precision=3)
        799.824

        # Calculate the mass of a peptide sequence with a charge of 2.
        >>> mass('PEPTIDE', charge=2, precision=3)
        801.375

        # Calculate the mass of a peptide sequence with a charge of -2.
        >>> mass('PEPTIDE', charge=-2, precision=3)
        797.345

        # Calcualte the mass of a modified peptide sequence.
        >>> mass('PE[3.14]PTIDE[Acetyl]', charge=2, precision=3)
        846.525

        # Calculate the mass of a peptide sequence with a charge of 2.
        >>> mass('PEPT[10][10]IDE', charge=2, precision=3)
        821.375

        # Calculate the mass of a peptide sequence with ambiguity
        >>> mass('(PEPT)[10]IDE', charge=2, precision=3)
        811.375

        >>> mass('(?DQ)NGTWEMESNENFEGYMK', precision=3)
        2307.905

        >>> mass('EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK', precision=3)
        1360.511

        >>> mass('{100}PEPTIDE', charge=0, precision=3, ion_type='p')
        899.36

        # Can also use ProFormaAnnotation to specify the charge
        >>> mass('PEPTIDE/2', precision=3)
        801.375

        # Or specify the charge directly
        >>> mass('PEPTIDE/2[+2Na+,+H+]', precision=3)
        846.346

        >>> mass('B', ion_type='by')  # doctest: +ELLIPSIS
        Traceback (most recent call last):
            ...
        peptacular.errors.AmbiguousAminoAcidError: ...

    """

    annotation = get_annotation_input(sequence=sequence, copy=True)

    override_annotation_properties(
        annotation=annotation,
        charge=charge,
    )

    return annotation.mass(
        ion_type=ion_type,
        monoisotopic=monoisotopic,
        isotope=isotope,
        loss=loss,
        precision=precision,
        use_isotope_on_mods=use_isotope_on_mods,
    )


def mz(
    sequence: str | ProFormaAnnotation,
    charge: int | None = None,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    precision: int | None = None,
    use_isotope_on_mods: bool = False,
) -> float:
    """
    Calculate the m/z (mass-to-charge ratio) of an amino acid 'sequence'.

    :param sequence: A sequence or ProFormaAnnotation.
    :type sequence: str | ProFormaAnnotation
    :param charge: The charge state, default is None.
    :type charge: int | None
    :param ion_type: The ion type. Default is 'p'.
    :type ion_type: str
    :param monoisotopic: If True, use monoisotopic mass else use average mass. Default is True.
    :type monoisotopic: bool
    :param isotope: The number of Neutrons to add/subtract from the final mass. Default is 0.
    :type isotope: int
    :param loss: The loss to add/subtract to the final mass. Default is 0.0.
    :type loss: float
    :param charge_adducts: The charge adducts. Default is None.
    :type charge_adducts: int | None
    :param isotope_mods: The isotope modifications. Default is None.
    :type isotope_mods: List[Mod] | None
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None

    :raise ValueError: If the ion type is not supported.

    :return: The Mass to Charge ratio (m/z) of the sequence.
    :rtype: float

    .. code-block:: python

        # Calculate the m/z of a peptide sequence.
        >>> mz('PEPTIDE', charge = 1, precision = 3)
        800.367

        # Calculate the b-ion m/z of a peptide sequence.
        >>> mz('PEPTIDE', charge = 1, ion_type='b', precision = 3)
        782.357

        # Calulate the average m/z of a peptide sequence.
        >>> mz('PEPTIDE', charge = 1, monoisotopic=False, precision = 3)
        800.831

        # Calculate the m/z of a peptide sequence with a charge of 2.
        >>> mz('PEPTIDE', charge=2, precision = 3)
        400.688

        # Calcualte the m/z of a modified peptide sequence.
        >>> mz('PE[3.14]PTIDE-[80]', charge=2, precision = 3)
        442.257

    """

    annotation = get_annotation_input(sequence=sequence, copy=True)

    override_annotation_properties(
        annotation=annotation,
        charge=charge,
    )

    return annotation.mz(
        ion_type=ion_type,
        monoisotopic=monoisotopic,
        isotope=isotope,
        loss=loss,
        precision=precision,
        use_isotope_on_mods=use_isotope_on_mods,
    )


def comp(
    sequence: str | ProFormaAnnotation,
    ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
    estimate_delta: bool = False,
    charge: int | None = None,
    isotope: int = 0,
    use_isotope_on_mods: bool = False,
) -> dict[str, int | float]:
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
    )

    return annotation.comp(
        ion_type=ion_type,
        isotope=isotope,
        use_isotope_on_mods=use_isotope_on_mods,
    )


def condense_to_mass_mods(
    sequence: str | ProFormaAnnotation,
    use_isotope_on_mods: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
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

    return annotation.condense_to_delta_mass(
        use_isotope_on_mods=use_isotope_on_mods,
        include_plus=include_plus,
        inplace=True,
    ).serialize(include_plus=include_plus, precision=precision)
