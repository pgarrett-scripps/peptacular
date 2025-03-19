"""
mass_calc.py is a simple module for computing the m/z and mass of an amino acid sequence.
"""

import warnings
from typing import Union, Optional, Tuple, List

from peptacular.constants import PROTON_MASS, AVERAGE_ATOMIC_MASSES, ELECTRON_MASS, ISOTOPIC_ATOMIC_MASSES, NEUTRON_MASS
from peptacular.chem.chem_calc import parse_chem_formula, _sequence_comp, _parse_mod_delta_mass_only, estimate_comp
from peptacular.chem.chem_util import chem_mass
from peptacular.mods.mod_db_setup import MONOSACCHARIDES_DB
from peptacular.chem.chem_constants import MONOISOTOPIC_AA_MASSES, AVERAGE_AA_MASSES, AVERAGE_FRAGMENT_ADJUSTMENTS, \
    MONOISOTOPIC_FRAGMENT_ADJUSTMENTS, MONOISOTOPIC_FRAGMENT_ION_ADJUSTMENTS, AVERAGE_FRAGMENT_ION_ADJUSTMENTS
from peptacular.proforma.proforma_parser import parse_static_mods, ProFormaAnnotation, Mod, parse_ion_elements
from peptacular.util import convert_type
from peptacular.errors import InvalidDeltaMassError, InvalidModificationMassError, UnknownAminoAcidError, \
    InvalidGlycanFormulaError, InvalidChemFormulaError, AmbiguousAminoAcidError
from peptacular.glycan import parse_glycan_formula
from peptacular.mods.mod_db import parse_psi_mass, parse_unimod_mass, is_unimod_str, is_psi_mod_str, parse_xlmod_mass, \
    parse_resid_mass, is_gno_str, parse_gno_mass, is_xlmod_str, is_resid_str
from peptacular.sequence.sequence_funcs import sequence_to_annotation
from peptacular.types import ChemComposition
from peptacular.proforma.input_convert import ModValue


def comp_mass(sequence: Union[str, ProFormaAnnotation],
              ion_type: str = 'p',
              charge: Optional[int] = None,
              isotope: int = 0,
              charge_adducts: Optional[str] = None,
              isotope_mods: Optional[List[ModValue]] = None,
              use_isotope_on_mods: bool = False) -> Tuple[ChemComposition, float]:
    """
    Get the chemical composition of a peptide sequence and its delta mass.

    :param sequence: A sequence or ProFormaAnnotation.
    :type sequence: str | ProFormaAnnotation
    :param ion_type: The ion type. Default is 'p'.
    :type ion_type: str
    :param charge: The charge state of the ion. Default is None.
    :type charge: int | None
    :param isotope: The number of Neutrons to add/subtract from the final mass. Default is 0.
    :type isotope: int
    :param charge_adducts: The charge adducts. Default is None.
    :type charge_adducts: str | None
    :param isotope_mods: The isotope modifications. Default is None.
    :type isotope_mods: List[ModValue] | None
    :param use_isotope_on_mods: If True, use the isotope on the modifications. Default is False.
    :type use_isotope_on_mods: bool

    :return: A tuple containing the chemical composition and delta mass.
    :rtype: Tuple[ChemComposition, float]

    .. code-block:: python

        # Unmodified Peptide
        >>> comp_mass('PEPTIDE')
        ({'C': 34, 'H': 53, 'N': 7, 'O': 15}, 0.0)

        >>> comp_mass('PEPTIDE', ion_type='b', charge=2)
        ({'C': 34, 'H': 53, 'N': 7, 'O': 14, 'e': -2}, 0.0)

        # Modified Peptide
        >>> comp_mass('PEPTIDE[1.0]')
        ({'C': 34, 'H': 53, 'N': 7, 'O': 15}, 1.0)

        # Charge Adducts
        >>> comp_mass('PEPTIDE', charge_adducts='+H+', charge=1)
        ({'C': 34, 'H': 54, 'N': 7, 'O': 15, 'e': -1}, 0.0)

        # Charge adducts in ProForma sequence
        >>> comp_mass('PEPTIDE/2[+2Na+]')
        ({'C': 34, 'H': 53, 'N': 7, 'O': 15, 'Na': 2, 'e': -2}, 0.0)

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence.copy()

    if charge is not None:
        annotation.charge = charge

    if charge_adducts is not None:
        annotation.charge_adducts = charge_adducts

    if isotope_mods is not None:
        annotation.isotope_mods = isotope_mods

    annotation.condense_static_mods(inplace=True)

    delta_mass = _pop_delta_mass_mods(annotation)  # sum delta mass mods and remove them from the annotation
    peptide_composition = _sequence_comp(annotation, ion_type, isotope, use_isotope_on_mods)

    if use_isotope_on_mods is True and delta_mass != 0:
        warnings.warn(
            'use_isotope_on_mods=True and delta_mass != 0. Cannot apply isotopic modifications to the delta mass.')

    return peptide_composition, delta_mass


def comp(sequence: Union[str, ProFormaAnnotation],
         ion_type: str = 'p',
         estimate_delta: bool = False,
         charge: Optional[int] = None,
         isotope: int = 0,
         charge_adducts: Optional[str] = None,
         isotope_mods: Optional[List[ModValue]] = None,
         use_isotope_on_mods: bool = False) -> ChemComposition:
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

        >>> comp('PEPTIDE[1.0]', estimate_delta=False)
        Traceback (most recent call last):
        ValueError: Non-zero delta mass (1.0) encountered without estimation enabled for sequence 'PEPTIDE[1.0]'.

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    composition, delta_mass = comp_mass(annotation, ion_type, charge, isotope, charge_adducts, isotope_mods,
                                        use_isotope_on_mods)

    if delta_mass != 0:

        if estimate_delta is False:
            raise ValueError(f"Non-zero delta mass ({delta_mass}) encountered without estimation enabled for "
                             f"sequence '{sequence}'.")

        if use_isotope_on_mods is True:
            delta_mass_comp = estimate_comp(delta_mass, annotation.isotope_mods)
            warnings.warn(
                'Applying isotopic modifications to the predicted composition. This will not be accurate.'
            )
        else:
            delta_mass_comp = estimate_comp(delta_mass, None)

        # Combine the compositions of the sequence and the delta mass
        for element in delta_mass_comp:
            composition.setdefault(element, 0)
            composition[element] += delta_mass_comp[element]

    return composition


def adjust_mass(base_mass: float,
                charge: int,
                ion_type: str = 'p',
                monoisotopic: bool = True,
                isotope: int = 0,
                loss: float = 0.0,
                charge_adducts: Optional[str] = None,
                precision: Optional[int] = None) -> float:
    """
    Adjust the mass.

    :param base_mass: The base mass of the sequence.
    :type base_mass: float
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
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None

    :return: The mass of the sequence.
    :rtype: float
    """

    if charge is None:
        charge = 0

    m = base_mass

    if charge_adducts is None:
        if ion_type in ('p', 'n'):
            charge_adduct_mass = PROTON_MASS * charge
        else:
            frag_ion_offset = MONOISOTOPIC_FRAGMENT_ION_ADJUSTMENTS[ion_type] if monoisotopic is True \
                else AVERAGE_FRAGMENT_ION_ADJUSTMENTS[ion_type]
            charge_adduct_mass = PROTON_MASS * (charge - 1) + frag_ion_offset
    else:
        charge_adduct_mass = _parse_charge_adducts_mass(charge_adducts, monoisotopic=monoisotopic)

    m += charge_adduct_mass

    # Base mass
    m += MONOISOTOPIC_FRAGMENT_ADJUSTMENTS[ion_type] if monoisotopic else AVERAGE_FRAGMENT_ADJUSTMENTS[ion_type]

    # m += (charge * charge_adduct_mass)
    m += isotope * NEUTRON_MASS + loss  # Add isotope and loss

    if precision is not None:
        m = round(m, precision)

    return m


def adjust_mz(base_mass: float,
              charge: int,
              precision: Optional[int] = None) -> float:
    """
    Adjust the mass to charge ratio (m/z).

    :param base_mass: The base mass of the sequence.
    :type base_mass: float
    :param charge: The charge state of the ion.
    :type charge: int
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None

    :return: The adjusted m/z of the sequence.
    :rtype: float
    """

    if charge is None:
        charge = 0

    m = base_mass

    m = m if charge == 0 else m / charge

    if precision is not None:
        m = round(m, precision)

    return m


def mass(sequence: Union[str, ProFormaAnnotation],
         charge: Optional[int] = None,
         ion_type: str = 'p',
         monoisotopic: bool = True,
         isotope: int = 0,
         loss: float = 0.0,
         charge_adducts: Optional[str] = None,
         isotope_mods: Optional[List[Mod]] = None,
         use_isotope_on_mods: bool = False,
         precision: Optional[int] = None) -> float:
    """
    Calculate the mass of an amino acid 'sequence'.

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
        peptacular.errors.AmbiguousAminoAcidError: Ambiguous amino acid: B! ...

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    if annotation.charge is not None and charge is None:
        charge = annotation.charge

    if annotation.charge_adducts is not None and charge_adducts is None:
        charge_adducts = annotation.charge_adducts[0]

    if annotation.isotope_mods is not None and isotope_mods is None:
        isotope_mods = annotation.isotope_mods

    if 'B' in annotation.sequence:
        raise AmbiguousAminoAcidError('B', 'Cannot determine the mass of a sequence with an ambiguous amino acid.')

    if 'Z' in annotation.sequence:
        raise AmbiguousAminoAcidError('Z', 'Cannot determine the mass of a sequence with an ambiguous amino acid.')

    if isotope_mods is not None and len(isotope_mods) > 0:
        peptide_composition, delta_mass = comp_mass(annotation, ion_type, charge, isotope, charge_adducts, isotope_mods,
                                                    use_isotope_on_mods)
        return chem_mass(peptide_composition, monoisotopic=monoisotopic, precision=precision) + delta_mass

    if ion_type not in ('p', 'n'):
        if charge == 0:
            warnings.warn('Calculating the mass of a fragment ion with charge state 0. Fragment ions should have a '
                          'charge state greater than 0 since the neutral mass doesnt exist.')

    m = 0.0
    if annotation.has_static_mods():
        static_map = parse_static_mods(annotation.static_mods)

        n_term_mod = static_map.get('N-Term')
        if n_term_mod is not None:
            m += sum(mod_mass(m, monoisotopic, precision=None) for m in n_term_mod)

        c_term_mod = static_map.get('C-Term')
        if c_term_mod is not None:
            m += sum(mod_mass(m, monoisotopic, precision=None) for m in c_term_mod)

        for aa, mod in static_map.items():

            if aa in ['N-Term', 'C-Term']:
                continue

            aa_count = annotation.sequence.count(aa)
            m += sum(mod_mass(m, monoisotopic, precision=None) for m in mod) * aa_count

    try:
        m += sum(MONOISOTOPIC_AA_MASSES[aa] if monoisotopic
                 else AVERAGE_AA_MASSES[aa] for aa in annotation.sequence)
    except KeyError as err:
        raise UnknownAminoAcidError(err) from err

    # Apply labile mods
    if annotation.has_labile_mods() and ion_type == 'p':
        for mod in annotation.labile_mods:
            m += mod_mass(mod)

    # Apply Unknown mods
    if annotation.has_unknown_mods():
        for mod in annotation.unknown_mods:
            m += mod_mass(mod)

    # Apply N-term mods
    if annotation.has_nterm_mods():
        for mod in annotation.nterm_mods:
            m += mod_mass(mod)

    # Apply intervals
    if annotation.has_intervals():
        for interval in annotation.intervals:
            if interval.mods is not None:
                for mod in interval.mods:
                    m += mod_mass(mod)

    # apply internal mods
    if annotation.has_internal_mods():
        for _, mods in annotation.internal_mods.items():
            for mod in mods:
                m += mod_mass(mod)

    # apply C-term mods
    if annotation.has_cterm_mods():
        for mod in annotation.cterm_mods:
            m += mod_mass(mod)

    return adjust_mass(m, charge, ion_type, monoisotopic, isotope, loss, charge_adducts, precision)


def mz(sequence: Union[str, ProFormaAnnotation],
       charge: Optional[int] = None,
       ion_type: str = 'p',
       monoisotopic: bool = True,
       isotope: int = 0,
       loss: float = 0.0,
       charge_adducts: Optional[str] = None,
       isotope_mods: Optional[List[Mod]] = None,
       precision: Optional[int] = None) -> float:
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
        400.687

        # Calcualte the m/z of a modified peptide sequence.
        >>> mz('PE[3.14]PTIDE-[80]', charge=2, precision = 3)
        442.257

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    if annotation.charge is not None and charge is None:
        charge = annotation.charge

    m = mass(sequence=annotation, charge=charge, ion_type=ion_type,
             monoisotopic=monoisotopic, isotope=isotope, loss=loss, charge_adducts=charge_adducts,
             isotope_mods=isotope_mods, precision=None)

    return adjust_mz(m, charge, precision)


def chem_mz(formula: Union[ChemComposition, str],
            charge: int = 1,
            monoisotopic: bool = True,
            precision: Optional[int] = None,
            sep: str = '') -> float:
    """
    Calculate the m/z of a chemical formula.
    """
    # TODO: Add charge adducts?
    m = chem_mass(formula, monoisotopic, precision, sep)
    return adjust_mz(m, charge, precision)


def glycan_mass(formula: Union[str, ChemComposition],
                monoisotopic: bool = True,
                precision: Optional[int] = None) -> float:
    """
    Calculate the mass of a glycan formula.

    :param formula: A glycan formula or ChemComposition.
    :type formula: str | Dict[str, int | float]
    :param monoisotopic: If True, use monoisotopic mass else use average mass. Default is True.
    :type monoisotopic: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None

    :raises UnknownGlycanError: If the glycan formula contains an unknown monosaccharide.

    :return: The mass of the glycan.
    :rtype: float

    .. code-block:: python

        # Calculate the mass of a glycan formula.
        >>> glycan_mass({'HexNAc': 2, 'Hex': 3, 'Neu': 1}, precision=3)
        1141.402

        >>> glycan_mass('HexNAc2Hex3Neu1', precision=3)
        1141.402

        >>> glycan_mass({'HexNAc': 2, 'Hex': 3, 'Neu': 1}, monoisotopic=False, precision=3)
        1142.027

        >>> glycan_mass('HexNAc2Hex3Neu1', monoisotopic=False, precision=3)
        1142.027

        # Example Error# show example error
        >>> glycan_mass('HexX', precision=3)
        Traceback (most recent call last):
        peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "HexX". Unknown glycan: "X"!

    """
    original_formula = formula
    if isinstance(formula, str):
        formula = parse_glycan_formula(formula)

    m = 0.0
    for monosaccharide, count in formula.items():

        if MONOSACCHARIDES_DB.contains_name(monosaccharide):
            entry = MONOSACCHARIDES_DB.get_entry_by_name(monosaccharide)
        elif MONOSACCHARIDES_DB.contains_synonym(monosaccharide):
            entry = MONOSACCHARIDES_DB.get_entry_by_synonym(monosaccharide)
        else:
            raise InvalidGlycanFormulaError(original_formula, f'Unknown monosaccharide: "{monosaccharide}"!')

        if monoisotopic:
            m += entry.mono_mass * count
        else:
            m += entry.avg_mass * count

    if precision is not None:
        m = round(m, precision)

    return m


def glycan_mz(formula: Union[str, ChemComposition],
              charge: int = 1,
              monoisotopic: bool = True,
              precision: Optional[int] = None) -> float:
    """
    Calculate the m/z of a glycan formula.
    """
    # TODO: Add charge adducts?
    m = glycan_mass(formula, monoisotopic, precision)
    return adjust_mz(m, charge, precision)


def mod_mass(mod: Union[str, Mod, List[Mod]], monoisotopic: bool = True, precision: Optional[int] = None) -> float:
    """
    Parse a modification string.

    :param mod: The modification string or Mod.
    :type mod: str | Mod
    :param monoisotopic: If True, use monoisotopic mass else use average mass. Default is True.
    :type monoisotopic: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None

    :raises UnknownModificationError: If the modification is unknown.
    :raises InvalidDeltaMassError: If the modification contains an invalid delta mass.
    :raises UnknownElementError: If the modification contains an unknown element.
    :raises UnknownGlycanError: If the modification contains an unknown glycan.
    :raises NotImplementedError: If the modification is not implemented.
    :raises InvalidModificationMassError: If the modification cannot be parsed.

    :return: The parsed modification mass.
    :rtype: float

    .. code-block:: python

        >>> mod_mass('Acetyl|INFO:newly discovered', precision=3)
        42.011

        >>> mod_mass(Mod('Acetyl|INFO:newly discovered', 2), precision=3)
        84.022

        >>> mod_mass('1', precision=3)
        1

        >>> mod_mass(['1', 1], precision=3)
        2

        >>> mod_mass('Acetyl|Obs:+42.010565', precision=3)
        42.011

        # example error
        >>> mod_mass('Acetdsyl|Obs:42d.010565', precision=3)
        Traceback (most recent call last):
        ...
        peptacular.errors.InvalidDeltaMassError: Invalid delta mass: 42d.010565

        # example invalid mass
        >>> mod_mass('info:HelloWorld', precision=3)
        Traceback (most recent call last):
        ...
        peptacular.errors.InvalidModificationMassError: Cannot determine mass for modification: "info:HelloWorld"

    """

    if isinstance(mod, list):
        return sum(mod_mass(m, monoisotopic, precision) for m in mod)

    if isinstance(mod, Mod):
        return mod_mass(mod.val, monoisotopic, precision) * mod.mult

    if isinstance(mod, int):
        return mod

    if isinstance(mod, float):
        return round(mod, precision) if precision is not None else mod

    mods = mod.split('|')
    for m in mods:
        m = _parse_mod_mass(m, monoisotopic, precision)
        if m is not None:
            return m

    raise InvalidModificationMassError(mod)


def _parse_obs_mass_from_proforma_str(obs_str: str, precision: Optional[int] = None) -> float:
    """
    Parse an observed mass string and return its mass.

    :param obs_str: The observed mass string to parse.
    :type obs_str: str
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None

    :raises InvalidDeltaMassError: If the observed mass string contains an invalid delta mass.

    :return: The mass of the observed mass string.
    :rtype: float


    .. code-block:: python

        # Parse an observed mass string.
        >>> _parse_obs_mass_from_proforma_str('42.0', precision=3)
        42.0

        >>> _parse_obs_mass_from_proforma_str('+42.0', precision=3)
        42.0

        >>> _parse_obs_mass_from_proforma_str('-42.0', precision=3)
        -42.0

        >>> _parse_obs_mass_from_proforma_str('d42.0', precision=3)
        Traceback (most recent call last):
        ...
        peptacular.errors.InvalidDeltaMassError: Invalid delta mass: d42.0

    """

    round_func = lambda x: round(x, precision) if precision is not None else x

    if obs_str.lower().startswith('obs:'):
        obs_str = ''.join(obs_str.split(':')[1:])

    # Try to Parse observed mass
    try:
        return round_func(float(obs_str))
    except ValueError as err:
        raise InvalidDeltaMassError(obs_str) from err


def _parse_glycan_mass_from_proforma_str(glycan_str: str, monoisotopic: bool, precision: Optional[int] = None) -> float:
    """
    Parse a glycan string and return its mass.

    :param glycan_str: The glycan string to parse.
    :type glycan_str: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None

    :raises InvalidGlycanFormulaError: If the glycan formula is invalid.

    :return: The mass of the glycan.
    :rtype: float

    .. code-block:: python

        # Monoisotopic Mass
        >>> _parse_glycan_mass_from_proforma_str('HexNAc2Hex3Neu1', monoisotopic=True, precision=3)
        1141.402

        # Average Mass
        >>> _parse_glycan_mass_from_proforma_str('HexNAc2Hex3Neu1', monoisotopic=False, precision=3)
        1142.027

        # Using a glycan name
        >>> _parse_glycan_mass_from_proforma_str('HexNAc', monoisotopic=False, precision=3)
        203.193

        # Using a glycan ID
        >>> _parse_glycan_mass_from_proforma_str('6BAAE1B1', monoisotopic=False, precision=3)
        72.063

        # Using a glycan ID
        >>> _parse_glycan_mass_from_proforma_str('Glycan:6BAAE1B1', monoisotopic=False, precision=3)
        72.063

        # Invalid glycan ID
        >>> _parse_glycan_mass_from_proforma_str('Glycan:XXX', monoisotopic=True, precision=3)
        Traceback (most recent call last):
        peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "XXX". Unknown glycan: "XXX"!

        # Invalid glycan str
        >>> _parse_glycan_mass_from_proforma_str('Glycan:HexXX', monoisotopic=True, precision=3)
        Traceback (most recent call last):
        peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "HexXX". Unknown glycan: "XX"!

    """

    round_func = lambda x: round(x, precision) if precision is not None else x

    if glycan_str.lower().startswith('glycan:'):
        glycan_str = ''.join(glycan_str.split(':')[1:])

    if MONOSACCHARIDES_DB.contains_id(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_id(glycan_str)
    elif MONOSACCHARIDES_DB.contains_name(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_name(glycan_str)
    elif MONOSACCHARIDES_DB.contains_synonym(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_synonym(glycan_str)
    else:
        entry = None

    if entry:
        if monoisotopic:
            return round_func(entry.mono_mass)
        return round_func(entry.avg_mass)

    else:  # Try to parse glycan formula
        try:
            return round_func(glycan_mass(glycan_str, monoisotopic, precision))
        except InvalidGlycanFormulaError as err:
            raise InvalidGlycanFormulaError(glycan_str, err.msg) from err


def _parse_chem_mass_from_proforma_str(chem_str: str, monoisotopic: bool, precision: Optional[int] = None) -> float:
    """
    Parse a chemical formula string and return its mass.

    :param chem_str: The chemical formula string to parse.
    :type chem_str: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None

    :raises InvalidChemFormulaError: If the chemical formula is invalid.

    :return: The mass of the chemical formula.
    :rtype: float

    .. code-block:: python

        # Monoisotopic Mass
        >>> _parse_chem_mass_from_proforma_str('C2H4O2', monoisotopic=True, precision=3)
        60.021

        # Average Mass
        >>> _parse_chem_mass_from_proforma_str('C2H4O2', monoisotopic=False, precision=3)
        60.052

        # Invalid chemical formula
        >>> _parse_chem_mass_from_proforma_str('C2H4O2X', monoisotopic=True, precision=3)
        Traceback (most recent call last):
        peptacular.errors.InvalidChemFormulaError: Error parsing chem formula: "C2H4O2X". Unknown element: "X"!

    """

    if chem_str.lower().startswith('formula:'):
        chem_str = ''.join(chem_str.split(':')[1:])

    comps = parse_chem_formula(chem_str)
    try:
        return chem_mass(comps, monoisotopic, precision)
    except InvalidChemFormulaError as err:
        raise InvalidChemFormulaError(chem_str, err.msg) from err


def _parse_mod_mass(mod: str, monoisotopic: bool = True, precision: Optional[int] = None) -> Union[float, None]:
    """
    Parse a modification and return its mass.

    :param mod: The modification to parse.
    :type mod: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None

    :raises UnknownModificationError: If the modification is unknown.
    :raises InvalidDeltaMassError: If the modification contains an invalid delta mass.
    :raises UnknownElementError: If the modification contains an unknown element.
    :raises UnknownGlycanError: If the modification contains an unknown glycan.
    :raises NotImplementedError: If the modification is not implemented.

    :return: The mass of the modification or None if the modification has no mass.
    :rtype:  float | None

    .. code-block:: python

        # Parse numeric modifications.
        >>> _parse_mod_mass('42.0')
        42.0

        >>> _parse_mod_mass('+42.0')
        42.0

        >>> _parse_mod_mass('-42.0')
        -42.0

        # Parse formula modifications.
        >>> _parse_mod_mass('Formula:C2')
        24.0

        # Parse observed modifications.
        >>> _parse_mod_mass('Obs:42.0')
        42.0

        # Parse UniMod modifications.
        >>> _parse_mod_mass('UniMod:1')
        42.010565
        >>> _parse_mod_mass('Acetyl')
        42.010565

        # Parse PSI-MOD modifications.
        >>> _parse_mod_mass('MOD:00231')
        521.395649
        >>> _parse_mod_mass('hexakis-L-cysteinyl hexairon hexasulfide')
        521.395649

        # Parse glycans
        >>> _parse_mod_mass('Glycan:HexNAc2Hex3Neu1')
        1141.4020671170001

        >>> _parse_mod_mass('U:+15.9949')
        15.9949

        >>> _parse_mod_mass('M:+15.9949')
        15.9949

        >>> _parse_mod_mass('#g1(0.1)')
        0.0

        >>> _parse_mod_mass('Acetyl#g1(0.1)')
        42.010565

        >>> _parse_mod_mass('Amidated')
        -0.984016

    """

    # localization fix
    if isinstance(mod, str) and '#' in mod:
        if mod.startswith('#'):  # for localized positions, return 0
            return 0.0

        mod = mod.split('#')[0]

    # Try to parse as a number first (Might cause issues if the mod is also unimod/psi ID)
    # Proforma2.0 standard requires that delta mass instances are always prefixed with a '+' or '-' but this would
    # require a lot of changes....and make user input more difficult
    mod = convert_type(mod)
    if isinstance(mod, int):
        return mod
    elif isinstance(mod, float):
        return round(mod, precision) if precision is not None else mod

    mod_lower = mod.lower()
    if mod_lower.startswith('glycan:'):
        # raises UnknownGlycanError
        return _parse_glycan_mass_from_proforma_str(mod, monoisotopic, precision)

    if is_gno_str(mod):
        return parse_gno_mass(mod, monoisotopic, precision)

    if is_xlmod_str(mod):
        # not implemented
        return parse_xlmod_mass(mod, monoisotopic, precision)

    if is_resid_str(mod):
        # not implemented
        return parse_resid_mass(mod, monoisotopic, precision)

    if mod_lower.startswith('info:'):  # Skip info modifications
        return None

    if is_psi_mod_str(mod):  # is a psi-mod modification
        # raises UnknownModificationError and InvalidDeltaMassError
        return parse_psi_mass(mod, monoisotopic, precision)

    if is_unimod_str(mod):  # is a unimod modification
        # raises UnknownModificationError and InvalidDeltaMassError
        return parse_unimod_mass(mod, monoisotopic, precision)

    # chemical formula
    if mod_lower.startswith('formula:'):
        # raises UnknownElementError
        return _parse_chem_mass_from_proforma_str(mod, monoisotopic, precision)

    # observed mass
    if mod_lower.startswith('obs:'):
        # raises InvalidDeltaMassError
        return _parse_obs_mass_from_proforma_str(mod, precision)

    return None


def _pop_delta_mass_mods(annotation: ProFormaAnnotation) -> float:
    """
    Pop the delta mass modifications from the modifications' dictionary. This leaves only modifications which
    can have their elemental composition calculated.

    :param annotation: The ProForma annotation.
    :type annotation: ProFormaAnnotation

    :return: The delta mass.
    :rtype: float

    .. code-block:: python

        >>> mods = ProFormaAnnotation(_sequence='', _nterm_mods = [42.0, -20.0])
        >>> _pop_delta_mass_mods(mods)
        22.0
        >>> mods
        ProFormaAnnotation(sequence=)

    """

    delta_mass = 0.0

    # Labile Mods
    if annotation.has_labile_mods():
        for i in range(len(annotation.labile_mods) - 1, -1, -1):
            mod = annotation.labile_mods[i]
            val = _parse_mod_delta_mass_only(mod)
            if val is not None:
                delta_mass += val
                annotation.labile_mods.pop(i)

    # Unknown Mods
    if annotation.has_unknown_mods():
        for i in range(len(annotation.unknown_mods) - 1, -1, -1):
            mod = annotation.unknown_mods[i]
            val = _parse_mod_delta_mass_only(mod)
            if val is not None:
                delta_mass += val
                annotation.unknown_mods.pop(i)

    # NTerm
    if annotation.has_nterm_mods():
        for i in range(len(annotation.nterm_mods) - 1, -1, -1):
            mod = annotation.nterm_mods[i]
            val = _parse_mod_delta_mass_only(mod)
            if val is not None:
                delta_mass += val
                annotation.nterm_mods.pop(i)

    # CTerm
    if annotation.has_cterm_mods():
        for i in range(len(annotation.cterm_mods) - 1, -1, -1):
            mod = annotation.cterm_mods[i]
            val = _parse_mod_delta_mass_only(mod)
            if val is not None:
                delta_mass += val
                annotation.cterm_mods.pop(i)

    # Intervals
    if annotation.has_intervals():
        for interval in annotation.intervals:
            if interval.has_mods():
                for i in range(len(interval.mods) - 1, -1, -1):
                    mod = interval.mods[i]
                    val = _parse_mod_delta_mass_only(mod)
                    if val is not None:
                        delta_mass += val
                        interval.mods.pop(i)

    # Internal mods
    if annotation.has_internal_mods():
        for k in annotation.internal_mods:
            for i in range(len(annotation.internal_mods[k]) - 1, -1, -1):
                mod = annotation.internal_mods[k][i]
                val = _parse_mod_delta_mass_only(mod)
                if val is not None:
                    delta_mass += val
                    annotation.internal_mods[k].pop(i)

    annotation.clear_empty_mods()

    return delta_mass


def _parse_adduct_mass(adduct: str,
                       precision: Optional[int] = None,
                       monoisotopic: bool = True) -> float:
    """
    Parse an adduct string and return its mass.

    :param adduct: The adduct string to parse.
    :type adduct: str
    :param precision: The precision of the mass. Default is None.
    :type precision: Optional[int]
    :param monoisotopic: Whether to use monoisotopic masses. Default is True.
    :type monoisotopic: bool

    :raises InvalidDeltaMassError: If the adduct contains an invalid delta mass.

    :return: The mass of the adduct.
    :rtype: float

    .. code-block:: python

        # Parse an adduct string.
        >>> _parse_adduct_mass('+Na+', precision=5)
        22.98922

        >>> _parse_adduct_mass('+2Na+', precision=5)
        45.97899

        >>> _parse_adduct_mass('+2Na-', precision=5)
        45.98009

        >>> _parse_adduct_mass('H+', precision=5)
        1.00728

        >>> _parse_adduct_mass('H-', precision=5)
        1.00837

    """

    m = 0.0
    element_count, element_symbol, element_charge = parse_ion_elements(adduct)

    if element_symbol == 'e':
        m += element_count * ELECTRON_MASS

    else:

        if monoisotopic is True:
            m += element_count * ISOTOPIC_ATOMIC_MASSES[element_symbol]
            m -= element_charge * ELECTRON_MASS
        else:
            m += element_count * AVERAGE_ATOMIC_MASSES[element_symbol]
            m -= element_charge * ELECTRON_MASS

        if precision is not None:
            m = round(m, precision)

    return m


def _parse_charge_adducts_mass(adducts: ModValue,
                               precision: Optional[int] = None,
                               monoisotopic: bool = True) -> float:
    """
    Parse the charge adducts and return their mass.

    :param adducts: The charge adducts to parse.
    :type adducts: ModValue
    :param precision: The precision of the mass. Default is None.
    :type precision: Optional[int]
    :param monoisotopic: Whether to use monoisotopic masses. Default is True.
    :type monoisotopic: bool

    :raises InvalidDeltaMassError: If the adduct contains an invalid delta mass.

    :return: The mass of the charge adducts.
    :rtype: float

    .. code-block:: python

        # Parse the charge adducts and return their mass.
        >>> _parse_charge_adducts_mass('+Na+,+H+', precision=5)
        23.9965

    """

    if isinstance(adducts, Mod):
        return _parse_charge_adducts_mass(adducts.val, precision, monoisotopic)

    if not isinstance(adducts, str):
        raise TypeError(f'Invalid type for adducts: {type(adducts)}! Must be a string.')

    if adducts == '+H+':
        return PROTON_MASS

    adducts = adducts.split(',')

    m = 0.0

    for adduct in adducts:
        m += _parse_adduct_mass(adduct, None, monoisotopic)

    if precision is not None:
        m = round(m, precision)

    return m


def ppm_error(theo: float, expt: float, precision: Optional[int] = None) -> float:
    """
    Calculate the parts per million (ppm) error between two values.

    :param theo: The theoretical value.
    :type theo: float
    :param expt: The experimental value.
    :type expt: float
    :param precision: The precision of the ppm error. Default is None.
    :type precision: int | None

    :return: The parts per million error.
    :rtype: float

    .. code-block:: python

        # Calculate the parts per million error between two values.
        >>> ppm_error(100.0, 100.1, 2)
        1000.0

    """

    ppm_err = ((expt - theo) / theo) * 1e6

    if precision is not None:
        return round(ppm_err, precision)

    return ppm_err


def dalton_error(theo: float, expt: float, precision: Optional[int] = None) -> float:
    """
    Calculate the Dalton error between two values.

    :param theo: The theoretical value.
    :type theo: float
    :param expt: The experimental value.
    :type expt: float
    :param precision: The precision of the Dalton error. Default is None.
    :type precision: int | None

    :return: The Dalton error.
    :rtype: float

    .. code-block:: python

        # Calculate the Dalton error between two values.
        >>> dalton_error(100.0, 100.1, 2)
        0.1

    """

    dalton_err = expt - theo

    if precision is not None:
        return round(dalton_err, precision)

    return dalton_err


def condense_to_mass_mods(sequence: Union[str, ProFormaAnnotation], include_plus: bool = False, precision: float = 6) \
        -> str:
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

        >>> condense_to_mass_mods('PEP[Phospho]TIDE')
        'PEP[79.966331]TIDE'

        >>> condense_to_mass_mods('PEP[Phospho]TIDE', include_plus=True)
        'PEP[+79.966331]TIDE'

        >>> condense_to_mass_mods('[Acetyl]-PEPTIDE')
        '[42.010565]-PEPTIDE'

        >>> condense_to_mass_mods('PEPTIDE-[Amidated]')
        'PEPTIDE-[-0.984016]'

        >>> condense_to_mass_mods('<13C>PEP[Phospho]TIDE')
        'P[5.016774]E[5.016774]P[84.983105]T[4.013419]I[6.020129]D[4.013419]E[5.016774]'

    """
    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    n_term_mods = annotation.pop_nterm_mods()
    c_term_mods = annotation.pop_cterm_mods()
    labile_mods = annotation.pop_labile_mods()

    # Split into segments and process each one
    segments = list(annotation.split())
    stripped_segments = [seg.strip() for seg in segments]

    new_annotation = annotation.strip()

    # Calculate mass differences
    for i, (segment, stripped) in enumerate(zip(segments, stripped_segments)):
        # Calculate mass difference
        mod_mass_val = mass(segment) - mass(stripped)
        if abs(mod_mass_val) > 1e-6:  # Only add if difference is significant
            new_annotation.add_internal_mod(i, round(mod_mass_val, precision))

    if n_term_mods is not None:
        n_term_mods_mass = sum(mod_mass(mod) for mod in n_term_mods)

        new_annotation.add_nterm_mods(round(n_term_mods_mass, precision))

    if c_term_mods is not None:
        c_term_mods_mass = sum(mod_mass(mod) for mod in c_term_mods)
        new_annotation.add_cterm_mods(round(c_term_mods_mass, precision))

    if labile_mods is not None:
        labile_mods_mass = sum(mod_mass(mod) for mod in labile_mods)
        new_annotation.add_labile_mods(round(labile_mods_mass, precision))

    return new_annotation.serialize(include_plus=include_plus)
