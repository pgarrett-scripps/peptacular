"""
mass.py is a simple module for computing the m/z and mass of an amino acid sequence.
"""
import warnings
from typing import Dict, Union, Optional

from peptacular.constants import PROTON_MASS
from peptacular.chem.chem import parse_chem_formula, _sequence_comp, _parse_mod_delta_mass_only, estimate_comp
from peptacular.chem.chem_util import chem_mass, _parse_charge_adducts_mass, NEUTRON_MASS
from peptacular.mods.mod_db_setup import MONOSACCHARIDES_DB
from peptacular.chem.chem_constants import MONOISOTOPIC_AA_MASSES, AVERAGE_AA_MASSES, MONOISOTOPIC_ION_ADJUSTMENTS, \
    AVERAGE_ION_ADJUSTMENTS
from peptacular.proforma.proforma import parse_static_mods, ProFormaAnnotation, Mod
from peptacular.util import convert_type
from peptacular.errors import InvalidDeltaMassError, InvalidModificationMassError, UnknownAminoAcidError, \
    InvalidGlycanFormulaError, InvalidChemFormulaError
from peptacular.glycan import parse_glycan_formula
from peptacular.mods.mod_db import parse_psi_mass, parse_unimod_mass, is_unimod_str, is_psi_mod_str, parse_xlmod_mass, \
    parse_resid_mass, is_gno_str, parse_gno_mass, is_xlmod_str, is_resid_str
from peptacular.sequence.sequence import sequence_to_annotation
from peptacular.types import ChemComposition


def comp_mass(sequence: Union[str, ProFormaAnnotation], ion_type: str = 'p',
              charge: Optional[int] = None, charge_adducts: Optional[str] = None) -> (Dict[str, int], float):
    """
    Get the chemical composition of a peptide sequence and its delta mass.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param ion_type: The ion type.
    :type ion_type: str
    :param charge: The charge state of the ion.
    :type charge: int
    :param charge_adducts: The charge adducts.
    :type charge_adducts: str

    :return: The chemical composition and delta mass of the peptide sequence.
    :rtype: tuple

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
        >>> comp_mass('PEPTIDE/2[+2Na+,+H+]')
        ({'C': 34, 'H': 54, 'N': 7, 'O': 15, 'e': -3, 'Na': 2}, 0.0)

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    annotation = annotation.condense_static_mods(inplace=False)
    delta_mass = _pop_delta_mass_mods(annotation)  # get delta mass and pop these mods
    peptide_composition = _sequence_comp(annotation, ion_type, charge, charge_adducts)
    return peptide_composition, delta_mass


def comp(sequence: Union[str, ProFormaAnnotation], ion_type: str = 'p', estimate_delta: bool = False,
         charge: Optional[int] = None, charge_adducts: Optional[str] = None) -> Dict[str, int]:
    """
    Calculates the elemental composition of a peptide sequence, including modifications,
    and optionally estimates the composition based on the delta mass from modifications.

    :param sequence: The peptide sequence to analyze.
    :type sequence: str
    :param ion_type: The type of ion ('p' for precursor, or other types indicating different ionization forms).
    :type ion_type: str
    :param estimate_delta: If True, estimate the composition based on the delta mass from modifications.
    :type estimate_delta: bool
    :param charge: The charge state of the ion.
    :type charge: int
    :param charge_adducts: The charge adducts.
    :type charge_adducts: str

    :raises ValueError: If delta_mass is nonzero and estimate_delta is False, indicating an unaccounted modification.

    :return: A dictionary with elements as keys and their counts as values.
    :raises ValueError: If delta_mass is nonzero and estimate_delta is False, indicating an unaccounted modification.

    .. code-block:: python

        # Calculate the elemental composition of a peptide sequence.
        >>> comp('PEPTIDE')
        {'C': 34, 'H': 53, 'N': 7, 'O': 15}

        >>> comp('PEPTIDE')
        {'C': 34, 'H': 53, 'N': 7, 'O': 15}

        >>> comp('PEPTIDE[1.0]', estimate_delta=True)
        {'C': 34.04446833455479, 'H': 53.06986041632441, 'N': 7.012225550345263, 'O': 15.01330250093913, 'S': 0.0003754919712730832}

        >>> comp('PEPTIDE[1.0]', estimate_delta=False)
        Traceback (most recent call last):
        ValueError: Non-zero delta mass (1.0) encountered without estimation enabled for sequence 'PEPTIDE[1.0]'.

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    composition, delta_mass = comp_mass(annotation, ion_type, charge, charge_adducts)

    if delta_mass != 0 and not estimate_delta:
        raise ValueError(f"Non-zero delta mass ({delta_mass}) encountered without estimation enabled for "
                         f"sequence '{sequence}'.")

    if delta_mass != 0:
        delta_mass_comp = estimate_comp(delta_mass, annotation.isotope_mods)

        # Combine the compositions of the sequence and the delta mass
        for element in delta_mass_comp:
            composition.setdefault(element, 0)
            composition[element] += delta_mass_comp[element]

    return composition


def mass(sequence: Union[str, ProFormaAnnotation],
         charge: Optional[int] = None,
         ion_type: str = 'p',
         monoisotopic: bool = True,
         isotope: int = 0,
         loss: float = 0.0,
         charge_adducts: Optional[str] = None,
         precision: Optional[int] = None) -> float:
    """
    Calculate the mass of an amino acid 'sequence'.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param charge: The charge state, default is [0].
    :type charge: int
    :param ion_type: The ion type, default is [y].
    :type ion_type: str
    :param monoisotopic: If True, use monoisotopic mass else use average mass, defaults to [True].
    :type monoisotopic: bool
    :param isotope: The isotope number, defaults to [0].
    :type isotope: int
    :param loss: The loss, defaults to [0.0].
    :type loss: float
    :param charge_adducts: The charge adducts, defaults to [None].
    :type charge_adducts: str
    :param precision: The precision of the mass, defaults to [None].
    :type precision: int

    :raise ValueError: If the ion type is not supported

    :return: Mass of the peptide sequence.
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
        2291.91

        >>> mass('EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK', precision=3)
        1360.511

        >>> mass('{100}PEPTIDE', charge=0, precision=3, ion_type='p')
        899.36

        # Can also use ProFormaAnnotation to specify the charge
        >>> mass('PEPTIDE/2', precision=3)
        801.375

        # Or specify the charge directly
        >>> mass('PEPTIDE/2[+2Na+,+H+]', precision=3)
        893.332

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    if annotation.charge is not None and charge is None:
        charge = annotation.charge
    if charge is None:
        charge = 0

    if annotation.charge_adducts is not None and charge_adducts is None:
        charge_adducts = annotation.charge_adducts[0]
    if charge_adducts is None:
        charge_adducts = '+H+'

    if 'B' in annotation.sequence or 'J' in annotation.sequence:
        raise ValueError('Invalid amino acid sequence: B and J are ambiguous')

    if annotation.isotope_mods is not None and len(annotation.isotope_mods) > 0:
        peptide_composition, delta_mass = comp_mass(annotation, ion_type)
        return chem_mass(peptide_composition, monoisotopic=monoisotopic, precision=precision) + delta_mass

    if ion_type != 'p':
        if charge == 0:
            warnings.warn('Calculating the mass of a fragment ion with charge state 0. Fragment ions should have a '
                          'charge state greater than 0 since the neutral mass doesnt exist.')
        _ = annotation.pop_labile_mods()

    m = 0.0
    if annotation.static_mods is not None:
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

        annotation.pop_static_mods()

    try:
        m += sum(MONOISOTOPIC_AA_MASSES[aa] if monoisotopic
                 else AVERAGE_AA_MASSES[aa] for aa in annotation.sequence)
    except KeyError as e:
        raise UnknownAminoAcidError(e)

    # Apply labile mods
    if annotation.labile_mods is not None:
        for mod in annotation.labile_mods:
            m += mod_mass(mod)

    # Apply Unknown mods
    if annotation.unknown_mods is not None:
        for mod in annotation.unknown_mods:
            m += mod_mass(mod)

    # Apply N-term mods
    if annotation.nterm_mods is not None:
        for mod in annotation.nterm_mods:
            m += mod_mass(mod)

    # Apply intervals
    if annotation.intervals is not None:
        for interval in annotation.intervals:
            if interval.mods is not None:
                for mod in interval.mods:
                    m += mod_mass(mod)

    # apply internal mods
    if annotation.internal_mods is not None:
        for k, mods in annotation.internal_mods.items():
            for mod in mods:
                m += mod_mass(mod)

    # apply C-term mods
    if annotation.cterm_mods is not None:
        for mod in annotation.cterm_mods:
            m += mod_mass(mod)

    charge_adduct_mass = _parse_charge_adducts_mass(charge_adducts)
    m += (charge * charge_adduct_mass)  # Add charge
    m += isotope * NEUTRON_MASS + loss  # Add isotope and loss
    m += MONOISOTOPIC_ION_ADJUSTMENTS[ion_type] if monoisotopic else AVERAGE_ION_ADJUSTMENTS[ion_type]  # +1 ion mass
    m -= PROTON_MASS  # Subtract proton mass

    if precision is not None:
        m = round(m, precision)

    return m


def mz(sequence: Union[str, ProFormaAnnotation],
       charge: Optional[int] = None,
       ion_type: str = 'p',
       monoisotopic: bool = True,
       isotope: int = 0,
       loss: float = 0.0,
       charge_adducts: Optional[str] = None,
       precision: Optional[int] = None) -> float:
    """
    Calculate the m/z (mass-to-charge ratio) of an amino acid 'sequence'.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param charge: The charge state, default is [0].
    :type charge: int
    :param ion_type: The ion type, default is [y].
    :type ion_type: str
    :param monoisotopic: If True, use monoisotopic mass else use average mass, defaults to [True].
    :type monoisotopic: bool
    :param isotope: The isotope number, defaults to [0].
    :type isotope: int
    :param loss: The loss, defaults to [0.0].
    :type loss: float
    :param charge_adducts: The charge adducts, defaults to +H+.
    :type charge_adducts: str
    :param precision: The precision of the m/z, defaults to [None].
    :type precision: int

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: m/z of the peptide sequence.
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

    if charge is None:
        charge = 0

    m = mass(sequence=annotation, charge=charge, ion_type=ion_type,
             monoisotopic=monoisotopic, isotope=isotope, loss=loss, charge_adducts=charge_adducts, precision=None)

    m = m if charge == 0 else m / charge

    if precision is not None:
        m = round(m, precision)

    return m


def glycan_mass(formula: Union[str, ChemComposition],
                monoisotopic: bool = True,
                precision: Optional[int] = None) -> float:
    """
    Calculate the mass of a glycan formula.

    :param formula: The glycan formula.
    :type formula: dict
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass.
    :type precision: int

    :raises UnknownGlycanError: If the glycan formula contains an unknown monosaccharide.

    :return: The mass of the glycan formula.
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


def mod_mass(mod: Union[str, Mod], monoisotopic: bool = True, precision: Optional[int] = None) -> float:
    """
    Parse a modification string.

    :param mod: The modification string.
    :type mod: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round to.
    :type precision: int

    :raises UnknownModificationError: If the modification is unknown.
    :raises InvalidDeltaMassError: If the modification contains an invalid delta mass.
    :raises UnknownElementError: If the modification contains an unknown element.
    :raises UnknownGlycanError: If the modification contains an unknown glycan.
    :raises NotImplementedError: If the modification is not implemented.
    :raises InvalidModificationMassError: If the modification cannot be parsed.

    :return: The parsed modification.
    :rtype: float

    .. code-block:: python

        >>> mod_mass('Acetyl|INFO:newly discovered', precision=3)
        42.011

        >>> mod_mass(Mod('Acetyl|INFO:newly discovered', 2), precision=3)
        84.022

        >>> mod_mass('1', precision=3)
        1

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
        peptacular.errors.InvalidModificationMassError: Cannot determine mass for modification: info:HelloWorld

    """

    if isinstance(mod, Mod):
        return mod_mass(mod.val, monoisotopic, precision)*mod.mult

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
    :param precision: The precision of the mass.
    :type precision: int

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
    except ValueError:
        raise InvalidDeltaMassError(obs_str)


def _parse_glycan_mass_from_proforma_str(glycan_str: str, monoisotopic: bool, precision: Optional[int] = None) -> float:
    """
    Parse a glycan string and return its mass.

    :param glycan_str: The glycan string to parse.
    :type glycan_str: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass.
    :type precision: int

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
        else:
            return round_func(entry.avg_mass)

    else:  # Try to parse glycan formula
        try:
            return round_func(glycan_mass(glycan_str, monoisotopic, precision))
        except InvalidGlycanFormulaError as e:
            raise InvalidGlycanFormulaError(glycan_str, e.msg) from e


def _parse_chem_mass_from_proforma_str(chem_str: str, monoisotopic: bool, precision: Optional[int] = None) -> float:
    """
    Parse a chemical formula string and return its mass.

    :param chem_str: The chemical formula string to parse.
    :type chem_str: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass.
    :type precision: int

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
    except InvalidChemFormulaError as e:
        raise InvalidChemFormulaError(chem_str, e.msg) from e


def _parse_mod_mass(mod: str, monoisotopic: bool = True, precision: Optional[int] = None) -> Union[float, None]:
    """
    Parse a modification and return its mass.

    :param mod: The modification to parse.
    :type mod: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass.
    :type precision: int

    :raises UnknownModificationError: If the modification is unknown.
    :raises InvalidDeltaMassError: If the modification contains an invalid delta mass.
    :raises UnknownElementError: If the modification contains an unknown element.
    :raises UnknownGlycanError: If the modification contains an unknown glycan.
    :raises NotImplementedError: If the modification is not implemented.

    :return: The mass of the modification.
    :rtype: float

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

    """

    # localization fix
    if isinstance(mod, str) and '#' in mod:
        if mod.startswith('#'):  # for localized positions, return 0
            return 0.0
        else:  # for only the original declaration consider the mass modification
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

    elif is_gno_str(mod):
        return parse_gno_mass(mod, monoisotopic, precision)

    elif is_xlmod_str(mod):
        # not implemented
        return parse_xlmod_mass(mod, monoisotopic, precision)

    elif is_resid_str(mod):
        # not implemented
        return parse_resid_mass(mod, monoisotopic, precision)

    elif mod_lower.startswith('info:'):  # Skip info modifications
        return None

    elif is_psi_mod_str(mod):  # is a psi-mod modification
        # raises UnknownModificationError and InvalidDeltaMassError
        return parse_psi_mass(mod, monoisotopic, precision)

    elif is_unimod_str(mod):  # is a unimod modification
        # raises UnknownModificationError and InvalidDeltaMassError
        return parse_unimod_mass(mod, monoisotopic, precision)

    # chemical formula
    elif mod_lower.startswith('formula:'):
        # raises UnknownElementError
        return _parse_chem_mass_from_proforma_str(mod, monoisotopic, precision)

    # observed mass
    elif mod_lower.startswith('obs:'):
        # raises InvalidDeltaMassError
        return _parse_obs_mass_from_proforma_str(mod, precision)


def _pop_delta_mass_mods(annotation: ProFormaAnnotation) -> float:
    """
    Pop the delta mass modifications from the modifications' dictionary. This leaves only modifications which
    can have their elemental composition calculated.

    :param mods: The modifications' dictionary.
    :type mods: dict

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

