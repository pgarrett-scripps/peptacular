"""
mass.py is a simple module for computing the m/z and mass of an amino acid sequence.
"""
from __future__ import annotations

from typing import Dict

from peptacular.chem import parse_chem_formula, _get_sequence_comp, _parse_mod_delta_mass_only, estimate_comp
from peptacular import constants
from peptacular.sequence.proforma import parse_static_mods, ProFormaAnnotation, Mod
from peptacular.util import convert_type
from peptacular.errors import UnknownElementError, UnknownGlycanError, InvalidDeltaMassError, \
    InvalidModificationMassError, UnknownAminoAcidError
from peptacular.glycan import parse_glycan_formula
from peptacular.mod_db import parse_psi_mass, parse_unimod_mass, is_unimod_str, is_psi_mod_str
from peptacular.sequence.sequence import parse_single_sequence
from peptacular.types import Chem_Composition


def comp_mass(sequence: str | ProFormaAnnotation, ion_type: str = 'p') -> (Dict[str, int], float):
    """
    Get the chemical composition of a peptide sequence and its delta mass.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param ion_type: The ion type.
    :type ion_type: str

    :return: The chemical composition and delta mass of the peptide sequence.
    :rtype: tuple

    .. code-block:: python

        # Get the chemical composition of a peptide sequence.
        >>> comp_mass('PEPTIDE')
        ({'C': 34, 'H': 53, 'N': 7, 'O': 15}, 0.0)

        >>> comp_mass('PEPTIDE[1.0]')
        ({'C': 34, 'H': 53, 'N': 7, 'O': 15}, 1.0)

    """

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    annotation = annotation.condense_static_mods(inplace=False)
    delta_mass = _pop_delta_mass_mods(annotation)  # get delta mass and pop these mods
    peptide_composition = _get_sequence_comp(annotation, ion_type)
    return peptide_composition, delta_mass


def comp(sequence: str | ProFormaAnnotation, ion_type: str = 'p', estimate_delta: bool = False) -> Dict[str, int]:
    """
    Calculates the elemental composition of a peptide sequence, including modifications,
    and optionally estimates the composition based on the delta mass from modifications.

    :param sequence: The peptide sequence to analyze.
    :param ion_type: The type of ion ('p' for precursor, or other types indicating different ionization forms).
    :param estimate_delta: If True, estimate the composition based on the delta mass from modifications.

    :raises ValueError: If delta_mass is nonzero and estimate_delta is False, indicating an unaccounted modification.

    :return: A dictionary with elements as keys and their counts as values.
    :raises ValueError: If delta_mass is nonzero and estimate_delta is False, indicating an unaccounted modification.

    .. code-block:: python

        # Calculate the elemental composition of a peptide sequence.
        >>> comp('PEPTIDE')
        {'C': 34, 'H': 53, 'N': 7, 'O': 15}

        >>> comp('PEPTIDE[1.0]', estimate_delta=True)
        {'C': 34.04446833455479, 'H': 53.06986041632441, 'N': 7.012225550345263, 'O': 15.01330250093913, 'S': 0.0003754919712730832}

        >>> comp('PEPTIDE[1.0]', estimate_delta=False)
        Traceback (most recent call last):
        ValueError: Non-zero delta mass (1.0) encountered without estimation enabled for sequence 'PEPTIDE[1.0]'.

    """

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    composition, delta_mass = comp_mass(annotation, ion_type)

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


def mass(sequence: str | ProFormaAnnotation,
         charge: int = 0,
         ion_type: str = 'p',
         monoisotopic: bool = True,
         isotope: int = 0,
         loss: float = 0.0,
         precision: int | None = None) -> float:
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

        >>> mass('PEPTIDEs', charge=0, precision=3, ion_type='p', isotope=1)
        Traceback (most recent call last):
        peptacular.errors.ProFormaFormatError: Invalid ProForma format: Expected either '[', '(', '?', '-', or '/' but got: s at index: 7

    """

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    if 'B' in annotation.sequence or 'J' in annotation.sequence:
        raise ValueError('Invalid amino acid sequence: B and J are ambiguous')

    if annotation.isotope_mods is not None and len(annotation.isotope_mods) > 0:
        peptide_composition, delta_mass = comp_mass(annotation, ion_type)
        return chem_mass(peptide_composition, monoisotopic=monoisotopic, precision=precision) + delta_mass

    if ion_type == 'p':
        ion_type = 'y'
    else:
        _ = annotation.pop_labile_mods()

    m = 0.0
    if annotation.static_mods is not None:
        static_map = parse_static_mods(annotation.static_mods)

        for aa, mod in static_map.items():
            aa_count = annotation.sequence.count(aa)
            mod_mass = sum(_parse_mod_mass(m, monoisotopic, precision=None) for m in mod) * aa_count
            m += mod_mass

        annotation.pop_static_mods()

    try:
        m += sum(constants.MONOISOTOPIC_AA_MASSES[aa] if monoisotopic
                 else constants.AVERAGE_AA_MASSES[aa] for aa in annotation.sequence)
    except KeyError as e:
        raise UnknownAminoAcidError(e)

    # Apply labile mods
    if annotation.labile_mods is not None:
        for mod in annotation.labile_mods:
            m += _parse_mod_mass(mod)

    # Apply Unknown mods
    if annotation.unknown_mods is not None:
        for mod in annotation.unknown_mods:
            m += _parse_mod_mass(mod)

    # Apply N-term mods
    if annotation.nterm_mods is not None:
        for mod in annotation.nterm_mods:
            m += _parse_mod_mass(mod)

    # Apply intervals
    if annotation.intervals is not None:
        for interval in annotation.intervals:
            if interval.mods is not None:
                for mod in interval.mods:
                    m += _parse_mod_mass(mod)

    # apply internal mods
    if annotation.internal_mods is not None:
        for k, mods in annotation.internal_mods.items():
            for mod in mods:
                m += _parse_mod_mass(mod)

    # apply C-term mods
    if annotation.cterm_mods is not None:
        for mod in annotation.cterm_mods:
            m += _parse_mod_mass(mod)

    m += (charge * constants.PROTON_MASS)  # Add charge
    m += isotope * constants.NEUTRON_MASS + loss  # Add isotope and loss
    m += constants.MONOISOTOPIC_ION_ADJUSTMENTS[ion_type] if monoisotopic \
        else constants.AVERAGE_ION_ADJUSTMENTS[ion_type]

    if precision is not None:
        m = round(m, precision)

    return m


def mz(sequence: str | ProFormaAnnotation,
       charge: int = 0,
       ion_type: str = 'p',
       monoisotopic: bool = True,
       isotope: int = 0,
       loss: float = 0.0,
       precision: int | None = None) -> float:
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

    m = mass(sequence=sequence, charge=charge, ion_type=ion_type,
             monoisotopic=monoisotopic, isotope=isotope, loss=loss, precision=None)

    m = m if charge == 0 else m / charge

    if precision is not None:
        m = round(m, precision)

    return m


def chem_mass(formula: Chem_Composition | str, monoisotopic: bool = True, precision: int | None = None) -> float:
    """
    Calculate the mass of a chemical formula.

    :param formula: The chemical formula.
    :type formula: dict
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to.
    :type precision: int

    :raises UnknownElementError: If the chemical formula contains an unknown element.

    :return: The mass of the chemical formula.
    :rtype: float

    .. code-block:: python

        # Calculate the mass of a chemical formula.
        >>> chem_mass({'C': 6, 'H': 12, 'O': 6}, precision=3)
        180.063

        >>> chem_mass({'13C': 6, 'H': 12, 'O': 6}, precision=3)
        186.084

        >>> chem_mass({'C': 6, 'H': 12, 'O': 6}, monoisotopic=False, precision=3)
        180.156

        # Use average masses for all elements except for iosotopes.
        >>> chem_mass({'13C': 6, 'H': 12, 'O': 6}, monoisotopic=False, precision=3)
        186.112

        # Example Error
        >>> chem_mass({'C': 6, 'H': 12, 'O': 6, 'X': 1}, precision=3)
        Traceback (most recent call last):
        peptacular.errors.UnknownElementError: Unknown element: X

    """

    if isinstance(formula, str):
        formula = parse_chem_formula(formula)

    m = 0.0
    for element, count in formula.items():

        if element not in constants.ISOTOPIC_ATOMIC_MASSES:
            raise UnknownElementError(element)

        if monoisotopic is True:
            m += constants.ISOTOPIC_ATOMIC_MASSES[element] * count
        else:
            if any(char.isdigit() for char in element):  # element is a isotope
                m += constants.ISOTOPIC_ATOMIC_MASSES[element] * count
            else:
                m += constants.AVERAGE_ATOMIC_MASSES[element] * count

    if precision is not None:
        m = round(m, precision)

    return m


def glycan_mass(formula: Chem_Composition | str,
                monoisotopic: bool = True,
                precision: int | None = None) -> float:
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
        >>> glycan_mass({'HesNAc':1}, precision=3)
        Traceback (most recent call last):
        peptacular.errors.UnknownGlycanError: Unknown glycan: HesNAc

    """

    if isinstance(formula, str):
        formula = parse_glycan_formula(formula)

    mass_table = constants.MONOSACCHARIDE_ID_TO_ISOTOPIC_MASSES if monoisotopic \
        else constants.MONOSACCHARIDE_ID_TO_AVERAGE_MASSES

    m = 0.0
    for monosaccharide, count in formula.items():

        if monosaccharide not in constants.MONOSACCHARIDE_NAME_TO_ID:
            raise UnknownGlycanError(monosaccharide)

        monosaccharide_id = constants.MONOSACCHARIDE_NAME_TO_ID[monosaccharide]
        m += mass_table[monosaccharide_id] * count

    if precision is not None:
        m = round(m, precision)

    return m


def parse_mod_mass(mod: str | Mod, monoisotopic: bool = True, precision: int | None = None) -> float:
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

        >>> parse_mod_mass('Acetyl|INFO:newly discovered', precision=3)
        42.011

        >>> parse_mod_mass('1', precision=3)
        1

        >>> parse_mod_mass('Acetyl|Obs:+42.010565', precision=3)
        42.011

        # example error
        >>> parse_mod_mass('Acetdsyl|Obs:42d.010565', precision=3)
        Traceback (most recent call last):
        ...
        peptacular.errors.InvalidDeltaMassError: Invalid delta mass: 42d.010565

        # example invalid mass
        >>> parse_mod_mass('info:HelloWorld', precision=3)
        Traceback (most recent call last):
        ...
        peptacular.errors.InvalidModificationMassError: Cannot determine mass for modification: info:HelloWorld

    """

    if isinstance(mod, Mod):
        return parse_mod_mass(mod.val, monoisotopic, precision)

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


def _parse_obs_mass_from_proforma_str(obs_str: str, precision: int | None = None) -> float:
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


def _parse_glycan_mass_from_proforma_str(glycan_str: str, monoisotopic: bool, precision: int | None = None) -> float:
    """
    Parse a glycan string and return its mass.

    :param glycan_str: The glycan string to parse.
    :type glycan_str: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass.
    :type precision: int

    :raises UnknownGlycanError: If the glycan string contains an unknown monosaccharide.

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
        >>> _parse_glycan_mass_from_proforma_str('6BAAXE121B1', monoisotopic=True, precision=3)
        Traceback (most recent call last):
        peptacular.errors.UnknownGlycanError: Unknown glycan: 6BAAXE121B1

        # Invalid glycan str
        >>> _parse_glycan_mass_from_proforma_str('HeSNAc2Hex3Neu1', monoisotopic=True, precision=3)
        Traceback (most recent call last):
        peptacular.errors.UnknownGlycanError: Unknown glycan: HeSNAc2Hex3Neu1

    """

    round_func = lambda x: round(x, precision) if precision is not None else x

    if glycan_str.lower().startswith('glycan:'):
        glycan_str = ''.join(glycan_str.split(':')[1:])

    glycan_id = constants.MONOSACCHARIDE_NAME_TO_ID.get(glycan_str, glycan_str)  # get id if possible
    if glycan_id in constants.MONOSACCHARIDE_ID_TO_ISOTOPIC_MASSES:  # Try to parse glycan ID

        if monoisotopic:
            return round_func(constants.MONOSACCHARIDE_ID_TO_ISOTOPIC_MASSES[glycan_id])
        else:
            return round_func(constants.MONOSACCHARIDE_ID_TO_AVERAGE_MASSES[glycan_id])
    else:  # Try to parse glycan formula
        try:
            return round_func(glycan_mass(glycan_str, monoisotopic, precision))
        except UnknownGlycanError as e:
            raise e


def _parse_chem_mass_from_proforma_str(chem_str: str, monoisotopic: bool, precision: int | None = None) -> float:
    """
    Parse a chemical formula string and return its mass.

    :param chem_str: The chemical formula string to parse.
    :type chem_str: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass.
    :type precision: int

    :raises UnknownElementError: If the chemical formula contains an unknown element.

    :return: The mass of the chemical formula.
    :rtype: float

    .. code-block:: python

        # Monoisotopic Mass
        >>> _parse_chem_mass_from_proforma_str('C2H4O2', monoisotopic=True, precision=3)
        60.021

        # Average Mass
        >>> _parse_chem_mass_from_proforma_str('C2H4O2', monoisotopic=False, precision=3)
        60.052

    """

    if chem_str.lower().startswith('formula:'):
        chem_str = ''.join(chem_str.split(':')[1:])

    comps = parse_chem_formula(chem_str)
    try:
        return chem_mass(comps, monoisotopic, precision)
    except UnknownElementError as e:
        raise e


def _parse_mod_mass(mod: str | Mod, monoisotopic: bool = True, precision: int | None = None) -> float | None:
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

    if isinstance(mod, Mod):
        return _parse_mod_mass(mod.val, monoisotopic, precision)

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

    elif mod_lower.startswith('gno:') or mod_lower.startswith('g:'):
        # not implemented
        raise NotImplementedError('GNO modifications are not implemented')

    elif mod_lower.startswith('xlmod:') or mod_lower.startswith('x:') or mod_lower.startswith('xl-mod:'):
        # not implemented
        raise NotImplementedError('XL-MOD modifications are not implemented')

    elif mod_lower.startswith('resid:') or mod_lower.startswith('r:'):
        # not implemented
        raise NotImplementedError('RESID modifications are not implemented')

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

    # TODO: Support for charge and adducts
    annotation.clear_empty_mods()

    return delta_mass
