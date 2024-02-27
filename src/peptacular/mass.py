"""
mass.py is a simple module for computing the m/z and mass of an amino acid sequence.
"""

from typing import Union, List, Dict

from peptacular.chem import parse_chem_formula, calculate_sequence_composition
from peptacular.constants import AVERAGE_AA_MASSES, MONOISOTOPIC_AA_MASSES, \
    NEUTRON_MASS, PROTON_MASS, ISOTOPIC_ATOMIC_MASSES, AVERAGE_ATOMIC_MASSES, \
    UNIMOD_ID_TO_MONO_MASSES, UNIMOD_ID_TO_AVERAGE_MASSES, MONOSACCHARIDE_ID_TO_ISOTOPIC_MASSES, \
    MONOSACCHARIDE_ID_TO_AVERAGE_MASSES, MONOSACCHARIDE_NAME_TO_ID, UNIMOD_NAME_TO_ID, PSI_MOD_NAME_TO_ID, \
    PSI_MOD_ID_TO_ISOTOPIC_MASSES, PSI_MOD_ID_TO_AVERAGE_MASSES
from peptacular.glycan import parse_glycan_formula
from peptacular.sequence.global_mods import parse_isotope_modifications, parse_static_modifications
from peptacular.sequence.sequence import get_modifications, strip_modifications
from peptacular.util import validate_ion_type


def calculate_chem_mass(formula: Union[Dict[str, int], str], monoisotopic: bool = True, precision: int = None) -> float:
    """
    Calculate the mass of a chemical formula.

    :param formula: The chemical formula.
    :type formula: dict
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to.
    :type precision: int

    :return: The mass of the chemical formula.
    :rtype: float

    .. code-block:: python

        # Calculate the mass of a chemical formula.
        >>> calculate_chem_mass({'C': 6, 'H': 12, 'O': 6}, precision=3)
        180.063

        >>> calculate_chem_mass({'13C': 6, 'H': 12, 'O': 6}, precision=3)
        186.084

        >>> calculate_chem_mass({'C': 6, 'H': 12, 'O': 6}, monoisotopic=False, precision=3)
        180.156

        >>> calculate_chem_mass({'13C': 6, 'H': 12, 'O': 6}, monoisotopic=False, precision=3)
        186.112

    """

    if isinstance(formula, str):
        formula = parse_chem_formula(formula)

    mass = 0.0
    for element, count in formula.items():

        if monoisotopic is True:
            mass += ISOTOPIC_ATOMIC_MASSES[element] * count
        else:
            if any(char.isdigit() for char in element):  # element is a isotope
                mass += ISOTOPIC_ATOMIC_MASSES[element] * count
            else:
                mass += AVERAGE_ATOMIC_MASSES[element] * count

    if precision is not None:
        mass = round(mass, precision)

    return mass


def calculate_glycan_mass(formula: Union[Dict[str, int], str], monoisotopic: bool = True,
                          precision: int = None) -> float:
    """
    Calculate the mass of a glycan formula.

    :param formula: The glycan formula.
    :type formula: dict
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass.
    :type precision: int

    :return: The mass of the glycan formula.
    :rtype: float

    .. code-block:: python

        # Calculate the mass of a glycan formula.
        >>> calculate_glycan_mass({'HexNAc': 2, 'Hex': 3, 'Neu': 1}, precision=3)
        1141.402

        >>> calculate_glycan_mass('HexNAc2Hex3Neu1', precision=3)
        1141.402

        >>> calculate_glycan_mass({'HexNAc': 2, 'Hex': 3, 'Neu': 1}, monoisotopic=False, precision=3)
        1142.027

        >>> calculate_glycan_mass('HexNAc2Hex3Neu1', monoisotopic=False, precision=3)
        1142.027

    """

    if isinstance(formula, str):
        formula = parse_glycan_formula(formula)

    mass_table = MONOSACCHARIDE_ID_TO_ISOTOPIC_MASSES if monoisotopic else MONOSACCHARIDE_ID_TO_AVERAGE_MASSES

    mass = 0.0
    for monosaccharide, count in formula.items():
        monosaccharide_id = MONOSACCHARIDE_NAME_TO_ID[monosaccharide]
        mass += mass_table[monosaccharide_id] * count

    if precision is not None:
        mass = round(mass, precision)

    return mass


def _parse_mod_mass(mod: str, monoisotopic: bool = True, precision: int = None) -> Union[float, None]:
    """
    Parse a modification and return its mass.

    :param mod: The modification to parse.
    :type mod: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The precision of the mass.
    :type precision: int

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
        >>> _parse_mod_mass('C2')
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
        >>> _parse_mod_mass('HexNAc2Hex3Neu1')
        1141.4020671170001

    """

    try:
        return round(float(mod), precision) if precision is not None else float(mod)
    except ValueError:
        pass

    parsed_mod = None
    mod_lower = mod.lower()
    if mod_lower.startswith('glycan:'):
        glycan = mod.split(':')[1]
        if glycan in MONOSACCHARIDE_ID_TO_ISOTOPIC_MASSES:
            if monoisotopic:
                parsed_mod = MONOSACCHARIDE_ID_TO_ISOTOPIC_MASSES[glycan]
            else:
                parsed_mod = MONOSACCHARIDE_ID_TO_AVERAGE_MASSES[glycan]
        else:
            parsed_mod = calculate_glycan_mass(glycan, monoisotopic)

    elif mod_lower.startswith('gno:') or mod_lower.startswith('g:'):
        # not implemented
        raise NotImplementedError('GNO modifications are not implemented')

    elif mod_lower.startswith('xlmod:') or mod_lower.startswith('x:') or mod_lower.startswith('xl-mod:'):
        # not implemented
        raise NotImplementedError('XL-MOD modifications are not implemented')

    elif mod_lower.startswith('psi-mod:') or mod_lower.startswith('m:') or mod_lower.startswith('mod:'):
        psi_mod_id = mod.split(':')[1]
        if psi_mod_id in PSI_MOD_NAME_TO_ID:
            psi_mod_id = PSI_MOD_NAME_TO_ID[mod]

        if psi_mod_id in PSI_MOD_ID_TO_ISOTOPIC_MASSES:
            if monoisotopic:
                parsed_mod = PSI_MOD_ID_TO_ISOTOPIC_MASSES[psi_mod_id]
            else:
                parsed_mod = PSI_MOD_ID_TO_AVERAGE_MASSES[psi_mod_id]
        else:
            raise ValueError(f'PSI-MOD id {psi_mod_id} not found')

    elif mod_lower.startswith('resid:') or mod_lower.startswith('r:'):
        # not implemented
        raise NotImplementedError('RESID modifications are not implemented')

    # unimod
    elif mod_lower.startswith('unimod:') or mod_lower.startswith('u:'):
        unimod_id = mod.split(':')[1]

        if unimod_id in UNIMOD_NAME_TO_ID:
            unimod_id = UNIMOD_NAME_TO_ID[unimod_id]

        if unimod_id in UNIMOD_ID_TO_MONO_MASSES:
            if monoisotopic:
                parsed_mod = UNIMOD_ID_TO_MONO_MASSES[unimod_id]
            else:
                parsed_mod = UNIMOD_ID_TO_AVERAGE_MASSES[unimod_id]
        else:
            raise ValueError(f'Unimod id {unimod_id} not found')

    # just info
    elif mod_lower.startswith('info:'):
        pass

    # chemical formula
    elif mod_lower.startswith('formula:'):
        formula = mod.split(':')[1]
        comps = parse_chem_formula(formula)
        parsed_mod = calculate_chem_mass(comps, monoisotopic, precision)

    # observed mass
    elif mod_lower.startswith('obs:'):
        try:
            parsed_mod = float(mod.split(':')[1])
        except ValueError:
            raise ValueError(f'Invalid observed mass {mod}')

    # Free floating UNIMOD modification
    elif mod in UNIMOD_NAME_TO_ID:
        unimod_id = UNIMOD_NAME_TO_ID[mod]

        if monoisotopic:
            parsed_mod = UNIMOD_ID_TO_MONO_MASSES[unimod_id]
        else:
            parsed_mod = UNIMOD_ID_TO_AVERAGE_MASSES[unimod_id]

    # Free floating PSI modification
    elif mod in PSI_MOD_NAME_TO_ID:
        psi_mod_id = PSI_MOD_NAME_TO_ID[mod]
        if monoisotopic:
            parsed_mod = PSI_MOD_ID_TO_ISOTOPIC_MASSES[psi_mod_id]
        else:
            parsed_mod = PSI_MOD_ID_TO_AVERAGE_MASSES[psi_mod_id]

    else:

        # try to parse as a chemical formula
        try:
            parsed_mod = calculate_chem_mass(mod, monoisotopic, precision)
        except:
            pass

        # try to parse glycan formual
        try:
            parsed_mod = calculate_glycan_mass(mod, monoisotopic, precision)
        except:
            pass

    if precision is not None:
        return round(parsed_mod, precision)

    return parsed_mod


def parse_mod(mod: str, monoisotopic: bool = True, precision: int = None) -> float:
    """
    Parse a modification string.

    :param mod: The modification string.
    :type mod: str
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round to.
    :type precision: int

    :raises ValueError: If the modification string is invalid.

    :return: The parsed modification.
    :rtype: float

    .. code-block:: python

        >>> parse_mod('Acetyl|INFO:newly discovered', precision=3)
        42.011

        >>> parse_mod('Acetyl|Obs:+42.010565', precision=3)
        42.011

    """

    if isinstance(mod, (int, float)):
        if precision is not None:
            return round(mod, precision)
        return mod

    mods = mod.split('|')
    for m in mods:
        m = _parse_mod_mass(m, monoisotopic, precision)
        if m is not None:
            return m

    raise ValueError(f'Invalid modification: {mod}')


def calculate_mass(sequence: str, charge: int = 0, ion_type: str = 'y', monoisotopic: bool = True,
                   isotope: int = 0, loss: float = 0.0, aa_masses: Dict = None, precision: int = None) -> float:
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
    :param aa_masses: A dictionary of amino acid masses, defaults to [None].
    :type aa_masses: Dict
    :param precision: The precision of the mass, defaults to [None].
    :type precision: int

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: Mass of the peptide sequence.
    :rtype: float

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> calculate_mass('PEPTIDE', precision=3)
        799.36

        >>> calculate_mass('<13C>PEPTIDE', precision=3)
        833.474

        >>> calculate_mass('<[10]@T>PEPTIDE', precision=3)
        809.36

        >>> calculate_mass('<[Formula:[13C6]H20]@T>PEPTIDE', precision=3)
        897.537

        # Calculate the b-ion mass of a peptide sequence.
        >>> calculate_mass('PEPTIDE', ion_type='b', precision=3)
        781.349

        # Calulate the average mass of a peptide sequence.
        >>> calculate_mass('PEPTIDE', monoisotopic=False, precision=3)
        799.824

        # Calculate the mass of a peptide sequence with a charge of 2.
        >>> calculate_mass('PEPTIDE', charge=2, precision=3)
        801.375

        # Calcualte the mass of a modified peptide sequence.
        >>> calculate_mass('PE[3.14]PTIDE[Acetyl]', charge=2, precision=3)
        846.525

        # Calculate the mass of a peptide sequence with a charge of 2.
        >>> calculate_mass('PEPT[10][10]IDE', charge=2, precision=3)
        821.375

    """

    validate_ion_type(ion_type=ion_type)

    # Parse modifications and strip them from sequence
    mods = get_modifications(sequence=sequence)

    labile_mods = mods.pop('l', None)

    stripped_sequence = strip_modifications(sequence=sequence)

    peptide_composition = calculate_sequence_composition(stripped_sequence, ion_type)

    if 'i' in mods:
        isotope_map = parse_isotope_modifications(mods['i'])

        for element, isotope_label in isotope_map.items():
            if element in peptide_composition:
                # Check if the isotope label already exists in the composition
                if isotope_label in peptide_composition:
                    # Add the count of the original element to the isotope label count
                    peptide_composition[isotope_label] += peptide_composition[element]
                else:
                    # If the isotope label doesn't exist, create it with the count of the original element
                    peptide_composition[isotope_label] = peptide_composition[element]

                del peptide_composition[element]

        mods.pop('i')

    mass = 0.0
    if 's' in mods:
        static_map = parse_static_modifications(mods['s'])

        for aa, mod in static_map.items():
            aa_count = stripped_sequence.count(aa)
            mod_mass = _parse_mod_mass(mod, monoisotopic, precision=None)*aa_count
            mass += mod_mass

        mods.pop('s')

    mass += calculate_chem_mass(peptide_composition, monoisotopic=monoisotopic, precision=None)

    # Calculate mass
    #mass += sum(aa_masses[aa] for aa in stripped_sequence)  # Add amino acids
    mass += sum(parse_mod(mod) for mods in mods.values() for mod in mods)  # Add modifications
    mass += (charge * PROTON_MASS)  # Add charge
    mass += isotope * NEUTRON_MASS + loss  # Add isotope and loss
    #mass += MONOISOTOPIC_ION_ADJUSTMENTS[ion_type] if monoisotopic else AVERAGE_ION_ADJUSTMENTS[ion_type]

    if precision is not None:
        mass = round(mass, precision)

    return mass


def calculate_mz(sequence: str, charge: int = 0, ion_type: str = 'y', monoisotopic: bool = True, isotope: int = 0,
                 loss: Union[List[float], float] = 0.0, aa_masses: Dict = None, precision: int = None) -> float:
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
    :param aa_masses: A dictionary of amino acid masses, defaults to [None].
    :type aa_masses: Dict
    :param precision: The precision of the m/z, defaults to [None].
    :type precision: int

    :raise ValueError: If the ion type is not one of 'a', 'b', 'c', 'x', 'y', or 'z'.

    :return: m/z of the peptide sequence.
    :rtype: float

    .. code-block:: python

        # Calculate the m/z of a peptide sequence.
        >>> calculate_mz('PEPTIDE', charge = 1, precision = 3)
        800.367

        # Calculate the b-ion m/z of a peptide sequence.
        >>> calculate_mz('PEPTIDE', charge = 1, ion_type='b', precision = 3)
        782.357

        # Calulate the average m/z of a peptide sequence.
        >>> calculate_mz('PEPTIDE', charge = 1, monoisotopic=False, precision = 3)
        800.831

        # Calculate the m/z of a peptide sequence with a charge of 2.
        >>> calculate_mz('PEPTIDE', charge=2, precision = 3)
        400.687

        # Calcualte the m/z of a modified peptide sequence.
        >>> calculate_mz('PE[3.14]PTIDE-[80]', charge=2, precision = 3)
        442.257

    """

    mass = calculate_mass(sequence=sequence, charge=charge, ion_type=ion_type,
                          monoisotopic=monoisotopic, isotope=isotope, loss=loss,
                          aa_masses=aa_masses, precision=None)

    mass = mass if charge == 0 else mass / charge

    if precision is not None:
        mass = round(mass, precision)

    return mass
