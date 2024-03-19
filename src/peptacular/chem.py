"""
chem.py contains functions for parsing and writing chemical formulas, and for calculating the composition of a sequence
with/and modifications.
"""
from __future__ import annotations

from typing import Dict, Union, List, Optional

from peptacular.chem_util import parse_chem_formula, _parse_chem_comp_from_proforma_str, write_chem_formula
from peptacular.mods.mod_setup import MONOSACCHARIDES_DB
from peptacular.sequence.sequence import sequence_to_annotation
from peptacular.constants import AA_COMPOSITIONS, ION_TYPE_COMPOSITION_ADJUSTMENTS, ISOTOPIC_AVERAGINE_MASS, AVERAGINE_RATIOS
from peptacular.errors import InvalidCompositionError, AmbiguousAminoAcidError, \
    UnknownAminoAcidError, DeltaMassCompositionError
from peptacular.glycan import glycan_comp
from peptacular.mods.mod_db import parse_unimod_comp, parse_psi_comp, is_psi_mod_str, is_unimod_str, parse_xlmod_comp, \
    parse_resid_comp, is_resid_str, is_xlmod_str, is_gno_str, parse_gno_comp
from peptacular.sequence.proforma import parse_static_mods, ProFormaAnnotation, Mod, parse_isotope_mods, \
    _parse_ion_elements
from peptacular.util import convert_type


def convert_glycan_formula_to_chem_formula(glycan: Dict[str, int] | str) -> str:
    """
    Converts a glycan dictionary to a chemical formula.

    :param glycan: A dictionary containing the glycan components and their counts, or a glycan formula string.
    :type glycan: dict | str

    :return: A chemical formula string.
    :rtype: str

    .. code-block:: python

            >>> convert_glycan_formula_to_chem_formula({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1})
            'C45H73N3O34'

    """

    return write_chem_formula(glycan_comp(glycan))


def mod_comp(mod: str | Mod) -> Dict[str, int]:
    """
    Parse a modification string.

    :param mod: The modification string.
    :type mod: str

    :raises InvalidCompositionError: If the modification string is invalid.

    :return: The parsed composition.
    :rtype: dict

    .. code-block:: python

        >>> mod_comp('Acetyl|INFO:newly discovered')
        {'H': 2, 'C': 2, 'O': 1}

        >>> mod_comp('Acetyl|Obs:+42.010565')
        {'H': 2, 'C': 2, 'O': 1}

    """

    if isinstance(mod, Mod):
        return mod_comp(mod.val)

    if isinstance(mod, (float, int)):
        raise InvalidCompositionError(mod)

    mods = mod.split('|')
    for m in mods:
        m = _parse_modification_composition(m)
        if m is not None:
            return m

    raise InvalidCompositionError(mod)


def estimate_comp(neutral_mass: float,
                  isotopic_mods: Optional[List[str] | List[Mod]] = None) -> Dict[str, float]:
    """
    Estimate the number of each element in a molecule based on its molecular mass using the averagine model.

    :param neutral_mass: The total neutral mass of the molecule.
    :type neutral_mass: float
    :param isotopic_mods: The isotopic modifications.
    :type isotopic_mods: list

    :return: The estimated number of each element in the molecule.
    :rtype: Dict[str, float]

    .. python::

        # Example usage
        >>> estimate_comp(1000)['C']
        44.468334554796016

        >>> estimate_comp(1000, ['13C'])['13C']
        44.468334554796016

    """

    composition = {atom: ratio * neutral_mass / ISOTOPIC_AVERAGINE_MASS for atom, ratio in
                   AVERAGINE_RATIOS.items()}

    if isotopic_mods:
        composition = _apply_isotope_mods_to_composition(composition, isotopic_mods)

    return composition


def _parse_glycan_comp_from_proforma_str(glycan_str: str) -> Dict[str, int]:
    """
    Parse a glycan string and return its mass.

    :param glycan_str: The glycan string to parse.
    :type glycan_str: str

    :raises UnknownGlycanError: If the glycan string contains an unknown monosaccharide.
    :raises InvalidFormulaError: If the formula is invalid.

    :return: The mass of the glycan.
    :rtype: float

    .. code-block:: python

        #  Get Composition
        >>> _parse_glycan_comp_from_proforma_str('HexNAc2Hex3Neu1')
        {'C': 43, 'H': 71, 'N': 3, 'O': 32}

        # Using a glycan name
        >>> _parse_glycan_comp_from_proforma_str('HexNAc')
        {'C': 8, 'H': 13, 'N': 1, 'O': 5}

        # Using a glycan ID
        >>> _parse_glycan_comp_from_proforma_str('6BAAE1B1')
        {'C': 3, 'H': 4, 'O': 2}

        # Using a glycan ID
        >>> _parse_glycan_comp_from_proforma_str('Glycan:6BAAE1B1')
        {'C': 3, 'H': 4, 'O': 2}

        # Invalid glycan ID
        >>> _parse_glycan_comp_from_proforma_str('6BAAXE121B1')
        Traceback (most recent call last):
        peptacular.errors.UnknownGlycanError: Unknown glycan: 6BAAXE121B1

        # Invalid glycan str
        >>> _parse_glycan_comp_from_proforma_str('HeSNAc2Hex3Neu1')
        Traceback (most recent call last):
        peptacular.errors.UnknownGlycanError: Unknown glycan: HeSNAc2Hex3Neu1

    """

    if glycan_str.lower().startswith('glycan:'):
        glycan_str = ''.join(glycan_str.split(':')[1:])

    if MONOSACCHARIDES_DB.contains_id(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_id(glycan_str)
    elif MONOSACCHARIDES_DB.contains_name(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_name(glycan_str)
    elif MONOSACCHARIDES_DB.contains_synonym(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_synonym(glycan_str)
    else:
        return parse_chem_formula(convert_glycan_formula_to_chem_formula(glycan_str))

    return parse_chem_formula(entry.composition)






def _parse_modification_composition(mod: str) -> Union[None, Dict[str, int]]:
    """
    Parses a modification composition and returns a dictionary with the element and their counts.

    :param mod: The modification composition.
    :type mod: str

    :raises InvalidCompositionError: If the modification composition string is invalid.
    :raises UnknownGlycanError: If the glycan formula contains an unknown glycan.
    :raises InvalidFormulaError: If the formula is invalid.
    :raises UnknownModificationError: If the modification is unknown.
    :raises NotImplementedError: If the modification is not implemented.

    :return: A dictionary with the element and their counts.
    :rtype: dict

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> _parse_modification_composition('U:2')
        {'H': 1, 'N': 1, 'O': -1}

        >>> _parse_modification_composition('Formula:[13C4]H12')
        {'13C': 4, 'H': 12}

        >>> _parse_modification_composition('Glycan:HexNAc2Hex3Neu5Gc1')
        {'C': 45, 'H': 73, 'N': 3, 'O': 34}

        >>> _parse_modification_composition('1')

    """

    if isinstance(mod, (int, float)):  # cannot get composition from a delta mass
        return None

    converted_mod = convert_type(mod)

    if isinstance(converted_mod, (int, float)):  # cannot get composition from a delta mass
        return None

    # localization fix
    if isinstance(mod, str) and '#' in mod:
        if mod.startswith('#'):  # for localized positions, return empty comp
            return {}
        else:  # for only the original declaration consider the mass modification
            mod = mod.split('#')[0]

    mod_lower = mod.lower()
    if mod_lower.startswith('glycan:'):
        return _parse_glycan_comp_from_proforma_str(mod)

    elif is_gno_str(mod):
        # not implemented
        return parse_chem_formula(parse_gno_comp(mod))

    elif is_xlmod_str(mod):
        # not implemented
        return parse_chem_formula(parse_xlmod_comp(mod))

    elif is_resid_str(mod):
        # not implemented
        return parse_chem_formula(parse_resid_comp(mod))

    elif mod_lower.startswith('info:'):  # cannot get comp for info
        pass

    elif mod_lower.startswith('obs:'):  # cannot get comp for observed mass
        pass

    elif is_psi_mod_str(mod):  # psi-mod
        return parse_chem_formula(parse_psi_comp(mod))

    elif is_unimod_str(mod):  # unimod
        return parse_chem_formula(parse_unimod_comp(mod))

    # chemical formula
    elif mod_lower.startswith('formula:'):
        return _parse_chem_comp_from_proforma_str(mod)


def _parse_delta_mass(delta_mass: str) -> float:
    """
    Parse the delta mass.

    :param delta_mass: The delta mass.
    :type delta_mass: str
    :return: The parsed delta mass.
    :rtype: float

    .. code-block:: python

        >>> _parse_delta_mass('42.0')
        42.0

        >>> _parse_delta_mass('Acetyl')

        >>> _parse_delta_mass('UniMod:1')

        >>> _parse_delta_mass('UniMod:+1')
        1.0

        >>> _parse_delta_mass('Obs:42.0')
        42.0
    """

    if isinstance(delta_mass, (int, float)):
        return delta_mass

    mass = None
    for mod_str in delta_mass.split('|'):

        # check if the mod contains a localization score
        if '#' in mod_str:
            mod_str = mod_str.split('#')[0]

        try:
            mass = float(mod_str)
            break
        except ValueError:
            pass

        mod_lower = mod_str.lower()
        if mod_lower.startswith('unimod:') or mod_lower.startswith('u:') or \
                mod_lower.startswith('psi-mod:') or mod_lower.startswith('mod:') or mod_lower.startswith('m:') or \
                mod_lower.startswith('resid:') or mod_lower.startswith('r:') or \
                mod_lower.startswith('xlmod:') or mod_lower.startswith('x:') or \
                mod_lower.startswith('gno:') or mod_lower.startswith('g:'):

            mod_str = mod_str.split(':')[1]

            if mod_str.startswith('+') or mod_str.startswith('-'):
                mass = float(mod_str)
                break

        if mod_str.lower().startswith('obs:'):
            mod_str = mod_str.split(':')[1]
            try:
                mass = float(mod_str)
                break
            except ValueError:
                pass

    return mass


def _get_sequence_comp(sequence: str | ProFormaAnnotation, ion_type: str, charge: int = None,
                       charge_adducts: str = None) -> Dict[str, int]:
    """
    Calculate the composition of a sequence.

    :param sequence: The sequence.
    :type sequence: str
    :param ion_type: The ion type.
    :type ion_type: str
    :param charge: The charge of the ion.
    :type charge: int
    :param charge_adducts: The charge adducts.
    :type charge_adducts: str

    :return: A dictionary with the element and their counts.
    :rtype: dict

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> _get_sequence_comp('PEPTIDE', 'y')
        {'C': 34, 'H': 53, 'N': 7, 'O': 15}

        >>> _get_sequence_comp('PEPTIDE', 'b')
        {'C': 34, 'H': 51, 'N': 7, 'O': 14}

        >>> _get_sequence_comp('<H>PEPTIDE', 'b')
        {'C': 34, 'H': 51, 'N': 7, 'O': 14}

        >>> _get_sequence_comp('{Unimod:2}PEPTIDE', 'p')
        {'C': 34, 'H': 54, 'N': 8, 'O': 14}

        >>> _get_sequence_comp('<13C>PEPTIDE', 'b')
        {'H': 51, 'N': 7, 'O': 14, '13C': 34}

        >>> _get_sequence_comp('PEPTIDE[Unimod:2]', 'y')
        {'C': 34, 'H': 54, 'N': 8, 'O': 14}

        >>> _get_sequence_comp('<[Unimod:2]@T>PEPTIDE', 'y')
        {'C': 34, 'H': 54, 'N': 8, 'O': 14}

        >>> _get_sequence_comp('<[Unimod:2]@N-Term>PEPTIDE', 'y')
        {'C': 34, 'H': 54, 'N': 8, 'O': 14}

        >>> _get_sequence_comp('PEPTIDE/2', 'y')
        {'C': 34, 'H': 55, 'N': 7, 'O': 15, 'e': -2}

        >>> _get_sequence_comp('<13C>PEPTIDE[Unimod:213413]', 'b')
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: Unimod:213413

        >>> _get_sequence_comp('I', 'by')
        {'C': 6, 'H': 12, 'N': 1, 'O': 1}

        # Ambiguous amino acid
        >>> _get_sequence_comp('B', 'by')
        Traceback (most recent call last):
        peptacular.errors.AmbiguousAminoAcidError: Ambiguous amino acid: B

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

    if ion_type == 'p':
        ion_type = 'y'
    else:
        _ = annotation.pop_labile_mods()

    if 'B' in annotation.sequence:
        raise AmbiguousAminoAcidError('B')

    if 'J' in annotation.sequence:
        raise AmbiguousAminoAcidError('J')

    # Get the composition of the base sequence
    composition = {}
    for aa in annotation.sequence:
        try:
            aa_comp = AA_COMPOSITIONS[aa]
        except KeyError:
            raise UnknownAminoAcidError(aa)
        for k, v in aa_comp.items():
            composition[k] = composition.get(k, 0) + v

    # Apply the adjustments for the ion type
    for k, v in ION_TYPE_COMPOSITION_ADJUSTMENTS[ion_type].items():
        composition[k] = composition.get(k, 0) + v

    # Pop isotopic mods
    isotopic_mods = annotation.pop_isotope_mods()
    static_mods = annotation.pop_static_mods()

    mod_compositions = [composition]

    # Loop over unknown mods and apply them
    if annotation.has_unknown_mods():
        for unknown_mod in annotation.unknown_mods:
            mod_compositions.append(mod_comp(unknown_mod))

    # loop over intervals and apply them
    if annotation.has_intervals():
        for interval in annotation.intervals:
            if interval.has_mods():
                for interval_mod in interval.mods:
                    mod_compositions.append(mod_comp(interval_mod))

    # apply labile mods
    if annotation.has_labile_mods():
        for labile_mod in annotation.labile_mods:
            mod_compositions.append(mod_comp(labile_mod))

    # apply N-term mods
    if annotation.has_nterm_mods():
        for nterm_mod in annotation.nterm_mods:
            mod_compositions.append(mod_comp(nterm_mod))

    # apply C-term mods
    if annotation.has_cterm_mods():
        for cterm_mod in annotation.cterm_mods:
            mod_compositions.append(mod_comp(cterm_mod))

    # apply internal mods
    if annotation.has_internal_mods():
        for k, internal_mods in annotation.internal_mods.items():
            for internal_mod in internal_mods:
                mod_compositions.append(mod_comp(internal_mod))

    # Apply static mods
    if static_mods:
        static_map = parse_static_mods(static_mods)

        n_term_mod = static_map.get('N-Term')
        if n_term_mod is not None:
            for m in n_term_mod:
                mod_compositions.append(mod_comp(m.val))

        c_term_mod = static_map.get('C-Term')
        if c_term_mod is not None:
            for m in c_term_mod:
                mod_compositions.append(mod_comp(m.val))

        for aa, mod in static_map.items():
            if aa in ['N-Term', 'C-Term']:
                continue
            aa_count = annotation.sequence.count(aa)
            for m in mod:
                mod_composition = mod_comp(m.val)

                mod_composition = {k: v * aa_count for k, v in mod_composition.items()}
                mod_compositions.append(mod_composition)

    # Add the charge adducts
    if charge != 0:
        adduct_comp = _parse_charge_adducts_comp(charge_adducts)
        adduct_comp = {k: v * charge for k, v in adduct_comp.items()}
        mod_compositions.append(adduct_comp)

    # Merge the compositions
    composition = {}
    for comp in mod_compositions:
        for k, v in comp.items():
            composition[k] = composition.get(k, 0) + v

    # Apply isotopic mods
    if isotopic_mods:
        composition = _apply_isotope_mods_to_composition(composition, isotopic_mods)

    return composition


def _parse_mod_delta_mass_only(mod: str | Mod) -> Union[float, None]:
    """
    Parse a modification string.

    :param mod: The modification string.
    :type mod: str

    :raises ValueError: If the modification string is invalid.

    :return: The parsed composition.
    :rtype: dict

    .. code-block:: python

        >>> _parse_mod_delta_mass_only('Acetyl|INFO:newly discovered')

        >>> _parse_mod_delta_mass_only('Acetyl|Obs:+42.010565')

        >>> _parse_mod_delta_mass_only('Obs:+42.010565')
        42.010565

        >>> _parse_mod_delta_mass_only('R:+1|INFO:Fantastic')
        1.0

    """

    if isinstance(mod, float) or isinstance(mod, int):
        return mod

    if isinstance(mod, Mod):
        return _parse_mod_delta_mass_only(mod.val)

    mods = mod.split('|')
    for m in mods:
        try:
            m = _parse_modification_composition(m)
        except DeltaMassCompositionError:
            continue
        if m is not None:
            return None

    for m in mods:
        m = _parse_delta_mass(m)
        if m is not None:
            return m

    raise ValueError(f'Invalid modification: {mod}')


def _apply_isotope_mods_to_composition(composition: Dict[str, int], isotopic_mods: List[Mod]) -> Dict[str, int]:
    """
    Apply isotopic modifications to a composition.

    :param composition: The composition.
    :type composition: dict
    :param isotopic_mods: The isotopic modifications.
    :type isotopic_mods: list

    :return: The modified composition.
    :rtype: dict

    .. code-block:: python

        # Example usage
        >>> _apply_isotope_mods_to_composition({'C': 6, 'H': 12, 'O': 6}, ['13C'])
        {'H': 12, 'O': 6, '13C': 6}

        >>> _apply_isotope_mods_to_composition({'C': 6, 'H': 12, 'O': 6}, ['13C', '15N'])
        {'H': 12, 'O': 6, '13C': 6}

    """

    isotope_map = parse_isotope_mods(isotopic_mods)

    for element, isotope_label in isotope_map.items():
        if element in composition:

            if element == isotope_label:
                continue

            # Check if the isotope label already exists in the composition
            if isotope_label in composition:
                # Add the count of the original element to the isotope label count
                composition[isotope_label] += composition[element]
            else:
                # If the isotope label doesn't exist, create it with the count of the original element
                composition[isotope_label] = composition[element]

            del composition[element]

    return composition


def _parse_adduct_comp(adduct: str) -> Dict[str, int]:
    """
    Parse an adduct string and return its mass.

    :param adduct: The adduct string to parse.
    :type adduct: str

    :raises InvalidDeltaMassError: If the adduct contains an invalid delta mass.

    :return: The mass of the adduct.
    :rtype: float

    .. code-block:: python

        # Parse an adduct string.
        >>> _parse_adduct_comp('+Na+')
        {'Na': 1, 'e': -1}

        >>> _parse_adduct_comp('+2Na+')
        {'Na': 2, 'e': -1}

        >>> _parse_adduct_comp('H+')
        {'H': 1, 'e': -1}

    """

    comp = {}
    element_count, element_symbol, element_charge = _parse_ion_elements(adduct)
    comp[element_symbol] = element_count
    comp['e'] = -1 * element_charge

    return comp


def _parse_charge_adducts_comp(adducts: Mod | str) -> Dict[str, int]:
    """
    Parse the charge adducts and return their mass.

    :param adducts: The charge adducts to parse.
    :type adducts: str

    :return: The mass of the charge adducts.
    :rtype: float

    .. code-block:: python

        # Parse the charge adducts and return their mass.
        >>> _parse_charge_adducts_comp('+Na+,+H+')
        {'Na': 1, 'e': -2, 'H': 1}

        >>> _parse_charge_adducts_comp('+H+')
        {'H': 1, 'e': -1}

    """

    if isinstance(adducts, Mod):
        adducts = adducts.val

    adducts = adducts.split(',')

    comps = []
    for adduct in adducts:
        comps.append(_parse_adduct_comp(adduct))

    # Merge the compositions
    composition = {}
    for comp in comps:
        for k, v in comp.items():
            composition[k] = composition.get(k, 0) + v

    return composition