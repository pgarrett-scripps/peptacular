"""
chem_calc.py contains functions for parsing and writing chemical formulas, and for calculating the composition of a
sequence with/and modifications.
"""
import warnings
from typing import Union, List, Optional

from peptacular.chem.chem_constants import ISOTOPIC_AVERAGINE_MASS
from peptacular.types import ChemComposition
from peptacular.proforma.input_convert import ModValue
from peptacular.chem.chem_util import parse_chem_formula, write_chem_formula
from peptacular.mods.mod_db_setup import MONOSACCHARIDES_DB
from peptacular.sequence.sequence_funcs import sequence_to_annotation
from peptacular.constants import (AA_COMPOSITIONS, AVERAGINE_RATIOS, NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS,
                                  FRAGMENT_ION_BASE_CHARGE_ADDUCTS)
from peptacular.errors import InvalidCompositionError, AmbiguousAminoAcidError, \
    UnknownAminoAcidError, DeltaMassCompositionError
from peptacular.glycan import glycan_comp
from peptacular.mods.mod_db import parse_unimod_comp, parse_psi_comp, is_psi_mod_str, is_unimod_str, parse_xlmod_comp, \
    parse_resid_comp, is_resid_str, is_xlmod_str, is_gno_str, parse_gno_comp
from peptacular.proforma.proforma_parser import parse_static_mods, ProFormaAnnotation, Mod, parse_isotope_mods, \
    parse_ion_elements
from peptacular.util import convert_type


def glycan_to_chem(glycan: Union[ChemComposition, str]) -> str:
    """
    Converts a glycan composition/formula to a chemical formula.

    :param glycan: A dictionary containing the glycan components and their counts, or a glycan formula string.
    :type glycan: Dict[str, int | float] | str

    :return: A chemical formula string.
    :rtype: str

    .. code-block:: python

            >>> glycan_to_chem({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1})
            'C45H73N3O34'

            >>> glycan_to_chem({'HexNAc': 2.3, 'Hex': 3.3, 'Neu5Gc': 1.1})
            'C50.3H81.6N3.4O37.9'

    """

    return write_chem_formula(glycan_comp(glycan))


def mod_comp(mod: ModValue) -> ChemComposition:
    """
    Parse a modification string.

    :param mod: The modification string.
    :type mod: str | int | float | Mod

    :raises InvalidCompositionError: If no composition can be retrieved from the modification string.

    :return: The composition of the modification.
    :rtype: Dict[str, int | float]

    .. code-block:: python

        >>> mod_comp('Acetyl|INFO:newly discovered')
        {'H': 2, 'C': 2, 'O': 1}

        >>> mod_comp('Acetyl|Obs:+42.010565')
        {'H': 2, 'C': 2, 'O': 1}

        >>> mod_comp(Mod('Acetyl|Obs:+42.010565', 2))
        {'H': 4, 'C': 4, 'O': 2}

    """

    if isinstance(mod, Mod):
        return {k: v * mod.mult for k, v in mod_comp(mod.val).items()}

    if isinstance(mod, (float, int)):
        raise InvalidCompositionError(mod)

    mods = mod.split('|')
    for m in mods:
        m = _parse_mod_comp(m)
        if m is not None:
            return m

    raise InvalidCompositionError(mod)


def estimate_comp(neutral_mass: float,
                  isotopic_mods: Optional[List[ModValue]] = None) -> ChemComposition:
    """
    Estimate the number of each element in a molecule based on its molecular mass using the averagine model.

    :param neutral_mass: The total neutral mass of the molecule.
    :type neutral_mass: float
    :param isotopic_mods: The isotopic modifications. Defaults to None.
    :type isotopic_mods: List[str | Mod] | None

    :return: The estimated composition.
    :rtype: Dict[str, int | float]

    . python::

        # Example usage
        >>> estimate_comp(1000)['C']
        44.468334554796016

        >>> estimate_comp(1000, ['13C'])['13C']
        44.468334554796016

    """

    composition = {atom: ratio * neutral_mass / ISOTOPIC_AVERAGINE_MASS for atom, ratio in
                   AVERAGINE_RATIOS.items()}

    if isotopic_mods is not None:
        composition = apply_isotope_mods_to_composition(composition, isotopic_mods)

    return composition


def _parse_glycan_comp(glycan_str: str) -> ChemComposition:
    """
    Parse a glycan string and return its mass.

    :param glycan_str: The glycan string to parse.
    :type glycan_str: str

    :raises UnknownGlycanError: If the glycan formula contains an unknown monosaccharide.
    :raises InvalidFormulaError: If the glycan formula is invalid.

    :return: The composition of the glycan.
    :rtype: Dict[str, int | float]

    .. code-block:: python

        #  Get Composition
        >>> _parse_glycan_comp('HexNAc2Hex3Neu1')
        {'C': 43, 'H': 71, 'N': 3, 'O': 32}

        # Float counts
        >>> _parse_glycan_comp('HexNAc2.3Hex3Neu1')
        {'C': 45.4, 'H': 74.9, 'N': 3.3, 'O': 33.5}

        # Using a glycan name
        >>> _parse_glycan_comp('HexNAc')
        {'C': 8, 'H': 13, 'N': 1, 'O': 5}

        # Using a glycan ID
        >>> _parse_glycan_comp('6BAAE1B1')
        {'C': 3, 'H': 4, 'O': 2}

        # Using a glycan ID
        >>> _parse_glycan_comp('Glycan:6BAAE1B1')
        {'C': 3, 'H': 4, 'O': 2}

        # Invalid glycan ID
        >>> _parse_glycan_comp('6BAAXE121B1') # doctest: +ELLIPSIS
        Traceback (most recent call last):
            ...
        peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "6BAAXE121B1". ...

        # Invalid glycan str
        _parse_glycan_comp('HeSNAc2Hex3Neu1') # doctest: +ELLIPSIS
        Traceback (most recent call last):
            ...
        peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "HeSNAc2Hex3Neu1". ...

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
        return parse_chem_formula(glycan_to_chem(glycan_str))

    return parse_chem_formula(entry.composition)


def _parse_mod_comp(mod: str) -> Union[ChemComposition, None]:
    """
    Parses a modification composition and returns a dictionary with the element and their counts.

    :param mod: The modification composition.
    :type mod: str

    :raises InvalidCompositionError: If the modification composition string is invalid.
    :raises UnknownGlycanError: If the glycan formula contains an unknown glycan.
    :raises InvalidFormulaError: If the formula is invalid.
    :raises UnknownModificationError: If the modification is unknown.
    :raises NotImplementedError: If the modification is not implemented.

    :return: The parsed composition, or None if the composition cannot be parsed.
    :rtype: Dict[str, int | float] | None

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> _parse_mod_comp('U:2')
        {'H': 1, 'N': 1, 'O': -1}

        >>> _parse_mod_comp('Formula:[13C4]H12')
        {'13C': 4, 'H': 12}

        >>> _parse_mod_comp('Glycan:HexNAc2Hex3Neu5Gc1')
        {'C': 45, 'H': 73, 'N': 3, 'O': 34}

        >>> _parse_mod_comp('1')

        >>> _parse_mod_comp('')

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
        mod = mod.split('#')[0]

    mod_lower = mod.lower()
    if mod_lower.startswith('glycan:'):
        return _parse_glycan_comp(mod)

    if is_gno_str(mod):
        # not implemented
        return parse_chem_formula(parse_gno_comp(mod))

    if is_xlmod_str(mod):
        # not implemented
        return parse_chem_formula(parse_xlmod_comp(mod))

    if is_resid_str(mod):
        # not implemented
        return parse_chem_formula(parse_resid_comp(mod))

    if mod_lower.startswith('info:'):  # cannot get comp for info
        return None

    if mod_lower.startswith('obs:'):  # cannot get comp for observed mass
        return None

    if is_psi_mod_str(mod):  # psi-mod
        return parse_chem_formula(parse_psi_comp(mod))

    if is_unimod_str(mod):  # unimod
        return parse_chem_formula(parse_unimod_comp(mod))

    # chemical formula
    if mod_lower.startswith('formula:'):
        return parse_chem_formula(mod.split(':')[1])

    return None


def _parse_mod_delta_mass(mod: str) -> Union[float, None]:
    """
    Parse the mod delta mass. If the mod cannot be parsed into a delta mass, it will return None.

    :param mod: The mod to parse.
    :type mod: str

    :return: The parsed delta mass, or None if the delta mass cannot be parsed.
    :rtype: float | None

    .. code-block:: python

        >>> _parse_mod_delta_mass('42.0')
        42.0

        >>> _parse_mod_delta_mass('Acetyl')

        >>> _parse_mod_delta_mass('UniMod:1')

        >>> _parse_mod_delta_mass('UniMod:+1')
        1.0

        >>> _parse_mod_delta_mass('Obs:42.0')
        42.0
    """

    if isinstance(mod, (int, float)):
        return mod

    mass = None
    for mod_str in mod.split('|'):

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


def _sequence_comp(annotation: Union[str, ProFormaAnnotation],
                   ion_type: str,
                   isotope: int = 0,
                   use_isotope_on_mods: bool = False) -> ChemComposition:
    """
    Calculate the composition of a sequence.

    :param annotation: The sequence or ProForma annotation.
    :type annotation: str | ProFormaAnnotation
    :param ion_type: The ion type.
    :type ion_type: str
    :param isotope: The number of Neutrons to add/subtract from the final mass. Default is 0.
    :type isotope: int
    :param use_isotope_on_mods: If True, the isotope modifications will be applied to the final composition.
    Default is False.
    :type use_isotope_on_mods: bool

    :raises UnknownModificationError: If the modification is unknown.
    :raises AmbiguousAminoAcidError: If the sequence contains an ambiguous amino acid.

    :return: The composition of the sequence.
    :rtype: Dict[str, int | float]

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> _sequence_comp('PEPTIDE/1', 'y')
        {'C': 34, 'H': 54, 'N': 7, 'O': 15, 'e': -1}

        >>> _sequence_comp('PEPTIDE/1', 'y', isotope=1)
        {'C': 34, 'H': 54, 'N': 7, 'O': 15, 'e': -1, 'n': 1}

        >>> _sequence_comp('G/1', 'i')
        {'C': 1, 'H': 4, 'N': 1, 'e': -1}

        >>> _sequence_comp('PEPTIDE/1', 'b')
        {'C': 34, 'H': 52, 'N': 7, 'O': 14, 'e': -1}

        >>> _sequence_comp('<H>PEPTIDE/1', 'b')
        {'C': 34, 'H': 52, 'N': 7, 'O': 14, 'e': -1}

        >>> _sequence_comp('{Unimod:2}PEPTIDE', 'p')
        {'C': 34, 'H': 54, 'N': 8, 'O': 14}

        >>> _sequence_comp('<13C>PEPTIDE/1', 'b')
        {'H': 52, 'N': 7, 'O': 14, 'e': -1, '13C': 34}

        >>> _sequence_comp('PEPTIDE[Unimod:2]/1', 'y')
        {'C': 34, 'H': 55, 'N': 8, 'O': 14, 'e': -1}

        >>> _sequence_comp('<[Unimod:2]@T>PEPTIDE/1', 'y')
        {'C': 34, 'H': 55, 'N': 8, 'O': 14, 'e': -1}

        >>> _sequence_comp('<[Unimod:2]@N-Term>PEPTIDE/1', 'y')
        {'C': 34, 'H': 55, 'N': 8, 'O': 14, 'e': -1}

        >>> _sequence_comp('PEPTIDE/2', 'p')
        {'C': 34, 'H': 55, 'N': 7, 'O': 15, 'e': -2}

        >>> _sequence_comp('PEPTIDE/2[+2Na+]', 'p')
        {'C': 34, 'H': 53, 'N': 7, 'O': 15, 'Na': 2, 'e': -2}

        >>> _sequence_comp('<13C>PEPTIDE[Formula:C10]', 'p')
        {'H': 53, 'N': 7, 'O': 15, '13C': 34, 'C': 10}

        >>> _sequence_comp('<13C>PEPTIDE[Formula:C10]', 'p', use_isotope_on_mods=True)
        {'H': 53, 'N': 7, 'O': 15, '13C': 44}

        >>> _sequence_comp('<13C>PEPTIDE[Unimod:213413]', 'b')
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: Unimod:213413

        >>> _sequence_comp('I', 'by')
        {'C': 6, 'H': 11, 'N': 1, 'O': 1}

        # Ambiguous amino acid
        >>> _sequence_comp('B', 'by') # doctest: +ELLIPSIS
        Traceback (most recent call last):
            ...
        peptacular.errors.AmbiguousAminoAcidError: Ambiguous amino acid: B! Cannot determine the composition ...

    """

    if isinstance(annotation, str):
        annotation = sequence_to_annotation(annotation)

    # If charge is not provided, set it to 0
    charge = annotation.charge
    if charge is None:
        charge = 0

    charge_adducts = annotation.charge_adducts
    if charge_adducts is not None:
        charge_adducts = annotation.charge_adducts[0]

    if charge_adducts is None:
        if ion_type in ('p', 'n'):
            charge_adducts = f'{charge}H+'
        else:
            charge_adducts = f'{charge-1}H+,{FRAGMENT_ION_BASE_CHARGE_ADDUCTS[ion_type]}'

    if ion_type not in ('p', 'n'):
        if charge == 0:
            warnings.warn('Calculating the comp of a fragment ion with charge state 0. Fragment ions should have a '
                          'charge state greater than 0 since the neutral form doesnt exist.')

    if 'B' in annotation.sequence:
        raise AmbiguousAminoAcidError('B',
                                      'Cannot determine the composition of a sequence with an ambiguous amino acid.')

    if 'Z' in annotation.sequence:
        raise AmbiguousAminoAcidError('Z',
                                      'Cannot determine the composition of a sequence with an ambiguous amino acid.')

    # Get the composition of the base sequence
    sequence_composition = {}
    for aa in annotation.sequence:
        try:
            aa_comp = AA_COMPOSITIONS[aa]
        except KeyError as err:
            raise UnknownAminoAcidError(aa) from err
        for k, v in aa_comp.items():
            sequence_composition[k] = sequence_composition.get(k, 0) + v

    # Apply the adjustments for the neutral fragment composition based on strictly the ion dissociation points.
    for k, v in NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS[ion_type].items():
        sequence_composition[k] = sequence_composition.get(k, 0) + v

    charge_adduct_comp = _parse_charge_adducts_comp(charge_adducts)

    for k, v in charge_adduct_comp.items():
        sequence_composition[k] = sequence_composition.get(k, 0) + v

    mod_composition = {}
    if annotation.has_unknown_mods():
        for unknown_mod in annotation.unknown_mods:
            for k, v in mod_comp(unknown_mod).items():
                mod_composition[k] = mod_composition.get(k, 0) + v

    if annotation.has_intervals():
        for interval in annotation.intervals:
            if interval.has_mods():
                for interval_mod in interval.mods:
                    for k, v in mod_comp(interval_mod).items():
                        mod_composition[k] = mod_composition.get(k, 0) + v

    if annotation.has_labile_mods() and ion_type == 'p':
        for labile_mod in annotation.labile_mods:
            for k, v in mod_comp(labile_mod).items():
                mod_composition[k] = mod_composition.get(k, 0) + v

    if annotation.has_nterm_mods():
        for nterm_mod in annotation.nterm_mods:
            for k, v in mod_comp(nterm_mod).items():
                mod_composition[k] = mod_composition.get(k, 0) + v

    if annotation.has_cterm_mods():
        for cterm_mod in annotation.cterm_mods:
            for k, v in mod_comp(cterm_mod).items():
                mod_composition[k] = mod_composition.get(k, 0) + v

    if annotation.has_internal_mods():
        for _, internal_mods in annotation.internal_mods.items():
            for internal_mod in internal_mods:
                for k, v in mod_comp(internal_mod).items():
                    mod_composition[k] = mod_composition.get(k, 0) + v

    if annotation.has_static_mods():
        static_map = parse_static_mods(annotation.static_mods)

        n_term_mod = static_map.get('N-Term')
        if n_term_mod is not None:
            for m in n_term_mod:
                for k, v in mod_comp(m.val).items():
                    mod_composition[k] = mod_composition.get(k, 0) + v

        c_term_mod = static_map.get('C-Term')
        if c_term_mod is not None:
            for m in c_term_mod:
                for k, v in mod_comp(m.val).items():
                    mod_composition[k] = mod_composition.get(k, 0) + v

        for aa, mod in static_map.items():
            if aa in ['N-Term', 'C-Term']:
                continue

            aa_count = annotation.sequence.count(aa)
            for m in mod:
                for k, v in mod_comp(m.val).items():
                    mod_composition[k] = mod_composition.get(k, 0) + v * aa_count

    mod_composition['n'] = mod_composition.get('n', 0) + isotope

    # Apply isotopic mods
    if annotation.has_isotope_mods():
        if use_isotope_on_mods:
            sequence_composition = apply_isotope_mods_to_composition(sequence_composition, annotation.isotope_mods)
            mod_composition = apply_isotope_mods_to_composition(mod_composition, annotation.isotope_mods)
        else:
            sequence_composition = apply_isotope_mods_to_composition(sequence_composition, annotation.isotope_mods)

    composition = {}
    for k, v in sequence_composition.items():
        composition[k] = composition.get(k, 0) + v

    for k, v in mod_composition.items():
        composition[k] = composition.get(k, 0) + v

    composition = {k: v for k, v in composition.items() if v != 0}

    return composition


def _parse_mod_delta_mass_only(mod: Union[str, Mod]) -> Union[float, None]:
    """
    Parse a modification string.

    :param mod: The modification string.
    :type mod: str | Mod

    :raises ValueError: If no delta mass or composition can be retrieved from the modification string.

    :return: The delta mass of the modification or None if the modification is not a delta mass.
    :rtype: float | None

    .. code-block:: python

        >>> _parse_mod_delta_mass_only('Acetyl|INFO:newly discovered')

        >>> _parse_mod_delta_mass_only('Acetyl|Obs:+42.010565')

        >>> _parse_mod_delta_mass_only('Obs:+42.010565')
        42.010565

        >>> _parse_mod_delta_mass_only('R:+1|INFO:Fantastic')
        1.0

    """

    if isinstance(mod, (int, float)):
        return mod

    if isinstance(mod, Mod):
        return _parse_mod_delta_mass_only(mod.val)

    mods = mod.split('|')
    for m in mods:
        try:
            m = _parse_mod_comp(m)
        except DeltaMassCompositionError:
            continue
        if m is not None:
            return None

    for m in mods:
        m = _parse_mod_delta_mass(m)
        if m is not None:
            return m

    raise ValueError(f'Invalid modification: {mod}')


def apply_isotope_mods_to_composition(composition: Union[ChemComposition, str],
                                      isotopic_mods: Optional[List[ModValue]]) -> ChemComposition:
    """
    Apply isotopic modifications to a composition.

    :param composition: The composition to apply the isotopic modifications to.
    :type composition: ChemComposition
    :param isotopic_mods: The isotopic modifications.
    :type isotopic_mods: List[str | Mod] | None

    :return: The modified composition.
    :rtype: Dict[str, int | float]

    .. code-block:: python

        # Apply isotopic modifications to a composition.
        >>> apply_isotope_mods_to_composition({'C': 6, 'H': 12, 'O': 6}, ['13C'])
        {'H': 12, 'O': 6, '13C': 6}

        # Apply multiple isotopic modifications to a composition.
        >>> apply_isotope_mods_to_composition({'C': 6, 'H': 12, 'O': 6}, ['13C', '15N'])
        {'H': 12, 'O': 6, '13C': 6}

        # Works with floats
        >>> apply_isotope_mods_to_composition({'C': 6.6, 'H': 12, 'O': 6}, ['13C', '15N'])
        {'H': 12, 'O': 6, '13C': 6.6}

    """
    if isinstance(composition, str):
        composition = parse_chem_formula(composition)
    else:
        composition = composition.copy()

    if isotopic_mods is None:
        return composition

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


def _parse_adduct_comp(adduct: str) -> ChemComposition:
    """
    Parse an adduct string and return its mass.

    :param adduct: The adduct string to parse.
    :type adduct: str

    :raises InvalidDeltaMassError: If the adduct contains an invalid delta mass.

    :return: The composition of the adduct.
    :rtype: Dict[str, int | float]

    .. code-block:: python

        # Parse an adduct string.
        >>> _parse_adduct_comp('+Na+')
        {'Na': 1, 'e': -1}

        >>> _parse_adduct_comp('+2Na+')
        {'Na': 2, 'e': -2}

        >>> _parse_adduct_comp('2H+')
        {'H': 2, 'e': -2}

    """

    comp = {}
    element_count, element_symbol, element_charge = parse_ion_elements(adduct)
    comp[element_symbol] = element_count
    comp['e'] = -1 * element_charge * element_count

    return comp


def _parse_charge_adducts_comp(adducts: ModValue) -> ChemComposition:
    """
    Parse the charge adducts and return their mass.

    :param adducts: The charge adducts to parse.
    :type adducts: ModValue

    :return: The composition of the charge adducts.
    :rtype: Dict[str, int | float]

    .. code-block:: python

        # Parse the charge adducts and return their mass.
        >>> _parse_charge_adducts_comp('+Na+,+H+')
        {'Na': 1, 'e': -2, 'H': 1}

        >>> _parse_charge_adducts_comp('2H+')
        {'H': 2, 'e': -2}

    """

    if isinstance(adducts, Mod):
        return _parse_charge_adducts_comp(adducts.val)

    if not isinstance(adducts, str):
        raise TypeError(f'Invalid type for adducts: {type(adducts)}! Must be a string.')

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
