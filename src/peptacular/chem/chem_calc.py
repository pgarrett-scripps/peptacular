"""
chem_calc.py contains functions for parsing and writing chemical formulas, and for calculating the composition of a
sequence with/and modifications.
"""

import copy
from typing import Counter, Iterable

from ..mod import Mod, MOD_VALUE_TYPES
from ..utils2 import parse_ion_elements
from ..util import parse_isotope_mods

from .chem_constants import ISOTOPIC_AVERAGINE_MASS
from .chem_util import parse_chem_formula, write_chem_formula
from ..mods.mod_db_setup import MONOSACCHARIDES_DB
from ..constants import (
    AVERAGINE_RATIOS,
)
from ..errors import (
    InvalidCompositionError,
    DeltaMassCompositionError,
)
from ..glycan import glycan_comp
from ..mods.mod_db import (
    parse_unimod_comp,
    parse_psi_comp,
    is_psi_mod_str,
    is_unimod_str,
    parse_xlmod_comp,
    parse_resid_comp,
    is_resid_str,
    is_xlmod_str,
    is_gno_str,
    parse_gno_comp,
)

from ..utils2 import convert_type


def glycan_to_chem(glycan: dict[str, int | float] | str) -> str:
    """
    Converts a glycan composition/formula to a chemical formula.

    :param glycan: A dictionary containing the glycan components and their counts, or a glycan formula string.
    :type glycan: dict[str, int | float] | str

    :return: A chemical formula string.
    :rtype: str

    .. code-block:: python

            >>> glycan_to_chem({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1})
            'C45H73N3O34'

            >>> glycan_to_chem({'HexNAc': 2.3, 'Hex': 3.3, 'Neu5Gc': 1.1})
            'C50.3H81.6N3.4O37.9'

    """

    return write_chem_formula(glycan_comp(glycan))


def mod_comp(mod: MOD_VALUE_TYPES | Mod) -> dict[str, int | float]:
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
        raise InvalidCompositionError(str(mod))

    mods = mod.split("|")
    for m in mods:
        m = _parse_mod_comp(m)
        if m is not None:
            return m

    raise InvalidCompositionError(mod)


def estimate_comp(
    neutral_mass: float, isotopic_mods: Iterable[MOD_VALUE_TYPES] | None = None
) -> dict[str, int | float]:
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

    composition = {
        atom: ratio * neutral_mass / ISOTOPIC_AVERAGINE_MASS
        for atom, ratio in AVERAGINE_RATIOS.items()
    }

    if isotopic_mods is not None:
        composition = apply_isotope_mods_to_composition(composition, isotopic_mods)

    return composition

def _parse_glycan_comp(glycan_str: str) -> dict[str, int | float]:
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

    if glycan_str.lower().startswith("glycan:"):
        glycan_str = "".join(glycan_str.split(":")[1:])

    if MONOSACCHARIDES_DB.contains_id(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_id(glycan_str)
    elif MONOSACCHARIDES_DB.contains_name(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_name(glycan_str)
    elif MONOSACCHARIDES_DB.contains_synonym(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_synonym(glycan_str)
    else:
        return parse_chem_formula(glycan_to_chem(glycan_str))

    if entry.composition is None:
        raise InvalidCompositionError(
            f"Error parsing glycan formula: '{glycan_str}'. "
            "Glycan does not have a valid composition."
        )

    return parse_chem_formula(entry.composition)


def _parse_mod_comp(mod: str) -> dict[str, int | float] | None:
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

    if isinstance(
        converted_mod, (int, float)
    ):  # cannot get composition from a delta mass
        return None

    # localization fix
    if isinstance(mod, str) and "#" in mod:  # type: ignore
        if mod.startswith("#"):  # for localized positions, return empty comp
            return {}
        mod = mod.split("#")[0]

    mod_lower = mod.lower()
    if mod_lower.startswith("glycan:"):
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

    if mod_lower.startswith("info:"):  # cannot get comp for info
        return None

    if mod_lower.startswith("obs:"):  # cannot get comp for observed mass
        return None

    if is_psi_mod_str(mod):  # psi-mod
        return parse_chem_formula(parse_psi_comp(mod))

    if is_unimod_str(mod):  # unimod
        return parse_chem_formula(parse_unimod_comp(mod))

    # chemical formula
    if mod_lower.startswith("formula:"):
        return parse_chem_formula(mod.split(":")[1])

    return None


def _parse_mod_delta_mass(mod: str) -> float | None:
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
    for mod_str in mod.split("|"):

        # check if the mod contains a localization score
        if "#" in mod_str:
            mod_str = mod_str.split("#")[0]

        try:
            mass = float(mod_str)
            break
        except ValueError:
            pass

        mod_lower = mod_str.lower()
        if (
            mod_lower.startswith("unimod:")
            or mod_lower.startswith("u:")
            or mod_lower.startswith("psi-mod:")
            or mod_lower.startswith("mod:")
            or mod_lower.startswith("m:")
            or mod_lower.startswith("resid:")
            or mod_lower.startswith("r:")
            or mod_lower.startswith("xlmod:")
            or mod_lower.startswith("x:")
            or mod_lower.startswith("gno:")
            or mod_lower.startswith("g:")
        ):

            mod_str = mod_str.split(":")[1]

            if mod_str.startswith("+") or mod_str.startswith("-"):
                mass = float(mod_str)
                break

        if mod_str.lower().startswith("obs:"):
            mod_str = mod_str.split(":")[1]
            try:
                mass = float(mod_str)
                break
            except ValueError:
                pass

    return mass


def parse_mod_delta_mass_only(mod: str | int | float | Mod) -> float | None:
    """
    Parse a modification string.

    :param mod: The modification string.
    :type mod: str | Mod

    :raises ValueError: If no delta mass or composition can be retrieved from the modification string.

    :return: The delta mass of the modification or None if the modification is not a delta mass.
    :rtype: float | None

    .. code-block:: python

        >>> parse_mod_delta_mass_only('Acetyl|INFO:newly discovered')

        >>> parse_mod_delta_mass_only('Acetyl|Obs:+42.010565')

        >>> parse_mod_delta_mass_only('Obs:+42.010565')
        42.010565

        >>> parse_mod_delta_mass_only('R:+1|INFO:Fantastic')
        1.0

    """

    if isinstance(mod, (int, float)):
        return mod

    if isinstance(mod, Mod):
        return parse_mod_delta_mass_only(mod.val)

    mods = mod.split("|")
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

    raise ValueError(f"Invalid modification: {mod}")


def apply_isotope_mods_to_composition(
    composition: dict[str, int | float] | str,
    isotopic_mods: Iterable[MOD_VALUE_TYPES] | None,
    inplace: bool = True,
) -> dict[str, int | float]:
    """
    Apply isotopic modifications to a composition.
    
    :param composition: The composition to apply the isotopic modifications to.
    :type composition: MutableMapping[str, int | float] | str
    :param isotopic_mods: The isotopic modifications.
    :type isotopic_mods: Iterable[MOD_VALUE_TYPES] | None
    :param inplace: Whether to modify the composition in place. Default is True.
    :type inplace: bool
    :return: The modified composition (same type as input, or dict if input was str).
    :rtype: Same as input type | dict[str, int | float]
    
    .. code-block:: python
    
        # Apply isotopic modifications to a composition.
        >>> apply_isotope_mods_to_composition({'C': 6, 'H': 12, 'O': 6}, ['13C'])
        {'H': 12, 'O': 6, '13C': 6}
        
        # Apply multiple isotopic modifications to a composition.
        >>> apply_isotope_mods_to_composition({'C': 6, 'H': 12, 'O': 6}, ['13C', '15N'])
        {'H': 12, 'O': 6, '13C': 6}
        
        # Works with Counter
        >>> counter = Counter({'C': 6, 'H': 12, 'O': 6})
        >>> result = apply_isotope_mods_to_composition(counter, ['13C'])
        >>> isinstance(result, Counter)  # True
    """
    if isinstance(composition, str):
        return apply_isotope_mods_to_composition(parse_chem_formula(composition), isotopic_mods, inplace=False)
        
    if not inplace:
        composition = copy.deepcopy(composition)

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



def get_isotope_composition_change(
    composition: dict[str, int | float] | str,
    isotopic_mods: Iterable[MOD_VALUE_TYPES] | None,
) -> dict[str, int | float]:
    if isinstance(composition, str):
        return get_isotope_composition_change(parse_chem_formula(composition), isotopic_mods)
    
    if isotopic_mods is None:
        return {}
    
    composition_change: dict[str, int | float] = {}
    isotope_map = parse_isotope_mods(isotopic_mods)
    
    for element, isotope_label in isotope_map.items():
        if element in composition:
            if element == isotope_label:
                continue
            
            # Record the removal of the original element
            composition_change[element] = composition_change.get(element, 0) - composition[element]
            
            # Record the addition of the isotope label
            composition_change[isotope_label] = composition_change.get(isotope_label, 0) + composition[element]
    
    return composition_change



def _parse_adduct_comp(adduct: str) -> dict[str, int | float]:
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

    comp: dict[str, int | float] = {}
    element_count, element_symbol, element_charge = parse_ion_elements(adduct)
    comp[element_symbol] = element_count
    comp["e"] = -1 * element_charge * element_count

    return comp


def parse_charge_adducts_comp(adducts: int | float | str | Mod) -> Counter[str]:
    """
    Parse the charge adducts and return their mass.

    :param adducts: The charge adducts to parse.
    :type adducts: MOD_VALUE_TYPES | Mod

    :return: The composition of the charge adducts.
    :rtype: Dict[str, int | float]

    .. code-block:: python

        # Parse the charge adducts and return their mass.
        >>> parse_charge_adducts_comp('+Na+,+H+')
        {'Na': 1, 'e': -2, 'H': 1}

        >>> parse_charge_adducts_comp('2H+')
        {'H': 2, 'e': -2}

    """

    if isinstance(adducts, Mod):
        return Counter(parse_charge_adducts_comp(adducts.val))

    if not isinstance(adducts, str):
        raise TypeError(f"Invalid type for adducts: {type(adducts)}! Must be a string.")

    adduct_comps: list[str] = adducts.split(",")

    comps: list[dict[str, int | float]] = []
    for adduct_comp in adduct_comps:
        comps.append(_parse_adduct_comp(adduct_comp))

    # Merge the compositions
    composition: Counter[str] = Counter()
    for comp in comps:
        composition.update(comp)

    return composition
