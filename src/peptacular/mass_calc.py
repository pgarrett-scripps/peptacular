"""
mass_calc.py is a simple module for computing the m/z and mass of an amino acid sequence.
"""

from .chem.chem_calc import (
    parse_chem_formula,
)
from .chem.chem_constants import (
    AVERAGE_FRAGMENT_ADJUSTMENTS,
    AVERAGE_FRAGMENT_ION_ADJUSTMENTS,
    MONOISOTOPIC_FRAGMENT_ADJUSTMENTS,
    MONOISOTOPIC_FRAGMENT_ION_ADJUSTMENTS,
)
from .chem.chem_util import chem_mass
from .constants import (
    AVERAGE_ATOMIC_MASSES,
    ELECTRON_MASS,
    ISOTOPIC_ATOMIC_MASSES,
    NEUTRON_MASS,
    PROTON_MASS,
    IonType,
)
from .errors import (
    InvalidChemFormulaError,
    InvalidDeltaMassError,
    InvalidGlycanFormulaError,
    InvalidModificationMassError,
)
from .funcs import convert_type, parse_ion_elements, round_to_precision
from .glycan import parse_glycan_formula
from .mod import MOD_VALUE_TYPES, Mod
from .mods.mod_db import (
    is_gno_str,
    is_psi_mod_str,
    is_resid_str,
    is_unimod_str,
    is_xlmod_str,
    parse_gno_mass,
    parse_psi_mass,
    parse_resid_mass,
    parse_unimod_mass,
    parse_xlmod_mass,
)
from .mods.mod_db_setup import MONOSACCHARIDES_DB


def _convert_adducts_to_str(
    charge_adducts: str | list[str] | list[Mod] | None = None,
) -> str | None:
    """
    Convert charge adducts to a string representation.

    :param charge_adducts: The charge adducts to convert.
    :type charge_adducts: str | List[str] | List[Mod] | None
    :return: A string representation of the charge adducts.
    :rtype: str

    :raises TypeError: If the charge adducts are not a string, list of strings, or Mod.

    .. code-block:: python
        # Convert charge adducts to a string representation.
        >>> _convert_adducts_to_str('+H+,+Na+')
        '+H+,+Na+'

        >>> _convert_adducts_to_str(['+H+', '+Na+'])
        '+H+,+Na+'

        >>> _convert_adducts_to_str(Mod('+H+', 2))
        '+H+,+H+'

        >>> _convert_adducts_to_str(None)
    """
    if isinstance(charge_adducts, Mod):
        return ",".join([str(charge_adducts.val)] * charge_adducts.mult)

    elif isinstance(charge_adducts, list):
        adduct_str: list[str] = []
        for adduct in charge_adducts:
            if isinstance(adduct, Mod):
                sub_list = [str(adduct.val)] * adduct.mult
                adduct_str.extend(sub_list)
            elif isinstance(adduct, str):  # type: ignore
                adduct_str.append(adduct)
            else:
                raise TypeError(
                    f"Invalid value for charge adducts: {adduct}! Must be a string or Mod, got: {type(adduct)}"
                )
        return ",".join(adduct_str)
    elif isinstance(charge_adducts, str):
        return charge_adducts
    elif charge_adducts is None:
        return None

    raise TypeError(
        f"Invalid value for charge adducts: {charge_adducts}! Must be a string, list of strings, or Mod, got: {type(charge_adducts)}"
    )


def adjust_mass(
    base_mass: float,
    charge: int | None,
    ion_type: str = "p",
    monoisotopic: bool = True,
    isotope: int = 0,
    loss: float = 0.0,
    charge_adducts: str | list[str] | list[Mod] | None = None,
    precision: int | None = None,
) -> float:
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

    adduct_str = _convert_adducts_to_str(charge_adducts)

    if adduct_str is None or adduct_str == "":
        if ion_type in (IonType.PRECURSOR, IonType.NEUTRAL):
            charge_adduct_mass = PROTON_MASS * charge
        else:
            frag_ion_offset = (
                MONOISOTOPIC_FRAGMENT_ION_ADJUSTMENTS[ion_type]
                if monoisotopic is True
                else AVERAGE_FRAGMENT_ION_ADJUSTMENTS[ion_type]
            )
            charge_adduct_mass = PROTON_MASS * (charge - 1) + frag_ion_offset
    else:
        charge_adduct_mass = _parse_charge_adducts_mass(
            adduct_str, monoisotopic=monoisotopic
        )

    m += charge_adduct_mass

    # Base mass
    m += (
        MONOISOTOPIC_FRAGMENT_ADJUSTMENTS[ion_type]
        if monoisotopic
        else AVERAGE_FRAGMENT_ADJUSTMENTS[ion_type]
    )

    # m += (charge * charge_adduct_mass)
    m += isotope * NEUTRON_MASS + loss  # Add isotope and loss

    if precision is not None:
        m = round(m, precision)

    return m


def adjust_mz(
    base_mass: float, charge: int | None, precision: int | None = None
) -> float:
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


def chem_mz(
    formula: dict[str, int | float] | str,
    charge: int = 1,
    monoisotopic: bool = True,
    precision: int | None = None,
    sep: str = "",
) -> float:
    """
    Calculate the m/z of a chemical formula.
    """
    # TODO: Add charge adducts?
    m = chem_mass(formula, monoisotopic, precision, sep)
    return adjust_mz(m, charge, precision)


def glycan_mass(
    formula: str | dict[str, int | float],
    monoisotopic: bool = True,
    precision: int | None = None,
) -> float:
    """
    Calculate the mass of a glycan formula.

    :param formula: A glycan formula or dict[str, int | float].
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
            raise InvalidGlycanFormulaError(
                str(original_formula), f'Unknown monosaccharide: "{monosaccharide}"!'
            )

        if monoisotopic:
            if entry.mono_mass is None:
                raise ValueError(
                    f'Cannot determine monoisotopic mass for: "{monosaccharide}"'
                )

            m += entry.mono_mass * count
        else:
            if entry.avg_mass is None:
                raise ValueError(
                    f'Cannot determine average mass for: "{monosaccharide}"'
                )

            m += entry.avg_mass * count

    return round_to_precision(m, precision)


def glycan_mz(
    formula: str | dict[str, int | float],
    charge: int = 1,
    monoisotopic: bool = True,
    precision: int | None = None,
) -> float:
    """
    Calculate the m/z of a glycan formula.
    """
    # TODO: Add charge adducts?
    m = glycan_mass(formula, monoisotopic, precision)
    return adjust_mz(m, charge, precision)

from functools import lru_cache

@lru_cache(maxsize=1024)
def _cached_mod_mass(mod_val: str | int | float, mod_mult: int, monoisotopic: bool, precision: int | None) -> float:
    """
    Cache mass calculations for Mod named tuples.
    
    :param mod_val: The modification value (string, int, or float).
    :type mod_val: str | int | float
    :param mod_mult: The modification multiplier.
    :type mod_mult: int
    :param monoisotopic: If True, use monoisotopic mass else use average mass.
    :type monoisotopic: bool
    :param precision: The precision of the mass.
    :type precision: int | None
    :return: The cached modification mass.
    :rtype: float
    """
    base_mass = _calculate_single_mod_mass(mod_val, monoisotopic, precision)
    return base_mass * mod_mult


def _calculate_single_mod_mass(mod_val: str | int | float, monoisotopic: bool, precision: int | None) -> float:
    """
    Calculate mass for a single modification value.
    
    :param mod_val: The modification value (string, int, or float).
    :type mod_val: str | int | float
    :param monoisotopic: If True, use monoisotopic mass else use average mass.
    :type monoisotopic: bool
    :param precision: The precision of the mass.
    :type precision: int | None
    :return: The modification mass.
    :rtype: float
    """
    if isinstance(mod_val, int):
        return mod_val
    
    if isinstance(mod_val, float):
        return round_to_precision(mod_val, precision)
    
    # String case
    mods = mod_val.split("|")
    for m in mods:
        mass = _parse_mod_mass(m, monoisotopic, precision)
        if mass is not None:
            return mass
    
    raise InvalidModificationMassError(mod_val)


def mod_mass(
    mod: str | int | float | Mod | list[Mod],
    monoisotopic: bool = True,
    precision: int | None = None,
) -> float:
    """
    Parse a modification string.
    
    [Original docstring remains the same...]
    """
    if isinstance(mod, list):
        return sum(mod_mass(m, monoisotopic, precision) for m in mod)
    
    if isinstance(mod, Mod):
        # Use cached version for Mod named tuples
        return _cached_mod_mass(mod.val, mod.mult, monoisotopic, precision)
    
    if isinstance(mod, int):
        return mod
    
    if isinstance(mod, float):
        return round_to_precision(mod, precision)
    
    # For string mods, calculate directly (not cached)
    mods = mod.split("|")
    for m in mods:
        mass = _parse_mod_mass(m, monoisotopic, precision)
        if mass is not None:
            return mass
    
    raise InvalidModificationMassError(mod)


def _parse_obs_mass_from_proforma_str(
    obs_str: str, precision: int | None = None
) -> float:
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

    if obs_str.lower().startswith("obs:"):
        obs_str = "".join(obs_str.split(":")[1:])

    obs_float = convert_type(obs_str)

    if not isinstance(obs_float, (int, float)):
        raise InvalidDeltaMassError(obs_str)

    return round_to_precision(float(obs_float), precision)


def _parse_glycan_mass_from_proforma_str(
    glycan_str: str, monoisotopic: bool, precision: int | None = None
) -> float:
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

    if glycan_str.lower().startswith("glycan:"):
        glycan_str = "".join(glycan_str.split(":")[1:])

    if MONOSACCHARIDES_DB.contains_id(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_id(glycan_str)
    elif MONOSACCHARIDES_DB.contains_name(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_name(glycan_str)
    elif MONOSACCHARIDES_DB.contains_synonym(glycan_str):
        entry = MONOSACCHARIDES_DB.get_entry_by_synonym(glycan_str)
    else:
        entry = None

    if entry:
        if monoisotopic is True:
            if entry.mono_mass is None:
                raise InvalidGlycanFormulaError(
                    glycan_str, f'Glycan "{glycan_str}" has no monoisotopic mass!'
                )
            return round_to_precision(entry.mono_mass, precision)
        if entry.avg_mass is None:
            raise InvalidGlycanFormulaError(
                glycan_str, f'Glycan "{glycan_str}" has no average mass!'
            )
        return round_to_precision(entry.avg_mass, precision)

    try:
        return round_to_precision(
            glycan_mass(glycan_str, monoisotopic, precision), precision
        )
    except InvalidGlycanFormulaError as err:
        raise InvalidGlycanFormulaError(glycan_str, err.msg) from err


def _parse_chem_mass_from_proforma_str(
    chem_str: str, monoisotopic: bool, precision: int | None = None
) -> float:
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

    if chem_str.lower().startswith("formula:"):
        chem_str = "".join(chem_str.split(":")[1:])

    comps = parse_chem_formula(chem_str)
    try:
        return chem_mass(comps, monoisotopic, precision)
    except InvalidChemFormulaError as err:
        raise InvalidChemFormulaError(chem_str, err.msg) from err


from functools import lru_cache


@lru_cache(maxsize=512)
def _parse_mod_mass_cached(
    mod: str | int | float,
    mod_type: str,  # Add type to differentiate 1.0 from 1
    monoisotopic: bool = True,
    precision: int | None = None,
) -> float | None:
    """
    Cached version of mod mass parsing.

    Args:
        mod: The modification value
        mod_type: Type name ("str", "int", "float") to avoid hash collisions
        monoisotopic: Whether to use monoisotopic masses
        precision: The precision of the mass

    Note: mod_type is only used for cache key differentiation since hash(1.0) == hash(1)
    """
    # Get type name for cache key to avoid hash(1.0) == hash(1) collision
    mod_type = type(mod).__name__

    # Only cache if input is hashable (str, int, float)
    if isinstance(mod, (str, int, float)):
        return _parse_mod_mass_cached(mod, mod_type, monoisotopic, precision)
    else:
        # Fallback for non-hashable types
        raise TypeError(f"Cannot parse mod mass for type {type(mod)}")


def _parse_mod_mass(
    mod: str, monoisotopic: bool = True, precision: int | None = None
) -> float | None:
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
    if isinstance(mod, str) and "#" in mod:  # type: ignore
        if mod.startswith("#"):  # for localized positions, return 0
            return 0.0

        mod = mod.split("#")[0]

    # Try to parse as a number first (Might cause issues if the mod is also unimod/psi ID)
    # Proforma2.0 standard requires that delta mass instances are always prefixed with a '+' or '-' but this would
    # require a lot of changes....and make user input more difficult

    mod = convert_type(mod)  # type: ignore

    if isinstance(mod, int):
        return mod

    if isinstance(mod, float):
        return round(mod, precision) if precision is not None else mod

    mod_lower = mod.lower()
    if mod_lower.startswith("glycan:"):
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

    if mod_lower.startswith("info:"):  # Skip info modifications
        return None

    if is_psi_mod_str(mod):  # is a psi-mod modification
        # raises UnknownModificationError and InvalidDeltaMassError
        return parse_psi_mass(mod, monoisotopic, precision)

    if is_unimod_str(mod):  # is a unimod modification
        # raises UnknownModificationError and InvalidDeltaMassError
        return parse_unimod_mass(mod, monoisotopic, precision)

    # chemical formula
    if mod_lower.startswith("formula:"):
        # raises UnknownElementError
        return _parse_chem_mass_from_proforma_str(mod, monoisotopic, precision)

    # observed mass
    if mod_lower.startswith("obs:"):
        # raises InvalidDeltaMassError
        return _parse_obs_mass_from_proforma_str(mod, precision)

    return None


def _parse_adduct_mass(
    adduct: str, precision: int | None = None, monoisotopic: bool = True
) -> float:
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

    if element_symbol == "e":
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


def _parse_charge_adducts_mass(
    adducts: MOD_VALUE_TYPES,
    precision: int | None = None,
    monoisotopic: bool = True,
) -> float:
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
        raise TypeError(
            f"Invalid value for charge adducts: {adducts}! Must be a string, got: {type(adducts)}"
        )

    if adducts == "+H+":
        return PROTON_MASS

    adducts_comps: list[str] = adducts.split(",")

    m = 0.0

    for adduct in adducts_comps:
        m += _parse_adduct_mass(adduct, None, monoisotopic)

    if precision is not None:
        m = round(m, precision)

    return m


def ppm_error(theo: float, expt: float, precision: int | None = None) -> float:
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


def dalton_error(theo: float, expt: float, precision: int | None = None) -> float:
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
