"""
mod_db.py - Functions for getting modification masses and compositions from unimod and psi-mod databases from
peptacular/data/psi and peptacular/data/unimod.
"""
from typing import Optional

from peptacular.errors import UnknownModificationError, InvalidDeltaMassError, InvalidCompositionError, \
    DeltaMassCompositionError, UnknownModificationMassError
from peptacular.mods.mod_db_setup import UNIMOD_DB, PSI_MOD_DB, XLMOD_DB, EntryDb, RESID_DB, GNO_DB


def _get_mass(db: EntryDb, mod_str: str, orig_str: str, monoisotopic: bool, precision: Optional[int]) -> float:
    """
    Helper function for getting the mass of a modification from a db.
    """
    round_func = lambda x: round(x, precision) if precision is not None else x

    if mod_str.startswith('+') or mod_str.startswith('-'):
        try:
            return round_func(float(mod_str))
        except ValueError as err:
            raise InvalidDeltaMassError(orig_str) from err

    if db.contains_id(mod_str):
        entry = db.get_entry_by_id(mod_str)
    elif db.contains_name(mod_str):
        entry = db.get_entry_by_name(mod_str)
    else:
        entry = None

    if not entry:
        raise UnknownModificationError(orig_str)

    if monoisotopic is True:
        m = entry.mono_mass
        if m is None:
            m = entry.calc_mono_mass
    else:
        m = entry.avg_mass
        if m is None:
            m = entry.calc_avg_mass

    if m is None:
        raise UnknownModificationMassError(orig_str)

    return round_func(m)



def _get_comp(db: EntryDb, mod_str: str, orig_str: str) -> str:
    """
    Helper function for getting the composition of a modification from a db.
    """
    if mod_str.startswith('+') or mod_str.startswith('-'):
        try:
            _ = float(mod_str)
        except ValueError as err:
            raise InvalidDeltaMassError(orig_str) from err

        raise DeltaMassCompositionError(orig_str)

    if db.contains_id(mod_str):
        entry = db.get_entry_by_id(mod_str)
    elif db.contains_name(mod_str):
        entry = db.get_entry_by_name(mod_str)
    else:
        entry = None

    if entry:
        comp = entry.composition
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    raise UnknownModificationError(orig_str)


def is_unimod_str(unimod_str: str) -> bool:
    """
    Check if a string is a valid unimod id or name.

    :param unimod_str: The unimod id or name.
    :type unimod_str: str

    :return: Whether the string is a valid unimod id or name.
    :rtype: bool

    .. code-block:: python

        >>> is_unimod_str('UniMod:1')
        True

        >>> is_unimod_str('U:1')
        True

        >>> is_unimod_str('1')
        True

        >>> is_unimod_str('Acetyl')
        True

        >>> is_unimod_str('13252454')
        False

    """
    unimod_str_lower = unimod_str.lower()
    return unimod_str_lower.startswith('unimod:') or unimod_str_lower.startswith(
        'u:') or UNIMOD_DB.contains_id(unimod_str) or UNIMOD_DB.contains_name(unimod_str)


def _strip_unimod_str(unimod_str: str) -> str:
    """
    Strip a unimod id or name to just the id.

    :param unimod_str: The unimod id or name.
    :type unimod_str: str

    :return: The unimod id.
    :rtype: str

    .. code-block:: python

        >>> _strip_unimod_str('UniMod:1')
        '1'

        >>> _strip_unimod_str('U:1')
        '1'

        >>> _strip_unimod_str('1')
        '1'

        >>> _strip_unimod_str('Acetyl')
        'Acetyl'

    """
    unimod_str_lower = unimod_str.lower()
    if unimod_str_lower.startswith('unimod:') or unimod_str_lower.startswith('u:'):
        return unimod_str.split(':')[1]
    return unimod_str


def parse_unimod_mass(mod_str: str, monoisotopic: bool, precision: Optional[int] = None) -> float:
    """
    Get the mass of an unimod modification.

    :param mod_str: The unimod id or name.
    :type mod_str: str
    :param monoisotopic: Whether to get the monoisotopic mass or the average mass.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to. Default is None.
    :type precision: Optional[int]

    :raises UnknownModificationError: If the mod entry cannot be found.
    :raises UnknownModificationMassError: If the entry has no mass.
    :raises InvalidDeltaMassError: If the mod is a delta mass and is not a valid number.

    :return: The mass of the modification.
    :rtype: float

    .. code-block:: python


        >>> parse_unimod_mass('UniMod:1', monoisotopic=True, precision=3)
        42.011

        >>> parse_unimod_mass('UniMod:1', monoisotopic=False, precision=3)
        42.037

        >>> parse_unimod_mass('U:1', monoisotopic=True, precision=3)
        42.011

        >>> parse_unimod_mass('1', monoisotopic=True, precision=3)
        42.011

        >>> parse_unimod_mass('Acetyl', monoisotopic=True, precision=3)
        42.011

        >>> parse_unimod_mass('+1', monoisotopic=True, precision=3)
        1.0

        >>> parse_unimod_mass('+1a', monoisotopic=True, precision=3)
        Traceback (most recent call last):
        peptacular.errors.InvalidDeltaMassError: Invalid delta mass: +1a

        >>> parse_unimod_mass('13252454', monoisotopic=False)
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: 13252454

    """

    return _get_mass(UNIMOD_DB, _strip_unimod_str(mod_str), mod_str, monoisotopic, precision)


def parse_unimod_comp(mod_str: str) -> str:
    """
    Get the composition of an UNIMOD modification.

    :param mod_str: The id or name of the mod.
    :type mod_str: str

    :raises UnknownModificationError: If the mod entry cannot be found.
    :raises InvalidCompositionError: If the composition is None.
    :raises DeltaMassCompositionError: If the mod is a delta mass and has no composition.
    :raises InvalidDeltaMassError: If the mod is a delta mass and is not a valid number.

    :return: The mods chemical formula (In Proforma2.0 Noation)
    :rtype: str

    .. code-block:: python

        >>> parse_unimod_comp('UniMod:1')
        'H2C2O1'

        >>> parse_unimod_comp('U:1')
        'H2C2O1'

        >>> parse_unimod_comp('1')
        'H2C2O1'

        >>> parse_unimod_comp('Acetyl')
        'H2C2O1'

        >>> parse_unimod_comp('U:+1')
        Traceback (most recent call last):
        peptacular.errors.DeltaMassCompositionError: Cannot retrieve composition for: U:+1

        >>> parse_unimod_comp('13252454')
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: 13252454

    """

    return _get_comp(UNIMOD_DB, _strip_unimod_str(mod_str), mod_str)


def is_psi_mod_str(psi_str: str) -> bool:
    """
    Check if a string is a PSI id or name.

    :param psi_str: The PSI id or name.
    :type psi_str: str

    :return: Whether the string is a PSI id or name.
    :rtype: bool

    .. code-block:: python

        >>> is_psi_mod_str('MOD:00046')
        True

        >>> is_psi_mod_str('M:00046')
        True

        >>> is_psi_mod_str('PSI-MOD:00046')
        True

        >>> is_psi_mod_str('00046')
        True

        >>> is_psi_mod_str('O-phospho-L-serine')
        True

        >>> is_psi_mod_str('13252454')
        False

    """
    psi_str_lower = psi_str.lower()
    return psi_str_lower.startswith('mod:') or psi_str_lower.startswith('m:') or psi_str_lower.startswith(
        'psi-mod:') or PSI_MOD_DB.contains_id(psi_str) or PSI_MOD_DB.contains_name(psi_str)


def _strip_psi_str(psi_str: str) -> str:
    """
    Strip a PSI id or name to just the id.

    :param psi_str: The PSI id or name.
    :type psi_str: str

    :return: The PSI id.
    :rtype: str

    .. code-block:: python

        >>> _strip_psi_str('MOD:00046')
        '00046'

        >>> _strip_psi_str('PSI-MOD:00046')
        '00046'

        >>> _strip_psi_str('M:00046')
        '00046'

        >>> _strip_psi_str('00046')
        '00046'

        >>> _strip_psi_str('O-phospho-L-serine')
        'O-phospho-L-serine'

    """
    psi_str_lower = psi_str.lower()
    if psi_str_lower.startswith('mod:') or psi_str_lower.startswith('m:') or psi_str_lower.startswith('psi-mod:'):
        return psi_str.split(':')[1]
    return psi_str


def parse_psi_mass(mod_str: str, monoisotopic: bool, precision: Optional[int] = None) -> float:
    """
    Get the mass of a PSI modification.

    :param mod_str: The unimod id or name.
    :type mod_str: str
    :param monoisotopic: Whether to get the monoisotopic mass or the average mass.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to. Default is None.
    :type precision: Optional[int]

    :raises UnknownModificationError: If the mod entry cannot be found.
    :raises UnknownModificationMassError: If the entry has no mass.
    :raises InvalidDeltaMassError: If the mod is a delta mass and is not a valid number.

    :return: The mass of the modification.
    :rtype: float

    .. code-block:: python


        >>> parse_psi_mass('MOD:00046', monoisotopic=True, precision=3)
        79.966

        >>> parse_psi_mass('MOD:00046', monoisotopic=False, precision=3)
        79.98

        >>> parse_psi_mass('PSI-MOD:00046', monoisotopic=False, precision=3)
        79.98

        >>> parse_psi_mass('M:00046', monoisotopic=False, precision=3)
        79.98

        >>> parse_psi_mass('00046', monoisotopic=False, precision=3)
        79.98

        >>> parse_psi_mass('O-phospho-L-serine', monoisotopic=False, precision=3)
        79.98

        >>> parse_psi_mass('M:+1', monoisotopic=False, precision=3)
        1.0

        >>> parse_psi_mass('M:+1a', monoisotopic=True, precision=3)
        Traceback (most recent call last):
        peptacular.errors.InvalidDeltaMassError: Invalid delta mass: M:+1a

        >>> parse_psi_mass('13252454', monoisotopic=False, precision=3)
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: 13252454

    """

    return _get_mass(PSI_MOD_DB, _strip_psi_str(mod_str), mod_str, monoisotopic, precision)


def parse_psi_comp(mod_str: str) -> str:
    """
    Get the composition of a PSI modification.

    :param mod_str: The id or name of the mod.
    :type mod_str: str

    :raises UnknownModificationError: If the mod entry cannot be found.
    :raises InvalidCompositionError: If the composition is None.
    :raises DeltaMassCompositionError: If the mod is a delta mass and has no composition.
    :raises InvalidDeltaMassError: If the mod is a delta mass and is not a valid number.

    :return: The mods chemical formula (In Proforma2.0 Noation)
    :rtype: str

    .. code-block:: python

        >>> parse_psi_comp('MOD:00046')
        'H1O3P1'

        >>> parse_psi_comp('M:00046')
        'H1O3P1'

        >>> parse_psi_comp('PSI-MOD:00046')
        'H1O3P1'

        >>> parse_psi_comp('00046')
        'H1O3P1'

        >>> parse_psi_comp('O-phospho-L-serine')
        'H1O3P1'

        # InvalidCompositionError when the modification is a delta mass offset
        >>> parse_psi_comp('M:+1')
        Traceback (most recent call last):
        peptacular.errors.DeltaMassCompositionError: Cannot retrieve composition for: M:+1

        # InvalidCompositionError when the composition is None
        >>> parse_psi_comp('M:modified D-asparagine residue')
        Traceback (most recent call last):
        peptacular.errors.InvalidCompositionError: Cannot retrieve composition for: M:modified D-asparagine residue

        # UnknownModificationError when the modification is not found
        >>> parse_psi_comp('13252454')
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: 13252454

    """

    return _get_comp(PSI_MOD_DB, _strip_psi_str(mod_str), mod_str)


def is_xlmod_str(xlmod_str: str) -> bool:
    """
    Check if a string is a xlmod id or name.

    :param xlmod_str:
    :return:
    """
    xlmod_str_lower = xlmod_str.lower()
    return xlmod_str_lower.startswith('xlmod:') or xlmod_str_lower.startswith('x:')


def _strip_xlmod_str(xlmod_str: str) -> str:
    """
    Strip a xlmod id or name to just the id.

    :param xlmod_str:
    :return:
    """
    xlmod_str_lower = xlmod_str.lower()
    if xlmod_str_lower.startswith('xlmod:') or xlmod_str_lower.startswith('x:'):
        return xlmod_str.split(':')[1]
    return xlmod_str


def parse_xlmod_mass(mod_str: str, monoisotopic: bool, precision: Optional[int] = None) -> float:
    """
    Get the mass of a xlmod modification.

    :param mod_str: The unimod id or name.
    :type mod_str: str
    :param monoisotopic: Whether to get the monoisotopic mass or the average mass.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to. Default is None.
    :type precision: Optional[int]

    :raises UnknownModificationError: If the mod entry cannot be found.
    :raises UnknownModificationMassError: If the entry has no mass.
    :raises InvalidDeltaMassError: If the mod is a delta mass and is not a valid number.

    :return: The mass of the modification.
    :rtype: float

    .. code-block:: python

        >>> parse_xlmod_mass('XLMOD:01000', monoisotopic=True, precision=3)
        156.079

        >>> parse_xlmod_mass('XLMOD:01000', monoisotopic=False, precision=3)
        156.179

        >>> parse_xlmod_mass('X:01000', monoisotopic=True, precision=3)
        156.079

        >>> parse_xlmod_mass('X:01000', monoisotopic=False, precision=3)
        156.179

        >>> parse_xlmod_mass('X:+1', monoisotopic=True, precision=3)
        1.0

        >>> parse_xlmod_mass('X:-3.1415', monoisotopic=False, precision=3)
        -3.142

    """

    return _get_mass(XLMOD_DB, _strip_xlmod_str(mod_str), mod_str, monoisotopic, precision)


def parse_xlmod_comp(mod_str: str) -> str:
    """
    Get the composition of a XLMOD modification.

    :param mod_str: The id or name of the mod.
    :type mod_str: str

    :raises UnknownModificationError: If the mod entry cannot be found.
    :raises InvalidCompositionError: If the composition is None.
    :raises DeltaMassCompositionError: If the mod is a delta mass and has no composition.
    :raises InvalidDeltaMassError: If the mod is a delta mass and is not a valid number.

    :return: The mods chemical formula (In Proforma2.0 Noation)
    :rtype: str

    .. code-block:: python

        >>> parse_xlmod_comp('XLMOD:01000')
        'C8H12O3'

        >>> parse_xlmod_comp('X:01000')
        'C8H12O3'

        >>> parse_xlmod_comp('X:+1')
        Traceback (most recent call last):
        peptacular.errors.DeltaMassCompositionError: Cannot retrieve composition for: X:+1

        >>> parse_xlmod_comp('X:-3.1415')
        Traceback (most recent call last):
        peptacular.errors.DeltaMassCompositionError: Cannot retrieve composition for: X:-3.1415

        >>> parse_xlmod_comp('13252454')
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: 13252454

    """

    return _get_comp(XLMOD_DB, _strip_xlmod_str(mod_str), mod_str)


def is_resid_str(resid_str: str) -> bool:
    """
    Check if a string is a RESID id or name.
    """
    resid_str_lower = resid_str.lower()
    return resid_str_lower.startswith('resid:') or resid_str_lower.startswith('r:')


def _strip_resid_str(resid_str: str) -> str:
    resid_str_lower = resid_str.lower()
    if resid_str_lower.startswith('resid:') or resid_str_lower.startswith('r:'):
        return resid_str.split(':')[1]
    return resid_str


def parse_resid_mass(mod_str: str, monoisotopic: bool, precision: Optional[int] = None) -> float:
    """
    Get the mass of a residue modification.

    :param mod_str: The unimod id or name.
    :type mod_str: str
    :param monoisotopic: Whether to get the monoisotopic mass or the average mass.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to. Default is None.
    :type precision: Optional[int]

    :raises UnknownModificationError: If the mod entry cannot be found.
    :raises UnknownModificationMassError: If the entry has no mass.
    :raises InvalidDeltaMassError: If the mod is a delta mass and is not a valid number.

    :return: The mass of the modification.
    :rtype: float

    .. code-block:: python

        >>> parse_resid_mass('RESID:AA0317', monoisotopic=True, precision=3)
        14.016

        >>> parse_resid_mass('RESID:AA0317', monoisotopic=False, precision=3)
        14.03

        >>> parse_resid_mass('R:AA0317', monoisotopic=True, precision=3)
        14.016

        >>> parse_resid_mass('R:+1', monoisotopic=True, precision=3)
        1.0

        >>> parse_resid_mass('R:-1', monoisotopic=True, precision=3)
        -1.0

    """
    return _get_mass(RESID_DB, _strip_resid_str(mod_str), mod_str, monoisotopic, precision)


def parse_resid_comp(mod_str: str) -> str:
    """
    Get the composition of a RESID modification.

    :param mod_str: The id or name of the mod.
    :type mod_str: str

    :raises UnknownModificationError: If the mod entry cannot be found.
    :raises InvalidCompositionError: If the composition is None.
    :raises DeltaMassCompositionError: If the mod is a delta mass and has no composition.
    :raises InvalidDeltaMassError: If the mod is a delta mass and is not a valid number.

    :return: The mods chemical formula (In Proforma2.0 Noation)
    :rtype: str

    .. code-block:: python

        >>> parse_resid_comp('RESID:AA0317')
        'C1H2'

        >>> parse_resid_comp('R:AA0317')
        'C1H2'

        >>> parse_resid_comp('R:+1')
        Traceback (most recent call last):
        peptacular.errors.DeltaMassCompositionError: Cannot retrieve composition for: R:+1

        >>> parse_resid_comp('R:-1')
        Traceback (most recent call last):
        peptacular.errors.DeltaMassCompositionError: Cannot retrieve composition for: R:-1

    """

    return _get_comp(RESID_DB, _strip_resid_str(mod_str), mod_str)


def is_gno_str(gno_str: str) -> bool:
    """
    Check if a string is a GNO id or name.

    :param gno_str: The GNO id or name.
    :type gno_str: str

    :return: Whether the string is a GNO id or name.
    :rtype: bool

    .. code-block:: python

        >>> is_gno_str('GNO:00000202')
        True

        >>> is_gno_str('G:00000202')
        True

        >>> is_gno_str('G:N-acetylneuraminic acid')
        True

        >>> is_gno_str('00000202')
        False

        >>> is_gno_str('N-acetylneuraminic acid')
        False

        >>> is_gno_str('13252454')
        False

    """

    gno_str_lower = gno_str.lower()
    return gno_str_lower.startswith('gno:') or gno_str_lower.startswith('g:')


def _strip_gno_str(gno_str: str) -> str:
    """
    Strip a GNO id or name to just the id.
    :param gno_str:
    :return:
    """

    gno_str_lower = gno_str.lower()
    if gno_str_lower.startswith('gno:') or gno_str_lower.startswith('g:'):
        return gno_str.split(':')[1]
    return gno_str


def parse_gno_mass(mod_str: str, monoisotopic: bool, precision: Optional[int] = None) -> float:
    """
    Get the mass of a GNO modification.

    :param mod_str: The unimod id or name.
    :type mod_str: str
    :param monoisotopic: Whether to get the monoisotopic mass or the average mass.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to. Default is None.
    :type precision: Optional[int]

    :raises UnknownModificationError: If the mod entry cannot be found.
    :raises UnknownModificationMassError: If the entry has no mass.
    :raises InvalidDeltaMassError: If the mod is a delta mass and is not a valid number.

    :return: The mass of the modification.
    :rtype: float

    .. code-block:: python

        >>> parse_gno_mass('GNO:G35503UV', monoisotopic=True, precision=3)
        1631.618

        >>> parse_gno_mass('GNO:G35503UV', monoisotopic=False, precision=3)
        1632.529

        >>> parse_gno_mass('G:G35503UV', monoisotopic=True, precision=3)
        1631.618

        >>> parse_gno_mass('G:G35503UV', monoisotopic=False, precision=3)
        1632.529

        >>> parse_gno_mass('G:+1', monoisotopic=True, precision=3)
        1.0

        >>> parse_gno_mass('G:-3.14', monoisotopic=True, precision=3)
        -3.14


    """
    return _get_mass(GNO_DB, _strip_gno_str(mod_str), mod_str, monoisotopic, precision)


def parse_gno_comp(mod_str: str) -> str:
    """
    Get the composition of a GNO modification.

    :param mod_str: The id or name of the mod.
    :type mod_str: str

    :raises UnknownModificationError: If the mod entry cannot be found.
    :raises InvalidCompositionError: If the composition is None.
    :raises DeltaMassCompositionError: If the mod is a delta mass and has no composition.
    :raises InvalidDeltaMassError: If the mod is a delta mass and is not a valid number.

    :return: The mods chemical formula (In Proforma2.0 Noation)
    :rtype: str

    .. code-block:: python

        >>> parse_gno_comp('GNO:G35503UV')
        'C64H105N5O43'

        >>> parse_gno_comp('G:G35503UV')
        'C64H105N5O43'

        >>> parse_gno_comp('G:+1')
        Traceback (most recent call last):
        peptacular.errors.DeltaMassCompositionError: Cannot retrieve composition for: G:+1

        >>> parse_gno_comp('G:-3.14')
        Traceback (most recent call last):
        peptacular.errors.DeltaMassCompositionError: Cannot retrieve composition for: G:-3.14

        >>> parse_gno_comp('G:13252454')
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: G:13252454

    """

    return _get_comp(GNO_DB, _strip_gno_str(mod_str), mod_str)
