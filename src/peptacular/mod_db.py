"""
mod_db.py - Functions for getting modification masses and compositions from unimod and psi-mod databases from
peptacular/data/psi and peptacular/data/unimod.
"""

from peptacular import constants
from peptacular.errors import UnknownModificationError, InvalidDeltaMassError, InvalidCompositionError


def is_unimod_str(unimod_str) -> bool:
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

    return unimod_str.lower().startswith('unimod:') or unimod_str.lower().startswith(
        'u:') or unimod_str in constants.UNIMOD_NAME_TO_ID or unimod_str in constants.UNIMOD_ID_TO_MONO_MASSES


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

    if unimod_str.lower().startswith('unimod:') or unimod_str.lower().startswith('u:'):
        return unimod_str.split(':')[1]
    else:
        return unimod_str


def parse_unimod_mass(unimod_str: str, monoisotopic: bool, precision: int = None) -> float:
    """
    Get the mass of an unimod modification.

    :param unimod_str: The unimod id or name.
    :type unimod_str: str
    :param monoisotopic: Whether to get the monoisotopic mass or the average mass.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to.
    :type precision: int

    :raises UnknownModificationError: If the unimod id/name is not found.
    :raises InvalidDeltaMassError: If the unimod id is a delta mass and is not a valid number.

    :return: The mass of the unimod modification.
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

    round_func = lambda x: round(x, precision) if precision is not None else x

    orig_str = unimod_str
    unimod_str = _strip_unimod_str(unimod_str)

    if unimod_str.startswith('+') or unimod_str.startswith('-'):
        try:
            return round_func(float(unimod_str))
        except ValueError:
            raise InvalidDeltaMassError(orig_str)

    elif unimod_str in constants.UNIMOD_NAME_TO_ID:  # is a name
        unimod_id = constants.UNIMOD_NAME_TO_ID[unimod_str]

        if monoisotopic is True:
            return round_func(constants.UNIMOD_ID_TO_MONO_MASSES[unimod_id])
        else:
            return round_func(constants.UNIMOD_ID_TO_AVERAGE_MASSES[unimod_id])

    elif unimod_str in constants.UNIMOD_ID_TO_MONO_MASSES:  # is an id

        if monoisotopic is True:
            return round_func(constants.UNIMOD_ID_TO_MONO_MASSES[unimod_str])
        else:
            return round_func(constants.UNIMOD_ID_TO_AVERAGE_MASSES[unimod_str])

    else:
        raise UnknownModificationError(orig_str)


def parse_unimod_comp(unimod_str: str) -> str:
    """
    Get the composition of an unimod modification.

    :param unimod_str: The unimod id or name.
    :type unimod_str: str

    :raises UnknownModificationError: If the unimod id/name is not found.
    :raises InvalidCompositionError: If the unimod id is a delta mass and has no composition.

    :return: The composition of the unimod modification.
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

        >>> parse_unimod_comp('+1')
        Traceback (most recent call last):
        peptacular.errors.InvalidCompositionError: Cannot retrieve composition for: +1

        >>> parse_unimod_comp('13252454')
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: 13252454

    """

    orig_str = unimod_str
    unimod_str = _strip_unimod_str(unimod_str)

    if unimod_str.startswith('+') or unimod_str.startswith('-'):
        raise InvalidCompositionError(orig_str)

    elif unimod_str in constants.UNIMOD_NAME_TO_ID:  # is a name
        unimod_id = constants.UNIMOD_NAME_TO_ID[unimod_str]
        comp = constants.UNIMOD_ID_TO_ISOTOPIC_COMPOSITIONS[unimod_id]
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    elif unimod_str in constants.UNIMOD_ID_TO_MONO_MASSES:  # is an id
        comp = constants.UNIMOD_ID_TO_ISOTOPIC_COMPOSITIONS[unimod_str]
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    else:
        raise UnknownModificationError(orig_str)


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

    return psi_str.lower().startswith('mod:') or psi_str.lower().startswith('m:') or psi_str.lower().startswith(
        'psi-mod:') or psi_str in constants.PSI_MOD_NAME_TO_ID or psi_str in constants.PSI_MOD_ID_TO_ISOTOPIC_MASSES


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

    if psi_str.lower().startswith('mod:') or psi_str.lower().startswith('m:') or psi_str.lower().startswith('psi-mod:'):
        return psi_str.split(':')[1]
    else:
        return psi_str


def parse_psi_mass(psi_str: str, monoisotopic: bool, precision: int = None) -> float:
    """
    Get the mass of a PSI modification.

    :param psi_str: The PSI id or name.
    :type psi_str: str
    :param monoisotopic: Whether to get the monoisotopic mass or the average mass.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to.
    :type precision: int

    :raises UnknownModificationError: If the PSI id/name is not found.
    :raises InvalidDeltaMassError: If the PSI id is a delta mass and is not a valid number.

    :return: The mass of the PSI modification.
    :rtype: float

    .. code-block:: python


        >>> parse_psi_mass('MOD:00046', monoisotopic=True, precision=3)
        79.966331

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

    round_func = lambda x: round(x, precision) if precision is not None else x

    orig_str = psi_str
    psi_str = _strip_psi_str(psi_str)

    if psi_str.startswith('+') or psi_str.startswith('-'):
        try:
            return round_func(float(psi_str))
        except ValueError:
            raise InvalidDeltaMassError(orig_str)

    elif psi_str in constants.PSI_MOD_NAME_TO_ID:  # is a name
        psi_id = constants.PSI_MOD_NAME_TO_ID[psi_str]
        if monoisotopic is True:
            return constants.PSI_MOD_ID_TO_ISOTOPIC_MASSES[psi_id]
        else:
            return constants.PSI_MOD_ID_TO_AVERAGE_MASSES[psi_id]

    elif psi_str in constants.PSI_MOD_ID_TO_ISOTOPIC_MASSES:  # is an id
        if monoisotopic is True:
            return constants.PSI_MOD_ID_TO_ISOTOPIC_MASSES[psi_str]
        else:
            return constants.PSI_MOD_ID_TO_AVERAGE_MASSES[psi_str]

    else:
        raise UnknownModificationError(orig_str)


def parse_psi_comp(psi_str: str) -> str:
    """
    Get the composition of a PSI modification.

    :param psi_str: The PSI id or name.
    :type psi_str: str

    :raises UnknownModificationError: If the PSI id/name is not found.
    :raises InvalidCompositionError: If the PSI id is a delta mass and has no composition.

    :return: The composition of the PSI modification.
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
        peptacular.errors.InvalidCompositionError: Cannot retrieve composition for: M:+1

        # InvalidCompositionError when the composition is None
        >>> parse_psi_comp('M:modified D-asparagine residue')
        Traceback (most recent call last):
        peptacular.errors.InvalidCompositionError: Cannot retrieve composition for: M:modified D-asparagine residue

        # UnknownModificationError when the modification is not found
        >>> parse_psi_comp('13252454')
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: 13252454

    """

    orig_str = psi_str

    if psi_str.lower().startswith('mod:') or psi_str.lower().startswith('m:') or psi_str.lower().startswith('psi-mod:'):
        psi_str = psi_str.split(':')[1]

    if psi_str.startswith('+') or psi_str.startswith('-'):
        raise InvalidCompositionError(orig_str)

    elif psi_str in constants.PSI_MOD_NAME_TO_ID:  # is a name
        psi_id = constants.PSI_MOD_NAME_TO_ID[psi_str]
        comp = constants.PSI_MOD_ID_TO_COMPOSITIONS[psi_id]
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    elif psi_str in constants.PSI_MOD_ID_TO_ISOTOPIC_MASSES:  # is an id
        comp = constants.PSI_MOD_ID_TO_COMPOSITIONS[psi_str]
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    else:
        raise UnknownModificationError(orig_str)
