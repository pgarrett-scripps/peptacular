"""
mod_db.py - Functions for getting modification masses and compositions from unimod and psi-mod databases from
peptacular/data/psi and peptacular/data/unimod.
"""

from peptacular import constants
from peptacular.errors import UnknownModificationError, InvalidDeltaMassError, InvalidCompositionError, \
    DeltaMassCompositionError


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
    unimod_str_lower = unimod_str.lower()
    return unimod_str_lower.startswith('unimod:') or unimod_str_lower.startswith(
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
    unimod_str_lower = unimod_str.lower()
    if unimod_str_lower.startswith('unimod:') or unimod_str_lower.startswith('u:'):
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
        raise DeltaMassCompositionError(orig_str)

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
    psi_str_lower = psi_str.lower()
    return psi_str_lower.startswith('mod:') or psi_str_lower.startswith('m:') or psi_str_lower.startswith(
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
    psi_str_lower = psi_str.lower()
    if psi_str_lower.startswith('mod:') or psi_str_lower.startswith('m:') or psi_str_lower.startswith('psi-mod:'):
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
        raise DeltaMassCompositionError(orig_str)

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
    else:
        return xlmod_str


def parse_xlmod_mass(xlmod_str: str, monoisotopic: bool, precision: int = None) -> float:
    """
    Get the mass of a xlmod modification.

    :param xlmod_str: The xlmod id or name.
    :type xlmod_str: str
    :param monoisotopic: Whether to get the monoisotopic mass or the average mass.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to.
    :type precision: int

    :raises UnknownModificationError: If the xlmod id/name is not found.
    :raises InvalidDeltaMassError: If the xlmod id is a delta mass and is not a valid number.

    :return: The mass of the xlmod modification.
    :rtype: float

    .. code-block:: python

        >>> parse_xlmod_mass('XLMOD:01000', monoisotopic=True, precision=3)
        156.07864431

        >>> parse_xlmod_mass('XLMOD:01000', monoisotopic=False, precision=3)
        156.17939099550614

        >>> parse_xlmod_mass('X:01000', monoisotopic=True, precision=3)
        156.07864431

        >>> parse_xlmod_mass('X:01000', monoisotopic=False, precision=3)
        156.17939099550614

        >>> parse_xlmod_mass('X:+1', monoisotopic=True, precision=3)
        1.0

        >>> parse_xlmod_mass('X:-3.1415', monoisotopic=False, precision=4)
        -3.1415

    """

    round_func = lambda x: round(x, precision) if precision is not None else x

    orig_str = xlmod_str
    xlmod_str = _strip_xlmod_str(xlmod_str)

    if xlmod_str.startswith('+') or xlmod_str.startswith('-'):
        try:
            return round_func(float(xlmod_str))
        except ValueError:
            raise InvalidDeltaMassError(orig_str)

    elif xlmod_str in constants.XLMOD_NAME_TO_ID:  # is a name
        xlmod_id = constants.XLMOD_NAME_TO_ID[xlmod_str]
        if monoisotopic is True:
            return constants.XLMOD_ID_TO_ISOTOPIC_MASSES[xlmod_id]
        else:
            return constants.XLMOD_ID_TO_AVERAGE_MASSES[xlmod_id]

    elif xlmod_str in constants.XLMOD_ID_TO_ISOTOPIC_MASSES:  # is an id
        if monoisotopic is True:
            return constants.XLMOD_ID_TO_ISOTOPIC_MASSES[xlmod_str]
        else:
            return constants.XLMOD_ID_TO_AVERAGE_MASSES[xlmod_str]

    else:
        raise UnknownModificationError(orig_str)


def parse_xlmod_comp(xlmod_str: str) -> str:
    """
    Get the composition of a xlmod modification.

    :param xlmod_str: The xlmod id or name.
    :type xlmod_str: str

    :raises UnknownModificationError: If the xlmod id/name is not found.
    :raises InvalidCompositionError: If the xlmod id is a delta mass and has no composition.

    :return: The composition of the xlmod modification.
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

    orig_str = xlmod_str

    if xlmod_str.lower().startswith('xlmod:') or xlmod_str.lower().startswith('x:'):
        xlmod_str = xlmod_str.split(':')[1]

    if xlmod_str.startswith('+') or xlmod_str.startswith('-'):
        raise DeltaMassCompositionError(orig_str)

    elif xlmod_str in constants.XLMOD_NAME_TO_ID:  # is a name
        xlmod_id = constants.XLMOD_NAME_TO_ID[xlmod_str]
        comp = constants.XLMOD_ID_TO_COMPOSITIONS[xlmod_id]
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    elif xlmod_str in constants.XLMOD_ID_TO_ISOTOPIC_MASSES:  # is an id
        comp = constants.XLMOD_ID_TO_COMPOSITIONS[xlmod_str]
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    else:
        raise UnknownModificationError(orig_str)


def is_resid_str(resid_str: str) -> bool:
    resid_str_lower = resid_str.lower()
    return resid_str_lower.startswith('resid:') or resid_str_lower.startswith('r:')


def _strip_resid_str(resid_str: str) -> str:
    resid_str_lower = resid_str.lower()
    if resid_str_lower.startswith('resid:') or resid_str_lower.startswith('r:'):
        return resid_str.split(':')[1]
    else:
        return resid_str


def parse_resid_mass(resid_str: str, monoisotopic: bool, precision: int = None) -> float:
    """
    Get the mass of a residue modification.
    :param resid_str:
    :param monoisotopic:
    :param precision:
    :return:

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
    round_func = lambda x: round(x, precision) if precision is not None else x

    orig_str = resid_str
    resid_str = _strip_resid_str(resid_str)

    if resid_str.startswith('+') or resid_str.startswith('-'):
        try:
            return round_func(float(resid_str))
        except ValueError:
            raise InvalidDeltaMassError(orig_str)

    elif resid_str in constants.RESID_NAME_TO_ID:  # is a name
        resid_id = constants.RESID_NAME_TO_ID[resid_str]
        if monoisotopic is True:
            return round_func(constants.RESID_ID_TO_ISOTOPIC_MASSES[resid_id])
        else:
            return round_func(constants.RESID_ID_TO_AVERAGE_MASSES[resid_id])

    elif resid_str in constants.RESID_ID_TO_ISOTOPIC_MASSES:  # is an id
        if monoisotopic is True:
            return round_func(constants.RESID_ID_TO_ISOTOPIC_MASSES[resid_str])
        else:
            return round_func(constants.RESID_ID_TO_AVERAGE_MASSES[resid_str])

    else:
        raise UnknownModificationError(orig_str)


def parse_resid_comp(resid_str: str) -> str:
    """
    Get the composition of a residue modification.
    :param resid_str:
    :return:

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

    orig_str = resid_str

    if resid_str.lower().startswith('resid:') or resid_str.lower().startswith('r:'):
        resid_str = resid_str.split(':')[1]

    if resid_str.startswith('+') or resid_str.startswith('-'):
        raise DeltaMassCompositionError(orig_str)

    elif resid_str in constants.RESID_NAME_TO_ID:  # is a name
        resid_id = constants.RESID_NAME_TO_ID[resid_str]
        comp = constants.RESID_ID_TO_COMPOSITIONS[resid_id]
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    elif resid_str in constants.RESID_ID_TO_ISOTOPIC_MASSES:  # is an id
        comp = constants.RESID_ID_TO_COMPOSITIONS[resid_str]
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    else:
        raise UnknownModificationError(orig_str)


# GNO

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

        >>> is_gno_str('00000202')
        True

        >>> is_gno_str('N-acetylneuraminic acid')
        True

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
    else:
        return gno_str


def parse_gno_mass(gno_str: str, monoisotopic: bool, precision: int = None) -> float:
    """
    Get the mass of a GNO modification.

    :param gno_str: The GNO id or name.
    :type gno_str: str
    :param monoisotopic: Whether to get the monoisotopic mass or the average mass.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to.
    :type precision: int

    :raises UnknownModificationError: If the GNO id/name is not found.
    :raises InvalidDeltaMassError: If the GNO id is a delta mass and is not a valid number.

    :return: The mass of the GNO modification.
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

    round_func = lambda x: round(x, precision) if precision is not None else x

    orig_str = gno_str
    gno_str = _strip_gno_str(gno_str)

    if gno_str.startswith('+') or gno_str.startswith('-'):
        try:
            return round_func(float(gno_str))
        except ValueError:
            raise InvalidDeltaMassError(orig_str)

    elif gno_str in constants.GNO_NAME_TO_ID:  # is a name
        gno_id = constants.GNO_NAME_TO_ID[gno_str]

        if monoisotopic is True:
            return round_func(constants.GNO_ID_TO_ISOTOPIC_MASSES[gno_id])
        else:
            return round_func(constants.GNO_ID_TO_AVERAGE_MASSES[gno_id])

    elif gno_str in constants.GNO_ID_TO_ISOTOPIC_MASSES:  # is an id

        if monoisotopic is True:
            return round_func(constants.GNO_ID_TO_ISOTOPIC_MASSES[gno_str])
        else:
            return round_func(constants.GNO_ID_TO_AVERAGE_MASSES[gno_str])

    else:
        raise UnknownModificationError(orig_str)


def parse_gno_comp(gno_str: str) -> str:
    """
    Get the composition of a GNO modification.

    :param gno_str: The GNO id or name.
    :type gno_str: str

    :raises UnknownModificationError: If the GNO id/name is not found.
    :raises InvalidCompositionError: If the GNO id is a delta mass and has no composition.

    :return: The composition of the GNO modification.
    :rtype: str

    .. code-block:: python

        >>> parse_gno_comp('GNO:G35503UV')
        'C64H105N5O43'

        >>> parse_gno_comp('G:G35503UV')
        'C64H105N5O43'

    """

    orig_str = gno_str

    if gno_str.lower().startswith('gno:') or gno_str.lower().startswith('g:'):
        gno_str = gno_str.split(':')[1]

    if gno_str.startswith('+') or gno_str.startswith('-'):
        raise DeltaMassCompositionError(orig_str)

    elif gno_str in constants.GNO_NAME_TO_ID:  # is a name
        gno_id = constants.GNO_NAME_TO_ID[gno_str]
        comp = constants.GNO_ID_TO_COMPOSITIONS[gno_id]
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    elif gno_str in constants.GNO_ID_TO_COMPOSITIONS:  # is an id
        comp = constants.GNO_ID_TO_COMPOSITIONS[gno_str]
        if comp is None:
            raise InvalidCompositionError(orig_str)
        return comp

    else:
        raise UnknownModificationError(orig_str)