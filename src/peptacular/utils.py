import sys
from typing import Iterable, Mapping, Protocol, runtime_checkable

from .constants import ModType, ModTypeLiteral
from .elements import ELEMENT_LOOKUP


@runtime_checkable
class SupportsStr(Protocol):
    """Protocol for any object that can be converted to string"""

    def __str__(self) -> str: ...


def handle_number_and_intern_mod(mod: SupportsStr | float | int) -> str:
    """Validate and intern a modification string"""
    if isinstance(mod, (float, int)):
        # ensure has +/- sign for positive/negative masses
        mod_str = f"{mod:+}"
    else:
        mod_str = str(mod).strip()
        if not mod_str:
            raise ValueError("Empty modification string is not allowed")
    return sys.intern(mod_str)


def get_mod_type(mod: ModTypeLiteral | ModType) -> ModType:
    # return ModType Enum for the given mod string
    if isinstance(mod, ModType):
        return mod

    if isinstance(mod, str):  # type: ignore
        for mod_type in ModType:
            if mod_type.value == mod:
                return mod_type
    else:
        raise ValueError(f"mod must be a string or ModType, got {type(mod)}")

    raise ValueError(f"Unknown mod type: {mod}")


def get_mods(
    mods: Iterable[ModTypeLiteral | ModType] | ModType | ModTypeLiteral | None,
) -> list[ModType]:
    """
    Get the list of modification types from the input.

    :param mods: Modification types as a ModType, iterable of ModTypes, or None.
    :type mods: None | ModType | Iterable[ModType]
    :return: List of ModType Enum values.
    :rtype: list[ModType]
    :raises ValueError: If mods is not None, ModType, or iterable of ModTypes.
    """

    if mods is None:
        return [mod_type for mod_type in ModType]
    elif isinstance(mods, (str, ModType)):
        # Single modification type
        return [get_mod_type(mods)]
    elif isinstance(mods, Iterable):  # type: ignore
        # List of modification types
        return [get_mod_type(mod) for mod in mods]

    raise ValueError(
        f"mods parameter must be str, list of str, or None, got {type(mods)}"
    )


def chem_mass(
    formula: Mapping[str, int | float],
    monoisotopic: bool = True,
) -> float:
    """
    Calculate the mass of a chemical formula or composition.

    :param formula: The chemical formula or composition.
    :type formula: Mapping[str, int | float] | str
    :param monoisotopic: Whether to use monoisotopic masses. Default is True.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to. Default is None.
    :type precision: int | None
    :param sep: The separator to use between element and counts. Default is ''.
    :type sep: str

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

        # Use average masses for all elements except for iosotopes.
        >>> chem_mass({'13C': 6, 'D': 12, 'O': 6}, monoisotopic=False, precision=3)
        198.186

        # Complex example using floats and particles
        >>> chem_mass('C4.45N8.22H59.99[13C34]S0.04e-16.33', monoisotopic=False, precision=3)
        672.437

        # Example Error
        >>> chem_mass("C6X2", precision=3)
        Traceback (most recent call last):
        peptacular.errors.InvalidChemFormulaError: Error parsing chem formula: "{'C': 6, 'X': 2}". Unknown element: "X"!

    """
    m: float = 0.0
    for element, count in formula.items():
        m += ELEMENT_LOOKUP[element].get_mass(monoisotopic=monoisotopic) * count
    return m
