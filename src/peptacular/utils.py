import sys
from typing import Protocol, runtime_checkable
from collections.abc import Iterable

from .constants import ModType, ModTypeLiteral


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
        return [get_mod_type(mods)]  # type: ignore
    elif isinstance(mods, Iterable):  # type: ignore
        # List of modification types
        return [get_mod_type(mod) for mod in mods]

    raise ValueError(
        f"mods parameter must be str, list of str, or None, got {type(mods)}"
    )



def ppm_error(theo: float, expt: float) -> float:
    """

    .. code-block:: python

        # Calculate the parts per million error between two values.
        >>> ppm_error(100.0, 100.1, 2)
        1000.0

    """
    return ((expt - theo) / theo) * 1e6


def dalton_error(theo: float, expt: float) -> float:
    """
    Calculate the Dalton error between two values.

    .. code-block:: python

        # Calculate the Dalton error between two values.
        >>> dalton_error(100.0, 100.1, 2)
        0.1

    """

    return expt - theo