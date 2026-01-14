from collections.abc import Iterable

from .data import ISOTOPES, Element
from .dclass import ElementInfo


def _handle_key_input(
    key: tuple[str | Element, int | None] | str | Element,
) -> tuple[Element, int | None]:
    """Helper to parse various key formats into (symbol, mass_number)"""
    if isinstance(key, tuple):
        if len(key) != 2:
            raise ValueError(f"Tuple key must have exactly 2 elements, got {len(key)}")
        symbol, mass_number = key

        if not isinstance(symbol, (str, Element)):
            raise TypeError(
                f"Symbol must be str or Element, got {type(symbol).__name__}"
            )
        if mass_number is not None and not isinstance(mass_number, int):
            raise TypeError(
                f"Mass number must be int or None, got {type(mass_number).__name__}"
            )

        if symbol == "D":
            symbol = "H"
            if mass_number is None:
                mass_number = 2
            if mass_number != 2:
                raise ValueError("Deuterium 'D' must have mass number 2")
        elif symbol == "T":
            symbol = "H"
            if mass_number is None:
                mass_number = 3
            if mass_number != 3:
                raise ValueError("Tritium 'T' must have mass number 3")

        if isinstance(symbol, str):
            try:
                symbol = Element(symbol)
            except ValueError:
                raise KeyError(f"Symbol '{symbol}' is not a valid Element")

        return (symbol, mass_number)
    elif isinstance(key, str):
        if not key:
            raise ValueError("Element key cannot be empty string")

        # Get digits prefix if present
        if key[0].isdigit():
            mass_number_str = ""
            i = 0
            while i < len(key) and key[i].isdigit():
                mass_number_str += key[i]
                i += 1
            if i >= len(key):
                raise ValueError(
                    f"Invalid isotope notation: '{key}' - no element symbol found"
                )
            mass_number = int(mass_number_str)
            symbol_str = key[i:]
            if not symbol_str[0].isupper():
                raise ValueError(
                    f"Invalid element symbol in '{key}' - must start with uppercase"
                )

            if symbol_str == "D":
                symbol_str = "H"
                if mass_number != 2:
                    raise ValueError("Deuterium 'D' must have mass number 2")
            elif symbol_str == "T":
                symbol_str = "H"
                if mass_number != 3:
                    raise ValueError("Tritium 'T' must have mass number 3")

            try:
                symbol = Element(symbol_str)
            except ValueError:
                raise KeyError(f"Symbol '{symbol_str}' is not a valid Element")

            return (symbol, mass_number)
        else:
            # No digit prefix - just element symbol
            if not key[0].isupper():
                raise ValueError(f"Element symbol must start with uppercase: '{key}'")
            mass_number = None
            symbol_str = key
            if symbol_str == "D":
                symbol_str = "H"
                mass_number = 2
            elif symbol_str == "T":
                symbol_str = "H"
                mass_number = 3

            try:
                symbol = Element(symbol_str)
            except ValueError:
                raise KeyError(f"Symbol '{symbol_str}' is not a valid Element")

            return (symbol, mass_number)

    elif isinstance(key, Element):
        return (key, None)

    else:
        raise TypeError(
            f"Key must be tuple[str|Element, int|None] or str or Element, got {type(key).__name__}"
        )


class ElementLookup:
    """
    Lookup class for element isotope data.

    Supports multiple lookup formats:
    - ('C', 12) -> Carbon-12
    - ('C', None) -> Most abundant carbon isotope (monoisotopic)
    - 'C' -> Most abundant carbon isotope
    - '13C' -> Carbon-13
    - 'D' -> Deuterium (2H)
    - '2H' -> Deuterium
    - 'T' -> Tritium (3H)

    If a specific isotope is not found, it will be automatically generated
    by adding/subtracting neutron masses from the monoisotopic isotope.

    The underlying data structure is a dict with keys: tuple[str, int | None]
    where the second element is the mass number, or None for monoisotopic.
    """

    # Neutron mass in Daltons
    NEUTRON_MASS = 1.00866491595

    def __init__(
        self, element_data: dict[tuple[Element, int | None], ElementInfo]
    ) -> None:
        """
        Initialize the element lookup.

        Args:
            element_data: Dictionary with keys (symbol, mass_number) where mass_number
                         can be None to indicate the monoisotopic (most abundant) isotope.
        """
        self.element_data: dict[tuple[Element, int | None], ElementInfo] = element_data

    def __getitem__(
        self, key: tuple[str | Element, int | None] | str | Element
    ) -> ElementInfo:
        """
        Get element info by various key formats.

        If a specific isotope is not found, it will be automatically generated.

        Args:
            key: Either:
                - tuple (symbol, mass_number): e.g., ('C', 12), ('C', None)
                - str with mass prefix: e.g., '13C', '2H'
                - str symbol only: e.g., 'C', 'D', 'T' (returns monoisotopic)

        Returns:
            ElementInfo for the requested isotope

        Raises:
            KeyError: If the element symbol doesn't exist at all
            ValueError: If the key format is invalid

        Examples:
            >>> lookup[('C', 12)]        # Carbon-12
            >>> lookup[('C', None)]      # Most abundant carbon (12C)
            >>> lookup['C']              # Most abundant carbon (12C)
            >>> lookup['13C']            # Carbon-13
            >>> lookup['14C']            # Carbon-14 (auto-generated if not in DB)
            >>> lookup['D']              # Deuterium (2H)
            >>> lookup['2H']             # Deuterium (same as 'D')
        """
        # Use lazy error handling - try to parse, only validate on exception
        symbol, mass_number = _handle_key_input(key)

        lookup_key = (symbol, mass_number)

        # Check if it exists
        if lookup_key in self.element_data:
            return self.element_data[lookup_key]

        raise KeyError(f"Isotope {symbol}-{mass_number} not found in ElementLookup")

    def _get_available_for_symbol(self, symbol: str | Element) -> list[int | None]:
        """Helper to get available mass numbers for an element symbol."""
        if isinstance(symbol, str):
            symbol = Element(symbol)

        available = [mass for (sym, mass) in self.element_data.keys() if sym == symbol]
        if not available:
            return []
        return sorted(
            available,
            key=lambda x: (
                x is None,
                x if x is not None else 0,
            ),  # None comes first, then sorted by mass
        )

    def __contains__(
        self, key: tuple[str | Element, int | None] | str | Element
    ) -> bool:
        """
        Check if an element/isotope exists in the lookup.

        Note: This only checks for existing entries, it does NOT trigger
        automatic isotope generation.

        Examples:
            >>> 'C' in lookup          # True
            >>> '13C' in lookup        # True if in DB
            >>> ('C', 12) in lookup    # True
            >>> 'Xx' in lookup         # False
        """
        try:
            symbol, mass_number = _handle_key_input(key)
            lookup_key = (symbol, mass_number)
            return lookup_key in self.element_data
        except (ValueError, TypeError, IndexError, KeyError):
            return False

    def __len__(self) -> int:
        """Return number of entries in the lookup."""
        return len(self.element_data)

    def __repr__(self) -> str:
        """String representation of the lookup."""
        n_elements = len({sym for sym, _ in self.element_data.keys()})
        return f"ElementLookup({len(self.element_data)} entries, {n_elements} elements)"

    def get_monoisotopic(self, symbol: str | Element) -> ElementInfo:
        """
        Get the most abundant (monoisotopic) isotope for an element.

        Args:
            symbol: Element symbol (e.g., 'C', 'H', 'N')

        Returns:
            ElementInfo for the most abundant isotope

        Examples:
            >>> lookup.get_monoisotopic('C')  # Returns 12C
            >>> lookup.get_monoisotopic('D')  # Returns 2H (deuterium)
        """

        for info in self.get_all_isotopes(symbol):
            if info.is_monoisotopic:
                return info

        raise KeyError(f"Monoisotopic isotope for '{symbol}' not found")

    def get_isotope(self, symbol: str | Element, mass_number: int) -> ElementInfo:
        """
        Get a specific isotope by symbol and mass number.

        Args:
            symbol: Element symbol (e.g., 'C', 'H')
            mass_number: Mass number (e.g., 13, 2)
            auto_generate: If True, generate missing isotopes automatically

        Returns:
            ElementInfo for the requested isotope

        Examples:
            >>> lookup.get_isotope('C', 13)  # Carbon-13
            >>> lookup.get_isotope('H', 2)   # Deuterium
            >>> lookup.get_isotope('C', 14)  # Carbon-14 (auto-generated)
        """
        if mass_number is None:
            raise ValueError("Mass number cannot be None for get_isotope()")
        try:
            element_symbol = Element(symbol) if isinstance(symbol, str) else symbol
            lookup_key: tuple[Element, int] = (element_symbol, mass_number)
        except ValueError as e:
            raise KeyError(f"Invalid element key: {symbol}-{mass_number}: {e}")
        if lookup_key not in self.element_data:
            raise KeyError(
                f"Isotope {symbol}-{mass_number} not found (auto_generate=False)"
            )
        return self.element_data[lookup_key]

    def get_all_isotopes(self, symbol: str | Element) -> list[ElementInfo]:
        """
        Get all isotopes for an element symbol (excluding the None entry).

        Args:
            symbol: Element symbol (e.g., 'C', 'H')
            include_generated: If True, include auto-generated isotopes (abundance=0)

        Returns:
            List of ElementInfo for all isotopes, sorted by mass number

        Examples:
            >>> carbon_isotopes = lookup.get_all_isotopes('C')
            >>> [iso.mass_number for iso in carbon_isotopes]
            [12, 13, 14]
        """
        if isinstance(symbol, str):
            try:
                symbol = Element(symbol)
            except ValueError:
                raise KeyError(f"Invalid element symbol '{symbol}'")

        isotopes: list[ElementInfo] = [
            elem
            for (sym, mass_number), elem in self.element_data.items()
            if sym == symbol and mass_number is not None
        ]

        if not isotopes:
            raise KeyError(f"No isotopes found for element '{symbol}'")

        return sorted(isotopes, key=lambda x: x.mass_number)

    def get_elements(self) -> list[str]:
        """
        Get list of all element symbols in the lookup.

        Returns:
            Sorted list of unique element symbols

        Examples:
            >>> lookup.get_elements()[:5]
            ['H', 'D', 'T', 'He', 'Li']
        """
        return sorted({sym for sym, _ in self.element_data.keys()})

    def mass(
        self,
        key: tuple[str | Element, int | None] | str | Element,
        monoisotopic: bool = True,
    ) -> float:
        """
        Get the mass for an element/isotope.

        IMPORTANT: If a specific isotope is provided (e.g., '13C', ('C', 13)),
        always returns the exact isotope mass regardless of monoisotopic parameter.
        The monoisotopic parameter only applies when requesting by symbol alone (e.g., 'C').

        Args:
            key: Element key (same formats as __getitem__)
            monoisotopic: Only applies when key is a symbol without mass number.
                         If True, return monoisotopic mass.
                         If False, return average mass.

        Returns:
            Mass in Daltons

        Examples:
            >>> lookup.mass('C')                      # 12.0 (monoisotopic)
            >>> lookup.mass('C', monoisotopic=False)  # 12.011 (average)
            >>> lookup.mass('13C')                    # 13.003... (always exact isotope mass)
            >>> lookup.mass('13C', monoisotopic=False) # 13.003... (still exact isotope mass!)
            >>> lookup.mass(('C', 13))                # 13.003... (always exact isotope mass)
        """
        # Parse the key to determine if specific isotope was requested
        _, mass_number = _handle_key_input(key)
        elem = self[key]

        # If a specific isotope was requested (mass_number is not None),
        # ALWAYS return the exact isotope mass
        if mass_number is not None:
            return elem.mass

        # Only symbol was provided - respect the monoisotopic parameter
        if monoisotopic:
            return elem.mass
        else:
            return elem.average_mass

    def __iter__(self) -> Iterable[ElementInfo]:
        """Iterator over all ElementInfo entries in the lookup."""
        return iter(self.element_data.values())

    def get_neutron_offsets_and_abundances(
        self, key: str | Element | ElementInfo
    ) -> list[tuple[int, float]]:
        # get the element info
        if isinstance(key, ElementInfo):
            key = key.symbol

        isotopes = self.get_all_isotopes(key)
        mono_isotope = self.get_monoisotopic(key)
        result: list[tuple[int, float]] = []
        for iso in isotopes:
            neutron_offset = iso.neutron_count - mono_isotope.neutron_count
            result.append((neutron_offset, iso.abundance))  # type: ignore
        return result

    def get_masses_and_abundances(
        self, key: str | Element | ElementInfo
    ) -> list[tuple[float, float]]:
        # get the element info

        if isinstance(key, ElementInfo):
            key = key.symbol

        isotopes = self.get_all_isotopes(key)
        result: list[tuple[float, float]] = []
        for iso in isotopes:
            result.append((iso.mass, iso.abundance))  # type: ignore
        return result

    def values(self) -> Iterable[ElementInfo]:
        """Get an iterable of all ElementInfo values in the lookup."""
        return self.element_data.values()

    def keys(self) -> Iterable[tuple[str, int | None]]:
        """Get an iterable of all keys in the lookup."""
        return self.element_data.keys()


# Create the global lookup instance
ELEMENT_LOOKUP = ElementLookup(ISOTOPES)


# {'C13': 10, 'H2': 5, 'O18': 8} -> {ElementInfo(...), 10, ElementInfo(...), 5, ElementInfo(...), 8}
def parse_composition(comp_dict: dict[str, int]) -> dict[ElementInfo, int]:
    """
    Parse a composition dictionary with string keys into ElementInfo keys.
    """
    parsed_comp: dict[ElementInfo, int] = {}

    for elem_key, count in comp_dict.items():
        elem_info = ELEMENT_LOOKUP[elem_key]
        parsed_comp[elem_info] = count

    return parsed_comp
