"""
Element_seup.py
"""

from copy import deepcopy
from dataclasses import dataclass
from typing import Callable, TypeVar
from collections import defaultdict

T = TypeVar("T")


@dataclass
class ElementInfo:
    """
    Class to store information about an element isotope
    """

    atomic_number: int
    atomic_symbol: str
    mass_number: int
    relative_atomic_mass: float
    isotopic_composition: float
    standard_atomic_weight: float | list[float] | None
    notes: str

    def __str__(self) -> str:
        return f"{self.mass_number}{self.atomic_symbol}"

    def average_mass_component(self) -> float:
        return self.relative_atomic_mass * self.isotopic_composition


def construct_element_info(
    atomic_number: int | None,
    atomic_symbol: str | None,
    mass_number: int | None,
    relative_atomic_mass: float | None,
    isotopic_composition: float | None,
    standard_atomic_weight: list[float] | float | None,
    notes: str | None,
) -> ElementInfo:

    if atomic_number is None:
        raise ValueError("Atomic number must be provided")
    if atomic_symbol is None:
        raise ValueError("Atomic symbol must be provided")
    if mass_number is None:
        raise ValueError("Mass number must be provided")
    if relative_atomic_mass is None:
        raise ValueError("Relative atomic mass must be provided")
    if isotopic_composition is None:
        raise ValueError("Isotopic composition must be provided")

    return ElementInfo(
        atomic_number=atomic_number,
        atomic_symbol=atomic_symbol,
        mass_number=mass_number,
        relative_atomic_mass=relative_atomic_mass,
        isotopic_composition=isotopic_composition,
        standard_atomic_weight=standard_atomic_weight,
        notes=notes if notes is not None else "",
    )

def get_element_info(chem_file_path: str) -> list[ElementInfo]:
    # From https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2

    element_infos: list[ElementInfo] = []
    with open(chem_file_path, "r") as file:
        # Initialize variables
        atomic_number = None
        atomic_symbol = None
        mass_number = None
        relative_atomic_mass = None
        isotopic_composition = None
        standard_atomic_weight = None
        notes = None

        # Process each line
        for line in file:
            line = line.strip()  # Remove leading and trailing whitespace
            if line == "":  # New block starts after an empty line

                e = construct_element_info(
                    atomic_number,
                    atomic_symbol,
                    mass_number,
                    relative_atomic_mass,
                    isotopic_composition,
                    standard_atomic_weight,
                    notes,
                )
                element_infos.append(e)

                # Reset variables
                atomic_number = None
                atomic_symbol = None
                mass_number = None
                relative_atomic_mass = None
                isotopic_composition = None
                standard_atomic_weight = None
                notes = None

            else:  # Extract key and value from the current line

                elems = line.split("=")
                key = elems[0].rstrip()
                value = elems[1].lstrip()

                if key == "Atomic Number":
                    atomic_number = int(value)
                elif key == "Atomic Symbol":
                    atomic_symbol = value
                elif key == "Mass Number":
                    mass_number = int(value)
                elif key == "Relative Atomic Mass":
                    relative_atomic_mass = float(value.split("(")[0])
                elif key == "Isotopic Composition":
                    if value == "":
                        isotopic_composition = 0.0
                    else:
                        isotopic_composition = float(value.split("(")[0])
                elif key == "Standard Atomic Weight":
                    if value == "":  # Can be empty
                        standard_atomic_weight = None
                    elif value.startswith("[") and value.endswith("]"):
                        standard_atomic_weight = list(
                            map(float, value[1:-1].split(","))
                        )
                    else:
                        standard_atomic_weight = float(value.split("(")[0])
                elif key == "Notes":
                    notes = value
                else:
                    raise ValueError(f"Unknown key: {key}")

        # Don't forget to add the last entry after exiting the loop
        if atomic_number is not None:
            e = construct_element_info(
                atomic_number,
                atomic_symbol,
                mass_number,
                relative_atomic_mass,
                isotopic_composition,
                standard_atomic_weight,
                notes,
            )
            element_infos.append(e)

    return element_infos


def _map_atomic_number_to_infos(
    infos: list[ElementInfo],
) -> dict[int, list[ElementInfo]]:
    d: dict[int, list[ElementInfo]] = defaultdict(list)
    for info in infos:
        d[info.atomic_number].append(info)
    return dict(d)


def _add_isotope_aliases(d: dict[str, T]) -> None:
    """Add common isotope aliases to the dictionary."""
    alias_mappings = [("T", "3T"), ("D", "2D"), ("3H", "3T"), ("2H", "2D")]

    for alias, target in alias_mappings:
        if target in d:
            d[alias] = d[target]


def _process_element_infos(
    elem_infos: list[ElementInfo],
    atomic_symbol_processor: Callable[[list[ElementInfo]], T],
    individual_isotope_processor: Callable[[ElementInfo], T],
    add_aliases: bool = True,
) -> dict[str, T]:
    """
    Generic function to process element infos with different processors.

    Args:
        elem_infos: List of ElementInfo objects
        atomic_symbol_processor: Function to process all infos for an atomic symbol
        individual_isotope_processor: Function to process individual isotope info
        add_aliases: Whether to add isotope aliases
    """
    # Map atomic number to all element infos
    aa_infos = _map_atomic_number_to_infos(elem_infos)

    d: dict[str, T] = {}

    for _, infos in aa_infos.items():
        infos.sort(key=lambda x: x.isotopic_composition, reverse=True)

        # First elem is monoisotopic
        monoisotopic_info = infos[0]

        # Process all isotopes for the atomic symbol
        d[monoisotopic_info.atomic_symbol] = atomic_symbol_processor(infos)

        # Add each individual isotope
        for info in infos:
            d[str(info)] = individual_isotope_processor(info)

    if add_aliases:
        _add_isotope_aliases(d)

    return d


def map_atomic_number_to_symbol(elem_infos: list[ElementInfo]) -> dict[int, str]:
    """Map atomic number to atomic symbol (using monoisotopic info)."""
    aa_infos = _map_atomic_number_to_infos(elem_infos)

    d: dict[int, str] = {}
    for _, infos in aa_infos.items():
        infos.sort(key=lambda x: x.isotopic_composition, reverse=True)
        monoisotopic_info = infos[0]
        d[monoisotopic_info.atomic_number] = monoisotopic_info.atomic_symbol

    return d


def map_atomic_symbol_to_average_mass(
    elem_infos: list[ElementInfo],
) -> dict[str, float]:
    """Map atomic symbol to average mass."""

    def process_atomic_symbol(infos: list[ElementInfo]) -> float:
        monoisotopic_info = infos[0]
        ave_mass = sum(info.average_mass_component() for info in infos)
        return monoisotopic_info.relative_atomic_mass if ave_mass == 0 else ave_mass

    def process_individual_isotope(info: ElementInfo) -> float:
        return info.relative_atomic_mass

    return _process_element_infos(
        elem_infos, process_atomic_symbol, process_individual_isotope, add_aliases=False
    )


def get_isotopic_atomic_compositions(elem_infos: list[ElementInfo]) -> dict[str, float]:
    """Get isotopic atomic compositions."""

    def process_atomic_symbol(infos: list[ElementInfo]) -> float:
        return infos[0].isotopic_composition  # monoisotopic

    def process_individual_isotope(info: ElementInfo) -> float:
        return info.isotopic_composition

    return _process_element_infos(
        elem_infos, process_atomic_symbol, process_individual_isotope
    )


def get_isotopic_atomic_masses(elem_infos: list[ElementInfo]) -> dict[str, float]:
    """Get isotopic atomic masses."""

    def process_atomic_symbol(infos: list[ElementInfo]) -> float:
        return infos[0].relative_atomic_mass  # monoisotopic

    def process_individual_isotope(info: ElementInfo) -> float:
        return info.relative_atomic_mass

    return _process_element_infos(
        elem_infos, process_atomic_symbol, process_individual_isotope
    )


def map_atomic_number_to_comp_neutron_offset(
    elem_infos: list[ElementInfo],
) -> dict[str, list[tuple[int, float]]]:
    """Map atomic number to composition with neutron offset."""

    def process_atomic_symbol(infos: list[ElementInfo]) -> list[tuple[int, float]]:
        monoisotopic_info = infos[0]
        result: list[tuple[int, float]] = []
        for info in infos:
            if info.isotopic_composition != 0.0:
                offset = info.mass_number - monoisotopic_info.mass_number
                result.append((offset, info.isotopic_composition))
        return result

    def process_individual_isotope(info: ElementInfo) -> list[tuple[int, float]]:
        return [(0, 1.0)]

    return _process_element_infos(
        elem_infos, process_atomic_symbol, process_individual_isotope
    )


def map_atomic_number_to_comp(
    elem_infos: list[ElementInfo],
) -> dict[str, list[tuple[float, float]]]:
    """Map atomic number to composition (mass, abundance pairs)."""

    def process_atomic_symbol(infos: list[ElementInfo]) -> list[tuple[float, float]]:
        result: list[tuple[float, float]] = []
        for info in infos:
            if info.isotopic_composition != 0.0:
                result.append((info.relative_atomic_mass, info.isotopic_composition))
        return result

    def process_individual_isotope(info: ElementInfo) -> list[tuple[float, float]]:
        return [(info.relative_atomic_mass, 1.0)]

    return _process_element_infos(
        elem_infos, process_atomic_symbol, process_individual_isotope
    )

def map_hill_order(elem_infos: list[ElementInfo]) -> dict[str, int]:
    # map atomic number to all element infos
    aa_infos = _map_atomic_number_to_infos(elem_infos)

    aa_infos = deepcopy(aa_infos)

    d: list[str] = []

    # make carbon first
    carbon_infos = aa_infos.pop(6)
    carbon_infos.sort(key=lambda x: x.isotopic_composition, reverse=True)

    monoisotopic_carbon = carbon_infos[0]
    d.append(monoisotopic_carbon.atomic_symbol)
    d.extend([str(info) for info in carbon_infos])

    # make hydrogen second
    hydrogen_infos = aa_infos.pop(1)
    hydrogen_infos.sort(key=lambda x: x.isotopic_composition, reverse=True)
    monoisotopic_hydrogen = hydrogen_infos[0]
    d.append(monoisotopic_hydrogen.atomic_symbol)
    d.extend([str(info) for info in hydrogen_infos])

    d.append("T")
    d.append("D")
    d.append("3H")
    d.append("2H")

    # add the rest alphabetically
    for _, infos in sorted(aa_infos.items(), key=lambda x: x[1][0].atomic_symbol):
        infos.sort(key=lambda x: x.isotopic_composition, reverse=True)
        d.append(infos[0].atomic_symbol)
        for info in infos:
            d.append(str(info))

    return {v: k for k, v in enumerate(d)}
