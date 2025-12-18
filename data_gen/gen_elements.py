
import peptacular as pt

"""

    3He
    ('He', 3): ElementInfo(
        number=2,
        symbol='He',
        mass_number=3,
        mass=3.0160293201,
        abundance=1.34e-06,
        average_mass=3.0160293201,
        is_monoisotopic=False,
    ),
    4He
    ('He', 4): ElementInfo(
        number=2,
        symbol='He',
        mass_number=4,
        mass=4.00260325413,
        abundance=0.99999866,
        average_mass=4.00260325413,
        is_monoisotopic=True,
    ),

    He
    ('He', None): ElementInfo(
        number=2,
        symbol='He',
        mass_number=4,
        mass=4.00260325413,
        abundance=0.99999866,
        average_mass=4.002601932120929,
        is_monoisotopic=None,
    ),


"""


def construct_element_info(
    number: int | None,
    symbol: str | None,
    mass_number: int | None,
    mass: float | None,
    abundance: float | None,
    average_mass: float | None = None,
) -> pt.ElementInfo:
    if number is None:
        raise ValueError("Atomic number must be provided")
    if symbol is None:
        raise ValueError("Atomic symbol must be provided")
    if mass_number is None:
        raise ValueError("Mass number must be provided")
    if mass is None:
        raise ValueError("Relative atomic mass must be provided")
    if abundance is None:
        raise ValueError("Isotopic composition must be provided")
    if average_mass is None:
        # If not provided, use the isotope mass as fallback
        average_mass = mass

    if symbol == 'D' or symbol == 'T':
        symbol = 'H'  # Deuterium

    return pt.ElementInfo(
        number=number,
        symbol=symbol,
        mass_number=mass_number,
        mass=mass,
        abundance=abundance,
        average_mass=average_mass,
        is_monoisotopic=False
    )


def get_element_info(chem_file_path: str) -> list[pt.ElementInfo]:
    # From https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2

    element_infos: list[pt.ElementInfo] = []
    with open(chem_file_path, "r") as file:
        # Initialize variables
        atomic_number = None
        atomic_symbol = None
        mass_number = None
        relative_atomic_mass = None
        isotopic_composition = None

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
                )
                element_infos.append(e)

                # Reset variables
                atomic_number = None
                atomic_symbol = None
                mass_number = None
                relative_atomic_mass = None
                isotopic_composition = None

            else:  # Extract key and value from the current line

                try:
                    elems = line.split("=")
                    key = elems[0].rstrip()
                    value = elems[1].lstrip()
                except IndexError as e:
                    print(f"Could not parse line: {line}")
                    raise e

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
                    pass
                elif key == "Notes":
                    pass
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
            )
            element_infos.append(e)

    return element_infos


def calculate_average_masses(elements: list[pt.ElementInfo]) -> dict[int, float]:
    """
    Calculate the average mass for each element based on natural isotopic abundances.
    
    Args:
        elements: List of all isotopes
    
    Returns:
        Dictionary mapping atomic_number -> average_mass
    """
    from collections import defaultdict
    
    # Group isotopes by atomic number
    by_atomic_number: dict[int, list[pt.ElementInfo]] = defaultdict(list)
    for elem in elements:
        by_atomic_number[elem.number].append(elem)
    
    average_masses: dict[int, float] = {}
    
    for atomic_number, isotopes in by_atomic_number.items():
        # Calculate weighted average mass
        total_mass = 0.0
        total_abundance = 0.0
        
        for isotope in isotopes:
            # Only include naturally occurring isotopes (abundance > 0)
            if isotope.abundance > 0.0:
                total_mass += isotope.mass * isotope.abundance
                total_abundance += isotope.abundance
        
        if total_abundance > 0.0:
            # Normal case: weighted average
            average_mass = total_mass / total_abundance
        else:
            # All isotopes are radioactive (abundance = 0)
            # Use the monoisotopic (most stable) mass
            # Sort by mass number and take the most common one
            sorted_isotopes = sorted(isotopes, key=lambda x: x.mass_number)
            average_mass = sorted_isotopes[0].mass
        
        average_masses[atomic_number] = average_mass
    
    return average_masses


def add_average_masses(elements: list[pt.ElementInfo]) -> list[pt.ElementInfo]:
    """
    Create new pt.ElementInfo objects with average_mass populated.
    
    Args:
        elements: List of pt.ElementInfo without average_mass
    
    Returns:
        New list of pt.ElementInfo with average_mass added
    """
    average_masses = calculate_average_masses(elements)
    
    new_elements: list[pt.ElementInfo] = []
    for elem in elements:
        avg_mass = average_masses[elem.number]
        new_elem = pt.ElementInfo(
            number=elem.number,
            symbol=elem.symbol,
            mass_number=elem.mass_number,
            mass=elem.mass,
            abundance=elem.abundance,
            average_mass=avg_mass,
            is_monoisotopic=elem.is_monoisotopic
        )
        new_elements.append(new_elem)
    
    return new_elements


def add_explicit_isotope(elements: dict[tuple[str, int | None], pt.ElementInfo]) -> None:
    """Given a list of pt.ElementInfo, add explicit isotopes for each element"""
    # groupby atomic number
    # sort by abundance
    # for the most abundant isotope, use its atomic symbol without mass number and add to list
    # for example C, 12 -> C, None
    from collections import defaultdict
    grouped_elements: dict[int, list[pt.ElementInfo]] = defaultdict(list)
    for elem in elements.values():
        grouped_elements[elem.number].append(elem)

    for _, isotope_list in grouped_elements.items():
        # sort by isotopic composition
        sorted_isotopes = sorted(isotope_list, key=lambda x: x.abundance, reverse=True)
        most_abundant = sorted_isotopes[0]
        key = (most_abundant.symbol, None)  # Use None to indicate most common isotope
        elements[key] = most_abundant


def fix_hydrogen_isotopes(elements: dict[tuple[str, int | None], pt.ElementInfo]) -> None:
    """Fix Hydrogen isotopes so that H, None points to Protium"""
    # Default uses H, D, T
    # add H, None, D, None, T, None
    # add H, 2, H, 3
    protium = elements[('H', 1)]
    deuterium = elements[('D', 2)]
    tritium = elements[('T', 3)]

    elements[('H', None)] = protium
    elements[('D', None)] = deuterium
    elements[('T', None)] = tritium

    elements[('H', 2)] = deuterium
    elements[('H', 3)] = tritium


def update_monoisotopic_flags(elements: list[pt.ElementInfo]) -> list[pt.ElementInfo]:
    """Update is_monoisotopic flags based on most abundant isotope per element"""
    from collections import defaultdict

    # Group isotopes by atomic number
    by_atomic_number: dict[int, list[pt.ElementInfo]] = defaultdict(list)
    for elem in elements:
        by_atomic_number[elem.number].append(elem)

    new_elements: list[pt.ElementInfo] = []

    for _, isotopes in by_atomic_number.items():
        # Find the most abundant isotope
        most_abundant = max(isotopes, key=lambda x: x.abundance)

        for isotope in isotopes:
            is_monoisotopic = (isotope == most_abundant)
            new_elem = pt.ElementInfo(
                number=isotope.number,
                symbol=isotope.symbol,
                mass_number=isotope.mass_number,
                mass=isotope.mass,
                abundance=isotope.abundance,
                average_mass=isotope.average_mass,
                is_monoisotopic=is_monoisotopic
            )
            new_elements.append(new_elem)

    return new_elements

def add_unspecified_isotope(elements: list[pt.ElementInfo]) -> list[pt.ElementInfo]:
    # For every element group, add one which has an unspecified isotope
    new_elems: list[pt.ElementInfo] = []

    for elem in elements:
        if elem.is_monoisotopic:
            unspecified = pt.ElementInfo(
                number=elem.number,
                symbol=elem.symbol,
                mass_number=None,
                mass=elem.mass,
                abundance=None,
                average_mass=elem.average_mass,
                is_monoisotopic=None
            )
            new_elems.append(unspecified)

        new_elems.append(elem.update(average_mass=elem.mass))

    return new_elems

def gen():
    """Generate the element_data.py file with hardcoded element data"""
    
    print("\n" + "="*60)
    print("GENERATING ELEMENT DATA")
    print("="*60)
    
    data_path = 'data_gen/data/elements.txt'
    print(f"  📖 Reading from: {data_path}")
    
    # Load and process element data
    elements = get_element_info(data_path)
    print(f"  ✓ Loaded {len(elements)} isotopes")

    # Update monoisotopic flags
    elements = update_monoisotopic_flags(elements)
    print(f"  ✓ Updated monoisotopic flags")
    
    # Add average masses
    elements = add_average_masses(elements)
    print(f"  ✓ Calculated average masses")

    # add unspecified isotopes and fix hydrogen
    elements = add_unspecified_isotope(elements)
    
    # Build lookup dictionary
    element_lookup: dict[tuple[str, int | None], pt.ElementInfo] = {}
    for elem in elements:
        key = (elem.symbol, elem.mass_number)
        element_lookup[key] = elem
    
    # Add explicit isotopes and fix hydrogen
    #add_explicit_isotope(element_lookup)
    #fix_hydrogen_isotopes(element_lookup)
    print(f"Total entries after processing: {len(element_lookup)}")
    
    # Sort by atomic number, then mass number
    element_lookup = dict(sorted(
        element_lookup.items(), 
        key=lambda x: (x[1].number, x[1].mass_number if x[1].mass_number is not None else 0)
    ))
    
    # Generate the output file
    output_file = 'src/peptacular/elements/data.py'
    #output_file = 'temp_elements_data.py'  # for testing
    print(f"\n  📝 Writing to: {output_file}")
    
    # Build element entries
    entries: list[str] = []
    for key, elem in element_lookup.items():
        symbol, mass_num = key
        key_str = f"(Element.{symbol}, None)" if mass_num is None else f"(Element.{symbol}, {mass_num})"
        
        entry = f'''    {key_str}: ElementInfo(
        number={elem.number},
        symbol=Element.{symbol},
        mass_number={elem.mass_number},
        mass={elem.mass},
        abundance={elem.abundance},
        average_mass={elem.average_mass},
        is_monoisotopic={elem.is_monoisotopic},
    ),'''
        entries.append(entry)
    
    entries_str = '\n'.join(entries)
    
    # Build Element StrEnum entries using unique symbols
    unique_symbols = sorted(set(elem.symbol for elem in element_lookup.values()))
    enum_entries: list[str] = []
    for symbol in unique_symbols:
        enum_entries.append(f'    {symbol} = "{symbol}"')
    
    enum_str = '\n'.join(enum_entries)
    
    # Write the complete file
    content = f'''"""Auto-generated element isotope data"""
# DO NOT EDIT - generated by construct_data.py

from enum import StrEnum

from .dclass import ElementInfo
import warnings


class Element(StrEnum):
    """Enumeration of element symbols."""
{enum_str}

    @classmethod
    def from_str(cls, symbol: str) -> "Element":
        """Get Element enum from string"""
        return cls(symbol)


try:
    ISOTOPES: dict[tuple[Element, int | None], ElementInfo] = {{
{entries_str}
    }}

except Exception as e:
    warnings.warn(
        f"Exception in element_data: {{e}}. Using empty dictionaries.",
        UserWarning,
        stacklevel=2
    )
    ISOTOPES: dict[tuple[Element, int | None], ElementInfo] = {{}} # type: ignore
'''
    
    with open(output_file, 'w') as f:
        f.write(content)
    
    print(f"✅ Successfully generated {output_file}")
    print(f"   Total entries: {len(element_lookup)}")
    
    # Print some statistics
    natural_isotopes = sum(1 for e in element_lookup.values() if e.abundance is not None and e.abundance > 0)
    monoisotopic_entries = sum(1 for k in element_lookup.keys() if k[1] is None)
    print(f"   Natural isotopes: {natural_isotopes}")
    print(f"   Monoisotopic entries: {monoisotopic_entries}")
    print(f"   Radioactive isotopes: {len(element_lookup) - natural_isotopes - monoisotopic_entries}")


if __name__ == "__main__":
    gen()

