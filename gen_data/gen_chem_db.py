import json
from copy import deepcopy
from dataclasses import dataclass
from typing import Union, List, Dict

NEUTRON_MASS = 1.00866491597


def gen_chem_db():
    # From https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2
    ELEMENTAL_INFO_FILE = "data/chem.txt"

    @dataclass
    class ElementInfo:
        atomic_number: int
        atomic_symbol: str
        mass_number: int
        relative_atomic_mass: float
        isotopic_composition: float
        standard_atomic_weight: Union[List[float], float, None]
        notes: str

        def __str__(self) -> str:
            return f"{self.mass_number}{self.atomic_symbol}"

        def average_mass_component(self) -> float:
            return self.relative_atomic_mass * self.isotopic_composition

    element_infos: List[ElementInfo] = []
    with open(ELEMENTAL_INFO_FILE, "r") as file:
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
                e = ElementInfo(atomic_number,
                                atomic_symbol,
                                mass_number,
                                relative_atomic_mass,
                                isotopic_composition,
                                standard_atomic_weight,
                                notes)
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

                elems = line.split('=')
                key = elems[0].rstrip()
                value = elems[1].lstrip()

                if key == "Atomic Number":
                    atomic_number = int(value)
                elif key == "Atomic Symbol":
                    atomic_symbol = value
                elif key == "Mass Number":
                    mass_number = int(value)
                elif key == "Relative Atomic Mass":
                    relative_atomic_mass = float(value.split('(')[0])
                elif key == "Isotopic Composition":
                    if value == '':
                        isotopic_composition = 0.0
                    else:
                        isotopic_composition = float(value.split('(')[0])
                elif key == "Standard Atomic Weight":
                    if value == '':  # Can be empty
                        standard_atomic_weight = None
                    elif value.startswith('[') and value.endswith(']'):
                        standard_atomic_weight = list(map(float, value[1:-1].split(',')))
                    else:
                        standard_atomic_weight = float(value.split('(')[0])
                elif key == "Notes":
                    notes = value
                else:
                    raise ValueError(f"Unknown key: {key}")

        # Don't forget to add the last entry after exiting the loop
        if atomic_number is not None:
            e = ElementInfo(atomic_number,
                            atomic_symbol,
                            mass_number,
                            relative_atomic_mass,
                            isotopic_composition,
                            standard_atomic_weight,
                            notes)
            element_infos.append(e)

    def map_atomic_number_to_infos(infos: List[ElementInfo]) -> Dict[int, List[ElementInfo]]:
        d = {}
        for info in infos:
            if info.atomic_number not in d:
                d[info.atomic_number] = [info]
            else:
                d[info.atomic_number].append(info)
        return d

    # map atomic number to all element infos
    atomic_number_to_infos = map_atomic_number_to_infos(element_infos)

    def map_atomic_number_to_symbol(aa_infos: Dict[int, List[ElementInfo]]) -> Dict[int, str]:
        d = {}
        for _, infos in aa_infos.items():
            infos.sort(key=lambda x: x.isotopic_composition, reverse=True)

            monoisotopic_info = infos[0]
            d[monoisotopic_info.atomic_number] = monoisotopic_info.atomic_symbol

        return d

    # map atomic number to symbol
    atomic_number_to_symbol = map_atomic_number_to_symbol(atomic_number_to_infos)

    def map_atomic_symbol_to_average_mass(aa_infos: Dict[int, List[ElementInfo]]) -> Dict[str, float]:
        d = {}
        for _, infos in aa_infos.items():
            infos.sort(key=lambda x: x.isotopic_composition, reverse=True)

            # first elem is monoisotopic
            monoisotopic_info = infos[0]
            d[monoisotopic_info.atomic_symbol] = sum([info.average_mass_component() for info in infos])

        return d

    # calculate average atomic mass for each element
    average_atomic_masses = map_atomic_symbol_to_average_mass(atomic_number_to_infos)

    def get_isotopic_atomic_compositions(aa_infos: Dict[int, List[ElementInfo]]) -> Dict[str, float]:
        d = {}
        for _, infos in aa_infos.items():
            infos.sort(key=lambda x: x.isotopic_composition, reverse=True)

            # first elem is monoisotopic
            monoisotopic_info = infos[0]
            d[monoisotopic_info.atomic_symbol] = monoisotopic_info.isotopic_composition

            for info in infos:
                d[str(info)] = info.isotopic_composition

        d['T'] = d['3T']
        d['D'] = d['2D']
        d['3H'] = d['3T']
        d['2H'] = d['2D']

        return d

    def get_isotopic_atomic_masses(aa_infos: Dict[int, List[ElementInfo]]) -> Dict[str, float]:
        d = {}
        for _, infos in aa_infos.items():
            infos.sort(key=lambda x: x.isotopic_composition, reverse=True)

            # first elem is monoisotopic
            monoisotopic_info = infos[0]
            d[monoisotopic_info.atomic_symbol] = monoisotopic_info.relative_atomic_mass

            for info in infos:
                d[str(info)] = info.relative_atomic_mass

        d['T'] = d['3T']
        d['D'] = d['2D']
        d['3H'] = d['3T']
        d['2H'] = d['2D']

        return d

    def map_atomic_number_to_comp_neutron_offset(aa_infos: Dict[int, List[ElementInfo]]) -> Dict[str, float]:
        d = {}
        for _, infos in aa_infos.items():
            infos.sort(key=lambda x: x.isotopic_composition, reverse=True)

            # first elem is monoisotopic
            monoisotopic_info = infos[0]
            d[monoisotopic_info.atomic_symbol] = []

            for info in infos: # Add all isotopes for each atomic symbol
                if info.isotopic_composition == 0.0:
                    continue
                d[monoisotopic_info.atomic_symbol].append(((info.mass_number - monoisotopic_info.mass_number), info.isotopic_composition))

            for info in infos: # Add each isotopic composition for each isotope
                d[str(info)] = [(0.0, 1.0)]

            d['T'] = d['3T']
            d['D'] = d['2D']
            d['3H'] = d['3T']
            d['2H'] = d['2D']

        return d

    def map_atomic_number_to_comp(aa_infos: Dict[int, List[ElementInfo]]) -> Dict[str, float]:
        d = {}
        for _, infos in aa_infos.items():
            infos.sort(key=lambda x: x.isotopic_composition, reverse=True)

            # first elem is monoisotopic
            monoisotopic_info = infos[0]
            d[monoisotopic_info.atomic_symbol] = []

            for info in infos: # Add all isotopes for each atomic symbol
                if info.isotopic_composition == 0.0:
                    continue
                d[monoisotopic_info.atomic_symbol].append((info.relative_atomic_mass, info.isotopic_composition))

            for info in infos: # Add each isotopic composition for each isotope
                d[str(info)] = [(info.relative_atomic_mass, 1.0)]

            d['T'] = d['3T']
            d['D'] = d['2D']
            d['3H'] = d['3T']
            d['2H'] = d['2D']

        return d

    def map_hill_order(aa_infos: Dict[int, List[ElementInfo]]) -> Dict[str, int]:

        aa_infos = deepcopy(aa_infos)

        d = []

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

        d.append('T')
        d.append('D')
        d.append('3H')
        d.append('2H')

        # add the rest alphabetically
        for atomic_number, infos in sorted(aa_infos.items(), key=lambda x: x[1][0].atomic_symbol):
            infos.sort(key=lambda x: x.isotopic_composition, reverse=True)
            d.append(infos[0].atomic_symbol)
            for info in infos:
                d.append(str(info))

        # map to index
        d = {v: k for k, v in enumerate(d)}

        return d


    atomic_number_to_comp_neutron_offset = map_atomic_number_to_comp_neutron_offset(atomic_number_to_infos)
    atomic_symbol_to_compositions = map_atomic_number_to_comp(atomic_number_to_infos)
    isotopic_atomic_masses = get_isotopic_atomic_masses(atomic_number_to_infos)

    hill_ordering = map_hill_order(atomic_number_to_infos)

    print(hill_ordering)
    print(atomic_number_to_comp_neutron_offset)
    print(atomic_symbol_to_compositions)
    print(isotopic_atomic_masses)
    print(atomic_number_to_symbol)
    print(average_atomic_masses)

    if not all([atomic_number_to_comp_neutron_offset, atomic_symbol_to_compositions, isotopic_atomic_masses, atomic_number_to_symbol, average_atomic_masses]):
        raise ValueError('Error parsing atomic data. Check the source file.')

    with open('../src/peptacular/data/element/hill_order.json', 'w') as f:
        json.dump(hill_ordering, f)

    # save to json files
    with open('../src/peptacular/data/element/isotopic_atomic_masses.json', 'w') as f:
        json.dump(isotopic_atomic_masses, f)

    with open('../src/peptacular/data/element/atomic_symbol_compositions.json', 'w') as f:
        json.dump(atomic_symbol_to_compositions, f)

    with open('../src/peptacular/data/element/atomic_symbol_neutron_offset_compositions.json', 'w') as f:
        json.dump(atomic_number_to_comp_neutron_offset, f)

    with open('../src/peptacular/data/element/atomic_number_to_symbol.json', 'w') as f:
        json.dump(atomic_number_to_symbol, f)

    with open('../src/peptacular/data/element/average_atomic_masses.json', 'w') as f:
        json.dump(average_atomic_masses, f)


if __name__ == '__main__':
    gen_chem_db()
