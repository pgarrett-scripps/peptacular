import json
from collections import Counter
from typing import List
from obo import read_obo
from peptacular import chem


def generate_psi_mod_db():
    UNIMOD_OBO = 'data/unimod.obo'
    with open(UNIMOD_OBO, 'r') as file:
        terms = read_obo(file)

    name_to_id = {}
    id_to_composition = {}
    id_to_isotopic_mass = {}
    id_to_average_mass = {}

    for term in terms:

        term_id = term.get('id')[0].split(':')[1]  # remove UNIMOD:
        term_name = term.get('name')[0]

        if term_name == 'unimod root node':
            continue

        property_values = {}
        for val in term.get('xref', []):
            elems = val.split('"')
            k = elems[0].rstrip()
            v = elems[1].strip()
            property_values[k] = v

        delta_composition = property_values.get('delta_composition')
        mono = property_values.get('delta_mono_mass')
        ave = property_values.get('delta_avge_mass')

        synonyms = term.get('synonym', [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        formula = None
        if delta_composition:
            delta_composition = delta_composition.split(' ')
            elem_counter = Counter()
            for comp in delta_composition:
                if '(' in comp:  # get number between ()
                    num = comp.split('(')[1].split(')')[0]
                    comp = comp.split('(')[0]
                    elem_counter[comp] += int(num)
                else:
                    elem_counter[comp] += 1

            formula = chem.write_chem_formula(elem_counter)

        """
        try:
            calc_mono_mass = mass.calculate_chem_mass(formula)
            calc_ave_mass = mass.calculate_chem_mass(formula, monoisotopic=False)
        except KeyError as e:
            try:
                calc_mono_mass = mass.calculate_glycan_mass(formula)
                calc_ave_mass = mass.calculate_glycan_mass(formula, monoisotopic=False)
            except ValueError as e:
                print(f'Error parsing {term_id} {term_name} {formula}, {e}')
                continue
        """

        if term_name in name_to_id:
            print(f'Warning: {term_name} already exists in id_to_composition')
            continue

        if term_id in id_to_composition:
            print(f'Warning: {term_id} already exists in id_to_composition')
            continue

        name_to_id[term_name] = term_id

        try:
            mono = float(mono)
        except ValueError:
            mono = None
        except TypeError:
            mono = None

        try:
            ave = float(ave)
        except ValueError:
            ave = None
        except TypeError:
            ave = None

        id_to_composition[term_id] = formula
        id_to_isotopic_mass[term_id] = mono
        id_to_average_mass[term_id] = ave

    print(name_to_id)
    print(id_to_composition)
    print(id_to_isotopic_mass)
    print(id_to_average_mass)

    # count number of bad entries
    bad_comps = 0
    for k, v in id_to_composition.items():
        if not v:
            bad_comps += 1

    bad_mono = 0
    for k, v in id_to_isotopic_mass.items():
        if not v:
            bad_mono += 1

    bad_ave = 0
    for k, v in id_to_average_mass.items():
        if not v:
            bad_ave += 1

    print(f'Number of bad composition entries: {bad_comps}')
    print(f'Number of bad isotopic mass entries: {bad_mono}')
    print(f'Number of bad average mass entries: {bad_ave}')

    if not all([name_to_id, id_to_composition, id_to_isotopic_mass, id_to_average_mass]):
        raise ValueError('Error parsing UNIMOD OBO file. Check the file format.')

    with open('../src/peptacular/data/unimod/id_to_chem_formula.json', 'w') as f:
        json.dump(id_to_composition, f)

    with open('../src/peptacular/data/unimod/id_to_isotopic_mass.json', 'w') as f:
        json.dump(id_to_isotopic_mass, f)

    with open('../src/peptacular/data/unimod/id_to_average_mass.json', 'w') as f:
        json.dump(id_to_average_mass, f)

    with open('../src/peptacular/data/unimod/name_to_id.json', 'w') as f:
        json.dump(name_to_id, f)


if __name__ == '__main__':
    generate_psi_mod_db()
