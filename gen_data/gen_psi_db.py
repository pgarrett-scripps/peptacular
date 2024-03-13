import json
from collections import Counter
from typing import List

from obo import read_obo
from peptacular import chem


def generate_psi_mod_db():
    UNIMOD_OBO = 'data/psi-mod.obo'
    with open(UNIMOD_OBO, 'r') as file:
        terms = read_obo(file)

    name_to_id = {}
    id_to_composition = {}
    id_to_isotopic_mass = {}
    id_to_average_mass = {}

    for term in terms:

        term_id = term.get('id')[0]
        term_name = term.get('name')[0]

        property_values = {}
        for val in term.get('xref', []):
            elems = val.split('"')
            if len(elems) < 2:
                continue

            k = elems[0].rstrip().replace(':', '')
            v = elems[1].strip()
            property_values[k] = v

        delta_composition = property_values.get('DiffFormula')
        mono = property_values.get('DiffMono')
        ave = property_values.get('DiffAvg')

        synonyms = term.get('synonym', [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        term_id = term_id.split(':', 1)[1]

        formula = None
        if delta_composition:
            delta_composition = delta_composition.split(' ')

            elem_counter = Counter()
            for comp, num in zip(delta_composition[::2], delta_composition[1::2]):
                elem_counter[comp.replace('(', '').replace(')', '')] = int(num)

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

    if not all([name_to_id, id_to_composition, id_to_isotopic_mass, id_to_average_mass]):
        raise ValueError('Error parsing PSI-MOD OBO file.')

    with open('../src/peptacular/data/psi/id_to_chem_formula.json', 'w') as f:
        json.dump(id_to_composition, f)

    with open('../src/peptacular/data/psi/id_to_isotopic_mass.json', 'w') as f:
        json.dump(id_to_isotopic_mass, f)

    with open('../src/peptacular/data/psi/id_to_average_mass.json', 'w') as f:
        json.dump(id_to_average_mass, f)

    with open('../src/peptacular/data/psi/name_to_id.json', 'w') as f:
        json.dump(name_to_id, f)


if __name__ == '__main__':
    generate_psi_mod_db()
