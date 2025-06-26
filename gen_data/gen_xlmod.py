import json
from collections import Counter
from typing import List

from obo import read_obo
from peptacular import chem
from peptacular.chem import _parse_isotope_component
from peptacular import chem_mass


def generate_xlmod_db():
    OBO = 'data/xlmod.obo'
    with open(OBO, 'r') as file:
        terms = read_obo(file)

    name_to_id = {}
    id_to_composition = {}
    id_to_isotopic_mass = {}
    id_to_average_mass = {}

    for term in terms:

        term_id = term.get('id')[0]
        term_name = term.get('name')[0]

        is_obsolete = term.get('is_obsolete', ['false'])[0]
        if is_obsolete in ['true', 'false']:
            is_obsolete = bool(is_obsolete=='true')
        else:
            raise ValueError(f'Invalid is_obsolete value: {is_obsolete}')

        if is_obsolete:
            continue

        property_values = {}
        for val in term.get('property_value', []):
            elems = val.split('"')
            if len(elems) < 2:
                continue

            k = elems[0].rstrip().replace(':', '')
            v = elems[1].strip()
            property_values[k] = v

        delta_composition = property_values.get('bridgeFormula', property_values.get('deadEndFormula', None))
        mono = property_values.get('monoIsotopicMass')

        synonyms = term.get('synonym', [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        term_id = term_id.split(':', 1)[1]
        ave = None
        formula = None
        if delta_composition:
            delta_composition = delta_composition.split(' ')

            elem_counter = Counter()
            for comp in delta_composition:
                if comp.startswith('-'):
                    comp = comp[1:]
                    comps = _parse_isotope_component(comp)
                    for elem, count in comps.items():
                        elem_counter[elem] -= count
                else:
                    comps = _parse_isotope_component(comp)
                    for elem, count in comps.items():
                        elem_counter[elem] += count

            formula = chem.write_chem_formula(elem_counter)
            ave = chem_mass(formula, monoisotopic=False)

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

    print('total entries:', len(id_to_composition))
    print(f'Number of bad composition entries: {bad_comps}')
    print(f'Number of bad isotopic mass entries: {bad_mono}')
    print(f'Number of bad average mass entries: {bad_ave}')

    if not all([name_to_id, id_to_composition, id_to_isotopic_mass, id_to_average_mass]):
        raise ValueError('Error parsing PSI-MOD OBO file.')

    with open('../src/peptacular/data/xlmod/id_to_chem_formula.json', 'w') as f:
        json.dump(id_to_composition, f)

    with open('../src/peptacular/data/xlmod/id_to_isotopic_mass.json', 'w') as f:
        json.dump(id_to_isotopic_mass, f)

    with open('../src/peptacular/data/xlmod/id_to_average_mass.json', 'w') as f:
        json.dump(id_to_average_mass, f)

    with open('../src/peptacular/data/xlmod/name_to_id.json', 'w') as f:
        json.dump(name_to_id, f)


if __name__ == '__main__':
    generate_xlmod_db()
