import json
from peptacular.mass import chem_mass
from typing import List
from obo import read_obo


def generate_monosaccharide_db():
    UNIMOD_OBO = 'data/monosaccharides.obo'
    with open(UNIMOD_OBO, 'r') as file:
        terms = read_obo(file)

    name_to_id = {}
    id_to_composition = {}
    id_to_isotopic_mass = {}
    id_to_average_mass = {}

    for term in terms:

        term_id = term.get('id')[0].split(':')[1]  # remove MONO:
        term_name = term.get('name')[0]

        property_values = {}
        for val in term.get('property_value', []):
            elems = val.split('"')
            k = elems[0].rstrip()
            v = elems[1].strip()
            property_values[k] = v

        delta_composition = property_values.get('has_chemical_formula')
        mono = property_values.get('has_monoisotopic_mass')

        synonyms = term.get('synonym', [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        # try to parse formula
        try:
            mono_mass = chem_mass(delta_composition, monoisotopic=True)
            ave_mass = chem_mass(delta_composition, monoisotopic=False)
        except:
            raise ValueError(f'Error parsing {term_id} {term_name} {mono} {delta_composition}')

        name_to_id[term_name] = term_id

        id_to_composition[term_id] = delta_composition
        id_to_isotopic_mass[term_id] = float(mono)
        id_to_average_mass[term_id] = float(ave_mass)


    print(name_to_id)
    print(id_to_composition)
    print(id_to_isotopic_mass)
    print(id_to_average_mass)

    if not all([name_to_id, id_to_composition, id_to_isotopic_mass, id_to_average_mass]):
        raise ValueError('Error parsing monosaccharides')

    with open('../src/peptacular/data/monosaccharide/id_to_isotopic_compositions.json', 'w') as f:
        json.dump(id_to_composition, f)

    with open('../src/peptacular/data/monosaccharide/id_to_isotopic_mass.json', 'w') as f:
        json.dump(id_to_isotopic_mass, f)

    with open('../src/peptacular/data/monosaccharide/id_to_average_mass.json', 'w') as f:
        json.dump(id_to_average_mass, f)

    with open('../src/peptacular/data/monosaccharide/name_to_id.json', 'w') as f:
        json.dump(name_to_id, f)


if __name__ == '__main__':
    generate_monosaccharide_db()
