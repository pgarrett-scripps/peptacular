import json
import regex as re
from collections import Counter
from typing import List
from obo import read_obo
from peptacular import chem
from peptacular.chem import glycan_comp, chem_mass


def generate_gno_db():
    OBO = "data/gno.obo"
    with open(OBO, "r") as file:
        terms = read_obo(file)

    name_to_id = {}
    id_to_composition = {}
    id_to_isotopic_mass = {}
    id_to_average_mass = {}

    for term in terms:
        term_id = term.get("id")[0].split(":")[1]  # remove GNO:
        term_name = term.get("name")[0]

        property_values = {}
        for val in term.get("property_value", []):
            try:
                elems = val.split('"')
                k = elems[0].rstrip()
                v = elems[1].strip()
                property_values[k] = v
            except IndexError:
                continue

        mono = None
        ave = None

        synonyms = term.get("synonym", [])
        if not isinstance(synonyms, List):
            synonyms = [synonyms]
        synonyms = [syn.split('"')[1] for syn in synonyms]

        formula = None

        if val := property_values.get("GNO:00000202"):
            composition = Counter()

            tokens = re.findall(r"([A-Za-z0-9]+)\((\d+)\)", val)
            for symbol, count in tokens:
                try:
                    comp = glycan_comp(symbol)
                except ValueError as e:
                    print(f"Error parsing1 {term_id} {term_name} {val}, {e}")
                    composition = None
                    break

                for elem, c in comp.items():
                    composition[elem] += c * int(count)

            if composition is not None:
                formula = chem.write_chem_formula(composition)

        elif val := property_values.get("GNO:00000101"):
            continue  # This property isnt correct

            try:
                comp = glycan_comp(val)
            except ValueError as e:
                # print(f'Error parsing2 {term_id} {term_name} {val}, {e}')
                continue

            formula = chem.write_chem_formula(comp)

        mono = None
        ave = None
        if formula is not None:
            mono = chem_mass(formula, monoisotopic=True)
            ave = chem_mass(formula, monoisotopic=False)
            # print(formula, mono, ave)
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
            print(f"Warning: {term_name} already exists in id_to_composition")
            continue

        if term_id in id_to_composition:
            print(f"Warning: {term_id} already exists in id_to_composition")
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

    # print(name_to_id)
    # print(id_to_composition)
    # print(id_to_isotopic_mass)
    # print(id_to_average_mass)

    # count number of bad entries
    bad_comps = 0
    for k, v in id_to_composition.items():
        if v == None:
            bad_comps += 1

    bad_mono = 0
    for k, v in id_to_isotopic_mass.items():
        if v == None:
            bad_mono += 1

    bad_ave = 0
    for k, v in id_to_average_mass.items():
        if v == None:
            bad_ave += 1

    print("Total number of entries:", len(terms))
    print(f"Number of bad composition entries: {bad_comps}")
    print(f"Number of bad isotopic mass entries: {bad_mono}")
    print(f"Number of bad average mass entries: {bad_ave}")

    if not all(
        [name_to_id, id_to_composition, id_to_isotopic_mass, id_to_average_mass]
    ):
        raise ValueError("Error parsing UNIMOD OBO file. Check the file format.")

    with open("../src/peptacular/data/gno/id_to_chem_formula.json", "w") as f:
        json.dump(id_to_composition, f)

    with open("../src/peptacular/data/gno/id_to_isotopic_mass.json", "w") as f:
        json.dump(id_to_isotopic_mass, f)

    with open("../src/peptacular/data/gno/id_to_average_mass.json", "w") as f:
        json.dump(id_to_average_mass, f)

    with open("../src/peptacular/data/gno/name_to_id.json", "w") as f:
        json.dump(name_to_id, f)


if __name__ == "__main__":
    generate_gno_db()
