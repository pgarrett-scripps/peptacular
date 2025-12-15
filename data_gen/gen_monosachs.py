from typing import Any, Generator
import warnings

import peptacular as pt

from utils import read_obo, get_id_and_name, is_obsolete


def _get_monosaccharide_entries(
    terms: list[dict[str, Any]],
) -> Generator[pt.MonosaccharideInfo, None, None]:
    for term in terms:
        term_id, term_name = get_id_and_name(term)
        term_id = term_id.replace("MONO:", "")

        if is_obsolete(term):
            continue

        property_values: dict[str, list[str]] = {}
        for val in term.get("property_value", []):
            elems = val.split('"')
            k = elems[0].rstrip()
            v = elems[1].strip()

            property_values.setdefault(k, []).append(v)

        delta_formula = property_values.get("has_chemical_formula", [None])
        delta_monoisotopic_mass = property_values.get("has_monoisotopic_mass", [None])
        delta_average_mass = property_values.get("has_average_mass", [None])

        if len(delta_formula) > 1:
            warnings.warn(
                f"[MONOSACCHARIDES] Multiple delta compositions for {term_id} {term_name} {delta_formula}"
            )
            delta_formula = delta_formula[0]
        elif len(delta_formula) == 1:
            delta_formula = delta_formula[0]

        if len(delta_monoisotopic_mass) > 1:
            warnings.warn(
                f"[MONOSACCHARIDES] Multiple delta mono masses for {term_id} {term_name} "
                f"{delta_monoisotopic_mass}"
            )
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]
        elif len(delta_monoisotopic_mass) == 1:
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]

        if len(delta_average_mass) > 1:
            warnings.warn(
                f"[MONOSACCHARIDES] Multiple delta average masses for {term_id} {term_name} "
                f"{delta_average_mass}"
            )
            delta_average_mass = delta_average_mass[0]
        elif len(delta_average_mass) == 1:
            delta_average_mass = delta_average_mass[0]


        formula = pt.ChargedFormula.from_string(f"Formula:{delta_formula}", allow_zero=True) if delta_formula else None

        # remove any elements with zero count from the composition

        comp_mass: float | None = None
        comp_avg_mass: float | None = None

        if formula:
            elements = tuple(elem for elem in formula.formula if elem.occurance != 0)

            # if charge is not None... raise an error
            if formula.charge is not None:
                raise ValueError(
                    f"Monosaccharide {term_id} {term_name} has non-zero charge in formula {delta_formula}"
                )

            formula = pt.ChargedFormula(formula=elements, charge=formula.charge)
            comp_mass = formula.get_mass(monoisotopic=True)
            comp_avg_mass = formula.get_mass(monoisotopic=False)


        delta_monoisotopic_mass = float(delta_monoisotopic_mass) if delta_monoisotopic_mass else None
        delta_average_mass = float(delta_average_mass) if delta_average_mass else None

        if delta_monoisotopic_mass is not None and comp_mass is not None:
            # assert that they are equal within 0.01 Da
            if abs(float(delta_monoisotopic_mass) - comp_mass) > 0.01:
                warnings.warn(
                    f"[MONOSACCHARIDES] Monoisotopic mass mismatch for {term_id} {term_name}: "
                    f"calculated {comp_mass}, reported {delta_monoisotopic_mass}"
                )
        if delta_average_mass is not None and comp_avg_mass is not None:
            # assert that they are equal within 0.01 Da
            if abs(float(delta_average_mass) - comp_avg_mass) > 0.01:
                warnings.warn(
                    f"[MONOSACCHARIDES] Average mass mismatch for {term_id} {term_name}: "
                    f"calculated {comp_avg_mass}, reported {delta_average_mass}"
                )

        if delta_monoisotopic_mass is None and comp_mass is not None:
            delta_monoisotopic_mass = comp_mass
        if delta_average_mass is None and comp_avg_mass is not None:
            delta_average_mass = comp_avg_mass

        yield pt.MonosaccharideInfo(
            id=term_id,
            name=term_name,
            formula=delta_formula,
            monoisotopic_mass=float(delta_monoisotopic_mass) if delta_monoisotopic_mass else None,
            average_mass=float(delta_average_mass) if delta_average_mass else None,
            dict_composition=formula.get_dict_composition() if formula else None
        )
        
def run():
    with open("data_gen/data/monosaccharides.obo", "r") as f:
        data = read_obo(f)

    with open("data_gen/data/additional_monosaccharides.obo", "r") as f:
        additional_data = read_obo(f)

    data.extend(additional_data)
    
    monosaccharides = list(_get_monosaccharide_entries(data))
    print(f"Found {len(monosaccharides)} monosaccharides")
    
    output_file = 'src/peptacular/mods/monosaccharide/data.py'
    print(f"Writing to {output_file}...")
    
    # Generate the monosaccharide entries
    entries: list[str] = []
    for ms in monosaccharides:
        # Format formula properly - None without quotes, strings with quotes
        formula_str = f'"{ms.formula}"' if ms.formula is not None else 'None'
        
        # Format composition properly - None without parse_composition, dict with it
        if ms.dict_composition is not None:
            composition_str = f'parse_composition({ms.dict_composition})'
        else:
            composition_str = 'None'
        
        entry = f'''    "{ms.name}": MonosaccharideInfo(
        id="{ms.id}",
        name="{ms.name}",
        formula={formula_str},
        monoisotopic_mass={ms.monoisotopic_mass},
        average_mass={ms.average_mass},
        dict_composition={ms.dict_composition},
    ),'''
        entries.append(entry)
    
    entries_str = '\n'.join(entries)
    
    # Write the complete file
    content = f'''"""Auto-generated monosaccharide data"""
# DO NOT EDIT - generated by gen_monosachs.py

from ...elements import parse_composition
from ..dclass import MonosaccharideInfo
import warnings


try:
    MONOSACCHARIDES: dict[str, MonosaccharideInfo] = {{
{entries_str}
    }}

except Exception as e:
    warnings.warn(
        f"Exception in monosaccharides_data: {{e}}. Using empty dictionaries.",
        UserWarning,
        stacklevel=2
    )
    MONOSACCHARIDES: dict[str, MonosaccharideInfo] = {{}} # type: ignore
'''
    
    with open(output_file, 'w') as f:
        f.write(content)
    
    print(f"✅ Successfully generated {output_file}")
    print(f"   Total entries: {len(monosaccharides)}")


if __name__ == "__main__":
    run()