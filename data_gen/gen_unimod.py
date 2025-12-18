from collections import defaultdict
from typing import Any, Generator
import warnings

import peptacular as pt

from utils import read_obo, get_id_and_name, is_obsolete


def _get_unimod_entries(
    terms: list[dict[str, Any]],
) -> Generator[pt.UnimodInfo, None, None]:
    for term in terms:
        term_id, term_name = get_id_and_name(term)
        term_id = term_id.replace("UNIMOD:", "")
        
        if is_obsolete(term):
            continue
        
        if term_name == "unimod root node":
            continue

        property_values: dict[str, list[str]] = {}
        for val in term.get("xref", []):
            elems = val.split('"')
            k = elems[0].rstrip()
            v = elems[1].strip()

            property_values.setdefault(k, []).append(v)

        delta_composition = property_values.get("delta_composition", [None])
        delta_monoisotopic_mass = property_values.get("delta_mono_mass", [None])
        delta_average_mass = property_values.get("delta_avge_mass", [None])

        if len(delta_composition) > 1:
            warnings.warn(
                f"[UNIMOD] Multiple delta compositions for {term_id} {term_name} {delta_composition}"
            )
            delta_composition = delta_composition[0]
        elif len(delta_composition) == 1:
            delta_composition = delta_composition[0]

        if len(delta_monoisotopic_mass) > 1:
            warnings.warn(
                f"[UNIMOD] Multiple delta mono masses for {term_id} {term_name} {delta_monoisotopic_mass}"
            )
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]
        elif len(delta_monoisotopic_mass) == 1:
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]

        if len(delta_average_mass) > 1:
            warnings.warn(
                f"[UNIMOD] Multiple delta average masses for {term_id} {term_name} {delta_average_mass}"
            )
            delta_average_mass = delta_average_mass[0]
        elif len(delta_average_mass) == 1:
            delta_average_mass = delta_average_mass[0]

        # Parse composition if available
        formula = None
        comp_mass: float | None = None
        comp_avg_mass: float | None = None

        if delta_composition:
            components = delta_composition.split()
            formuals_counts: list[tuple[str, int]] = []
            for comp in components:
                if '(' in comp and ')' in comp:
                    key = comp.strip().split('(')[0]
                    count = int(comp.strip().split('(')[1].replace(')', ''))
                else:
                    key = comp.strip()
                    count = 1
                formuals_counts.append((key, count))

            # replace glycan with formuals
            for i, (key, count) in enumerate(formuals_counts):

                if key == 'HexA':
                    key = 'aHex'

                if key == 'Pent':
                    key = 'Pen'

                if key == 'Su':
                    key = 'Sug'

                if key == 'Sulf':
                    key = 'sulfate'
                
                if key == 'Ac':
                    formuals_counts[i] = ('C2H2O', count)
                    continue
                    
                if key == 'Kdn':
                    formuals_counts[i] = ('C9H14O8', count)
                    continue

                if key == 'Me':
                    formuals_counts[i] = ('CH2', count)
                    continue


                monossacharide = pt.MONOSACCHARIDE_LOOKUP.get(key)
                if monossacharide is not None:
                    glycan_formula = monossacharide.formula
                    formuals_counts[i] = (glycan_formula, count)

            composition: dict[pt.ElementInfo, int] = defaultdict(int)
            try:
                for formula_part, count in formuals_counts:

                    if count == 0:
                        continue  # skip zero counts

                    # if formula_part starts with digit enclose in parentheses
                    if formula_part[0].isdigit():
                        formula_part = f"[{formula_part}]"    
                
                    chem_formula = pt.ChargedFormula.from_string(f"Formula:{formula_part}", allow_zero=True)
                    
                    # Remove any elements with zero count from the composition
                    elements = tuple(elem for elem in chem_formula.formula if elem.occurance != 0)
                    
                    # Check for non-zero charge
                    if chem_formula.charge is not None:
                        raise ValueError(
                            f"Unimod {term_id} {term_name} has non-zero charge in formula {delta_composition}"
                        )
            
                    comp_comp = pt.ChargedFormula(formula=elements, charge=chem_formula.charge).get_composition()
                    for elem, elem_count in comp_comp.items():
                        composition[elem] += elem_count * count
            except Exception as e:
                raise ValueError(
                    f"Error parsing composition for Unimod {term_id} {term_name} with formula {delta_composition}: {e}"
                ) from e

            #print({term_name: {str(k): v for k, v in composition.items()}})
            merged_formula = pt.ChargedFormula.from_composition(composition)
            #print(f"  Merged formula: {merged_formula}")

            calc_average_mass = merged_formula.get_mass(monoisotopic=False)
            #print(f"  Calculated average mass: {calc_average_mass}")
            calc_mono_mass = merged_formula.get_mass(monoisotopic=True)
            #print(f"  Calculated monoisotopic mass: {calc_mono_mass}")

            comp_mass = calc_mono_mass
            comp_avg_mass = calc_average_mass
            formula = merged_formula
        
        delta_monoisotopic_mass = float(delta_monoisotopic_mass) if delta_monoisotopic_mass else None
        delta_average_mass = float(delta_average_mass) if delta_average_mass else None

        # Validate masses if both calculated and reported are available
        if delta_monoisotopic_mass is not None and comp_mass is not None:
            if abs(float(delta_monoisotopic_mass) - comp_mass) > 0.01:
                symbol = '🔴' if abs(float(delta_monoisotopic_mass) - comp_mass) > 1.0 else '⚠️'
                warnings.warn(
                    f"\n  {symbol}  UNIMOD MASS MISMATCH [{term_id}] {term_name}\n"
                    f"      Monoisotopic: calculated={comp_mass:.6f}, reported={delta_monoisotopic_mass}\n"
                    f"      Formula: {delta_composition}"
                )
        if delta_average_mass is not None and comp_avg_mass is not None:
            if abs(float(delta_average_mass) - comp_avg_mass) > 0.2:
                symbol = '⚠️⚠️' if abs(float(delta_average_mass) - comp_avg_mass) > 1.0 else '⚠️'
                warnings.warn(
                    f"\n  {symbol}  UNIMOD MASS MISMATCH [{term_id}] {term_name}\n"
                    f"      Average: calculated={comp_avg_mass:.6f}, reported={delta_average_mass}\n"
                    f"      Formula: {delta_composition}"
                )

        # Use calculated masses if reported masses are missing
        if delta_monoisotopic_mass is None and comp_mass is not None:
            delta_monoisotopic_mass = comp_mass
        if delta_average_mass is None and comp_avg_mass is not None:
            delta_average_mass = comp_avg_mass

        yield pt.UnimodInfo(
            id=term_id,
            name=term_name,
            formula=str(formula).lstrip('Formula:'),
            monoisotopic_mass=float(delta_monoisotopic_mass) if delta_monoisotopic_mass else None,
            average_mass=float(delta_average_mass) if delta_average_mass else None,
            dict_composition=formula.get_dict_composition() if formula else None
        )


def run():
    print("\n" + "="*60)
    print("GENERATING UNIMOD DATA")
    print("="*60)
    
    print("  📖 Reading from: data_gen/data/unimod.obo")
    with open("data_gen/data/unimod.obo", "r") as f:
        data = read_obo(f)
    
    unimod_entries = list(_get_unimod_entries(data))
    print(f"  ✓ Parsed {len(unimod_entries)} Unimod entries")

    # print stats on number of entries missing mono avg or formula
    missing_mono = sum(1 for mod in unimod_entries if mod.monoisotopic_mass is None)
    missing_avg = sum(1 for mod in unimod_entries if mod.average_mass is None)
    missing_formula = sum(1 for mod in unimod_entries if mod.formula is None)
    if missing_mono > 0 or missing_avg > 0 or missing_formula > 0:
        print(f"\n  ⚠️  Data Completeness:")
        if missing_mono > 0:
            print(f"      Missing monoisotopic mass: {missing_mono}")
        if missing_avg > 0:
            print(f"      Missing average mass: {missing_avg}")
        if missing_formula > 0:
            print(f"      Missing formula: {missing_formula}")

    
    output_file = 'src/peptacular/mods/unimod/data.py'
    print(f"\n  📝 Writing to: {output_file}")
    
    # Generate the Unimod entries
    entries: list[str] = []
    for mod in unimod_entries:
        # Format formula properly - None without quotes, strings with quotes
        formula_str = f'"{mod.formula}"' if mod.formula is not None else 'None'
        
        # Format composition properly - None without parse_composition, dict with it
        if mod.dict_composition is not None:
            composition_str = f'parse_composition({mod.dict_composition})'
        else:
            composition_str = 'None'
        
        entry = f'''    "{mod.id}": UnimodInfo(
        id="{mod.id}",
        name="{mod.name}",
        formula={formula_str},
        monoisotopic_mass={mod.monoisotopic_mass},
        average_mass={mod.average_mass},
        dict_composition={mod.dict_composition},
    ),'''
        entries.append(entry)
    
    entries_str = '\n'.join(entries)
    
    # Write the complete file
    content = f'''"""Auto-generated Unimod data"""
# DO NOT EDIT - generated by gen_unimod.py

from ..dclass import UnimodInfo
import warnings


try:
    UNIMOD_MODIFICATIONS: dict[str, UnimodInfo] = {{
{entries_str}
    }}

    UNIMOD_NAME_TO_ID: dict[str, str] = {{
        mod.name: mod.id
        for mod in UNIMOD_MODIFICATIONS.values()
    }}
except Exception as e:
    warnings.warn(
        f"Exception in unimod_data: {{e}}. Using empty dictionaries.",
        UserWarning,
        stacklevel=2
    )
    UNIMOD_MODIFICATIONS: dict[str, UnimodInfo] = {{}} # type: ignore
    UNIMOD_NAME_TO_ID: dict[str, str] = {{}} # type: ignore
'''
    
    with open(output_file, 'w') as f:
        f.write(content)
    
    print(f"✅ Successfully generated {output_file}")
    print(f"   Total entries: {len(unimod_entries)}")


if __name__ == "__main__":
    run()