from typing import Any, Generator
import warnings

import peptacular as pt

from utils import read_obo, get_id_and_name, is_obsolete


def _get_psimod_entries(
    terms: list[dict[str, Any]],
) -> Generator[pt.PsimodInfo, None, None]:
    for term in terms:
        term_id, term_name = get_id_and_name(term)
        term_id = term_id.replace("MOD:", "")
        
        if is_obsolete(term):
            continue

        property_values: dict[str, list[str]] = {}
        for val in term.get("xref", []):
            elems = val.split('"')
            if len(elems) < 2:
                continue
            
            k = elems[0].rstrip().replace(":", "")
            v = elems[1].strip()

            property_values.setdefault(k, []).append(v)

        delta_composition = property_values.get("DiffFormula", [None])
        delta_monoisotopic_mass = property_values.get("DiffMono", [None])
        delta_average_mass = property_values.get("DiffAvg", [None])

        if len(delta_composition) > 1:
            warnings.warn(
                f"[PSI-MOD] Multiple delta compositions for {term_id} {term_name} {delta_composition}"
            )
            delta_composition = delta_composition[0]
        elif len(delta_composition) == 1:
            delta_composition = delta_composition[0]

        if len(delta_monoisotopic_mass) > 1:
            warnings.warn(
                f"[PSI-MOD] Multiple delta mono masses for {term_id} {term_name} {delta_monoisotopic_mass}"
            )
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]
        elif len(delta_monoisotopic_mass) == 1:
            delta_monoisotopic_mass = delta_monoisotopic_mass[0]

        if len(delta_average_mass) > 1:
            warnings.warn(
                f"[PSI-MOD] Multiple delta average masses for {term_id} {term_name} {delta_average_mass}"
            )
            delta_average_mass = delta_average_mass[0]
        elif len(delta_average_mass) == 1:
            delta_average_mass = delta_average_mass[0]

        # Check if values are 'none'
        if delta_monoisotopic_mass == "none":
            delta_monoisotopic_mass = None

        if delta_average_mass == "none":
            delta_average_mass = None

        if delta_composition == "none":
            delta_composition = None

        # Parse composition if available
        formula = None
        composition = None
        parsed_formula = None

        if isinstance(delta_composition, str):
            delta_composition_parts = delta_composition.split()
            
            # PSI-MOD format: "(12)C 7 H 12 (14)N 2 (16)O 1"
            # Parse element/isotope and count pairs
            formula_parts: list[str] = []
            i = 0
            while i < len(delta_composition_parts):
                elem_part = delta_composition_parts[i]
                
                # Get the count (next element in list)
                if i + 1 < len(delta_composition_parts):
                    count_str = delta_composition_parts[i + 1]
                    try:
                        count = int(count_str)
                    except ValueError:
                        # If next part is not a number, it's another element
                        count = 1
                        i += 1
                        continue
                    i += 2
                else:
                    count = 1
                    i += 1
                
                # Parse element with optional isotope: "(12)C" or "H" or "(13)C"
                if '(' in elem_part and ')' in elem_part:
                    # Has isotope notation like "(12)C"
                    isotope = elem_part[elem_part.index('(') + 1:elem_part.index(')')]
                    element = elem_part[elem_part.index(')') + 1:]
                    
                    # Format as "[13C6]" for isotope notation, with count inside brackets
                    if count == 1:
                        formula_parts.append(f"[{isotope}{element}]")
                    else:
                        # All counts (positive, negative, -1) go inside brackets
                        formula_parts.append(f"[{isotope}{element}{count}]")
                else:
                    # Regular element like "H" or "O"
                    element = elem_part
                    if count == 1:
                        formula_parts.append(element)
                    elif count == -1:
                        formula_parts.append(f"{element}-1")
                    elif count < 0:
                        formula_parts.append(f"{element}{count}")
                    else:
                        formula_parts.append(f"{element}{count}")
            
            formula_str = ''.join(formula_parts)
            try:
                parsed_formula = pt.ChargedFormula.from_string(f"Formula:{formula_str}", allow_zero=True)
                formula = formula_str
                composition = parsed_formula.get_dict_composition()
            except Exception as e:
                warnings.warn(
                    f"[PSI-MOD] Error parsing formula for {term_id} {term_name}: {delta_composition} -> {formula_str}, {e}"
                )
                formula = None
                composition = None
        
        # Convert mass strings to floats
        mono_mass: float | None = None
        avg_mass: float | None = None
        
        if delta_monoisotopic_mass is not None and isinstance(delta_monoisotopic_mass, str):
            mono_mass = float(delta_monoisotopic_mass)
        elif isinstance(delta_monoisotopic_mass, float):
            mono_mass = delta_monoisotopic_mass

        if delta_average_mass is not None and isinstance(delta_average_mass, str):
            avg_mass = float(delta_average_mass)
        elif isinstance(delta_average_mass, float):
            avg_mass = delta_average_mass
        
        # Validate formula masses against provided masses
        if parsed_formula is not None:
            calc_mono = parsed_formula.get_mass(monoisotopic=True)
            calc_avg = parsed_formula.get_mass(monoisotopic=False)
            
            if mono_mass is not None and abs(calc_mono - mono_mass) > 0.01:
                # red sign if > 1.0 Da difference
                symbol = '🔴' if abs(calc_mono - mono_mass) > 1.0 else '⚠️'
                warnings.warn(
                    f"\n  {symbol}  PSI-MOD MASS MISMATCH [{term_id}] {term_name}\n"
                    f"      Monoisotopic: calculated={calc_mono:.6f}, reported={mono_mass:.6f}\n"
                    f"      Formula: {str(formula).lstrip('Formula:')}"
                )
            
            if avg_mass is not None and abs(calc_avg - avg_mass) > 0.2:
                symbol = '⚠️⚠️' if abs(calc_avg - avg_mass) > 1.0 else '⚠️'
                warnings.warn(
                    f"\n  {symbol}  PSI-MOD MASS MISMATCH [{term_id}] {term_name}\n"
                    f"      Average: calculated={calc_avg:.6f}, reported={avg_mass:.6f}\n"
                    f"      Formula: {str(formula).lstrip('Formula:')}"
                )

        yield pt.PsimodInfo(
            id=term_id,
            name=term_name,
            formula=formula,
            monoisotopic_mass=mono_mass,
            average_mass=avg_mass,
            dict_composition=composition
        )


def run():
    print("\n" + "="*60)
    print("GENERATING PSI-MOD DATA")
    print("="*60)
    
    print("  📖 Reading from: data_gen/data/PSI-MOD.obo")
    with open("data_gen/data/PSI-MOD.obo", "r") as f:
        data = read_obo(f)
    
    unimod_entries = list(_get_psimod_entries(data))
    print(f"  ✓ Parsed {len(unimod_entries)} PSI-MOD entries")

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

    # write a file of missing entries
    missing_entries_file = 'data_gen/data/psimod_missing_entries.txt'
    print(f"\n  📄 Writing missing entries to: {missing_entries_file}")
    with open(missing_entries_file, 'w') as f:
        for mod in unimod_entries:
            if mod.monoisotopic_mass is None or mod.average_mass is None or mod.formula is None:
                f.write(f"{mod.id}\t{mod.name}\n")
    
    output_file = 'src/peptacular/mods/psimod/data.py'
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
        
        entry = f'''    "{mod.id}": PsimodInfo(
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
    content = f'''"""Auto-generated PSI-MOD data"""
# DO NOT EDIT - generated by gen_psimod.py

from ..dclass import PsimodInfo
import warnings


try:
    PSI_MODIFICATIONS: dict[str, PsimodInfo] = {{
{entries_str}
    }}

    PSI_NAME_TO_ID: dict[str, str] = {{
        mod.name: mod.id
        for mod in PSI_MODIFICATIONS.values()
    }}
except Exception as e:
    warnings.warn(
        f"Exception in psimod_data: {{e}}. Using empty dictionaries.",
        UserWarning,
        stacklevel=2
    )
    PSI_MODIFICATIONS: dict[str, PsimodInfo] = {{}} # type: ignore
    PSI_NAME_TO_ID: dict[str, str] = {{}} # type: ignore
'''
    
    with open(output_file, 'w') as f:
        f.write(content)
    
    print(f"✅ Successfully generated {output_file}")
    print(f"   Total entries: {len(unimod_entries)}")


if __name__ == "__main__":
    run()