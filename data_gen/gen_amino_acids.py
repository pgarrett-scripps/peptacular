"""Generate amino acid data with calculated masses"""

from utils import calculate_mass, format_composition_string

# Raw amino acid composition data
AA_COMPOSITIONS: dict[str, dict[str, int]] = {
    "G": {"C": 2, "H": 3, "N": 1, "O": 1},  # Glycine
    "A": {"C": 3, "H": 5, "N": 1, "O": 1},  # Alanine
    "S": {"C": 3, "H": 5, "N": 1, "O": 2},  # Serine
    "P": {"C": 5, "H": 7, "N": 1, "O": 1},  # Proline
    "V": {"C": 5, "H": 9, "N": 1, "O": 1},  # Valine
    "T": {"C": 4, "H": 7, "N": 1, "O": 2},  # Threonine
    "C": {"C": 3, "H": 5, "N": 1, "O": 1, "S": 1},  # Cysteine
    "I": {"C": 6, "H": 11, "N": 1, "O": 1},  # Isoleucine
    "L": {"C": 6, "H": 11, "N": 1, "O": 1},  # Leucine
    "J": {"C": 6, "H": 11, "N": 1, "O": 1},  # Leucine or Isoleucine (ambiguous)
    "N": {"C": 4, "H": 6, "N": 2, "O": 2},  # Asparagine
    "D": {"C": 4, "H": 5, "N": 1, "O": 3},  # Aspartic acid
    "Q": {"C": 5, "H": 8, "N": 2, "O": 2},  # Glutamine
    "K": {"C": 6, "H": 12, "N": 2, "O": 1},  # Lysine
    "E": {"C": 5, "H": 7, "N": 1, "O": 3},  # Glutamic acid
    "M": {"C": 5, "H": 9, "N": 1, "O": 1, "S": 1},  # Methionine
    "H": {"C": 6, "H": 7, "N": 3, "O": 1},  # Histidine
    "F": {"C": 9, "H": 9, "N": 1, "O": 1},  # Phenylalanine
    "R": {"C": 6, "H": 12, "N": 4, "O": 1},  # Arginine
    "Y": {"C": 9, "H": 9, "N": 1, "O": 2},  # Tyrosine
    "W": {"C": 11, "H": 10, "N": 2, "O": 1},  # Tryptophan
    "U": {"C": 3, "H": 5, "N": 1, "O": 1, "Se": 1},  # Selenocysteine
    "O": {"C": 12, "H": 19, "N": 3, "O": 2},  # Pyrrolysine
    "X": {},  # Unknown amino acid (treated as 0 mass)
}

# Amino acid names
AA_NAMES: dict[str, str] = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic acid",
    "C": "Cysteine",
    "E": "Glutamic acid",
    "Q": "Glutamine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
    "B": "Asparagine or Aspartic acid",
    "Z": "Glutamine or Glutamic acid",
    "X": "Any amino acid",
    "U": "Selenocysteine",
    "O": "Pyrrolysine",
    "J": "Leucine or Isoleucine",
}

# Three letter codes
AA_THREE_LETTER: dict[str, str] = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "E": "Glu",
    "Q": "Gln",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "B": "Asx",
    "Z": "Glx",
    "X": "Xaa",
    "U": "Sec",
    "O": "Pyl",
    "J": "Xle",
}

def gen():
    """Generate the amino_acid_data.py file with calculated masses"""
    
    print("\n" + "="*60)
    print("GENERATING AMINO ACID DATA")
    print("="*60)
    
    # Calculate masses for each amino acid
    aa_data: dict[str, dict] = {}
    
    for aa_code, composition in AA_COMPOSITIONS.items():
        name = AA_NAMES[aa_code]
        three_letter = AA_THREE_LETTER[aa_code]
        
        # Determine if this is an ambiguous or mass-ambiguous amino acid
        is_ambiguous = aa_code in ["B", "Z", "J"]
        is_mass_ambiguous = aa_code in ["B", "Z"]
        
        monoisotopic_mass = calculate_mass(composition, monoisotopic=True)
        average_mass = calculate_mass(composition, monoisotopic=False)
        formula = format_composition_string(composition)
        
        aa_data[aa_code] = {
            "name": name,
            "three_letter": three_letter,
            "formula": formula,
            "composition": composition,
            "monoisotopic_mass": monoisotopic_mass,
            "average_mass": average_mass,
            "is_ambiguous": is_ambiguous,
            "is_mass_ambiguous": is_mass_ambiguous,
        }
    
    # Add B and Z which don't have compositions but have names
    for aa_code in ["B", "Z"]:
        if aa_code not in aa_data:
            aa_data[aa_code] = {
                "name": AA_NAMES[aa_code],
                "three_letter": AA_THREE_LETTER[aa_code],
                "formula": None,
                "composition": None,
                "monoisotopic_mass": None,
                "average_mass": None,
                "is_ambiguous": True,
                "is_mass_ambiguous": True,
            }
    
    # Sort by amino acid code for consistent output
    sorted_aa_codes = sorted(aa_data.keys())
    
    # Generate the output file
    output_file = 'src/peptacular/amino_acids/data.py'
    print(f"\n  📝 Writing to: {output_file}")
    
    # Build AminoAcidInfo entries
    info_entries: list[str] = []
    for aa_code in sorted_aa_codes:
        data = aa_data[aa_code]
        
        # Format optional fields
        formula = f'"{data["formula"]}"' if data["formula"] is not None else "None"
        mono_mass = f'{data["monoisotopic_mass"]:.10f}' if data["monoisotopic_mass"] is not None else "None"
        avg_mass = f'{data["average_mass"]:.10f}' if data["average_mass"] is not None else "None"
        dict_composition = data["composition"] if data["composition"] is not None else "None"
        
        entry = f'''    AminoAcid.{aa_code}: AminoAcidInfo(
        id=AminoAcid.{aa_code},
        name="{data["name"]}",
        three_letter_code="{data["three_letter"]}",
        formula={formula},
        monoisotopic_mass={mono_mass},
        average_mass={avg_mass},
        dict_composition={dict_composition},
        is_mass_ambiguous={data["is_mass_ambiguous"]},
        is_ambiguous={data["is_ambiguous"]},
    ),'''
        info_entries.append(entry)
    
    infos_str = '\n'.join(info_entries)
    
    # Build StrEnum entries using single-letter codes
    enum_entries: list[str] = []
    for aa_code in sorted_aa_codes:
        enum_entries.append(f'    {aa_code} = "{aa_code}"')
    
    enum_str = '\n'.join(enum_entries)
    
    # Write the complete file
    content = f'''"""Auto-generated amino acid data with calculated masses"""
# DO NOT EDIT - generated by gen_amino_acids.py

from enum import StrEnum

from .dclass import AminoAcidInfo


class AminoAcid(StrEnum):
    """Enumeration of amino acid single-letter codes."""
{enum_str}

    @classmethod
    def from_str(cls, aa: str) -> "AminoAcid":
        """Get AminoAcid enum from string"""
        return cls(aa.upper())


AMINO_ACID_INFOS: dict[AminoAcid, AminoAcidInfo] = {{
{infos_str}
}}
'''
    
    with open(output_file, 'w') as f:
        f.write(content)
    
    print(f"\n✅ Successfully generated {output_file}")
    print(f"   Total amino acids: {len(aa_data)}")
    
    # Print some statistics
    with_mass = sum(1 for d in aa_data.values() if d["monoisotopic_mass"] is not None)
    without_mass = len(aa_data) - with_mass
    print(f"   With calculated mass: {with_mass}")
    print(f"   Without mass (ambiguous): {without_mass}")
    
    # Print a few examples
    print("\n  📊 Example masses:")
    for aa_code in ["G", "A", "W", "U"]:
        if aa_code in aa_data:
            data = aa_data[aa_code]
            if data["monoisotopic_mass"]:
                print(f"   {aa_code} ({data['name']}): {data['monoisotopic_mass']:.6f} Da")


if __name__ == "__main__":
    gen()
