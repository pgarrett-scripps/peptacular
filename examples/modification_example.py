import peptacular as pt

# Apply static modifications
peptide = pt.apply_static_mods("PEPTIDE-[2]", {"P": ["phospho"], "(?<=P)E": [1]})
assert peptide == "P[phospho]E[1]P[phospho]TIDE-[2]"

# Apply variable modifications
peptides = pt.apply_variable_mods(
    "PEPTIDE-[2]", {"P": [["phospho"]], "(?<=P)E": [[1]]}, max_mods=2
)
print(peptides)
assert peptides == [
    "P[phospho]E[1]PTIDE-[2]",
    "P[phospho]EP[phospho]TIDE-[2]",
    "P[phospho]EPTIDE-[2]",
    "PE[1]P[phospho]TIDE-[2]",
    "PE[1]PTIDE-[2]",
    "PEP[phospho]TIDE-[2]",
    "PEPTIDE-[2]",
]
