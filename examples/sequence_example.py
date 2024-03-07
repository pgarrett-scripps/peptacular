import peptacular as pt

# Create modified peptide sequence, sequence can be str, float or int, 'PEPTIDE' -> 'P(1.2345)EPTID(1)E[Amide]'
peptide = pt.add_mods('PEPTIDE', {0: [1.2345], 5: [1], 'c': ['Amide']})
assert peptide == 'P[1.2345]EPTID[1]E-[Amide]'

# Parse modifications from peptide sequence, 'P[1.2345]EPTID[1]E-[Amide]' -> {0: [1.2345], 5: [1], 'c': ['Amide']}
parsed_modifications = pt.get_mods(peptide)
assert parsed_modifications == {0: [1.2345], 5: [1], 'c': ['Amide']}

# Strip modifications from peptide sequence, 'P[1.2345]EPTID[1]E-[Amide]' -> 'PEPTIDE'
stripped_peptide = pt.strip_mods(peptide)
assert stripped_peptide == 'PEPTIDE'