from peptacular import sequence

# Create modified peptide sequence, mods can be str, float or int
peptide = sequence.add_modifications('PEPTIDE', {0: 1.2345, 5: 1, 7: 'Amide'})
assert peptide == 'P(1.2345)EPTID(1)E[Amide]'

# Parse modifications from peptide sequence
parsed_modifications = sequence.get_modifications(peptide)
assert parsed_modifications == {0: 1.2345, 5: 1, 7: 'Amide'}

# Strip modifications from peptide sequence
stripped_peptide = sequence.strip_modifications(peptide)
assert stripped_peptide == 'PEPTIDE'