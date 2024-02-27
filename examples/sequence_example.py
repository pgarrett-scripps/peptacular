from peptacular import sequence

# Create modified peptide sequence, sequence can be str, float or int, 'PEPTIDE' -> 'P(1.2345)EPTID(1)E[Amide]'
peptide = sequence.add_modifications('PEPTIDE', {0: 1.2345, 5: 1, 7: 'Amide'})
assert peptide == 'P(1.2345)EPTID(1)E[Amide]'

# Parse modifications from peptide sequence, 'P(1.2345)EPTID(1)E[Amide]' -> {0:1.2345, 5:1, 7:'Amide'}
parsed_modifications = sequence.get_modifications(peptide)
assert parsed_modifications == {0: 1.2345, 5: 1, 7: 'Amide'}

# Strip modifications from peptide sequence, 'P(1.2345)EPTID(1)E[Amide]' -> 'PEPTIDE'
stripped_peptide = sequence.strip_modifications(peptide)
assert stripped_peptide == 'PEPTIDE'
