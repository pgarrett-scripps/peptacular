from peptacular.mass import calculate_mass, calculate_mz

# Calculate mass (monoisotopic)
peptide_mass = calculate_mass('PEP(1.0)TIDE[2.0]', charge=2)
assert peptide_mass == 804.3745169602499

# Calculate m/z (monoisotopic)
peptide_mz = calculate_mz('PEP(1.0)TIDE[2.0]', charge=2)
assert peptide_mz == 402.18725848012497

# Calculate m/z (average)
peptide_mz = calculate_mz('PEP(1.0)TIDE[2.0]', charge=2, monoisotopic=False)
assert peptide_mz == 402.42130880862