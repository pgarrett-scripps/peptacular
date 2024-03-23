import peptacular as pt

# Calculate mass (monoisotopic)
peptide_mass = pt.mass('PEP[1.0]TIDE-[2.0]', charge=2, precision=3)
assert peptide_mass == 804.375

# Calculate m/z (monoisotopic)
peptide_mz = pt.mz('PEP[1.0]TIDE-[2.0]', charge=2, precision=3)
assert peptide_mz == 402.187

# Calculate m/z (average)
peptide_mz = pt.mz('PEP[1.0]TIDE-[2.0]', charge=2, monoisotopic=False, precision=3)
assert peptide_mz == 402.419

# For a given ion type
peptide_mz = pt.mz('PEP[1.0]TIDE-[2.0]', ion_type='y', charge=2, precision=3)
assert peptide_mz == 402.419