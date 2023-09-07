from peptacular.digest import digest

# Can use enzyme name from Constants.PROTEASES
peptides = digest('TIDERTIDEKT(1)IDE[2]', enzyme_regex='trypsin/P', missed_cleavages=2)
assert set(peptides) == {'TIDER', 'TIDERTIDEK', 'TIDERTIDEKT(1)IDE[2]', 'TIDEK',
                    'TIDEKT(1)IDE[2]', 'T(1)IDE[2]'}

# or use a custom regex
peptides = digest('TIDERTIDEKT(1)IDE[2]', enzyme_regex='([KR])', missed_cleavages=2)
assert set(peptides) == {'TIDER', 'TIDERTIDEK', 'TIDERTIDEKT(1)IDE[2]', 'TIDEK',
                    'TIDEKT(1)IDE[2]', 'T(1)IDE[2]'}

# Also supports semi-enzymatic digestion (might want to use min/max len to filter)
peptides = digest('TIDERTIDEKT(1)IDE[2]', enzyme_regex='trypsin/P', missed_cleavages=2,
                  semi=True, min_len=10)
assert set(peptides) == {'TIDERTIDEK', 'TIDERTIDEKT(1)IDE[2]', 'TIDERTIDEKT(1)ID', 'TIDERTIDEKT(1)I',
                    'TIDERTIDEKT(1)', 'IDERTIDEKT(1)IDE[2]', 'DERTIDEKT(1)IDE[2]',
                    'ERTIDEKT(1)IDE[2]', 'RTIDEKT(1)IDE[2]'}

# Non-specific digestion
peptides = digest('PEPTIDE', enzyme_regex='non-specific')
print(peptides)