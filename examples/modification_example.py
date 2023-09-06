from peptacular.sequence import apply_static_modifications, apply_variable_modifications

# Apply static modifications (returns single sequence)
peptide = apply_static_modifications('PEPTIDE[2]', {'P': 'phospho', '(?<=P)E': 1})
assert peptide == 'P(phospho)E(1)P(phospho)TIDE[2]'

# Apply variable modifications (returns all possible combinations)
peptides = apply_variable_modifications('PEPTIDE[2]', {'P': 'phospho', '(?<=P)E': 1}, max_mods=2)
assert peptides == ['P(phospho)E(1)PTIDE[2]', 'P(phospho)EP(phospho)TIDE[2]', 'P(phospho)EPTIDE[2]',
                    'PE(1)P(phospho)TIDE[2]', 'PE(1)PTIDE[2]', 'PEP(phospho)TIDE[2]', 'PEPTIDE[2]']