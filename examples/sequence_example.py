import peptacular as pt

# Add mods to sequence
sequence = pt.add_mods('P[1.2345]EPTIDE', {5: 1, 'cterm': 'Amide'})
assert sequence == 'P[1.2345]EPTID[1]E-[Amide]'

# Get mods: mods are returned as a list of Mod objects Mod(value, multiplier)
mods = pt.get_mods(sequence)
assert mods == {'cterm': [pt.Mod('Amide', 1)], 0: [pt.Mod(1.2345, 1)], 5: [pt.Mod(1, 1)]}

# Strip mods
stripped_sequence = pt.strip_mods(sequence)
assert stripped_sequence == 'PEPTIDE'

# Pop mods
stripped_sequence, mods = pt.pop_mods(sequence)
assert stripped_sequence == 'PEPTIDE'
assert mods == {'cterm': [pt.Mod('Amide', 1)], 0: [pt.Mod(1.2345, 1)], 5: [pt.Mod(1, 1)]}

# Reverse sequence
reverse_sequence = pt.reverse(sequence, swap_terms=False)
assert reverse_sequence == 'ED[1]ITPEP[1.2345]-[Amide]'
