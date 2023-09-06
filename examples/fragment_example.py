from peptacular.fragment import fragment

# Calculate the m/z values for the y+ fragments. Order [y5+, y4+, ..., y1+]
fragments = fragment('P(1.0)TIDE[2.0]', ion_types='y', charges=1, monoisotopic=True)
assert fragments == [577.27188355666, 479.21911970781, 378.17144123940005,
                     265.08737726227, 150.06043423844]

# Or for multiple ion types and charges. Order: [y3+, y2+, y1+, y3++, y2++, y1++, b3+, ...,  b1++]
fragments = fragment('P(1.0)EP', ion_types=['y', 'b'], charges=[1, 2], monoisotopic=True)
assert fragments == [343.16596193614004, 245.11319808729002, 116.07060499932,
                     172.08661920145502, 123.06023727703, 58.538940733045,
                     325.15539725244, 228.10263340359, 99.06004031562,
                     163.08133685960502, 114.55495493517999, 50.033658391195004]


# Using average mass
fragment('PEPTIDE', ion_types=['y', 'b'], charges=[1, 2], monoisotopic=False)

# Calculate fragment masses instead of m/z values
fragment('PEPTIDE', ion_types=['y', 'b'], charges=[1, 2], monoisotopic=True, mz=False)
