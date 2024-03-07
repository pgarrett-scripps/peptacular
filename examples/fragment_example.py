import peptacular as pt

# Calculate the m/z values for the y+ fragments.
fragments = pt.fragment('P[1.0]TIDE-[2.0]', ion_types='y', charges=1, monoisotopic=True)

# Or for multiple ion types and charges.
fragments = pt.fragment('P[1.0]EP', ion_types=['y', 'b'], charges=[1, 2], monoisotopic=True)

# Or for internal ions
fragments = pt.fragment('P[1.0]EP', ion_types=['by'], charges=[1, 2], monoisotopic=True)
