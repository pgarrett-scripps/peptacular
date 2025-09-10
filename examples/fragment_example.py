import peptacular as pt

# Calculate the m/z values for the y+ fragments.
pt.fragment("P[1.0]TIDE-[2.0]", ion_types="y", charges=1, monoisotopic=True)

# Or for multiple ion types and charges.
pt.fragment("P[1.0]EP", ion_types=["y", "b"], charges=[1, 2], monoisotopic=True)

# Or for internal ions
pt.fragment("P[1.0]EP", ion_types="by", charges=[1, 2], monoisotopic=True)

# Immonium ions
pt.fragment("P[1.0]EP", ion_types="i", charges=1, monoisotopic=True)

# Can also return M/Z values rather than Fragment objects
pt.fragment("P[1.0]EP", ion_types="y", charges=1, monoisotopic=True, return_type="mz")
