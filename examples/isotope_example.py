import peptacular as pt

# 1) Get the isotopic distribution for a chemical composition:
formula = {'C': 12, 'H': 6, 'N': 3}
isotopes = pt.isotopic_distribution(chemical_formula=formula, max_isotopes=3)
assert isotopes == [(192.05617, 1.0), (193.05321, 0.010959894014211729), (193.05952, 0.1297887395127868)]

# 2) Get the isotopic distribution for a peptide sequence:
sequence = 'PEPTIDE'
composition = pt.comp(sequence)
isotopes = pt.isotopic_distribution(chemical_formula=composition, max_isotopes=3)
assert isotopes ==  [(799.35997, 1.0), (800.36332, 0.3677347619528959), (801.36668, 0.06562576793973895)]


# 3) Get the estimated isotopic distribution for a given mass value
mass = 1000.0
composition = pt.estimate_comp(mass)
isotopes = pt.isotopic_distribution(chemical_formula=composition, max_isotopes=3)
assert isotopes == [(1000.0000038305802, 1.0), (1001.0033538305802, 0.47589204488021836),
                     (1002.0067138305802, 0.11066305966308926)]

# 4) By default, the isotopic_distribution function uses the masses of the elements, but it is also possible to use
#    neutron offsets from the monoisotopic peak.
sequence = 'PEPTIDE'
composition = pt.comp(sequence)
isotopes = pt.isotopic_distribution(chemical_formula=composition, max_isotopes=3, use_neutron_count=True)
assert isotopes == [(0, 1.0), (1, 0.4051174337315902), (2, 0.11084816056223826)]
