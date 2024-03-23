
# Peptacular

A spectacularly simple package for working with peptide sequences.

## ReadTheDocs
https://peptacular.readthedocs.io/en/latest/index.html

## Installation

```bash
pip install peptacular
```

### Proforma Notation:
- https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.1c00771/suppl_file/pr1c00771_si_001.pdf

### Modification Types:
- Modifications can be represented as str, int, or float.
- During parsing, the module automatically identifies the modification type based on its representation.

## Ion Types:
 - Terminal: a, b, c, x, y, z
 - Internal: ax, ay, az, bx, by, bz, cx, cy, cz
 - Immonium: i
 - special: p, n (precursor, none)

## Working with Sequences
```python
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
```

## Calculating mass and m/z

```python
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
```

## Building Fragment Ions

```python
import peptacular as pt

# Calculate the m/z values for the y+ fragments.
pt.fragment('P[1.0]TIDE-[2.0]', ion_types='y', charges=1, monoisotopic=True)

# Or for multiple ion types and charges.
pt.fragment('P[1.0]EP', ion_types=['y', 'b'], charges=[1, 2], monoisotopic=True)

# Or for internal ions
pt.fragment('P[1.0]EP', ion_types='by', charges=[1, 2], monoisotopic=True)

# Immonium ions
pt.fragment('P[1.0]EP', ion_types='i', charges=1, monoisotopic=True)

# Can also return M/Z values rather than Fragment objects
pt.fragment('P[1.0]EP', ion_types='y', charges=1, monoisotopic=True, return_type='mz')
```

## Digesting Sequences
```python
import peptacular as pt

# Can use enzyme name from Constants.PROTEASES
peptides = pt.digest('TIDERTIDEKT[1]IDE-[2]', enzyme_regex='trypsin/P', missed_cleavages=2)
assert set(peptides) == {'TIDER', 'TIDERTIDEK', 'TIDERTIDEKT[1]IDE-[2]', 'TIDEK',
                    'TIDEKT[1]IDE-[2]', 'T[1]IDE-[2]'}

# or use a custom regex
peptides = pt.digest('TIDERTIDEKT[1]IDE-[2]', enzyme_regex='([KR])', missed_cleavages=2)
assert set(peptides) == {'TIDER', 'TIDERTIDEK', 'TIDERTIDEKT[1]IDE-[2]', 'TIDEK',
                    'TIDEKT[1]IDE-[2]', 'T[1]IDE-[2]'}

# Also supports semi-enzymatic digestion (might want to use min/max len to filter)
peptides = pt.digest('TIDERTIDEKT[1]IDE-[2]', enzyme_regex='trypsin/P', missed_cleavages=2,
                  semi=True, min_len=10)
assert set(peptides) == {'TIDERTIDEK', 'TIDERTIDEKT[1]IDE-[2]', 'TIDERTIDEKT[1]ID', 'TIDERTIDEKT[1]I',
                    'TIDERTIDEKT[1]', 'IDERTIDEKT[1]IDE-[2]', 'DERTIDEKT[1]IDE-[2]',
                    'ERTIDEKT[1]IDE-[2]', 'RTIDEKT[1]IDE-[2]'}
```

## Static and Variable Modifications

```python
import peptacular as pt

# Apply static modifications
peptide = pt.apply_static_mods('PEPTIDE-[2]', {'P': ['phospho'], '(?<=P)E': [1]})
assert peptide == 'P[phospho]E[1]P[phospho]TIDE-[2]'

# Apply variable modifications
peptides = pt.apply_variable_mods('PEPTIDE-[2]', {'P': [['phospho']], '(?<=P)E': [[1]]}, max_mods=2)
print(peptides)
assert peptides == ['P[phospho]E[1]PTIDE-[2]', 'P[phospho]EP[phospho]TIDE-[2]', 'P[phospho]EPTIDE-[2]',
                    'PE[1]P[phospho]TIDE-[2]', 'PE[1]PTIDE-[2]', 'PEP[phospho]TIDE-[2]', 'PEPTIDE-[2]']

```

## Isotopic Distribution
```python
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
```
