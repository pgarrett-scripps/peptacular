
# Peptacular

A spectacularly simple package for working with peptide sequences.

Now level 2-ProForma compliant.

Note: This package is still in development and new released will have breaking changes. I'm hoping for v2.0.0 to be
a stable release. 

## ReadTheDocs
https://peptacular.readthedocs.io/en/latest/index.html

## Installation

```bash
pip install peptacular
```

### Proforma Notation:
- https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.1c00771/suppl_file/pr1c00771_si_001.pdf

### Modification Types:
- Modifications can be represented as strings, integers, or floats.
- During parsing, the module automatically identifies the modification type based on its representation.

## Working with Sequences

```python
import peptacular as pt

# Create modified peptide sequence, sequence can be str, float or int, 'PEPTIDE' -> 'P(1.2345)EPTID(1)E[Amide]'
peptide = pt.add_mods('PEPTIDE', {0: [1.2345], 5: [1], 'c': ['Amide']})
assert peptide == 'P[1.2345]EPTID[1]E-[Amide]'

# Parse modifications from peptide sequence, 'P[1.2345]EPTID[1]E-[Amide]' -> {0: [1.2345], 5: [1], 'c': ['Amide']}
parsed_modifications = pt.get_mods(peptide)
assert parsed_modifications == {0: [1.2345], 5: [1], 'c': ['Amide']}

# Strip modifications from peptide sequence, 'P[1.2345]EPTID[1]E-[Amide]' -> 'PEPTIDE'
stripped_peptide = pt.strip_mods(peptide)
assert stripped_peptide == 'PEPTIDE'
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
```

## Building Fragment Ions

```python
import peptacular as pt

# Calculate the m/z values for the y+ fragments.
fragments = pt.fragment('P[1.0]TIDE-[2.0]', ion_types='y', charges=1, monoisotopic=True)

# Or for multiple ion types and charges.
fragments = pt.fragment('P[1.0]EP', ion_types=['y', 'b'], charges=[1, 2], monoisotopic=True)

# Or for internal ions
fragments = pt.fragment('P[1.0]EP', ion_types=['by'], charges=[1, 2], monoisotopic=True)
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
