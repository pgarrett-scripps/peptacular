
# Peptacular: A Peptide Toolkit

Peptacular is a comprehensive toolkit designed for the manipulation, interpretation, and analysis of peptide sequences.
It is mainly for researchers and scientists working in the field of proteomics, Peptacular provides functionalities 
that make it easier to handle peptide sequences, especially those with modifications commonly encountered in mass 
spectrometry-based proteomics.

## Note

I am currently refactoring the code very often. Since I don't think anyone else is
currently using this but me I am not concerned with breaking changes. I will
update this note when I feel the code is stable enough to be used by others.

## Installation

```bash
pip install peptacular
```

## Working with Sequences
```python
from peptacular import sequence

# create modified peptide sequence, mods can be str, float or int
peptide = sequence.add_modifications('PEPTIDE', {0: 1.2345, 5: 1, 7: 'Amide'})

assert peptide == 'P(1.2345)EPTID(1)E[Amide]'

parsed_modifications = sequence.parse_modifications(peptide) 
assert parsed_modifications == {0: 1.2345, 5: 1, 7: 'Amide'}

stripped_peptide = sequence.strip_modifications(peptide)
assert stripped_peptide == 'PEPTIDE'
```

## Calculating Mass and M/Z
```python
from peptacular.mass import calculate_mass, calculate_mz

peptide_mass = calculate_mass('PEP(1.0)TIDE[2.0]', charge=0, monoisotopic=True)
assert peptide_mass == 802.3599640267099

peptide_mz = calculate_mz('PEP(1.0)TIDE[2.0]', charge=2, monoisotopic=True)
assert peptide_mz == 402.18725848012497
```

## Building fragment ions
```python
from peptacular.fragment import fragment

# Calculate single fragment ion series
fragments = fragment('P(1.0)TIDE[2.0]', ion_types='y', charges=1, monoisotopic=True)
assert list(fragments) == [577.27188355666, 479.21911970781, 378.17144123940005, 
                           265.08737726227, 150.06043423844]

# Or multiple ion series
fragments = fragment('P(1.0)EP', ion_types=['y', 'b'], charges=[1,2], monoisotopic=True)
assert list(fragments) == [343.16596193614004, 245.11319808729002, 116.07060499932, 
                           172.08661920145502, 123.06023727703, 58.538940733045, 
                           325.15539725244, 228.10263340359, 99.06004031562, 
                           163.08133685960502, 114.55495493517999, 50.033658391195004]
```

## Digesting Proteins
```python
from peptacular.digest import digest

# Can use enzyme name from Constants.PROTEASES
peptides = digest('TIDERTIDEKT(1)IDE[2]', enzyme_regex='trypsin/P', missed_cleavages=2)
assert peptides == ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKT(1)IDE[2]', 'TIDEK', 
                    'TIDEKT(1)IDE[2]', 'T(1)IDE[2]']

# or use a custom regex
peptides = digest('TIDERTIDEKT(1)IDE[2]', enzyme_regex='([KR])', missed_cleavages=2)
assert peptides == ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKT(1)IDE[2]', 'TIDEK', 
                    'TIDEKT(1)IDE[2]', 'T(1)IDE[2]']
```

## Apply Modifications
```python
from peptacular.sequence import apply_static_modifications, apply_variable_modifications

# Apply static modifications, key can be a regex string
peptide = apply_static_modifications('PEPTIDE[2]', {'P': 'phospho', '(?<=P)E': 1})
assert peptide == 'P(phospho)E(1)P(phospho)TIDE[2]'

# Apply variable modifications
peptides = apply_variable_modifications('PEPTIDE[2]', {'P': 'phospho'}, max_mods=2)
assert peptides == ['P(phospho)EP(phospho)TIDE[2]', 'P(phospho)EPTIDE[2]', 
                    'PEP(phospho)TIDE[2]', 'PEPTIDE[2]']
```