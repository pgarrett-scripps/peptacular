
# Peptacular

A spectacularly simple package for working with peptide sequences. 

## ReadTheDocs
https://peptacular.readthedocs.io/en/latest/index.html

## Installation

```bash
pip install peptacular
```

### Modification Notation:
- Term modifications (N-terminus and C-terminus) are specified using square brackets: `[]`.
- Residue modifications within the sequence are specified using parentheses: `()`.

### Modification Types:
- Modifications can be represented as strings, integers, or floats.
- During parsing, the module automatically identifies the modification type based on its representation.

### Example Sequences:
- `PEPT(1)IDE` - T modification of 1 (index 3)
- `PEPTIDE(1.234)[2]` - C-Term modification of 2 (index 7) and E modification of 1.234 (index 6)
- `[Acetyl]PEPTIDE(1.234)` - N-Term modification of 'Acetyl' (index -1) and E modification of 1.234 (index 6)

## Working with Sequences

```python
from peptacular import sequence

# Create modified peptide sequence, mods can be str, float or int
peptide = sequence.add_modifications('PEPTIDE', {0: 1.2345, 5: 1, 7: 'Amide'})
assert peptide == 'P(1.2345)EPTID(1)E[Amide]'

# Parse modifications from peptide sequence
parsed_modifications = sequence.get_modifications(peptide)
assert parsed_modifications == {0: 1.2345, 5: 1, 7: 'Amide'}

# Strip modifications from peptide sequence
stripped_peptide = sequence.strip_modifications(peptide)
assert stripped_peptide == 'PEPTIDE'
```

## Calculating Mass and M/Z
```python
from peptacular.mass import calculate_mass, calculate_mz

# Calculate mass (using monoisotopic mass)
peptide_mass = calculate_mass('PEP(1.0)TIDE[2.0]', charge=2)
assert peptide_mass == 804.3745169602499

# Calculate m/z (using average mass)
peptide_mz = calculate_mz('PEP(1.0)TIDE[2.0]', charge=2, monoisotopic=False)
assert peptide_mz == 402.42130880862
```

## Building fragment ions
```python
from peptacular.fragment import fragment

# Calculate m/z values for y+ fragments
fragments = fragment('P(1.0)TIDE[2.0]', ion_types='y', charges=1, monoisotopic=True)
assert list(fragments) == [577.27188355666, 479.21911970781, 378.17144123940005, 
                           265.08737726227, 150.06043423844]

# Can also specify multiple ion types and charges
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

# Also supports semi-enzymatic digestion (might want to use min/max len to filter)
peptides = digest('TIDERTIDEKT(1)IDE[2]', enzyme_regex='trypsin/P', missed_cleavages=2, 
                  semi=True, min_len=10)
assert peptides == ['TIDERTIDEK', 'TIDERTIDEKT(1)IDE[2]', 'TIDERTIDEKT(1)ID', 'TIDERTIDEKT(1)I', 
                    'TIDERTIDEKT(1)', 'IDERTIDEKT(1)IDE[2]', 'DERTIDEKT(1)IDE[2]', 
                    'ERTIDEKT(1)IDE[2]', 'RTIDEKT(1)IDE[2]']
```

## Apply Modifications
```python
from peptacular.sequence import apply_static_modifications, apply_variable_modifications

# Apply static modifications
peptide = apply_static_modifications('PEPTIDE[2]', {'P': 'phospho', '(?<=P)E': 1})
assert peptide == 'P(phospho)E(1)P(phospho)TIDE[2]'

# Apply variable modifications
peptides = apply_variable_modifications('PEPTIDE[2]', {'P': 'phospho', '(?<=P)E': 1}, max_mods=2)
assert peptides == ['P(phospho)E(1)PTIDE[2]', 'P(phospho)EP(phospho)TIDE[2]', 'P(phospho)EPTIDE[2]',
                    'PE(1)P(phospho)TIDE[2]', 'PE(1)PTIDE[2]', 'PEP(phospho)TIDE[2]', 'PEPTIDE[2]']
```
