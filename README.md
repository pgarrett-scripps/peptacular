
# Peptacular: A Peptide Toolkit

Peptacular is a comprehensive toolkit designed for the manipulation, interpretation, and analysis of peptide sequences. It is mainly for researchers and scientists working in the field of proteomics, Peptacular provides functionalities that make it easier to handle peptide sequences, especially those with modifications commonly encountered in mass spectrometry-based proteomics.

## Installation

```bash
pip install peptacular
```

## Usage

Peptide sequences are represented as strings. Modifications are represented as a dictionary with the position of the modification as the key and the modification as the value. the sequcne module allows for easy conversion between the string and dictionary for of a modified peptide.

```python
from peptacular import sequence, mass, fragment

peptide = 'PEPTIDE'
modifications = {0: '1.2345', 5: 'Oxidation'}

modified_peptide = sequence.add_modifications(peptide, modifications)

assert modified_peptide == 'P(1.2345)EPTID(Oxidation)E'

parsed_modifications = sequence.parse_modifications(modified_peptide)  # {0: '1.2345', 4: 'Oxidation'}
stripped_peptide = sequence.strip_modifications(modified_peptide)  # 'PEPTIDE'

assert parsed_modifications == modifications
assert stripped_peptide == peptide

# calculate mass
peptide_mass = mass.calculate_mass(modified_peptide, charge=0, monoisotopic=True)

# calculate fragments
fragments = fragment.calculate_fragment_mz_series(modified_peptide, ion_type='y', charge=1, monoisotopic=True)
```


