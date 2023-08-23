
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

Digest.py utilizes the _spans module, which represents digested peptides as simple span objects of the original protein 
sequence. A span is simply a tuple of 3 integers: [Int, Int, Int], with the first 2 positions being used
for the start and end of the peptide within the protein sequence. The last int is used
to specify the number of missed cleavages of the peptide sequence. While you really only, need to know 
the start and end positions to recreate the peptide sequence, having the missed cleavage information can
be used to greatly speed up semi enzymatic digestion since many duplicate semi enzymatic peptide sequences
are generated when digesting peptide with the same start/stop index. To more efficiently apply semi enzymatic digestion
rules to peptides with a variable number of missed cleavages, it is most efficeint to first
group peptides by start position, then generate all semi enzymatic peptides (with that start index) based on the 
peptide with the most missed cleavages of the group. Next, you group the peptides based on
end index and do the same. This can greatly increase the efficiency of complex in silico digestion.

For Example:

```python
sequence = 'PEPKTIDEKS'
enzymatic_peptides = 'PEP', 'PEPKTIDE', 'PEPKTIDEKS', 'KTIDE', 'KTIDEKS', 'S'
spans = [0, 3, 0], []
```