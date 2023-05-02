# Peptacular

Peptacular is a Python library for simulating the enzymatic digestion of protein 
sequences. It allows users to define custom enzymes and cleavage rules using 
regular expressions and to generate peptide products with varying levels of 
missed cleavages. The library also supports semi-enzymatic and non-enzymatic digestion.

## Features
- Custom enzyme definition using regular expressions
- Control over the number of missed cleavages
- Supports semi-enzymatic and non-enzymatic digestion
- Filters peptides based on minimum and maximum length
- Converts peptide results to a pandas DataFrame for easy analysis

- ## Installation
To install Peptacular, run:

```bash
pip install peptacular
```

## Usage
Here is a basic example of using Peptacular to digest a protein sequence with a custom enzyme:

```python
from peptacular.protein import digest_protein, peptides_to_df

protein_sequence = 'PEPKTIDEPERPTIDE'

enzyme_regexes = (
    [('([KR])([^P])', 1)],

    [])

missed_cleavages = 1
min_len = 3
max_len = 20
non_enzymatic = False
semi_enzymatic = False

peptides = digest_protein(
    protein_sequence,
    enzyme_regexes,
    missed_cleavages,
    min_len,
    max_len,
    non_enzymatic,
    semi_enzymatic
)

peptide_df = peptides_to_df(peptides)
print(peptide_df)
```