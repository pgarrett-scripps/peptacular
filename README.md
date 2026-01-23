# Peptacular

<div align="center">
  <img src="peptacular_logo.png" alt="Peptacular Logo" width="400"/>
  
  # Peptacular
  
  A Python package for peptide sequence analysis built around **ProForma 2.1 notation**. Calculate masses, generate fragments, predict isotopic patterns, and more.
  
  [![PyPI version](https://badge.fury.io/py/peptacular.svg)](https://badge.fury.io/py/peptacular)
  [![License: MIT](https://img.shields.io/badge/License-MIT-g.svg)](https://opensource.org/licenses/MIT)
  
</div>

## Features

- **Nearly Complete ProForma 2.1 Parsing**
- **Mass/MZ/Composition Calculations**
- **Predicted Isotopic Distributions**
- **Enymatic Protein Digestion** 
- **Fragment Ion Generation** 
- **Physiochemical Property Calculations** 
- **Built in Parallel Processing** 

## Installation

```bash
pip install peptacular
```

## Quick Start

See docs for more detail.

```python
import peptacular as pt

# Parse a ProForma sequence
peptide = pt.parse("PEM[Oxidation]TIDE")

# Calculate mass and m/z
mass = peptide.mass()  # 931.3854 Da
mz = peptide.mz(charge=2)      # 466.1963 (charge +2)

# Get elemental composition
comp = peptide.comp()

# Generate isotopic distribution
iso_dist = peptide.isotopic_distribution()
for iso in iso_dist[:3]:
    print(f"m/z: {iso.mass:.3f}, abundance: {iso.abundance:.3f}")

# Fragment the peptide
for frag in peptide.fragment(ion_types=['b', 'y'], charges=[1, 2]):
    print(f"{frag.ion_type}{frag.position}+{frag.charge}: {frag.mz:.3f}")

# Digest
protein = pt.parse("PEM[Oxidation]TRPEPTIDEKPEPTIDEIDE")
for span in protein.digest(pt.Proteases.TRYPSIN):
    print(f"  {protein[span].serialize()}")
```

## ProForma 2.1 Compliance

See [PROFORMA_COMPLIANCE.md](PROFORMA_COMPLIANCE.md) for detailed compliance status.

## Contributing

Contributions welcome! Check the examples directory for code style and documentation patterns.

## License

MIT

## Citation

Working on a JOSS submitions, but in the meantime use:

https://doi.org/10.5281/zenodo.15054278

