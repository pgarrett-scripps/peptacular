# Peptacular

<div align="center">
  <img src="peptacular_logo.png" alt="Peptacular Logo" width="400"/>
  
  # Peptacular
  
  A Python package for peptide sequence analysis built around **ProForma 2.0 notation**. Calculate masses, generate fragments, predict isotopic patterns, and analyze physicochemical properties.
  
  [![PyPI version](https://badge.fury.io/py/peptacular.svg)](https://badge.fury.io/py/peptacular)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  
</div>

## Features

- **ProForma 2.0 Parsing** - Parse and serialize ProForma sequences with modifications, charge states, and isotope labels
- **Mass Calculations** - Monoisotopic and average mass with modification, isotope, and charge state support
- **Isotopic Distributions** - Generate isotopic patterns with configurable resolution and abundance filtering
- **Protein Digestion** - Simulate enzymatic digestion with missed cleavages and semi-specific options
- **Fragment Generation** - Create b/y/a/x/c/z ions with neutral losses, isotopes, and charge states
- **Property Calculations** - Hydrophobicity, pI, aromaticity, secondary structure predictions
- **Format Conversion** - Import sequences from IP2, DIANN, Casanovo, and MS2PIP
- **Parallel Processing** - Automatic multiprocessing for batch operations

Peptacular uses lazy loading of data wherever possible. Internally, all mods are represented as strings and are only parsed when needed (for example, when calculating mass or composition). If you need strict validation, set validate=True, e.g. pt.parse('PEPTIDE', validate=True). Similarly, accessing the mods via their property (like annotation.nterm_mods) returns a Mods object that contains a reference to the modification strings; this is also lazily loaded and only parses the mods when needed. Values are cached where appropriate, so Peptacular should become faster with repeated calculations.

## Installation

```bash
pip install peptacular
```

## Quick Start

```python
import peptacular as pt

# Parse a ProForma sequence
peptide = pt.parse("PEM[Oxidation]TIDE/2")

# Calculate mass and m/z
mass = peptide.mass()  # 931.3854 Da
mz = peptide.mz()      # 466.1963 (charge +2)

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
protein = pt.parse("PEM[Oxidation]TRPEPTIDEKPEPTIDEIDE/2")
for span in protein.digest(pt.Proteases.TRYPSIN):
    print(f"  {protein[span].serialize()}")
```

## ProForma 2.0 Compliance

Peptacular implements **Base-ProForma** and most of **Level 2-ProForma**:

### Supported
- Unmodified sequences: `PEPTIDE`
- Named modifications (Unimod/PSI-MOD): `PEM[Oxidation]TIDE`, `PEM[MOD:00425]TIDE`
- Accession numbers: `PEM[UNIMOD:35]TIDE`
- Delta masses: `PEM[+15.995]TIDE`
- Formula modifications: `PEM[Formula:O]TIDE`
- Terminal modifications: `[Acetyl]-PEPTIDE-[Amidated]`
- Labile modifications: `{Glycan:Hex}PEPTIDE`
- Multiple modifications: `PEM[Oxidation][Dioxidation]TIDE`
- Charge states: `PEPTIDE/2`
- Charge adducts: `PEPTIDE/[Na:z+1]`, `PEPTIDE/[Na:z+1^2,H:z+1]`
- Global isotope labels: `<13C>PEPTIDE`
- Global fixed modifications: `<[Oxidation]@M>PEPTIDE`
- Ambiguous positions: `[Oxidation]?PEPTIDE`
- Position ranges: `PRT(ESFRMS)[+19.0523]ISK`
- Glycan compositions: `N[Glycan:Hex5HexNAc4]K`

### Not Yet Supported
- Cross-linking notation (XL-MOD, branch points)
- Position sets with group notation: `PEP[Oxidation#1]M[#1]AT`
- Localization scores: `PEP[Oxidation#1(0.95)]M[#1(0.05)]AT`
- RESID and GNO modifications
- Chimeric spectra notation

See [PROFORMA_COMPLIANCE.md](PROFORMA_COMPLIANCE.md) for detailed compliance status.


## Contributing

Contributions welcome! Check the examples directory for code style and documentation patterns.

## License

MIT

## Citation

https://doi.org/10.5281/zenodo.15054278

## See Also

- **ProForma** - Proteoform notation standard
- **Pyteomics** - Python proteomics toolkit
- **Biopython** - Biological computation library