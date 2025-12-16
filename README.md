# Peptacular

A Python package for peptide sequence analysis built around **ProForma 2.0 notation**. Calculate masses, generate fragments, predict isotopic patterns, and analyze physicochemical properties.

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
```

## Documentation by Example

The examples directory contains working code for all major features:

### Core Functionality

- **[proforma_notation.py](examples/proforma_notation.py)** - ProForma 2.0 notation parsing and serialization
- **[annotation.py](examples/annotation.py)** - Creating and manipulating peptide annotations
- **[mass_mz_comp.py](examples/mass_mz_comp.py)** - Mass, m/z, and elemental composition

### Analysis Tools

- **[isotope_examples.py](examples/isotope_examples.py)** - Isotopic distribution calculations
- **[fragment.py](examples/fragment.py)** - Fragment ion generation
- **[digest.py](examples/digest.py)** - Protein digestion simulation
- **[physiochemical_properties.py](examples/physiochemical_properties.py)** - Peptide property calculations

### Utilities

- **[converters.py](examples/converters.py)** - Format conversion between different tools
- **[parrallel.py](examples/parrallel.py)** - Batch processing with parallelization

Run all examples:
```bash
python examples/run_all_examples.py
```

## Common Use Cases

### Parse and Modify Sequences

```python
import peptacular as pt

# Parse ProForma sequences
simple = pt.parse("PEPTIDE")
modified = pt.parse("[Acetyl]-PEM[Oxidation]TIDE/2")
labeled = pt.parse("<13C>PEPTIDE")

# Programmatic modification
annot = pt.ProFormaAnnotation(sequence="PEPTIDE")
annot.set_nterm_mods({"Acetyl": 1})
annot.set_internal_mods_at_index(2, {"Oxidation": 1})
annot.set_charge(2)
print(annot.serialize())  # [Acetyl]-PEM[Oxidation]TIDE/2
```

### Mass / MZ / Composition

```python
# Calculate masses with different parameters
precursor_mass = peptide.mass(ion_type="p")
neutral_mass = peptide.mass(ion_type="n")
mz_2plus = peptide.mz(charge=2)

# With isotopes and losses
mass_c13 = peptide.mass(isotopes=1)
mass_h2o_loss = peptide.mass(losses={pt.NeutralDelta.WATER: 1})

# Adduct masses
mass_na = peptide.mass(charge='Na:z+1')
```

### Isotopic Patterns

```python
# Generate isotopic distribution
dist = peptide.isotopic_distribution(
    max_isotopes=5,
    min_abundance_threshold=0.01,
    distribution_resolution=3
)

# For fragments with modifications
dist_fragment = peptide.isotopic_distribution(
    ion_type=pt.IonType.Y,
    charge=2,
    isotopes=1,
    losses={pt.NeutralDelta.WATER: 1}
)
```

### Digestion

```python
protein = pt.parse("PEPTIDEKPEPTIDERPEPTIDER")

# Tryptic digestion
for span in protein.digest(pt.Proteases.TRYPSIN):
    peptide = protein[span]
    print(f"{peptide.serialize()} - missed cleavages: {span.missed_cleavages}")

# With missed cleavages and length filtering
for span in protein.digest(
    pt.Proteases.TRYPSIN,
    missed_cleavages=2,
    min_len=7,
    max_len=30
):
    print(protein[span].serialize())

# Semi-enzymatic
for span in protein.digest(pt.Proteases.TRYPSIN, semi=True):
    print(protein[span].serialize())
```

### Fragment Ion Generation

```python
# Generate b and y ions
for frag in peptide.fragment(ion_types=['b', 'y']):
    print(f"{frag.annotation}: {frag.mz:.3f}")

# With multiple charge states and losses
for frag in peptide.fragment(
    ion_types=['b', 'y'],
    charges=[1, 2],
    losses=[pt.NeutralDelta.WATER, pt.NeutralDelta.AMMONIA]
):
    print(f"{frag.annotation}: {frag.mz:.3f}")

# Internal fragments
for frag in peptide.fragment(ion_types=['ax'], min_len=3, max_len=5):
    print(f"{frag.annotation}: {frag.mz:.3f}")
```

### Physicochemical Properties

```python
# Simple properties
hydro = peptide.hydrophobicity
pi = peptide.pi
arom = peptide.aromaticity

# Secondary structure prediction
ss = peptide.secondary_structure()
print(f"Alpha helix: {ss['alpha_helix']:.1f}%")
print(f"Beta sheet: {ss['beta_sheet']:.1f}%")

# Custom property calculations
custom = peptide.calc_property(
    scale=pt.HydrophobicityScale.KYTE_DOOLITTLE,
    aggregation_method='avg'
)
```

### Format Conversion

```python
# From other formats
ip2_seq = pt.convert_ip2_sequence('K.PEPTIDE.K')
diann_seq = pt.convert_diann_sequence('_PEM(ox)TIDE_')
casanovo_seq = pt.convert_casanovo_sequence('+43.006PEPTIDE')

# To MS2PIP format
unmod_seq, mod_str = peptide.to_ms2_pip()
```

### Batch Processing

```python
# Process multiple sequences in parallel
sequences = ['PEPTIDE', 'PROTEIN', 'FRAGMENT', 'SEQUENCE']

# Automatic parallelization for lists
masses = pt.mass(sequences)
compositions = pt.comp(sequences)
annotations = pt.parse(sequences)

# Properties
hydrophobicity = pt.hydrophobicity(sequences)
pi_values = pt.pi(sequences)
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

## Performance

Key optimizations:

- Binary exponentiation for fast isotopic distribution calculations
- Efficient composition tracking with Counter objects
- Automatic parallelization for batch operations
- Frozen dataclasses minimize memory overhead

## Requirements

- Python 3.10+
- No external dependencies for core functionality

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