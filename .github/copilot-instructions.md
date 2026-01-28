


# ProForma Notation - Basic Summary


1 - Never make summary documentation unles specifically asked.
2 - check makfile for commands

## Documentation & Comments

### Docstring Format

Use **Google-style docstrings** but keep them minimal - type hints handle the rest.

**Simple function:**
```python
def calculate_mass(sequence: str, charge: int = 1) -> float:
    """Calculate the mass-to-charge ratio of a peptide."""
```

**When you need more detail:**
```python
def find_isotopes(mz: float, tolerance: float = 0.01) -> list[Peak]:
    """Find isotopic peaks within the tolerance window.
    
    Uses a greedy algorithm to identify the most intense peaks first,
    then searches for their isotopic patterns.
    """
```

**Classes:**
```python
class Peptide:
    """Represents a peptide sequence with ProForma modifications."""
```

### What to Document

- **One-line summary** for all public functions/classes
- **Additional details** only when the implementation is non-obvious
- **Don't repeat** what's already in type hints
- **Private functions** (`_name`) can skip docstrings if obvious

### Building Docs
```bash
cd docs
make html
# View at docs/_build/html/index.html
```

see **proforma.schema.json** for the full ProForma 2.0 json object specification.

## What is ProForma?

ProForma is a **standardized text notation for representing peptides and proteins with modifications**. It's designed to be both human-readable and machine-parsable, allowing scientists to precisely describe modified peptide sequences in mass spectrometry data.

## Core Concept

Think of it as a way to write: **"amino acid sequence + where modifications are located + what those modifications are"**

## Basic Examples

### 1. Simple Unmodified Peptide
```
PEPTIDE
```
Just amino acids using standard one-letter codes (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)

### 2. Peptide with Modification
```
PEM[Oxidation]TIDE
```
- Methionine (M) is oxidized
- Modifications go in square brackets `[]` right after the modified amino acid

### 3. Multiple Modifications
```
PEM[Oxidation]TIS[Phospho]DE
```
- M is oxidized
- S is phosphorylated

### 4. Terminal Modifications
```
[Acetyl]-PEPTIDE
[iTRAQ4plex]-PEPTIDE-[Amidated]
```
- N-terminal modifications: `[mod]-` before sequence
- C-terminal modifications: `-[mod]` after sequence

## Ways to Specify Modifications

ProForma supports multiple ways to describe the same modification:

```
EM[Oxidation]TIDE              # By name (Unimod)
EM[UNIMOD:35]TIDE              # By accession number
EM[+15.995]TIDE                # By mass change
EM[Formula:O]TIDE              # By chemical formula
```

## Key Advanced Features

### Ambiguous Modification Position
When you know a modification exists but not exactly where:
```
[Phospho]?PEPTIDE              # Phospho is somewhere, location unknown
```

### Multiple Possible Sites
```
PEP[Phospho#g1]TIS[#g1]DE     # Phospho is on either T or S
```

### Labile Modifications
Modifications that fall off during fragmentation:
```
{Glycan:Hex}PEPTIDE            # Glycan present but lost in MS2
```

### Cross-linked Peptides
This is somewhat handled at the parsing level but will not will not be implmented in the codebased. Dont worry about this too much.
```
PEPTK[#XL1]IDE//SEQK[#XL1]    # Two peptides linked together
```

### Chimeric Spectra
Multiple peptides in same spectrum:
This is somewhat handled at the parsing level but will not will not be implmented int eh codebased. Dont worry about this too much.
```
PEPTIDE+SEQUENCE               # Two co-eluting peptides
```

### Charge States
```
PEPTIDE/2                      # Charge state +2
```

### Charge Adducts
```
PEPTIDE/[Na+:z+1]              # Sodium adduct with +1 charge
PEPTIDE/[Na+:z+1^2]            # added 2 times (total charge: +2)
EPT[Formula:Zn:z+2]IDE/[Na:z+1^2] # total +4

```

both charge and charge adduct cannot occur simultaneously.

```
PEPTIDE/[Na+:z+1^2]              # Sodium adduct with +1 charge 2 times
```



## Compliance Levels

ProForma has different levels of complexity:

1. **Base-ProForma** - Simple sequences with basic modifications
2. **Level 2-ProForma** - Adds ambiguity, formulas, delta masses
3. **Extensions** - Specialized features for:
   - Top-down proteomics
   - Cross-linking
   - Glycoproteomics
   - Advanced complexity

## Common Use Cases

### Bottom-up Proteomics
```
[Acetyl]-EM[Oxidation]EVTSES[Phospho]PEK
```
Typical tryptic peptide with PTMs

### Top-down Proteomics
```
<[Oxidation]@M>FULLPROTEINSEQUENCE...
```
Full protein with fixed modifications

### Glycopeptide
```
NEEYN[Glycan:Hex5HexNAc4]K
```
N-glycosylation site

### Cross-linking
```
PEPTK[XLMOD:02001#XL1]IDE//SEQK[#XL1]
```
DSS cross-link between two lysines

## Why ProForma?

**Before ProForma:** Everyone used different formats to describe modified peptides
- Hard to share data
- Hard to write software that works with different tools
- Ambiguous representations

**With ProForma:** Standard notation means:
- Data can be easily exchanged between labs
- Software tools can interoperate
- Unambiguous communication of results
- Integration with databases (Unimod, PSI-MOD, etc.)

## Key Design Principles

1. **Human readable** - Scientists can read and understand it
2. **Machine parsable** - Software can reliably parse it
3. **Extensible** - Can add new features as needs evolve
4. **Precise** - Captures uncertainty and ambiguity when present
5. **Standards-based** - Uses controlled vocabularies (Unimod, PSI-MOD, etc.)


