# Peptacular AI Coding Agent Instructions

## Project Overview
**Peptacular** is a lightweight Python package for parsing and working with **ProForma 2.0 compliant** peptide & protein sequences. It's designed for mass spectrometry proteomics with emphasis on performance and minimal dependencies (only `regex` required for core functionality).

Current branch: `mem-opt` - focuses on memory optimization improvements.

## Architecture & Core Concepts

### ProForma Annotation System
The **central data structure** is `ProFormaAnnotation` (in `src/peptacular/proforma/annotation.py`):
- Uses `__slots__` for memory efficiency (no `__dict__`)
- Immutable modification objects (`Mod`) that are shared/cached across instances
- Factory pattern: methods return new ProFormaAnnotation objects for chaining
- Supports all ProForma 2.0 features: modifications, charge states, intervals, ambiguity, glycans

Example:
```python
ann = ProFormaAnnotation.parse("[Acetyl]-PEPTIDE[Phospho]")
ann.add_nterm_mods("Formyl").mass()  # Methods chain
```

### Key Module Organization
- `proforma/` - ProForma parsing, annotation manipulation, combinatorics
  - `annotation.py` - Main ProFormaAnnotation class with mixins
  - `parser.py` - ProFormaParser for parsing ProForma strings
  - `dclasses/` - Data classes (Mod, ModList, ModDict, Interval, etc.)
- `sequence/` - Functional API for sequence operations (works with strings or annotations)
- `mods/` - Modification databases (Unimod, PSI-MOD, XLMOD, GNO) loaded from JSON
- `fragment/` - Peptide fragmentation (b/y/a/x/c/z ions, internal fragments)
- `digestion/` - Enzymatic and non-enzymatic peptide digestion
- `chem/` - Chemical composition and mass calculations

### Performance Patterns
1. **Caching**: `@lru_cache` extensively used (e.g., `mass_calc.py`, `mod.py`, `dclasses/modlist.py`)
2. **Memory optimization**: 
   - `__slots__` on core classes
   - Shared immutable `Mod` objects (factory pattern in `mod.py`)
   - Lazy evaluation with `@cached_property`
3. **Parallel processing**: Auto-multiprocessing in functional API (`sequence/parrallel.py`)
   - Detects free-threaded Python for thread-based parallelism
   - Falls back to process-based for standard Python
4. **Return type flexibility**: Many functions accept `return_type` parameter (e.g., `'str'`, `'annotation'`, `'span'`)
   - Use `'str'` for large datasets to avoid ProFormaAnnotation overhead

## Development Workflow

### Build & Test Commands (via `uv`)
```bash
make install-dev     # Install all dependencies including dev tools
make test            # Run pytest
make test-cov        # Run tests with HTML coverage report (htmlcov/)
make lint            # Run ruff linter
make format          # Auto-format with ruff (imports + formatting)
make docs            # Build Sphinx docs
make serve-docs      # Serve docs at localhost:8000
make run-gen         # Regenerate modification databases from OBO files
```

**Critical**: Uses `uv` for dependency management (not pip). `pyproject.toml` defines dependency groups.

### Testing Philosophy (see `tests/RULES.md`)
- Tests should be **basic** - one feature per test, relatively short
- Limit asserts (ideally only final sequence via `serialize()`)
- Use simple peptide sequences, avoid complexity
- Compare serialized output: `assert serialize(result) == "EXPECTED[Mod]"`

### Code Quality Standards
- `ruff` for linting/formatting (configured in `pyproject.toml`)
- Pylint max score: 9.0+
- Type annotations required
- Good variable names defined in `pyproject.toml` (e.g., `mz`, `aa`, `db` are acceptable short names)

## Modification Databases

### Database Generation (`generate_mod_dbs.py`)
Parses OBO files (`data/*.obo`) to generate JSON databases:
- **Unimod, PSI-MOD, XLMOD, GNO, Monosaccharides**
- Stores: ID→composition, ID→isotopic mass, ID→average mass, name→ID mappings
- Run `make run-gen` to regenerate from source OBO files

### Database Usage
```python
from peptacular.mods import get_unimod_mass, get_psimod_comp
mass = get_unimod_mass("Phospho")  # Returns monoisotopic mass
comp = get_psimod_comp("MOD:00696")  # Returns chemical composition dict
```

Modification lookups support:
- Database IDs (e.g., `"UNIMOD:21"`)
- Common names (e.g., `"Phospho"`)
- Delta masses (e.g., `"+79.966"`, `"-17.026"`)

## Common Patterns & Conventions

### Dual API: Functional vs. Object-Oriented
```python
# Functional API (works with strings or annotations)
import peptacular as pt
mass = pt.mass("PEPTIDE[Phospho]")
frags = pt.fragment("PEPTIDE", ion_types=['b', 'y'])

# OOP API (via ProFormaAnnotation)
ann = pt.ProFormaAnnotation.parse("PEPTIDE[Phospho]")
mass = ann.mass()
frags = ann.fragment(ion_types=['b', 'y'])
```

### Modification Representation
- `Mod(value, multiplier)` - e.g., `Mod("Phospho", 2)` for 2x phosphorylation
- Internal mods: `{position: [Mod, ...]}` - 0-indexed on unmodified sequence
- Terminal mods: `"nterm"`, `"cterm"` keys in mod dictionaries
- Special mods: `"isotope"`, `"static"`, `"labile"`, `"unknown"`, `"charge_adducts"`

### ProForma Notation Examples
```python
"[Acetyl]-PEPTIDE-[Amide]"              # Terminal mods
"PEPT[+79.966]IDE"                       # Delta mass
"PEPTIDE[UNIMOD:21]"                     # Unimod reference
"<13C>[Acetyl]-PEPTIDE"                  # Isotope label + mod
"PEPTIDE/2"                              # Charge state
"PEP{Hex}TIDE"                           # Glycan
"PEPTIDE[Phospho]^2"                     # Multiplier
```

### Fragment Ion Types
- Forward: `a`, `b`, `c`
- Reverse: `x`, `y`, `z`
- Internal: `ax`, `ay`, `az`, `bx`, `by`, `bz`, `cx`, `cy`, `cz`
- Precursor: `p` (precursor ion)
- Neutral: `n` (neutral loss)

## Critical Files to Understand

- `src/peptacular/__init__.py` - Public API exports (use `import peptacular as pt`)
- `src/peptacular/proforma/annotation.py` - Core ProFormaAnnotation class
- `src/peptacular/proforma/parser.py` - ProForma parsing logic
- `src/peptacular/sequence/mass.py` - Mass calculation with overloads
- `src/peptacular/fragment/mixin.py` - Fragment generation via mixin
- `src/peptacular/constants.py` - Amino acid compositions, ion types, element data
- `examples/` - Working examples for each major feature

## Known Issues & Limitations (from README/CHANGELOG)
- Internal fragment ion mass calculation may not be fully accurate (ay, by, cy likely correct)
- GNO and RESID mods currently disabled (parsing issues)
- Project under active development - API may change

## ProForma Compliance
Fully supports **ProForma 2.0 specification** (see https://pubs.acs.org/doi/10.1021/acs.jproteome.1c00771):
- Ambiguity groups `(A|B|C)`
- Intervals `[1-5:Phospho]`
- Global modifications `<[Acetyl]@K>`
- Chimeric sequences with `+` delimiter
- Connected sequences with `//` delimiter

## Dependencies & Optional Features
- **Core**: Only `regex` required
- **Optional groups** (in `pyproject.toml`):
  - `plotting`: plotly, pandas
  - `benchmark`: lxml, pyteomics
  - `dev`: pytest, ruff, sphinx, pympler

When adding features requiring new deps, use dependency groups appropriately.

## Debugging Tips
- Coverage reports in `htmlcov/` after `make test-cov`
- Profiling scripts in root: `try_*.py` files for testing specific features
- Use `serialize()` to inspect ProFormaAnnotation state
- Check `examples/` for canonical usage patterns
