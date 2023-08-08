# Changelog

All notable changes to this project will be documented in this file.

## [0.3.0]

### Added:
- score.py - a module for scoring peptide matches

### Changes:
- made Fragment dataclass frozen

### Removed:
- removed indexing code from mass.py


## [0.2.0]

### Added:
- Fragmenter.py - an easier to use peptide fragmenter (returns dataclasses)
- Internal fragment ion support 

## Changes:
- fixed bug with C-Term PTMs

## [0.1.0]

### Added:
- spans.py to handle generating peptide spans from cleavage site info
- masses.py to handle generating peptide masses from sequences

### Removed:
- refseq module - was a dumb idea

## [0.0.6]

### Added:
- fragment ion
- peptide mass and mz functions

### Changes
- converted peptide.py to sequence.py (makes more sense since proteins are also sequences of amino acids)

## [0.0.5]

### Added:
- license
- changelog
- documentation

### Changes
- src-based layout
- pyproject.toml

### Removed
- removed fasta module -> independent FastaFrames package
- calculate_peptide_mass(
- moved streamlit components to root

## [0.0.3]

### Added
- added semi, and non-enzymatic digestion