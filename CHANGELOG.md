# Changelog

All notable changes to this project will be documented in this file.

## Next release?
- make fragment faster by using numpys cumsum function and caching single modified residue masses

## [1.2.0] 
### Added:
- Added support for custom aa masses to mass.py and fragment.py

## [1.1.1] 
### Added:
- Added isotopes and loss to fragment.py

## [1.1.0] 
### Added:
- speed_test.py to examples and more examples!
- more term functions! pop_c_term & pop_n_term & _get_c_term_index & _get_n_term_index
- split_sequence function to sequence.py, which splits a sequence into a list of modified single residues
- pop_modifications function to sequence.py, since it's useful!

### Changed
- fragment.fragment now returns a list rather than a generator
- identify_cleavage_sites no longer returns cleavage sites at the beginning or end of the sequence (x)PEPTIDE(x)
- removed the x's appended onto sequence when calculating sites
- fixed bug with non-specific digestion
- made term functions more robust, using _get_c_term_index & _get_n_term_index when possible
- split term.py into term.modification and term.residue (it just got too large)

## [1.0.1] 

### Added:
- term.py - a module for handling terminal modifications
- readthedocs 'sphynx-style' documentation
- examples folder

### Changed:
- split up sequence module into sequence and digest
- split mass module into fragment and mass
- updated term modification notation to use square brackets
- renamed most functions... again
- updated documentation to be compatible with sphynx

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