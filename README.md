

# Peptacular

[![DOI](https://zenodo.org/badge/591504879.svg)](https://doi.org/10.5281/zenodo.15054278)[![Python package](https://github.com/pgarrett-scripps/peptacular/actions/workflows/python-package.yml/badge.svg)](https://github.com/pgarrett-scripps/peptacular/actions/workflows/python-package.yml)[![Pylint](https://github.com/pgarrett-scripps/peptacular/actions/workflows/pylint.yml/badge.svg)](https://github.com/pgarrett-scripps/peptacular/actions/workflows/pylint.yml)


**peptacular** is an extremely lightweight package with only one dependency: ``regex``.

It contains functions for parsing and working with Proforma2.0 compliant peptide & protein sequences.

If you use **peptacular** in your research, please cite: https://doi.org/10.5281/zenodo.15054278

## Documentation
https://peptacular.readthedocs.io/en/latest/index.html


## Installation

```bash
pip install peptacular
```

## Proforma Notation:
- ProForma_v2.pdf
- https://pubs.acs.org/doi/10.1021/acs.jproteome.1c00771


## Warnings

- The internal fragment ion mass calculation may not be accurate. Fairly certain that the ay, by, and cy internal fragments are correct.
- GNO and RESID mods are disabled for now since very few were able to be parsed.
- Project is still under development. I will be adding more features and fixing bugs as I find them.