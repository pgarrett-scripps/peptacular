
# Peptacular

**peptacular** is an extremely lightweight package with only one dependency: ``regex``.

It contains functions for parsing and working with Proforma2.0 compliant peptide & protein sequences.

## Documentation
https://peptacular.readthedocs.io/en/latest/


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