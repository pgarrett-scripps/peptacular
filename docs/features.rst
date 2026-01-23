Features
========

Peptacular is a Python package for peptide sequence analysis built around ProForma 2.0 notation.

Core Features
-------------

ProForma 2.0 Parsing
~~~~~~~~~~~~~~~~~~~~

Parse and serialize ProForma sequences with modifications, charge states, and isotope labels.

- Full support for Base-ProForma and Level 2-ProForma and select other features
- Named modifications from Unimod and PSI-MOD databases
- Delta mass and formula-based modifications
- Terminal, labile, and ambiguous modifications
- Charge states and adducts
- Global isotope labels and fixed modifications

Mass Calculations
~~~~~~~~~~~~~~~~~

Calculate monoisotopic and average masses with comprehensive modification support.

- Monoisotopic and average mass calculations
- Support for isotope labeling and neutral losses
- Charge state and adduct handling
- Precursor and fragment ion masses
- Elemental composition tracking

Isotopic Distributions
~~~~~~~~~~~~~~~~~~~~~~

Generate theoretical isotopic patterns with configurable resolution and abundance filtering.

- High-performance isotopic pattern generation
- Configurable resolution and abundance thresholds
- Support for fragments with modifications
- Binary exponentiation algorithm for efficiency

Protein Digestion
~~~~~~~~~~~~~~~~~

Simulate enzymatic digestion with various proteases and cleavage rules.

- Built-in protease database (trypsin, chymotrypsin, pepsin, etc.)
- Missed cleavage support
- Semi-specific and non-specific digestion
- Length filtering
- Custom protease definitions

Fragment Generation
~~~~~~~~~~~~~~~~~~~

Create theoretical fragment ions for MS/MS analysis.

- b, y, a, x, c, z ion types
- Internal fragments (ax, by, etc.) & Immonium ions
- Multiple charge states / adducts
- Neutral losses (H2O, NH3, custom)
- Isotopes

Property Calculations
~~~~~~~~~~~~~~~~~~~~~

Calculate various physicochemical properties of peptides.

- Hydrophobicity (multiple scales: Kyte-Doolittle, Eisenberg, etc.)
- Isoelectric point (pI)
- Aromaticity
- Aliphatic index
- Instability index
- GRAVY score
- Secondary structure predictions
- Custom property scales

Format Conversion
~~~~~~~~~~~~~~~~~

Import and export sequences from popular proteomics tools.

**Supported Formats:**

- IP2 (flanking amino acids notation)
- DIANN (modification notation)
- Casanovo (de novo sequencing output)
- MS2PIP (fragment prediction input/output)

Parallel Processing
~~~~~~~~~~~~~~~~~~~

Automatic multiprocessing for batch operations.

- Automatic parallelization for list inputs
- No code changes needed
- Applies to all major functions (mass, digest, fragment, etc.)
- Efficient memory usage

Design Principles
-----------------

Lazy Loading / Caching
~~~~~~~~~~~~~~~~~~~~~~~

Modifications are stored as strings and only parsed when needed for calculations. 
This minimizes overhead and improves performance for operations that don't require 
mass and composition data.

.. testcode::

   import peptacular as pt
   
   # Modifications are not parsed until needed
   peptide = pt.parse("PEM[Oxidation]TIDE")
   
   # First call parses the modifications and caches the return values of 'Oxidation'
   mass1 = peptide.mass()
   
   # Second call uses cached values directly (No re-parsing)
   mass2 = peptide.mass()


**Validation**

By default, parsing does not validate inputs. This allows for users to potenitally create annotations
with invalid modifications/sequences. If you need strict validation:

.. testcode::

   import peptacular as pt
   
   # Strict validation
   a = pt.parse('PEPTIDE', validate=True)
   a.static_mods = 'INVALID_MOD'  # Raises ValueError

   a.validate = False  # Disable validation
   a.static_mods = 'INVALID_MOD'  # Works without error

   # get validation status
   print(a.validate)

