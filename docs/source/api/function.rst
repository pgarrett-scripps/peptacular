Core Functions
==============

peptacular provides convenient top-level functions that work with both strings and ProFormaAnnotation objects.

Parsing and Serialization
--------------------------

.. autofunction:: peptacular.parse
.. autofunction:: peptacular.serialize

Mass and Composition
---------------------

.. autofunction:: peptacular.mass
.. autofunction:: peptacular.mz  
.. autofunction:: peptacular.comp

Fragmentation
-------------

.. autofunction:: peptacular.fragment
.. autofunction:: peptacular.internal_fragments
.. autofunction:: peptacular.terminal_fragments

Digestion
---------

.. autofunction:: peptacular.digest
.. autofunction:: peptacular.cleave_after
.. autofunction:: peptacular.cleave_before

Modification Functions
----------------------

.. autofunction:: peptacular.add_mods
.. autofunction:: peptacular.pop_mods
.. autofunction:: peptacular.strip_mods
.. autofunction:: peptacular.condense_static_mods

Sequence Manipulation
---------------------

.. autofunction:: peptacular.slice_sequence
.. autofunction:: peptacular.reverse_sequence
.. autofunction:: peptacular.shuffle_sequence

Isotopic Distributions
----------------------

.. autofunction:: peptacular.isotope_pattern
.. autofunction:: peptacular.centroid_pattern

Scoring Functions
-----------------

.. autofunction:: peptacular.score_spectrum
.. autofunction:: peptacular.match_peaks

Examples
--------

Basic Usage
^^^^^^^^^^^

.. code-block:: python

    import peptacular as pt

    # Parse and analyze
    mass = pt.mass('PEPTIDE')
    mz = pt.mz('PEPTIDE/2')
    
    # Fragment
    fragments = pt.fragment('PEPTIDE', ion_types=['b', 'y'])
    
    # Digest
    peptides = pt.digest('PROTEIN', enzyme='trypsin')

Advanced Examples
^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Work with modifications
    modified = pt.add_mods('PEPTIDE', {'nterm': ['+42']})
    stripped = pt.strip_mods('PEPT[+16]IDE')
    
    # Isotope patterns
    pattern = pt.isotope_pattern('PEPTIDE', charge=1)
    
    # Scoring
    score = pt.score_spectrum(fragments, experimental_mz, experimental_intensity)
