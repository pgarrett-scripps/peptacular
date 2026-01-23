API Reference
=============

.. automodule:: peptacular
   :members:
   :undoc-members:
   :show-inheritance:

Core
----

Peptacular contains a functional and object-oriented API for working with peptides and proteins. Everything 
can be accessed through the peptacular namespace (```import peptacular as pt```), but for clarity the API is broken down into sections below.


Sequence
~~~~~~~~

Processing
^^^^^^^^^^

.. automodule:: peptacular.sequence.digestion
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.transformations
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.mod_builder
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.subseqs
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.basic
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.combinatoric
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.converters
   :members:
   :undoc-members:
   :show-inheritance:

Mass/Comp/Isotope/Fragment
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: peptacular.sequence.isotope
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.mass_funcs
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: peptacular.sequence.fragmentation
   :members:
   :undoc-members:
   :show-inheritance:

Property
^^^^^^^^

.. automodule:: peptacular.sequence.properties
   :members:
   :undoc-members:
   :show-inheritance:

Data Classes
~~~~~~~~~~~~

Some of the core data classes used throughout Peptacular.

.. autoclass:: peptacular.annotation.mod.Interval
.. autoclass:: peptacular.annotation.mod.Mod
.. autoclass:: peptacular.annotation.mod.Mods
.. autoclass:: peptacular.spans.Span
.. autoclass:: peptacular.isotope.IsotopicData
.. autoclass:: peptacular.isotope.IsotopeLookup

Data Modules
------------

These modules are unsed internally to provide data on the various entities Peptacular works with. 
They are provided by the tacular package, but are accessable from the peptacular namespace.

Amino Acids
~~~~~~~~~~~

.. automodule:: peptacular.amino_acids
   :members:
   :undoc-members:
   :show-inheritance:

Elements
~~~~~~~~

.. automodule:: peptacular.elements
   :members:
   :undoc-members:
   :show-inheritance:

Psimods
~~~~~~~

.. automodule:: peptacular.psimod
   :members:
   :undoc-members:
   :show-inheritance:

Unimods
~~~~~~~

.. automodule:: peptacular.unimod
   :members:
   :undoc-members:
   :show-inheritance:

Monosaccharides
~~~~~~~~~~~~~~~

.. automodule:: peptacular.monosaccharides
   :members:
   :undoc-members:
   :show-inheritance:

Ion Types
~~~~~~~~~

.. automodule:: peptacular.ion_types
   :members:
   :undoc-members:
   :show-inheritance:

Neutral Deltas
~~~~~~~~~~~~~~

.. automodule:: peptacular.neutral_deltas
   :members:
   :undoc-members:
   :show-inheritance:

Reference Molecules
~~~~~~~~~~~~~~~~~~~

.. automodule:: peptacular.refmol
   :members:
   :undoc-members:
   :show-inheritance: