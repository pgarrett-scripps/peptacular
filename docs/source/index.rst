.. peptacular documentation master file, created by
   sphinx-quickstart on Wed Mar 19 12:53:11 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to peptacular's documentation!
======================================

Peptacular is a Python library for peptide-related bioinformatics. It provides tools for:

* Peptide mass calculations
* Enzymatic digestion
* Spectral matching and scoring
* Glycan and modification handling
* Sequence manipulation and more

The library is designed to be easy to use, with all functions available directly from the main package:

.. code-block:: python

   import peptacular as pt
   
   # Calculate the mass of a peptide
   mass = pt.mass("PEPTIDE")

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules/getting_started
   modules/examples

Module Reference
=================

.. toctree::
   :maxdepth: 4
   :caption: Module Documentation:

   api/peptacular

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
