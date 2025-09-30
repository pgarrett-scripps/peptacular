ProFormaAnnotation Class
========================

The ``ProFormaAnnotation`` class is the core component of peptacular, providing comprehensive support for ProForma 2.0 compliant peptide sequence representation and manipulation.

Core Properties
---------------

Parsing and Serialization
--------------------------
.. automethod:: peptacular.ProFormaAnnotation.parse
.. automethod:: peptacular.ProFormaAnnotation.serialize


Mass and Composition Analysis
-----------------------------
.. automethod:: peptacular.ProFormaAnnotation.mass
.. automethod:: peptacular.ProFormaAnnotation.mz
.. automethod:: peptacular.ProFormaAnnotation.comp

Modification Management
-----------------------

The ProFormaAnnotation class has a ton of methods for managing modifications on the sequence. Valid modifications are of types:
int, float, str, or Mod. Internally modifications are kept within a list, and cast to type Mod when added. Internal mods
are stored in a dictionary mapping their position to a modlist. Intervals are stored within their own Interval object, which has
a start, end, ambiguity flag, and modification list.

Appending Modifications
^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.append_isotope_mod
.. automethod:: peptacular.ProFormaAnnotation.append_static_mod
.. automethod:: peptacular.ProFormaAnnotation.append_labile_mod
.. automethod:: peptacular.ProFormaAnnotation.append_unknown_mod
.. automethod:: peptacular.ProFormaAnnotation.append_nterm_mod
.. automethod:: peptacular.ProFormaAnnotation.append_cterm_mod
.. automethod:: peptacular.ProFormaAnnotation.append_charge_adduct
.. automethod:: peptacular.ProFormaAnnotation.append_internal_mod_at_index
.. automethod:: peptacular.ProFormaAnnotation.append_interval

Extending Modifications
^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.extend_isotope_mods
.. automethod:: peptacular.ProFormaAnnotation.extend_static_mods
.. automethod:: peptacular.ProFormaAnnotation.extend_labile_mods
.. automethod:: peptacular.ProFormaAnnotation.extend_unknown_mods
.. automethod:: peptacular.ProFormaAnnotation.extend_nterm_mods
.. automethod:: peptacular.ProFormaAnnotation.extend_cterm_mods
.. automethod:: peptacular.ProFormaAnnotation.extend_charge_adducts
.. automethod:: peptacular.ProFormaAnnotation.extend_internal_mods
.. automethod:: peptacular.ProFormaAnnotation.extend_internal_mods_at_index
.. automethod:: peptacular.ProFormaAnnotation.extend_intervals

Setting Modifications
^^^^^^^^^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.set_isotope_mods
.. automethod:: peptacular.ProFormaAnnotation.set_static_mods
.. automethod:: peptacular.ProFormaAnnotation.set_labile_mods
.. automethod:: peptacular.ProFormaAnnotation.set_unknown_mods
.. automethod:: peptacular.ProFormaAnnotation.set_nterm_mods
.. automethod:: peptacular.ProFormaAnnotation.set_cterm_mods
.. automethod:: peptacular.ProFormaAnnotation.set_charge_adducts
.. automethod:: peptacular.ProFormaAnnotation.set_internal_mods
.. automethod:: peptacular.ProFormaAnnotation.set_internal_mods_at_index
.. automethod:: peptacular.ProFormaAnnotation.set_intervals
.. automethod:: peptacular.ProFormaAnnotation.set_charge

Popping Modifications
^^^^^^^^^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.pop_isotope_mods
.. automethod:: peptacular.ProFormaAnnotation.pop_static_mods
.. automethod:: peptacular.ProFormaAnnotation.pop_labile_mods
.. automethod:: peptacular.ProFormaAnnotation.pop_unknown_mods
.. automethod:: peptacular.ProFormaAnnotation.pop_nterm_mods
.. automethod:: peptacular.ProFormaAnnotation.pop_cterm_mods
.. automethod:: peptacular.ProFormaAnnotation.pop_charge_adducts
.. automethod:: peptacular.ProFormaAnnotation.pop_internal_mods
.. automethod:: peptacular.ProFormaAnnotation.pop_intervals
.. automethod:: peptacular.ProFormaAnnotation.pop_charge

Removing Modifications
^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.remove_isotope_mods
.. automethod:: peptacular.ProFormaAnnotation.remove_static_mods
.. automethod:: peptacular.ProFormaAnnotation.remove_labile_mods
.. automethod:: peptacular.ProFormaAnnotation.remove_unknown_mods
.. automethod:: peptacular.ProFormaAnnotation.remove_nterm_mods
.. automethod:: peptacular.ProFormaAnnotation.remove_cterm_mods
.. automethod:: peptacular.ProFormaAnnotation.remove_charge_adducts
.. automethod:: peptacular.ProFormaAnnotation.remove_internal_mods
.. automethod:: peptacular.ProFormaAnnotation.remove_intervals
.. automethod:: peptacular.ProFormaAnnotation.remove_charge

Getting Modifications
^^^^^^^^^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.get_isotope_mods
.. automethod:: peptacular.ProFormaAnnotation.get_static_mods
.. automethod:: peptacular.ProFormaAnnotation.get_labile_mods
.. automethod:: peptacular.ProFormaAnnotation.get_unknown_mods
.. automethod:: peptacular.ProFormaAnnotation.get_nterm_mods
.. automethod:: peptacular.ProFormaAnnotation.get_cterm_mods
.. automethod:: peptacular.ProFormaAnnotation.get_charge_adducts
.. automethod:: peptacular.ProFormaAnnotation.get_internal_mods
.. automethod:: peptacular.ProFormaAnnotation.get_internal_mods_at_index
.. automethod:: peptacular.ProFormaAnnotation.get_intervals
.. automethod:: peptacular.ProFormaAnnotation.get_charge

Bulk Modifications
^^^^^^^^^^^^^^^^^^

Bulk operations apply to the whole proforma annotation. For instance the append_mods() method takes a dictionary of modifications to apply.
Where the keys are 'nterm', 'cterm', 'isotope', 'static', 'labile', 'unknown', 'interval', 'internal', 'charge', 'charge_adducts'

.. automethod:: peptacular.ProFormaAnnotation.append_mods
.. automethod:: peptacular.ProFormaAnnotation.extend_mods
.. automethod:: peptacular.ProFormaAnnotation.set_mods
.. automethod:: peptacular.ProFormaAnnotation.pop_mods
.. automethod:: peptacular.ProFormaAnnotation.remove_mods
.. automethod:: peptacular.ProFormaAnnotation.get_mods

Sequence Manipulation
---------------------
.. automethod:: peptacular.ProFormaAnnotation.slice
.. automethod:: peptacular.ProFormaAnnotation.split
.. automethod:: peptacular.ProFormaAnnotation.join
.. automethod:: peptacular.ProFormaAnnotation.reverse
.. automethod:: peptacular.ProFormaAnnotation.shift
.. automethod:: peptacular.ProFormaAnnotation.shuffle
.. automethod:: peptacular.ProFormaAnnotation.sort

Utility Methods
---------------

Query Methods
^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.contains_mod
.. automethod:: peptacular.ProFormaAnnotation.get_mods
.. automethod:: peptacular.ProFormaAnnotation.count_residues
.. automethod:: peptacular.ProFormaAnnotation.percent_residues

Comparison Methods
^^^^^^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.is_subsequence
.. automethod:: peptacular.ProFormaAnnotation.find_indices
.. automethod:: peptacular.ProFormaAnnotation.coverage
.. automethod:: peptacular.ProFormaAnnotation.percent_coverage


Advanced Features
-----------------

Ambiguity Handling
^^^^^^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.annotate_ambiguity

Combinatorics
^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.permutations
.. automethod:: peptacular.ProFormaAnnotation.combinations
.. automethod:: peptacular.ProFormaAnnotation.combinations_with_replacement
.. automethod:: peptacular.ProFormaAnnotation.product

Modification Building
^^^^^^^^^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.build_mods
.. automethod:: peptacular.ProFormaAnnotation.condense_static_mods
.. automethod:: peptacular.ProFormaAnnotation.condense_mods_to_intervals
.. automethod:: peptacular.ProFormaAnnotation.condense_to_delta_mass

Sliding Windows
^^^^^^^^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.sliding_windows

Digestion
^^^^^^^^
.. automethod:: peptacular.ProFormaAnnotation.get_left_semi_enzymatic_sequences
.. automethod:: peptacular.ProFormaAnnotation.get_right_semi_enzymatic_sequences
.. automethod:: peptacular.ProFormaAnnotation.get_semi_enzymatic_sequences
.. automethod:: peptacular.ProFormaAnnotation.get_non_enzymatic_sequences
.. automethod:: peptacular.ProFormaAnnotation.regex_digest
.. automethod:: peptacular.ProFormaAnnotation.digest

