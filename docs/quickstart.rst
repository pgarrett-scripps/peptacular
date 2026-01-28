Quick Start
===========

Installation
------------

Install Peptacular using pip:

.. code-block:: bash

   pip install peptacular

Basic Usage
-----------

.. testsetup:: *

   import peptacular as pt

Parse/Serialize
~~~~~~~~~~~~~~~~~~~

.. testcode::

   import peptacular as pt
   
   # Parse a ProForma sequence
   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   print(peptide.charge_state)
   print(peptide.serialize())

.. testoutput::

   2
   PEM[Oxidation]TIDE/2


Fragment Generation
~~~~~~~~~~~~~~~~~~~

.. testcode::

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   
   # Fragment the peptide
   for frag in peptide.fragment(ion_types=['b', 'y'], charges=[1, 2])[:5]:
       print(f"{frag.ion_type}{frag.position}+{frag.charge_state}: {frag.mz:.3f}")

.. testoutput::
   
   b1+1: 98.060
   b1+2: 49.534
   b2+1: 227.103
   b2+2: 114.055
   b3+1: 374.138

Protein Digestion
~~~~~~~~~~~~~~~~~

.. testcode::

   # Digest a protein
   protein = pt.parse("PEM[Oxidation]TRPEPTIDEKPEPTIDEIDE/2")
   for span in protein.digest(pt.Proteases.TRYPSIN):
       print(f"  {protein[span].serialize()}")

.. testoutput::
   :hide:

   ...

Parse and Modify Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. testcode::

   # Parse ProForma sequences
   simple = pt.parse("PEPTIDE")
   modified = pt.parse("[Acetyl]-PEP[Oxidation]TIDE/2")
   labeled = pt.parse("<13C>PEPTIDE")

   # Programmatic modification
   annot = pt.ProFormaAnnotation(sequence="PEPTIDE")
   annot.set_nterm_mods({"Acetyl": 1})
   annot.set_internal_mods_at_index(2, {"Oxidation": 1})
   annot.set_charge(2)
   print(annot.serialize())

.. testoutput::

   [Acetyl]-PEP[Oxidation]TIDE/2

Mass / MZ / Composition
~~~~~~~~~~~~~~~~~~~~~~~

.. testcode::

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   
   # Calculate masses with different parameters
   precursor_mass = peptide.mass(ion_type="p")
   neutral_mass = peptide.mass(ion_type="n")
   mz_2plus = peptide.mz(charge=2)

   # With isotopes and losses
   mass_c13 = peptide.mass(isotopes=1)
   mass_h2o_loss = peptide.mass(deltas={'H-2O-1': 1})

   # Adduct masses
   mass_na = peptide.mass(charge='Na:z+1')

Isotopic Patterns
~~~~~~~~~~~~~~~~~

.. testcode::

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   
   # Generate isotopic distribution
   dist = peptide.isotopic_distribution(
       max_isotopes=5,
       min_abundance_threshold=0.01,
       distribution_resolution=3
   )

   # For fragments with modifications
   dist_fragment = peptide.isotopic_distribution(
       ion_type='y',
       charge=2,
       isotopes=1,
       deltas={'H-2O-1': 1}
   )

Digestion
~~~~~~~~~

.. testcode::

   protein = pt.parse("PEPTIDEKPEPTIDERPEPTIDER")

   # Tryptic digestion
   for span in protein.digest('trypsin'):
       peptide = protein[span]
       print(f"{peptide.serialize()} - missed cleavages: {span.missed_cleavages}")

.. testoutput::
   :hide:

   ...

.. testcode::

   protein = pt.parse("PEPTIDEKPEPTIDERPEPTIDER")
   
   # With missed cleavages and length filtering
   for span in protein.digest(
       'trypsin',
       missed_cleavages=2,
       min_len=7,
       max_len=30
   ):
       print(protein[span].serialize())

.. testoutput::
   :hide:

   ...

.. testcode::

   protein = pt.parse("PEPTIDEKPEPTIDERPEPTIDER")
   
   # Semi-enzymatic
   for span in protein.digest('trypsin', semi=True):
       print(protein[span].serialize())

.. testoutput::
   :hide:

   ...


Physicochemical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. testcode::

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   
   # Simple properties
   hydro = peptide.prop.hydrophobicity
   pi = peptide.prop.pi
   arom = peptide.prop.aromaticity

   # Secondary structure prediction
   ss = peptide.prop.secondary_structure()
   print(f"Alpha helix: {ss['alpha_helix']:.1f}%")
   print(f"Beta sheet: {ss['beta_sheet']:.1f}%")

.. testoutput::
   :hide:

   ...

.. testcode::

   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   
   # Custom property calculations
   custom = peptide.prop.calc_property(
       scale=pt.HydrophobicityScale.KYTE_DOOLITTLE,
       aggregation_method='avg'
   )

Format Conversion
~~~~~~~~~~~~~~~~~

.. testcode::

   # From other formats
   ip2_seq = pt.convert_ip2_sequence('K.PEPTIDE.K')
   diann_seq = pt.convert_diann_sequence('_PEM[Carbamidomethyl]TIDE_')
   casanovo_seq = pt.convert_casanovo_sequence('+43.006PEPTIDE')

   # To MS2PIP format
   peptide = pt.parse("PEM[Oxidation]TIDE")
   unmod_seq, mod_str = peptide.to_ms2_pip()

Batch Processing
~~~~~~~~~~~~~~~~

.. testcode::

   # Process multiple sequences in parallel
   sequences = ['PEPTIDE', 'PEPTIDES', 'PEPTIDEST']

   # Automatic parallelization for lists
   masses = pt.mass(sequences)
   compositions = pt.comp(sequences)
   annotations = pt.parse(sequences)

   # Properties
   hydrophobicity = pt.hydrophobicity(sequences)
   pi_values = pt.pi(sequences)