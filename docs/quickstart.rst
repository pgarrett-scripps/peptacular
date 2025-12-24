Quick Start
===========

Installation
------------

Install Peptacular using pip:

.. code-block:: bash

   pip install peptacular

Basic Usage
-----------

Parse/Serialize
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import peptacular as pt

   # Parse a ProForma sequence
   peptide = pt.parse("PEM[Oxidation]TIDE/2")
   print(peptide.charge_state)  # 2
   print(peptide.serialize())  # PEM[Oxidation]TIDE/2


Fragment Generation
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Fragment the peptide
   for frag in peptide.fragment(ion_types=['b', 'y'], charges=[1, 2]):
       print(f"{frag.ion_type}{frag.position}+{frag.charge}: {frag.mz:.3f}")

Protein Digestion
~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Digest a protein
   protein = pt.parse("PEM[Oxidation]TRPEPTIDEKPEPTIDEIDE/2")
   for span in protein.digest(pt.Proteases.TRYPSIN):
       print(f"  {protein[span].serialize()}")

Parse and Modify Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Parse ProForma sequences
   simple = pt.parse("PEPTIDE")
   modified = pt.parse("[Acetyl]-PEM[Oxidation]TIDE/2")
   labeled = pt.parse("<13C>PEPTIDE")

   # Programmatic modification
   annot = pt.ProFormaAnnotation(sequence="PEPTIDE")
   annot.set_nterm_mods({"Acetyl": 1})
   annot.set_internal_mods_at_index(2, {"Oxidation": 1})
   annot.set_charge(2)
   print(annot.serialize())  # [Acetyl]-PEM[Oxidation]TIDE/2

Mass / MZ / Composition
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Calculate masses with different parameters
   precursor_mass = peptide.mass(ion_type="p")
   neutral_mass = peptide.mass(ion_type="n")
   mz_2plus = peptide.mz(charge=2)

   # With isotopes and losses
   mass_c13 = peptide.mass(isotopes=1)
   mass_h2o_loss = peptide.mass(losses={pt.NeutralDelta.WATER: 1})

   # Adduct masses
   mass_na = peptide.mass(charge='Na:z+1')

Isotopic Patterns
~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Generate isotopic distribution
   dist = peptide.isotopic_distribution(
       max_isotopes=5,
       min_abundance_threshold=0.01,
       distribution_resolution=3
   )

   # For fragments with modifications
   dist_fragment = peptide.isotopic_distribution(
       ion_type=pt.IonType.Y,
       charge=2,
       isotopes=1,
       losses={pt.NeutralDelta.WATER: 1}
   )

Digestion
~~~~~~~~~

.. code-block:: python

   protein = pt.parse("PEPTIDEKPEPTIDERPEPTIDER")

   # Tryptic digestion
   for span in protein.digest(pt.Proteases.TRYPSIN):
       peptide = protein[span]
       print(f"{peptide.serialize()} - missed cleavages: {span.missed_cleavages}")

   # With missed cleavages and length filtering
   for span in protein.digest(
       pt.Proteases.TRYPSIN,
       missed_cleavages=2,
       min_len=7,
       max_len=30
   ):
       print(protein[span].serialize())

   # Semi-enzymatic
   for span in protein.digest(pt.Proteases.TRYPSIN, semi=True):
       print(protein[span].serialize())


Physicochemical Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Simple properties
   hydro = peptide.hydrophobicity
   pi = peptide.pi
   arom = peptide.aromaticity

   # Secondary structure prediction
   ss = peptide.secondary_structure()
   print(f"Alpha helix: {ss['alpha_helix']:.1f}%")
   print(f"Beta sheet: {ss['beta_sheet']:.1f}%")

   # Custom property calculations
   custom = peptide.calc_property(
       scale=pt.HydrophobicityScale.KYTE_DOOLITTLE,
       aggregation_method='avg'
   )

Format Conversion
~~~~~~~~~~~~~~~~~

.. code-block:: python

   # From other formats
   ip2_seq = pt.convert_ip2_sequence('K.PEPTIDE.K')
   diann_seq = pt.convert_diann_sequence('_PEM(ox)TIDE_')
   casanovo_seq = pt.convert_casanovo_sequence('+43.006PEPTIDE')

   # To MS2PIP format
   unmod_seq, mod_str = peptide.to_ms2_pip()

Batch Processing
~~~~~~~~~~~~~~~~

.. code-block:: python

   # Process multiple sequences in parallel
   sequences = ['PEPTIDE', 'PROTEIN', 'FRAGMENT', 'SEQUENCE']

   # Automatic parallelization for lists
   masses = pt.mass(sequences)
   compositions = pt.comp(sequences)
   annotations = pt.parse(sequences)

   # Properties
   hydrophobicity = pt.hydrophobicity(sequences)
   pi_values = pt.pi(sequences)
