ProForma Compliance
=======================

Peptacular implements **Base-ProForma** and most of **Level 2-ProForma** features.

Supported Features
------------------

Basic Sequences
~~~~~~~~~~~~~~~

**Unmodified sequences**

.. code-block:: text

   PEPTIDE

Named Modifications
~~~~~~~~~~~~~~~~~~~

**Unimod/PSI-MOD modifications**

.. code-block:: text

   PEM[Oxidation]TIDE
   PEM[U:Oxidation]TIDE

**Accession numbers**

.. code-block:: text

   PEM[UNIMOD:35]TIDE

Mass Modifications
~~~~~~~~~~~~~~~~~~

**Delta masses**

.. code-block:: text

   PEM[+15.995]TIDE
   PEM[U:+15.995]TIDE

**Formula modifications**

.. code-block:: text

   PEM[Formula:O]TIDE

Terminal Modifications
~~~~~~~~~~~~~~~~~~~~~~

**N-terminal and C-terminal**

.. code-block:: text

   [Acetyl]-PEPTIDE
   PEPTIDE-[Amidated]
   [Acetyl]-PEPTIDE-[Amidated]

Labile Modifications
~~~~~~~~~~~~~~~~~~~~

**Modifications lost during fragmentation**

.. code-block:: text

   {Glycan:Hex}PEPTIDE

Multiple Modifications
~~~~~~~~~~~~~~~~~~~~~~

**Multiple modifications on same residue**

.. code-block:: text

   PEM[Oxidation][Dioxidation]TIDE

Charge States
~~~~~~~~~~~~~

**Simple charge states**

.. code-block:: text

   PEPTIDE/2
   PEPTIDE/3

Charge Adducts
~~~~~~~~~~~~~~

**Adduct ions**

.. code-block:: text

   PEPTIDE/[Na:z+1]
   PEPTIDE/[Na:z+1^2,H:z+1]

Global Labels
~~~~~~~~~~~~~

**Isotope labels**

.. code-block:: text

   <13C>PEPTIDE
   <15N>PEPTIDE

**Fixed modifications**

.. code-block:: text

   <[Oxidation]@M>PEPTIDE
   <[Carbamidomethyl]@C>PEPTIDE

Ambiguous Modifications
~~~~~~~~~~~~~~~~~~~~~~~~

**Unknown positions**

.. code-block:: text

   [Oxidation]?PEPTIDE

**Position ranges**

.. code-block:: text

   PRT(ESFRMS)[+19.0523]ISK

Ambiguous Sequences
~~~~~~~~~~~~~~~~~~~~~~~~

**Unknown residues**
.. code-block:: text

   PEXTIDEB

**Unknown Regions**

.. code-block:: text

   PE(?PTI)DE


Glycan Compositions
~~~~~~~~~~~~~~~~~~~

**N-glycosylation**

.. code-block:: text

   N[Glycan:Hex5HexNAc4]K
   N[Glycan:Hex5HexNAc4Neu5Ac2]K

Partially Supported
-------------------

The following features are partially implemented, and will not
result in a parsing error unless validation is enabled. Although, any 
method or function which requires calculating the mass/compisition of these 
features will raise an error.


Position Sets with Groups
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   PEP[Oxidation#1]M[#1]AT

Localization Scores
~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   PEP[Oxidation#1(0.95)]M[#1(0.05)]AT

Alternative Vocabularies
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   # RESID / GNO / XLMOD modifications not yet supported
   PEM[RESID:AA0581]TIDE
   N[GNO:G00001]K
   PEPTK[XLMOD:02001]IDE

Unsupported
-----------

The following features are Unsupported and will raise parsing errors.

Chimeric Spectra
~~~~~~~~~~~~~~~~  

.. code-block:: text

   # Multiple co-fragmented peptides
   PEPTIDE+SEQUENCE

Cross-linking
~~~~~~~~~~~~~

.. code-block:: text

   # XL-MOD notation, branch points
   PEPTK[XLMOD:02001#XL1]IDE//SEQK[#XL1]

For More Information
--------------------

See `PROFORMA_COMPLIANCE.md <https://github.com/pgarrett-scripps/peptacular/blob/main/PROFORMA_COMPLIANCE.md>`_ 
for detailed compliance status and updates.
