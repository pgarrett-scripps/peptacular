Getting Started
===============

**peptacular** is an extremely lightweight package with only one dependency: ``regex``.

It contains functions for parsing and working with Proforma2.0 compliant peptide & protein sequences.

Installation
------------

**From pip**:

.. code-block:: bash

   pip install peptacular

**From source**:

.. code-block:: bash

   git clone https://github.com/pgarrett-scripps/peptacular.git
   cd peptacular
   pip install .


Basic Usage
-----------

All modules and functions in **peptacular** are available under the ``peptacular`` namespace but it is recommended to import as follows:

.. code-block:: python

    import peptacular as pt

**ProForma** sequence parsing in **peptacular** is lazy, meaning only the required notation is
validated during parsing. Invalid modifications are not checked until they are explicitly needed,
such as when calculating mass or composition.

**For example:**

.. code-block:: python

    pt.pop_mods('PEP[INVALID]TIDE')  # Successfully runs
    pt.mass('PEP[INVALID]TIDE')  # Raises an error due to the invalid modification


Key Features
------------

**peptacular** is fully compliant with **ProForma 2.0** and includes functions for:

- **Digestion:** Perform single, multi, or sequential digestion.
- **Fragmentation:** Generate internal, terminal, immonium, and neutral loss fragments.
- **Mass and Composition:** Calculate mass, m/z, or elemental composition of peptides.
- **Modifications:** Apply or remove static and variable modifications.
- **Parsing and Serializing:** Handle ProForma 2.0-compliant sequence parsing and serialization.
- **Isotopic Distributions:** Simulate isotopic patterns.
- **Scoring:** Compare theoretical fragments against experimental spectra.

See the **Examples** section for more detailed use cases.


Sequence Handling
-----------------

All functions in ``pt.sequence`` accept peptide sequences as strings but internally convert them to
**ProformaAnnotation** objects. After processing, they are converted back to strings.

When applying multiple sequence operations on the same peptide, it is more efficient to first convert the
sequence to a **ProformaAnnotation** and use its methods directly.

**Example: Converting between `ProformaAnnotation` and `str`:**

.. code-block:: python

    import peptacular st pt

    annot = pt.parse('PEPTIDE')
    seq = pt.serialize(annot) # or annot.serialize()
    assert 'PEPTIDE' == seq


- This returns either a ``ProformaAnnotation`` or a ``MultiProformaAnnotation`` object.
- ``ProformaAnnotation`` is used for single, linear sequences (the most common use case).
- ``MultiProformaAnnotation`` handles crosslinked or multiple sequences.

**Crosslinked and Multiple Sequences:**

- **Crosslinked**: ``{sequence1}\\{sequence2}``
- **Disconnected**: ``{sequence1}+{sequence2}``

``MultiProformaAnnotation`` contains a list of individual ``ProformaAnnotation`` objects along with their
connection type.


ProForma Notation
----------------------

**ProForma 2.0** was introduced by the **Proteomics Standards Initiative (PSI)** to standardize the representation of peptide sequences, including modifications.

- 📄 **Reference Paper:** `ProForma 2.0 Specification <https://pubs.acs.org/doi/10.1021/acs.jproteome.1c00771>`_
- 📚 **Latest Specification:** `ProForma 2.0 GitHub <https://github.com/HUPO-PSI/ProForma/tree/master/SpecDocument>`_


**Basic Syntax Overview**

- **N-terminal:** ``[+100]-PEPTIDE``
- **C-terminal:** ``PEPTIDE-[+100]``
- **Internal:** ``PEPT[+100]IDE``
- **Global:** ``<[+100]@C>PEPTIDE`` or ``<[+100]@C,P>PEPTIDE``
- **Isotope:** ``<13C>PEPTIDE``  or ``<15N><13C>PEPTIDE``
- **Labile:** ``{+100}PEPTIDE``

Global, isotope, and labile mods are specified before N-terminal modification, or first residue if no terminal mod is present.

**Combined Example:**

.. code-block:: python

    pt.parse('<[+20]@C><13C>{+75}[-40]-PEPT[+50]IDE-[+200]')

    # Returns
    ProFormaAnnotation(
        sequence='PEPTIDE',
        isotope_mods=[Mod('13C', 1)],
        static_mods=[Mod('[+20]@C', 1)],
        labile_mods=[Mod(75, 1)],
        nterm_mods=[Mod(-40, 1)],
        cterm_mods=[Mod(200, 1)],
        internal_mods={3: [Mod(50, 1)]}
    )

**Specifying Proforma Modifications**

The ``Mod`` object contains:
- The **modification string**
- The **number of times** it is applied

You can apply **multiple modifications** at the same position by adding them sequentially:
- ``[+100][+30]`` → Two separate modifications
- ``[+100]^2`` → The same modification applied twice

**Modification Types**

- **Mass-based:** ``[+100]``, ``[100]``, ``[-100]``
- **Chemical formula:** ``[Formula:C12H20O2]``
- **UNIMOD:** ``[Oxidation]``, ``[UNIMOD:21]``, ``[U:21]``
- **PSI-MOD:** ``[L-methionine sulfoxide]``, ``[MOD:00046]``, ``[M:00046]``
- **RESID:** ``[R:L-methionine (R)-sulfoxide]``, ``[RESID:AA0037]``
- **GNO:** ``[GNO:G02815KT]``

While the prefix for unimod and psi-mods are not required (U: and M:), it is still reccommended to use them.
