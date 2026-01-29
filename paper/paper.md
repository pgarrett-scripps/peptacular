---
title: 'Peptacular: A Python Library for ProForma 2.0-Compliant Peptide Analysis'
tags:
  - Python
  - Proteomics
  - Mass Spectrometry
  - ProForma
  - Bioinformatics 
authors:
  - name: Patrick Tyler Garrett 
    orcid: 0000-0002-8434-9693 
    affiliation: 1
  - name: John R. Yates III
    orcid: 0000-0001-5267-1672 
    corresponding: true
    affiliation: 1
affiliations:
  - name: The Scripps Research Institute, United States
    index: 1
date: 12 December 2025
bibliography: paper.bib
---

# Summary

Mass spectrometry-based proteomics relies on computational methods to identify and characterize amino acid sequences [@angel-2012]. These sequences, which range from short peptides to complete proteins, exhibiting considerable chemical complexity. This complexity includes post-translational modifications, variable charge states, charge adducts, neutral losses, and isotopic patterns [@smith-2013]. ProForma notation provides a standardized, human-readable representation of this complexity, facilitating consistent communication of peptide and protein identifications across the proteomics community [@leduc-2018].

**Peptacular** is a pure Python library designed for working with **ProForma** notation in high-throughput proteomics applications. The library provides functionality for sequence manipulation, molecular mass and composition calculation, physicochemical property prediction, fragment ion generation, in silico enzymatic digestion, isotopic distribution computation, and format conversion between ProForma and vendor-specific notations.

# Statement of Need

The proteomics field has historically lacked standardization for peptide and protein sequences. Individual software tools have implemented proprietary notations, creating barriers to data integration and reanalysis. While the Proteomics Standards Initiative developed ProForma notation to address this issue, adoption remains limited due to its relative novelty and the scarcity of computational tools supporting the standard.

Modern proteomics experiments routinely identify hundreds of thousands of peptide-spectrum matches. Results are typically exported as tabular files (CSV, TSV, Parquet) compatible with Python data analysis libraries such as **pandas** [@team-2025] and **polars** [@vink-2025]. However, existing proteomics libraries lack efficient mechanisms for performing sequence-based calculations at this scale or integrating seamlessly with dataframe-centric workflows.

Peptacular addresses these limitations by providing two complementary APIs: an object-oriented interface, and a functional interface. The functional API accepts serialized ProForma sequences, applies transformations using parallel processing, and returns serialized results. This design enables direct integration with dataframes without sacrificing computational performance. The object-oriented interface allows for direct manipulation of the annotation objetcs without the need for serialiazation or parsing, improving performance for more complex operations.

# Features and Implementation

Peptacular's employs lazy evaluation to minimize memory overhead: modifications remain in serialized form until explicitly required for calculations. Since Python strings are immutable and proteomics datasets typically contain repeated modifications, Peptacular utilizes caching wherever possible, yielding performance and memory improvements.

The object-oriented API uses a factory pattern to construct Annotation objects, enabling sequence manipulation including modification placement, isotopic labeling, and custom fragment ion definitions. The functional API provides stateless functions that operate directly on serialized sequences, facilitating parallel execution across CPU cores.

Modification reference data, including masses, compositions, and identifiers from Unimod [@creasy-2004], and other databases, are embedded within the package as Python modules rather than external files. This design eliminates file I/O overhead for parsing the supported ontologies, improving performance. Additionally, Peptacular uses conditional initialization: modification-related data structures are only instantiated when the corresponding modification type is present in the input sequences, reducing startup overhead and memory consumption.

The library supports bidirectional conversion between ProForma notation and common vendor-specific formats, facilitating integration with existing proteomics workflows while promoting migration toward standardized representations.

# Related Software

Several existing tools address aspects of peptide manipulation and property calculation, including pyteomics [@goloborodko-2013] and biopython [@cock-2009]. However, these tools were developed prior to the establishment of ProForma 2.0 notation and consequently lack native support for many modern ProForma-specific features. Furthermore, these packages include extensive dependencies and functionality that extends well beyond the scope of peptide notation standardization. Peptacular addresses this gap by providing a focused implementation designed to serve as computational infrastructure connecting different software tools and packages within the proteomics ecosystem. 

# Proforma Compliance

## 5.3.1 Base-ProForma Compliance

- [x] **Amino acids (+UO)** - `AAHCFKUOT` (§6.1)
- [x] **Unimod names** - `PEM[Oxidation]AT` (§6.2.1)
- [x] **PSI-MOD names** - `PEM[monohydroxylated residue]AT` (§6.2.1)
- [x] **Unimod numbers** - `PEM[UNIMOD:35]AT` (§6.2.2)
- [x] **PSI-MOD numbers** - `PEM[MOD:00425]AT` (§6.2.2)
- [x] **Delta masses** - `PEM[+15.995]AT` (§6.2.3)
- [x] **N-terminal modifications** - `[Carbamyl]-QPEPTIDE` (§6.3)
- [x] **C-terminal modifications** - `PEPTIDEG-[Methyl]` (§6.3)
- [x] **Labile modifications** - `{Glycan:Hex}EM[U:Oxidation]EV` (§6.4)
- [x] **Multiple modifications** - `MPGNW[Oxidation][Carboxymethyl]PESQE` (§6.5)
- [x] **Information tag** - `ELV[INFO:AnyString]IS` (§6.6)

## 5.3.2 Level 2-ProForma Compliance

- [x] **Ambiguous amino acids** - `BZJX` (§7.1)
- [x] **Prefixed delta masses** - `PEM[U:+15.995]AT` (§7.2)
- [x] **Mass gap** - `PEX[+147.035]AT` (§7.3)
- [x] **Formulas** - `PEM[Formula:O]AT`, `PEM[Formula:[17O1]]AT` (§7.4)
- [x] **Mass with interpretation** - `PEM[+15.995|Oxidation]AT` (§7.5)
- [x] **Unknown mod position** - `[Oxidation]?PEMAT` (§7.6.1)
- [x] **Set of positions** - `PEP[Oxidation#1]M[#1]AT` (§7.6.2)
- [x] **Range of positions** - `PRT(ESFRMS)[+19.0523]ISK` (§7.6.3)
- [x] **Position scores** - `PEP[Oxidation#1(0.95)]M[#1(0.05)]AT` (§7.6.4)
- [x] **Range position scores** - `(PEP)[Oxidation#1(0.95)]M[#1(0.05)]AT` (§7.6.5)
- [x] **Amino acid ambiguity** - `(?VCH)AT` (§7.7)
- [x] **Modification prefixes** - `PEPM[U:Oxidation]AS[M:O-phospho-L-serine]` (§7.8)

## 5.3.3 Level 2-ProForma + Top-Down

- [x] **RESID modifications** - `EM[R:L-methionine sulfone]EM[RESID:AA0581]` (§8.1)
- [x] **Names** - `(>Heavy chain)EVQLVESG` (§8.2)

## 5.3.4 Level 2-ProForma + Cross-Linking

- [x] **XL-MOD modifications** - `EVTK[X:Aryl azide]LEK[XLMOD:00114]SEFD` (§9.1)
- [ ] **Cross-linkers (intrachain)** - `EVTK[X:Aryl azide#XL1]LEK[#XL1]SEFD` (§9.2.1)
- [ ] **Cross-linkers (interchain)** - `EVTK[X:Aryl azide#XL1]L//EK[#XL1]SEFD` (§9.2.2)
- [ ] **Branches** - `ED[MOD:00093#BRANCH]//D[#BRANCH]ATR` (§9.3) 

## 5.3.5 Level 2-ProForma + Glycans

- [x] **GNO modifications** - `NEEYN[GNO:G59626AS]K` (§10.1)
- [x] **Glycan compositions** - `NEEYN[Glycan:Hex5HexNAc4NeuAc1]K` (§10.2)

## 5.3.6 Level 2-ProForma + Advanced Complexity

- [x] **Charged formulas** - `SEQUEN[Formula:Zn1:z+2]CE` (§11.1) 
- [x] **Controlling placement** - `PTI(MERMERME)[+32|Position:E]PTIDE` (§11.2) 
- [x] **Global isotope** - `<13C>CARBON` (§11.3.1)
- [x] **Fixed modifications** - `<[Oxidation]@M>ATPEMILTCMGCLK` (§11.3.2)
- [x] **Chimeric spectra** - `NEEYN+SEQUEN` (§11.4) 
- [x] **Charges** - `SEQUEN/2`, `SEQUEN/[Na:z+1,H:z+1]` (§11.5)
- [x] **Ion notation** - `SEQUEN-[b-type-ion]` (§11.6)

# AI usage disclosure

Generative AI models were used to assist in developing this software's codebase, tests, and documentation, primarily Claude Sonnet 4.5, Gemini 2.0 Pro, and GitHub Copilot's autocomplete extension. These tools were accessed through the Copilot extension in VS Code and were used for code generation, debugging assistance, and documentation writing. No AI models were used in the creation of this manuscript. 
