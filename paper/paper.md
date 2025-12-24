```yaml
title: 'Peptacular: A Python Library for ProForma 2.0-Compliant Peptide Analysis'
tags:
  - Python
  - proteomics
  - mass spectrometry
  - ProForma
  - bioinformatics 
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
```

---

# Summary

Mass spectrometry-based proteomics relies on computational methods to identify and characterize amino acid sequences [@angel-2012]. These sequences, which range from short peptides to complete proteins, exhibit considerable chemical complexity. This complexity includes post-translational modifications, variable charge states, charge adducts, neutral losses, and isotopic patterns [@smith-2013]. ProForma notation provides a standardized, human-readable representation of this complexity, facilitating consistent communication of peptide and protein identifications across the proteomics community [@leduc-2018].

**Peptacular** is a pure Python library designed for working with **ProForma** notation in high-throughput proteomics applications. The library provides functionality for sequence manipulation, molecular mass and composition calculation, physicochemical property prediction, fragment ion generation, in silico enzymatic digestion, isotopic distribution computation, and format conversion between ProForma and vendor-specific notations.

# Statement of Need

The proteomics field has historically lacked standardization in peptide and protein sequence representation. Individual software tools have implemented proprietary notations, creating barriers to data integration and reanalysis. While the Proteomics Standards Initiative developed ProForma notation to address this issue, adoption remains limited due to its relative novelty and the scarcity of computational tools supporting the standard.

Modern proteomics experiments routinely identify hundreds of thousands of peptide-spectrum matches. Results are typically exported as tabular files (CSV, TSV, Parquet) compatible with Python data analysis libraries such as **pandas** [@&#x74;eam-2025] and **polars** [@vink-2025]. However, existing proteomics libraries lack efficient mechanisms for performing sequence-based calculations at this scale or integrating seamlessly with dataframe-centric workflows. 

Peptacular addresses these limitations by providing two complementary APIs: an object-oriented interface for complex sequence manipulation, and a functional interface designed for batch processing. The functional API accepts serialized ProForma sequences, applies transformations using parallel processing, and returns serialized results. This design enables direct integration with pandas and polars dataframes without sacrificing computational performance.

# Features and Implementation

Peptacular's architecture prioritizes computational efficiency through several design decisions. The library employs lazy evaluation to minimize memory overhead: modifications remain in serialized form until explicitly required for calculations. Since Python strings are immutable and proteomics datasets typically contain repeated modifications, Peptacular implements caching of modification parsing and serialization, yielding performance and memory improvements.

The object-oriented API uses a factory pattern to construct Annotation objects, enabling sequence manipulation including modification placement, isotopic labeling, and custom fragment ion definitions. The functional API provides stateless functions that operate directly on serialized sequences, facilitating parallel execution across CPU cores with near-linear scaling.

Modification reference data, including masses, compositions, and identifiers from Unimod [@creasy-2004], and other databases, are embedded within the package as Python modules rather than external files. This design eliminates file I/O overhead during multiprocessing initialization, improving parallel performance. Additionally, Peptacular uses conditional initialization: modification-related data structures are only instantiated when the corresponding modification type is present in the input sequences, reducing startup overhead and memory consumption.

The library supports bidirectional conversion between ProForma notation and common vendor-specific formats, facilitating integration with existing proteomics workflows while promoting migration toward standardized representations.

# Related Software

Several existing tools address aspects of peptide manipulation and property calculation, including pyteomics [@goloborodko-2013], rdkit [@goloborodko-2013], and biopython [@cock-2009]. However, these tools were developed prior to the establishment of ProForma 2.0 notation and consequently lack native support for many newer ProForma-specific features. Furthermore, these packages include extensive dependencies and functionality that extends well beyond the scope of peptide notation standardization. Peptacular addresses this gap by providing a focused implementation designed to serve as computational infrastructure connecting different software tools and packages within the proteomics ecosystem.