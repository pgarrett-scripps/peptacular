"""
This booster module will use xgboost to predict certain properties of a peptide
such as retention time, collisional cross section, and MS/MS spectra.

Input will be a feature vector containing the aa counts as well as the start and end modification.

Feature vector:

1) 22 dim (counts of each amino acids)
2) ?Dim encoding of first amino acid
3) ?Dim encoding of last amino acid
4) 22 dim (sum of modification mass for each amino acid type)
5) 1 dim (sum of nterm modification mass)
6) 1 dim (sum of cterm modification mass)
7) 2 dim (peptide length, peptide mass)
8) 1 dim (charge state) [optional]
9) 5 dim (sliding hydrophobicity window)
10) 5 dim (sliding other window)
"""

from typing import Sequence

from ...sequence.basic import serialize

from ...proforma.annotation import ProFormaAnnotation
from .adapters import MS2PIPAdapter, SageAdapter, PeptideData
from .aligner import LinearAligner
from .predictor import RtPredictor
from .utils import train_test_split


__all__ = [
    "MS2PIPAdapter",
    "SageAdapter",
    "PeptideData",
    "LinearAligner",
    "RtPredictor",
    "train_test_split",
]

