"""
chem_constants.py
"""

from typing import Dict

from peptacular.chem.chem_util import chem_mass
from peptacular.constants import NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS, ISOTOPIC_ATOMIC_MASSES, AA_COMPOSITIONS, \
    AVERAGE_ATOMIC_MASSES, AVERAGINE_RATIOS, FRAGMENT_ION_COMPOSITION_ADJUSTMENTS, FRAGMENT_ION_COMPOSITIONS

# NEUTRAL FRAGMENT MASSES
MONOISOTOPIC_FRAGMENT_ADJUSTMENTS: Dict[str, float] = \
    {aa: chem_mass(comp) for aa, comp in NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS.items()}
AVERAGE_FRAGMENT_ADJUSTMENTS: Dict[str, float] = \
    {aa: chem_mass(comp, monoisotopic=False) for aa, comp in NEUTRAL_FRAGMENT_COMPOSITION_ADJUSTMENTS.items()}

MONOISOTOPIC_FRAGMENT_ION_ADJUSTMENTS: Dict[str, float] = {aa: chem_mass(comp) for aa, comp in
                                                           FRAGMENT_ION_COMPOSITIONS.items()}
AVERAGE_FRAGMENT_ION_ADJUSTMENTS: Dict[str, float] = {aa: chem_mass(comp, monoisotopic=False) for aa, comp in
                                                      FRAGMENT_ION_COMPOSITIONS.items()}

# FRAGMENT ION MASSES (+1)
MONOISOTOPIC_ION_ADJUSTMENTS: Dict[str, float] = \
    {aa: chem_mass(comp) for aa, comp in FRAGMENT_ION_COMPOSITION_ADJUSTMENTS.items()}
AVERAGE_ION_ADJUSTMENTS: Dict[str, float] = \
    {aa: chem_mass(comp, monoisotopic=False) for aa, comp in FRAGMENT_ION_COMPOSITION_ADJUSTMENTS.items()}

# AA MASSES
MONOISOTOPIC_AA_MASSES: Dict[str, float] = {aa: chem_mass(comp) for aa, comp in AA_COMPOSITIONS.items()}
AVERAGE_AA_MASSES: Dict[str, float] = {aa: chem_mass(comp, monoisotopic=False) for aa, comp in AA_COMPOSITIONS.items()}
ISOTOPIC_AVERAGINE_MASS: float = sum(v * ISOTOPIC_ATOMIC_MASSES[k] for k, v in AVERAGINE_RATIOS.items())
AVERAGE_AVERAGINE_MASS: float = sum(v * AVERAGE_ATOMIC_MASSES[k] for k, v in AVERAGINE_RATIOS.items())
