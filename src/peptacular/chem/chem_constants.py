"""
chem_constants.py
"""

from typing import Dict

from ..chem.chem_util import chem_mass
from ..constants import (
    AA_COMPOSITIONS,
    AVERAGE_ATOMIC_MASSES,
    AVERAGINE_RATIOS,
    ISOTOPIC_ATOMIC_MASSES,
    NEUTRAL_FRAGMENT_ION_COMPOSITIONS,
)

MONOISOTOPIC_NEUTRAL_FRAGMENT_ION: Dict[str, float] = {
    aa: chem_mass(comp, monoisotopic=True)
    for aa, comp in NEUTRAL_FRAGMENT_ION_COMPOSITIONS.items()
}
AVERAGE_NEUTRAL_FRAGMENT_ION: Dict[str, float] = {
    aa: chem_mass(comp, monoisotopic=False)
    for aa, comp in NEUTRAL_FRAGMENT_ION_COMPOSITIONS.items()
}

# AA MASSES
MONOISOTOPIC_AA_MASSES: Dict[str, float] = {
    aa: chem_mass(comp) for aa, comp in AA_COMPOSITIONS.items()
}
AVERAGE_AA_MASSES: Dict[str, float] = {
    aa: chem_mass(comp, monoisotopic=False) for aa, comp in AA_COMPOSITIONS.items()
}
ISOTOPIC_AVERAGINE_MASS: float = sum(
    v * ISOTOPIC_ATOMIC_MASSES[k] for k, v in AVERAGINE_RATIOS.items()
)
AVERAGE_AVERAGINE_MASS: float = sum(
    v * AVERAGE_ATOMIC_MASSES[k] for k, v in AVERAGINE_RATIOS.items()
)
