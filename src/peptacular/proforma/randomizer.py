"""
Randomizer for ProForma annotations.
"""

from random import randint, choice, sample
from enum import Enum, auto
from typing import List, Union

from peptacular.proforma.proforma_dataclasses import Mod, Interval
from peptacular.proforma.proforma_parser import ProFormaAnnotation


_UNIMOD_LEVEL_BASE_MOD_VALS = ['Oxidation', 'UNIMOD:10']
_UNIMOD_LEVEL2_MOD_VALS = ['U:Oxidation', 'U:10', 'U:+1', 'U:-1', 'U:+3.1415', 'U:-3.1415']
_PSI_LEVEL_BASE_MOD_VALS = ['O-phospho-L-serine', 'MOD:00046']
_PSI_LEVEL2_MOD_VALS = ['M:O-phospho-L-serine', 'M:00046', 'M:+1', 'M:-1', 'M:+3.1415', 'M:-3.1415']
_DELTA_MASS_MOD_VALS = ['+1', '-1', '+3.1415', '-3.1415']
_MOD_INFO_VALS = ['INFO:Cool', 'INFO:Awesome', 'INFO:Radical', 'INFO:Amazing', 'INFO:Fantastic']
_CHEM_FORMULA_MOD_VALS = ['Formula:C12H22O11', 'Formula:[13C6]H12O6[12C-4]', 'Formula:CHO', 'Formula:C2H-5O']
_GLYCAN_MOD_VALS = ['Glycan:HexNAc2Hex3Neu1', 'Glycan:Hex', 'Glycan:6BAAE1B1']
_GNO_MOD_VALS = ['GNO:G59626AS', 'GNO:G62765YT', 'G:G59626AS', 'G:G62765YT', 'G:+1', 'G:-1', 'G:+3.1415', 'G:-3.1415']
_RESID_MOD_VALS = ['RESID:AA0581', 'RESID:AA0037', 'R:AA0581', 'R:AA0037', 'R:+1', 'R:-1', 'R:+3.1415', 'R:-3.1415']
_ISOTOPE_MOD_VALS = ['13C', '15N', '18O', '2H', 'T', 'D']
_STATIC_MOD_VALS = ['[Oxidation]@M', '[Oxidation]@M,C,D', '[+1]@C', '[-1]@C', '[+3.1415]@C', '[-3.1415]@C']
_XLMOD_VALS = ['XLMOD:02001', 'XLMOD:02010', 'XLMOD:02000', 'X:02001', 'X:02010', 'X:02000']
_CHARGE_ADDUCT_VALS = ['+H+', '+2Na+,-H+', '+2Na+,+H+', '2I-', '+e-']

_TOP_DOWN_MODS = _CHEM_FORMULA_MOD_VALS + _RESID_MOD_VALS
_CROSS_LINKING_MODS = _XLMOD_VALS
_GLYCAN_MODS = _GLYCAN_MOD_VALS + _GNO_MOD_VALS

_BASE_AMINO_ACIDS = "VWPSDCYTAIMHGQENFLKR"
_LEVEL2_AMINO_ACIDS = _BASE_AMINO_ACIDS + 'OUBZXJ'
_LEVEL2_AMINO_ACIDS_WITHOUT_AMBIGUITY = _BASE_AMINO_ACIDS + 'OU'

_BASE_MODS = _UNIMOD_LEVEL_BASE_MOD_VALS + _PSI_LEVEL_BASE_MOD_VALS + _DELTA_MASS_MOD_VALS
_LEVEL2_MODS = _UNIMOD_LEVEL2_MOD_VALS + _PSI_LEVEL2_MOD_VALS + _BASE_MODS


class ProformaComplianceLevel(Enum):
    """
    Enum for ProForma compliance levels.
    """
    BASE = auto()
    LEVEL2 = auto()
    TOP_DOWN = auto()
    CROSS_LINKING = auto()
    GLYCAN = auto()
    SPECTRUM = auto()


def _random_sequence(amino_acids: str, min_sequence_length: int, max_sequence_length: int) -> str:
    """
    Generate a random sequence of amino acids.
    """
    return ''.join(choice(amino_acids) for _ in range(randint(min_sequence_length, max_sequence_length)))


def random_sequence(level: ProformaComplianceLevel, min_sequence_length: int = 5, max_sequence_length: int = 50,
                    sequence_ambiguity: bool = True) -> str:
    """
    Generate a random sequence based on the compliance level.
    """
    if level == ProformaComplianceLevel.BASE:
        return _random_sequence(_BASE_AMINO_ACIDS, min_sequence_length, max_sequence_length)

    if level == ProformaComplianceLevel.LEVEL2:
        if sequence_ambiguity:
            return _random_sequence(_LEVEL2_AMINO_ACIDS, min_sequence_length, max_sequence_length)
        return _random_sequence(_LEVEL2_AMINO_ACIDS_WITHOUT_AMBIGUITY, min_sequence_length, max_sequence_length)

    raise ValueError("Invalid level")


def _random_mod(mods: List[str], count: int, info: bool) -> Mod:
    mod = choice(mods)
    if info:
        for _ in range(randint(1, 2)):
            mod += f"|{choice(_MOD_INFO_VALS)}"
    return Mod(mod, count)


def random_mod(level: ProformaComplianceLevel, count: int = 1, info: bool = False) -> Mod:
    """
    Generate a random modification based on the compliance level.
    """
    if level == ProformaComplianceLevel.BASE:
        return _random_mod(_BASE_MODS, count, info)
    if level == ProformaComplianceLevel.LEVEL2:
        return _random_mod(_LEVEL2_MODS, count, info)
    if level == ProformaComplianceLevel.TOP_DOWN:
        return _random_mod(_TOP_DOWN_MODS, count, info)
    if level == ProformaComplianceLevel.CROSS_LINKING:
        return _random_mod(_CROSS_LINKING_MODS, count, info)
    if level == ProformaComplianceLevel.GLYCAN:
        return _random_mod(_GLYCAN_MODS, count, info)
    if level == ProformaComplianceLevel.SPECTRUM:
        return _random_mod([], count, info)

    raise ValueError("Invalid level")


def random_interval(level: ProformaComplianceLevel, start_index: int, end_index: int) -> Interval:
    """
    Generate a random interval within a given range.
    """
    return _random_interval(level, start_index, end_index)


def _random_interval(level: ProformaComplianceLevel, start_index: int, end_index: int) -> Interval:
    """
    Generate a random interval within a given range.
    """
    # Ensure 'end' is strictly greater than 'start' to avoid zero-length intervals
    start = randint(start_index, end_index - 1)
    end = randint(start + 1, end_index)
    ambiguous = choice([True, False])
    mods = [random_mod(level, 1) for _ in range(randint(0, 2))]
    return Interval(start, end, ambiguous, mods)


def random_intervals(level: ProformaComplianceLevel, sequence: str, num_intervals: int) -> List[Interval]:
    """
    Generate a list of random intervals within a given sequence.
    """
    intervals = []
    if num_intervals == 0:
        return intervals
    segment_length = len(sequence) // num_intervals
    for i in range(num_intervals):
        start_index = segment_length * i
        end_index = start_index + segment_length - 1 if i < num_intervals - 1 else len(sequence)
        interval = _random_interval(level, start_index, end_index)
        intervals.append(interval)
    return intervals


def compliance_randomizer(level: Union[ProformaComplianceLevel],
                          min_sequence_length: int = 5,
                          max_sequence_length: int = 50,
                          mod_prob: float = 0.1,
                          sequence_ambiguity: bool = True) -> ProFormaAnnotation:
    """
    1) Base Level Support (Technical name: Base-ProForma Compliant)
    Represents the lowest level of compliance, this level involves providing support for:
    - Amino acid sequences
    - Protein modifications using two of the supported CVs/ontologies: Unimod and PSIMOD.
    - Protein modifications using delta masses (without prefixes)
    - N-terminal, C-terminal and labile modifications.
    - Ambiguity in the modification position, including support for localisation scores.
    - Ambiguity in the amino acid sequence.
    - INFO tag.

    2) Additional Separate Support (Technical name: level 2-ProForma compliant)
    These features are independent of each other:
    - Unusual amino acids (O and U).
    - Ambiguous amino acids (e.g. X, B, Z). This would include support for sequence tags of
    known mass (using the character X).
    - Protein modifications using delta masses (using prefixes for the different
    CVs/ontologies).
    - Use of prefixes for Unimod (U:) and PSI-MOD (M:) names.
    - Support for the joint representation of experimental data and its interpretation.
    """
    if isinstance(level, int):
        if level == 1:
            level = ProformaComplianceLevel.BASE
        elif level == 2:
            level = ProformaComplianceLevel.LEVEL2
        else:
            raise ValueError("Invalid level")

    if level not in (ProformaComplianceLevel.BASE, ProformaComplianceLevel.LEVEL2):
        raise ValueError("Invalid level")

    # Amino acid sequences
    sequence = random_sequence(level, min_sequence_length, max_sequence_length, sequence_ambiguity)
    annotation = ProFormaAnnotation(sequence)

    # Protein modifications using two of the supported CVs/ontologies: Unimod and PSIMOD.
    # Protein modifications using delta masses (without prefixes)
    num_internal_mods = int(len(sequence) * mod_prob)
    internal_mods = [random_mod(level=level, count=1, info=choice([True, False]))
                     for _ in range(num_internal_mods)]
    internal_mod_indices = sample(range(len(sequence)), num_internal_mods)

    for i, mod in zip(internal_mod_indices, internal_mods):
        annotation.add_internal_mod(i, mod)

    # N-terminal, C-terminal and labile modifications.
    for _ in range(randint(0, 3)):
        annotation.add_nterm_mods(random_mod(level=level, count=1, info=choice([True, False])))
    for _ in range(randint(0, 3)):
        annotation.add_cterm_mods(random_mod(level=level, count=1, info=choice([True, False])))
    for _ in range(randint(0, 3)):
        annotation.add_labile_mods(random_mod(level=level, count=1, info=choice([True, False])))
    for _ in range(randint(0, 3)):
        annotation.add_unknown_mods(random_mod(level=level, count=randint(1, 3), info=choice([True, False])))

    # Ambiguity in the modification position, including support for localisation scores.
    # Ambiguity in the amino acid sequence.
    intervals = random_intervals(level, sequence, randint(0, 2))
    annotation.add_intervals(intervals)

    # Localization score
    if choice([True, False]):
        score_mod = random_mod(level=level, count=1, info=False)
        score_mod = Mod(str(score_mod.val) + '#g1(0.9)', 1)
        additional_mods = [Mod('#g1(0.1)', 1) for _ in range(randint(1, 3))]
        score_mods = [score_mod] + additional_mods
        score_mods_indices = sample(range(len(sequence)), len(score_mods))

        for i, mod in zip(score_mods_indices, score_mods):
            # annotation.pop_internal_mod(i)
            annotation.add_internal_mod(i, mod)

    return annotation


def top_down_randomizer(annotation: ProFormaAnnotation, mod_prob: float = 0.1):
    """
    3) Top-Down Extensions (Technical name: level 2-ProForma + top-down compliant)
    - Additional CV/ontologies for protein modifications: RESID (the prefix R MUST be used
    for RESID CV/ontology term names)
    - Chemical formulas (this feature occurs in two places in this list).
    """

    # Protein modifications using two of the supported CVs/ontologies: Unimod and PSIMOD.
    # Protein modifications using delta masses (without prefixes)
    num_internal_mods = int(len(annotation) * mod_prob)
    internal_mods = [random_mod(level=ProformaComplianceLevel.TOP_DOWN, count=1, info=choice([True, False]))
                     for _ in range(num_internal_mods)]
    internal_mod_indices = sample(range(len(annotation)), num_internal_mods)

    for i, mod in zip(internal_mod_indices, internal_mods):
        annotation.add_internal_mod(i, mod)


def cross_linking_randomizer(annotation: ProFormaAnnotation):
    """
    4) Cross-Linking Extensions (Technical name: level 2-ProForma + cross-linking
    compliant)
    - Cross-linked peptides (using the XL-MOD CV/ontology, the prefix X MUST be used for
    XL-MOD CV/ontology term names).
    """

    # Cross-linking
    if choice([True, False]):
        cross_link_mod = random_mod(level=ProformaComplianceLevel.CROSS_LINKING, count=1, info=False)
        cross_link_mod = Mod(str(cross_link_mod.val) + '#XL1', 1)
        additional_mods = [Mod('#XL1', 1) for _ in range(randint(1, 3))]
        cross_link_mods = [cross_link_mod] + additional_mods
        score_mods_indices = sample(range(len(annotation)), len(cross_link_mods))

        for i, mod in zip(score_mods_indices, cross_link_mods):
            # annotation.pop_internal_mod(i)
            annotation.add_internal_mod(i, mod)


def glycan_randomizer(annotation: ProFormaAnnotation, mod_prob: float = 0.1):
    """
    5) Glycan Extensions (Technical name: level 2-ProForma + glycans compliant)
    - Additional CV/ontologies for protein modifications: GNO (the prefix G MUST be used
    for GNO CV/ontology term names)
    - Glycan composition.
    - Chemical formulas (this feature occurs in two places in this list).
    """

    # Protein modifications using two of the supported CVs/ontologies: Unimod and PSIMOD.
    # Protein modifications using delta masses (without prefixes)
    num_internal_mods = int(len(annotation) * mod_prob)
    internal_mods = [random_mod(level=ProformaComplianceLevel.GLYCAN, count=1, info=choice([True, False]))
                     for _ in range(num_internal_mods)]
    internal_mod_indices = sample(range(len(annotation)), num_internal_mods)

    for i, mod in zip(internal_mod_indices, internal_mods):
        annotation.add_internal_mod(i, mod)


def spectrum_randomizer(annotation: ProFormaAnnotation):
    """
    6) Spectral Support (Technical name: level 2-ProForma + mass spectrum compliant)
    - Charge and chimeric spectra are special cases (see Appendix II).
    - Global modifications (e.g., every C is C13).
    """

    # Add charge and isotope type
    if choice([True, False]):
        annotation.add_charge(choice([1, 2, 3]))

    for _ in range(randint(0, 3)):
        mod = _random_mod(_ISOTOPE_MOD_VALS, 1, False)
        annotation.add_isotope_mods(mod)

    for _ in range(randint(0, 3)):
        mod = _random_mod(_STATIC_MOD_VALS, 1, False)
        annotation.add_static_mods(mod)

    if annotation.has_charge():

        # Add adducts
        if choice([True, False]):
            annotation.add_charge_adducts(Mod(choice(_CHARGE_ADDUCT_VALS), 1))
