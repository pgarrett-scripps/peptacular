from enum import StrEnum

from .data import *  # type: ignore (Too many imports to care)


class SecondaryStructureMethod(StrEnum):
    """Enum for secondary structure prediction methods."""

    DELEAGE_ROUX = "DeleageRoux"
    LEVITT = "Levitt"
    CHOU_FASMAN = "ChouFasman"


class SecondaryStructureType(StrEnum):
    """Enum for secondary structure types."""

    ALPHA_HELIX = "alpha_helix"
    BETA_SHEET = "beta_sheet"
    BETA_TURN = "beta_turn"
    COIL = "coil"


class PropertyScale(StrEnum):
    """Base enum for all amino acid property scales."""

    pass


class HydrophobicityScale(PropertyScale):
    """Hydrophobicity scales."""

    KYTE_DOOLITTLE = "hphob_kyte_doolittle"
    ADOBERIN = "hphob_adoberin"
    ABRAHAM_LEO = "hphob_abraham_leo"
    AGROS = "hphob_agros"
    RAO_ARGOS = "hphob_rao_argos"
    BLACK_MOULD = "hphob_black_mould"
    BULL_BREESE = "hphob_bull_breese"
    CASARI_SIPPL = "hphob_casari_sippl"
    CID = "hphob_cid"
    COWAN_3_4 = "hphob_cowan_3_4"
    COWAN_7_5 = "hphob_cowan_7_5"
    EISENBERG = "hphob_eisenberg"
    ENGELMAN = "hphob_engelman"
    FASMAN = "hphob_fasman"
    FAUCHERE = "hphob_fauchere"
    GOLDSACK = "hphob_goldsack"
    GUY = "hphob_guy"
    JONES = "hphob_jones"
    JURETIC = "hphob_juretic"
    KIDERA = "hphob_kidera"
    MIYAZAWA = "hphob_miyazawa"
    PARKER = "hphob_parker"
    PONNUSWAMY = "hphob_ponnuswamy"
    MANAVALAN = "hphob_manavalan"
    ROSE = "hphob_rose"
    ROSEMAN = "hphob_roseman"
    SWEET = "hphob_sweet"
    TANFORD = "hphob_tanford"
    WILSON = "hphob_wilson"
    ZIMMERMAN = "hphob_zimmerman"
    CHOTHIA = "hphob_chothia"
    JANIN = "hphob_janin"
    WOLFENDEN = "hphob_wolfenden"
    WELLING = "hphob_welling"


class SecondaryStructureScale(PropertyScale):
    """Secondary structure prediction scales."""

    DELEAGE_ROUX_ALPHA_HELIX = "deleage_roux_alpha_helix"
    DELEAGE_ROUX_BETA_SHEET = "deleage_roux_beta_sheet"
    DELEAGE_ROUX_BETA_TURN = "deleage_roux_beta_turn"
    DELEAGE_ROUX_COIL = "deleage_roux_coil"
    LEVITT_ALPHA_HELIX = "levitt_alpha_helix"
    LEVITT_BETA_SHEET = "levitt_beta_sheet"
    LEVITT_BETA_TURN = "levitt_beta_turn"
    CHOU_FASMAN_ALPHA_HELIX = "chou_fasman_alpha_helix"
    CHOU_FASMAN_BETA_SHEET = "chou_fasman_beta_sheet"
    CHOU_FASMAN_BETA_TURN = "chou_fasman_beta_turn"


class SurfaceAccessibilityScale(PropertyScale):
    """Surface accessibility and area scales."""

    VERGOTEN = "surface_accessibility_vergoten"
    JANIN = "surface_accessiblility_janin"  # Note: keeping original typo
    ACCESSIBLE_RESIDUES = "accessible_residues"
    AVERAGE_BURIED_AREA = "average_buried_area"


class ChargeScale(PropertyScale):
    """Charge and pK related scales."""

    PK_NTERMINAL = "pk_nterminal"
    PK_CTERMINAL = "pk_cterminal"
    PK_SIDECHAIN = "pk_sidechain"


class PolarityScale(PropertyScale):
    """Polarity scales."""

    GRANTHAM = "polarity_grantham"
    ZIMMERMAN = "polarity_zimmerman"


class HPLCScale(PropertyScale):
    """HPLC retention time scales."""

    MEEK_2_1 = "hplc_meek_2_1"
    BROWNE = "hplc_browne"
    MEEK_7_4 = "hplc_meek_7_4"
    BROWNE_TFA = "hplc_browne_tfa"


class BetaStrandScale(PropertyScale):
    """Beta strand propensity scales."""

    PARALLEL = "beta_strand_parallel"
    ANTIPARALLEL = "beta_strand_antiparallel"
    TOTAL = "beta_strand_total"


class PhysicalPropertyScale(PropertyScale):
    """Physical and chemical property scales."""

    MOLECULAR_WEIGHTS = "molecular_weights"
    BULKINESS = "bulkiness"
    REFRACTIVITY = "refractivity"
    FLEXIBILITY_VIHINEN = "flexibility_vihinen"
    HYDROPHILICITY_HOP_WOOD = "hydrophilicity_hop_wood"
    RATIOSIDE = "ratioside"
    MUTABILITY = "mutability"
    CODONS = "codons"
    RECOGNITION_FACTORS = "recognition_factors"
    TRANSMEMBRANE_TENDENCY = "transmembrane_tendency"


class CompositionScale(PropertyScale):
    """Amino acid composition scales."""

    MCCALDRON = "aa_composition_mccaldron"
    SWISSPROT = "aa_composition_swissprot"


# Updated dictionaries using the enums
secondary_structure_scales: dict[str, dict[str, float]] = {
    SecondaryStructureScale.DELEAGE_ROUX_ALPHA_HELIX: deleage_roux_alpha_helix,
    SecondaryStructureScale.DELEAGE_ROUX_BETA_SHEET: deleage_roux_beta_sheet,
    SecondaryStructureScale.DELEAGE_ROUX_BETA_TURN: deleage_roux_beta_turn,
    SecondaryStructureScale.DELEAGE_ROUX_COIL: deleage_roux_coil,
    SecondaryStructureScale.LEVITT_ALPHA_HELIX: levitt_alpha_helix,
    SecondaryStructureScale.LEVITT_BETA_SHEET: levitt_beta_sheet,
    SecondaryStructureScale.LEVITT_BETA_TURN: levitt_beta_turn,
    SecondaryStructureScale.CHOU_FASMAN_ALPHA_HELIX: chou_fasman_alpha_helix,
    SecondaryStructureScale.CHOU_FASMAN_BETA_SHEET: chou_fasman_beta_sheet,
    SecondaryStructureScale.CHOU_FASMAN_BETA_TURN: chou_fasman_beta_turn,
}

secondary_structure_scales_by_name: dict[str, dict[str, dict[str, float]]] = {
    SecondaryStructureMethod.DELEAGE_ROUX: {
        SecondaryStructureType.ALPHA_HELIX: deleage_roux_alpha_helix,
        SecondaryStructureType.BETA_SHEET: deleage_roux_beta_sheet,
        SecondaryStructureType.BETA_TURN: deleage_roux_beta_turn,
        SecondaryStructureType.COIL: deleage_roux_coil,
    },
    SecondaryStructureMethod.LEVITT: {
        SecondaryStructureType.ALPHA_HELIX: levitt_alpha_helix,
        SecondaryStructureType.BETA_SHEET: levitt_beta_sheet,
        SecondaryStructureType.BETA_TURN: levitt_beta_turn,
    },
    SecondaryStructureMethod.CHOU_FASMAN: {
        SecondaryStructureType.ALPHA_HELIX: chou_fasman_alpha_helix,
        SecondaryStructureType.BETA_SHEET: chou_fasman_beta_sheet,
        SecondaryStructureType.BETA_TURN: chou_fasman_beta_turn,
    },
}


all_property_scales: dict[str, dict[str, float]] = {
    # Hydrophobicity scales
    HydrophobicityScale.KYTE_DOOLITTLE: hphob_kyte_doolittle,
    HydrophobicityScale.ADOBERIN: hphob_adoberin,
    HydrophobicityScale.ABRAHAM_LEO: hphob_abraham_leo,
    HydrophobicityScale.AGROS: hphob_agros,
    HydrophobicityScale.RAO_ARGOS: hphob_rao_argos,
    HydrophobicityScale.BLACK_MOULD: hphob_black_mould,
    HydrophobicityScale.BULL_BREESE: hphob_bull_breese,
    HydrophobicityScale.CASARI_SIPPL: hphob_casari_sippl,
    HydrophobicityScale.CID: hphob_cid,
    HydrophobicityScale.COWAN_3_4: hphob_cowan_3_4,
    HydrophobicityScale.COWAN_7_5: hphob_cowan_7_5,
    HydrophobicityScale.EISENBERG: hphob_eisenberg,
    HydrophobicityScale.ENGELMAN: hphob_engelman,
    HydrophobicityScale.FASMAN: hphob_fasman,
    HydrophobicityScale.FAUCHERE: hphob_fauchere,
    HydrophobicityScale.GOLDSACK: hphob_goldsack,
    HydrophobicityScale.GUY: hphob_guy,
    HydrophobicityScale.JONES: hphob_jones,
    HydrophobicityScale.JURETIC: hphob_juretic,
    HydrophobicityScale.KIDERA: hphob_kidera,
    HydrophobicityScale.MIYAZAWA: hphob_miyazawa,
    HydrophobicityScale.PARKER: hphob_parker,
    HydrophobicityScale.PONNUSWAMY: hphob_ponnuswamy,
    HydrophobicityScale.MANAVALAN: hphob_manavalan,
    HydrophobicityScale.ROSE: hphob_rose,
    HydrophobicityScale.ROSEMAN: hphob_roseman,
    HydrophobicityScale.SWEET: hphob_sweet,
    HydrophobicityScale.TANFORD: hphob_tanford,
    HydrophobicityScale.WILSON: hphob_wilson,
    HydrophobicityScale.ZIMMERMAN: hphob_zimmerman,
    HydrophobicityScale.CHOTHIA: hphob_chothia,
    HydrophobicityScale.JANIN: hphob_janin,
    HydrophobicityScale.WOLFENDEN: hphob_wolfenden,
    HydrophobicityScale.WELLING: hphob_welling,
    # Secondary structure scales
    SecondaryStructureScale.DELEAGE_ROUX_ALPHA_HELIX: deleage_roux_alpha_helix,
    SecondaryStructureScale.DELEAGE_ROUX_BETA_SHEET: deleage_roux_beta_sheet,
    SecondaryStructureScale.DELEAGE_ROUX_BETA_TURN: deleage_roux_beta_turn,
    SecondaryStructureScale.DELEAGE_ROUX_COIL: deleage_roux_coil,
    SecondaryStructureScale.LEVITT_ALPHA_HELIX: levitt_alpha_helix,
    SecondaryStructureScale.LEVITT_BETA_SHEET: levitt_beta_sheet,
    SecondaryStructureScale.LEVITT_BETA_TURN: levitt_beta_turn,
    SecondaryStructureScale.CHOU_FASMAN_ALPHA_HELIX: chou_fasman_alpha_helix,
    SecondaryStructureScale.CHOU_FASMAN_BETA_SHEET: chou_fasman_beta_sheet,
    SecondaryStructureScale.CHOU_FASMAN_BETA_TURN: chou_fasman_beta_turn,
    # Surface accessibility scales
    SurfaceAccessibilityScale.VERGOTEN: surface_accessibility_vergoten,
    SurfaceAccessibilityScale.JANIN: surface_accessiblility_janin,
    SurfaceAccessibilityScale.ACCESSIBLE_RESIDUES: accessible_residues,
    SurfaceAccessibilityScale.AVERAGE_BURIED_AREA: average_buried_area,
    # Charge scales
    ChargeScale.PK_NTERMINAL: pk_nterminal,
    ChargeScale.PK_CTERMINAL: pk_cterminal,
    ChargeScale.PK_SIDECHAIN: pk_sidechain,
    # Polarity scales
    PolarityScale.GRANTHAM: polarity_grantham,
    PolarityScale.ZIMMERMAN: polarity_zimmerman,
    # HPLC scales
    HPLCScale.MEEK_2_1: hplc_meek_2_1,
    HPLCScale.BROWNE: hplc_browne,
    HPLCScale.MEEK_7_4: hplc_meek_7_4,
    HPLCScale.BROWNE_TFA: hplc_browne_tfa,
    # Beta strand scales
    BetaStrandScale.PARALLEL: beta_strand_parallel,
    BetaStrandScale.ANTIPARALLEL: beta_strand_antiparallel,
    BetaStrandScale.TOTAL: beta_strand_total,
    # Physical property scales
    PhysicalPropertyScale.MOLECULAR_WEIGHTS: molecular_weights,
    PhysicalPropertyScale.BULKINESS: bulkiness,
    PhysicalPropertyScale.REFRACTIVITY: refractivity,
    PhysicalPropertyScale.FLEXIBILITY_VIHINEN: flexibility_vihinen,
    PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD: hydrophilicity_hopp_wood,
    PhysicalPropertyScale.RATIOSIDE: ratioside,
    PhysicalPropertyScale.MUTABILITY: mutability,
    PhysicalPropertyScale.CODONS: codons,
    PhysicalPropertyScale.RECOGNITION_FACTORS: recognition_factors,
    PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY: transmembrane_tendency,
    # Composition scales
    CompositionScale.MCCALDRON: aa_composition_mccaldron,
    CompositionScale.SWISSPROT: aa_composition_swissprot,
}

hydrophobicity_scales: dict[str, dict[str, float]] = {
    HydrophobicityScale.KYTE_DOOLITTLE: hphob_kyte_doolittle,
    HydrophobicityScale.ADOBERIN: hphob_adoberin,
    HydrophobicityScale.ABRAHAM_LEO: hphob_abraham_leo,
    HydrophobicityScale.AGROS: hphob_agros,
    HydrophobicityScale.BLACK_MOULD: hphob_black_mould,
    HydrophobicityScale.BULL_BREESE: hphob_bull_breese,
    HydrophobicityScale.CASARI_SIPPL: hphob_casari_sippl,
    HydrophobicityScale.CID: hphob_cid,
    HydrophobicityScale.COWAN_3_4: hphob_cowan_3_4,
    HydrophobicityScale.COWAN_7_5: hphob_cowan_7_5,
    HydrophobicityScale.EISENBERG: hphob_eisenberg,
    HydrophobicityScale.ENGELMAN: hphob_engelman,
    HydrophobicityScale.FASMAN: hphob_fasman,
    HydrophobicityScale.FAUCHERE: hphob_fauchere,
    HydrophobicityScale.GOLDSACK: hphob_goldsack,
    HydrophobicityScale.GUY: hphob_guy,
    HydrophobicityScale.JONES: hphob_jones,
    HydrophobicityScale.JURETIC: hphob_juretic,
    HydrophobicityScale.KIDERA: hphob_kidera,
    HydrophobicityScale.MIYAZAWA: hphob_miyazawa,
    HydrophobicityScale.PARKER: hphob_parker,
    HydrophobicityScale.PONNUSWAMY: hphob_ponnuswamy,
    HydrophobicityScale.ROSE: hphob_rose,
    HydrophobicityScale.ROSEMAN: hphob_roseman,
    HydrophobicityScale.SWEET: hphob_sweet,
    HydrophobicityScale.TANFORD: hphob_tanford,
    HydrophobicityScale.WILSON: hphob_wilson,
    HydrophobicityScale.ZIMMERMAN: hphob_zimmerman,
    HydrophobicityScale.RAO_ARGOS: hphob_rao_argos,
    HydrophobicityScale.MANAVALAN: hphob_manavalan,
    HydrophobicityScale.CHOTHIA: hphob_chothia,
    HydrophobicityScale.JANIN: hphob_janin,
    HydrophobicityScale.WOLFENDEN: hphob_wolfenden,
    HydrophobicityScale.WELLING: hphob_welling,
}


# Surface accessibility scales
surface_accessibility_scales: dict[str, dict[str, float]] = {
    SurfaceAccessibilityScale.VERGOTEN: surface_accessibility_vergoten,
    SurfaceAccessibilityScale.JANIN: surface_accessiblility_janin,
    SurfaceAccessibilityScale.ACCESSIBLE_RESIDUES: accessible_residues,
    SurfaceAccessibilityScale.AVERAGE_BURIED_AREA: average_buried_area,
}


hplc_scales: dict[str, dict[str, float]] = {
    HPLCScale.MEEK_2_1: hplc_meek_2_1,
    HPLCScale.BROWNE: hplc_browne,
    HPLCScale.MEEK_7_4: hplc_meek_7_4,
    HPLCScale.BROWNE_TFA: hplc_browne_tfa,
}


# Hydrophilicity scales
hydrophilicity_scales: dict[str, dict[str, float]] = {
    PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD: hydrophilicity_hopp_wood
}

# Flexibility scales
flexibility_scales: dict[str, dict[str, float]] = {
    PhysicalPropertyScale.FLEXIBILITY_VIHINEN: flexibility_vihinen
}

# Polarity scales (enum-keyed)
polarity_scales: dict[str, dict[str, float]] = {
    PolarityScale.GRANTHAM: polarity_grantham,
    PolarityScale.ZIMMERMAN: polarity_zimmerman,
}

# Composition scales (enum-keyed)
composition_scales: dict[str, dict[str, float]] = {
    CompositionScale.MCCALDRON: aa_composition_mccaldron,
    CompositionScale.SWISSPROT: aa_composition_swissprot,
}

# Physical property scales — exhaustive mapping using PhysicalPropertyScale
physical_property_scales: dict[str, dict[str, float]] = {
    PhysicalPropertyScale.MOLECULAR_WEIGHTS: molecular_weights,
    PhysicalPropertyScale.BULKINESS: bulkiness,
    PhysicalPropertyScale.REFRACTIVITY: refractivity,
    PhysicalPropertyScale.FLEXIBILITY_VIHINEN: flexibility_vihinen,
    PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD: hydrophilicity_hopp_wood,
    PhysicalPropertyScale.RATIOSIDE: ratioside,
    PhysicalPropertyScale.MUTABILITY: mutability,
    PhysicalPropertyScale.CODONS: codons,
    PhysicalPropertyScale.RECOGNITION_FACTORS: recognition_factors,
    PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY: transmembrane_tendency,
}
