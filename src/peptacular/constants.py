from enum import StrEnum
from typing import Final, Literal

# ============================================================================
# Enums
# ============================================================================

PROTON_MASS: Final[float] = 1.00727646688
ELECTRON_MASS: Final[float] = 0.00054857990946
NEUTRON_MASS: Final[float] = 1.00866491597
C13_NEUTRON_MASS: Final[float] = 1.003350
PEPTIDE_AVERAGINE_NEUTRON_MASS: Final[float] = 1.002856


class AminoAcid(StrEnum):
    """An amino acid including unconventional ones (U/O) and ambiguous ones (X/J/B/Z)"""

    A = "A"
    B = "B"
    C = "C"
    D = "D"
    E = "E"
    F = "F"
    G = "G"
    H = "H"
    I = "I"
    J = "J"
    K = "K"
    L = "L"
    M = "M"
    N = "N"
    O = "O"
    P = "P"
    Q = "Q"
    R = "R"
    S = "S"
    T = "T"
    U = "U"
    V = "V"
    W = "W"
    X = "X"
    Y = "Y"
    Z = "Z"

    @classmethod
    def from_str(cls, aa: str) -> "AminoAcid":
        """Get AminoAcid enum from string"""
        return cls(aa.upper())


VALID_AMINO_ACIDS = {aa.value for aa in AminoAcid}


class CV(StrEnum):
    """One of the five supported controlled vocabularies"""

    UNIMOD = "UNIMOD"
    PSI_MOD = "MOD"
    RESID = "RESID"
    GNOME = "GNO"
    XL_MOD = "XLMOD"
    CUSTOM = "CUSTOM"
    OBSERVED = "OBSERVED"


CV_TO_NAME_PREFIX = {
    CV.UNIMOD: "",
    CV.PSI_MOD: "",
    CV.RESID: "R:",
    CV.GNOME: "G:",
    CV.XL_MOD: "X:",
    CV.CUSTOM: "C:",
}

CV_TO_ACCESSION_PREFIX = {
    CV.UNIMOD: "UNIMOD:",
    CV.PSI_MOD: "MOD:",
    CV.RESID: "RESID:",
    CV.GNOME: "GNO:",
    CV.XL_MOD: "XLMOD:",
}

CV_TO_MASS_PREFIX = {
    CV.UNIMOD: "UNIMOD:",
    CV.PSI_MOD: "MOD:",
    CV.RESID: "RESID:",
    CV.GNOME: "GNO:",
    CV.XL_MOD: "XLMOD:",
    CV.OBSERVED: "Obs:",
}


class Element(StrEnum):
    """All elements"""

    He = "He"
    Li = "Li"
    Be = "Be"
    Ne = "Ne"
    Na = "Na"
    Mg = "Mg"
    Al = "Al"
    Si = "Si"
    Cl = "Cl"
    Ar = "Ar"
    Ca = "Ca"
    Sc = "Sc"
    Ti = "Ti"
    Cr = "Cr"
    Mn = "Mn"
    Fe = "Fe"
    Co = "Co"
    Ni = "Ni"
    Cu = "Cu"
    Zn = "Zn"
    Ga = "Ga"
    Ge = "Ge"
    As = "As"
    Se = "Se"
    Br = "Br"
    Kr = "Kr"
    Rb = "Rb"
    Sr = "Sr"
    Zr = "Zr"
    Nb = "Nb"
    Mo = "Mo"
    Tc = "Tc"
    Ru = "Ru"
    Rh = "Rh"
    Pd = "Pd"
    Ag = "Ag"
    Cd = "Cd"
    In = "In"
    Sn = "Sn"
    Sb = "Sb"
    Te = "Te"
    Xe = "Xe"
    Cs = "Cs"
    Ba = "Ba"
    La = "La"
    Ce = "Ce"
    Pr = "Pr"
    Nd = "Nd"
    Pm = "Pm"
    Sm = "Sm"
    Eu = "Eu"
    Gd = "Gd"
    Tb = "Tb"
    Dy = "Dy"
    Ho = "Ho"
    Er = "Er"
    Tm = "Tm"
    Yb = "Yb"
    Lu = "Lu"
    Hf = "Hf"
    Ta = "Ta"
    Re = "Re"
    Os = "Os"
    Ir = "Ir"
    Pt = "Pt"
    Au = "Au"
    Hg = "Hg"
    Tl = "Tl"
    Pb = "Pb"
    Bi = "Bi"
    Po = "Po"
    At = "At"
    Rn = "Rn"
    Fr = "Fr"
    Ra = "Ra"
    Ac = "Ac"
    Th = "Th"
    Pa = "Pa"
    Np = "Np"
    Pu = "Pu"
    Am = "Am"
    Cm = "Cm"
    Bk = "Bk"
    Cf = "Cf"
    Es = "Es"
    Fm = "Fm"
    Md = "Md"
    No = "No"
    Lr = "Lr"
    Rf = "Rf"
    Db = "Db"
    Sg = "Sg"
    Bh = "Bh"
    Hs = "Hs"
    Mt = "Mt"
    Ds = "Ds"
    Rg = "Rg"
    Cn = "Cn"
    Nh = "Nh"
    Fl = "Fl"
    Mc = "Mc"
    Lv = "Lv"
    Ts = "Ts"
    Og = "Og"
    U = "U"
    W = "W"
    I = "I"
    Y = "Y"
    V = "V"
    K = "K"
    S = "S"
    B = "B"
    C = "C"
    N = "N"
    O = "O"
    F = "F"
    H = "H"
    P = "P"


class MonosaccharideName(StrEnum):
    """A named monosaccharide"""

    aHex = "aHex"
    Dec = "Dec"
    dHex = "dHex"
    en_aHex = "en,aHex"
    Fuc = "Fuc"
    Hep = "Hep"
    Hex = "Hex"
    HexN = "HexN"
    HexNAc = "HexNAc"
    HexNAcS = "HexNAcS"
    HexNS = "HexNS"
    HexP = "HexP"
    HexS = "HexS"
    Neu = "Neu"
    NeuAc = "NeuAc"
    NeuGc = "NeuGc"
    Non = "Non"
    Oct = "Oct"
    Pen = "Pen"
    Phosphate = "Phosphate"
    S = "S"
    Sug = "Sug"
    Sulfate = "Sulfate"
    Tet = "Tet"
    Tri = "Tri"


class Terminal(StrEnum):
    """Terminal position specification"""

    ANYWHERE = "Anywhere"
    N_TERM = "N-term"
    C_TERM = "C-term"

    @classmethod
    def from_str(cls, term: str) -> "Terminal":
        """Get Terminal enum from string"""
        term_upper = term.upper()
        if term_upper == "N-TERM":
            return cls.N_TERM
        elif term_upper == "C-TERM":
            return cls.C_TERM
        raise ValueError(f"Unknown terminal type: {term}")


# str enum
class ModType(StrEnum):
    NTERM = "nterm"
    CTERM = "cterm"
    ISOTOPE = "isotope"
    STATIC = "static"
    LABILE = "labile"
    UNKNOWN = "unknown"
    INTERVAL = "interval"
    INTERNAL = "internal"
    CHARGE = "charge"


ModTypeLiteral = Literal[
    "nterm",
    "cterm",
    "isotope",
    "static",
    "labile",
    "unknown",
    "interval",
    "internal",
    "charge",
]


class ParrallelMethod(StrEnum):
    PROCESS = "process"
    THREAD = "thread"
    SEQUENTIAL = "sequential"


ParrallelMethodLiteral = Literal["process", "thread", "sequential"]
