# Copyright 2003 Yair Benita.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Indices to be used with ProtParam."""

# Turn black code style off
# fmt: off

# Hydrophobicity

# Kyte & Doolittle index of hydrophobicity
# J. Mol. Biol. 157:105-132(1982).
# "KyteDoolittle"
hphob_kyte_doolittle = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
      "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
      "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

# Aboderin hydrophobicity index
# International J. of Biochemistry, 2(11), 537-544.
# "Aboderin"
hphob_adoberin = {"A": 5.1, "R": 2.0, "N": 0.6, "D": 0.7, "C": 0.0,
      "Q": 1.4, "E": 1.8, "G": 4.1, "H": 1.6, "I": 9.3,
      "L": 10.0, "K": 1.3, "M": 8.7, "F": 9.6, "P": 4.9,
      "S": 3.1, "T": 3.5, "W": 9.2, "Y": 8.0, "V": 8.5}

# Abraham & Leo hydrophobicity index
# Proteins: Structure, Function and Genetics 2:130-152(1987).
# "AbrahamLeo"
hphob_abraham_leo = {"A": 0.44, "R": -2.42, "N": -1.32, "D": -0.31, "C": 0.58,
      "Q": -0.71, "E": -0.34, "G": 0.0, "H": -0.01, "I": 2.46,
      "L": 2.46, "K": -2.45, "M": 1.1, "F": 2.54, "P": 1.29,
      "S": -0.84, "T": -0.41, "W": 2.56, "Y": 1.63, "V": 1.73}

# Argos hydrophobicity index
# European Journal of Biochemistry, 128(2-3), 565-575.
# "Argos"
hphob_agros = {"A": 0.61, "R": 0.6, "N": 0.06, "D": 0.46, "C": 1.07,
      "Q": 0.0, "E": 0.47, "G": 0.07, "H": 0.61, "I": 2.22,
      "L": 1.53, "K": 1.15, "M": 1.18, "F": 2.02, "P": 1.95,
      "S": 0.05, "T": 0.05, "W": 2.65, "Y": 1.88, "V": 1.32}

"""
Amino acid scale: Membrane buried helix parameter.
Author(s):
Rao M.J.K., Argos P.
Reference: Biochim. Biophys. Acta 869:197-214(1986).
https://web.expasy.org/protscale/pscale/Hphob.Argos.html
"""
hphob_rao_argos = {"A": 1.360, "R": 0.150, "N": 0.330, "D": 0.110, "C": 1.270,
      "Q": 0.330, "E": 0.250, "G": 1.090, "H": 0.680, "I": 1.440,
      "L": 1.470, "K": 0.090, "M": 1.420, "F": 1.570, "P": 0.540,
      "S": 0.970, "T": 1.080, "W": 1.000, "Y": 0.830, "V": 1.370}

# Black & Mould hydrophobicity index
# Anal. Biochem. 193:72-82(1991).
# "BlackMould"
hphob_black_mould = {"A": 0.616, "R": 0.0, "N": 0.236, "D": 0.028, "C": 0.68,
      "Q": 0.251, "E": 0.043, "G": 0.501, "H": 0.165, "I": 0.943,
      "L": 0.943, "K": 0.283, "M": 0.738, "F": 1.0, "P": 0.711,
      "S": 0.359, "T": 0.45, "W": 0.878, "Y": 0.88, "V": 0.825}

# Bull & Breese hydrophobicity index
# Arch. Biochem. Biophys. 161:665-670(1974)
# "BullBreese"
hphob_bull_breese = {"A": 0.61, "R": 0.69, "N": 0.89, "D": 0.61, "C": 0.36,
      "Q": 0.97, "E": 0.51, "G": 0.81, "H": 0.69, "I": -1.45,
      "L": -1.65, "K": 0.46, "M": -0.66, "F": -1.52, "P": -0.17,
      "S": 0.42, "T": 0.29, "W": -1.2, "Y": -1.43, "V": -0.75}

# Casari & Sippl hydrophobic potential
# Journal of molecular biology, 224(3), 725-732.
# "Casari"
hphob_casari_sippl = {"A": 0.2, "R": -0.7, "N": -0.5, "D": -1.4, "C": 1.9,
      "Q": -1.1, "E": -1.3, "G": -0.1, "H": 0.4, "I": 1.4,
      "L": 0.5, "K": -1.6, "M": 0.5, "F": 1.0, "P": -1.0,
      "S": -0.7, "T": -0.4, "W": 1.6, "Y": 0.5, "V": 0.7}

# Cid hydrophobicity index
# Protein engineering, 5(5), 373-375.
# "Cid"
hphob_cid = {"A": 0.02, "R": -0.42, "N": -0.77, "D": -1.04, "C": 0.77,
      "Q": -1.1, "E": -1.14, "G": -0.8, "H": 0.26, "I": 1.81,
      "L": 1.14, "K": -0.41, "M": 1.0, "F": 1.35, "P": -0.09,
      "S": -0.97, "T": -0.77, "W": 1.71, "Y": 1.11, "V": 1.13}

# Cowan hydrophobicity indices at ph 3.4 and 7.5
# Peptide Research 3:75-80(1990).
# "Cowan3.4" "Conwan7.5"
_hphob_cowan = {3.4 : {"A": 0.42, "R": -1.56, "N": -1.03, "D": -0.51, "C": 0.84,
             "Q": -0.96, "E": -0.37, "G": 0.0, "H": -2.28, "I": 1.81,
             "L": 1.8, "K": -2.03, "M": 1.18, "F": 1.74, "P": 0.86,
             "S": -0.64, "T": -0.26, "W": 1.46, "Y": 0.51, "V": 1.34},
      7.5 : {"A": 0.35, "R": -1.5, "N": -0.99, "D": -2.15, "C": 0.76,
             "Q": -0.93, "E": -1.95, "G": 0.0, "H": -0.65, "I": 1.83,
             "L": 1.8, "K": -1.54, "M": 1.1, "F": 1.69, "P": 0.84,
             "S": -0.63, "T": -0.27, "W": 1.35, "Y": 0.39, "V": 1.32}
      }

hphob_cowan_3_4 = _hphob_cowan[3.4]
hphob_cowan_7_5 = _hphob_cowan[7.5]

# Eisenberg Normalized consensus hydrophobicity scale
# J. Mol. Biol. 179:125-142(1984)
# "Eisenberg"
hphob_eisenberg = {"A": 0.62, "R": -2.53, "N": -0.78, "D": -0.9, "C": 0.29,
      "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.4, "I": 1.38,
      "L": 1.06, "K": -1.5, "M": 0.64, "F": 1.19, "P": 0.12,
      "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08}

# Engelman Hydrophobic Transfer Free Energies
# Annual review of biophysics and biophysical chemistry, 15(1), 321-353.
# "Engelman"
hphob_engelman = {"A": -1.6, "R": 12.3, "N": 4.8, "D": 9.2, "C": -2,
      "Q": 4.1, "E": 8.2, "G": -1, "H": 3, "I": -3.1,
      "L": -2.8, "K": 8.8, "M": -3.4, "F": -3.7, "P": 0.2,
      "S": -0.6, "T": -1.2, "W": -1.9, "Y": 0.7, "V": -2.6}

# Fasman hydrophobicity index
# (1989). Prediction of protein structure and the principles of protein conformation. Springer.
# "Fasman"
hphob_fasman = {"A": -0.21, "R": 2.11, "N": 0.96, "D": 1.36, "C": -6.04,
      "Q": 1.52, "E": 2.3, "G": 0, "H": -1.23, "I": -4.81,
      "L": -4.68, "K": 3.88, "M": -3.66, "F": -4.65, "P": 0.75,
      "S": 1.74, "T": 0.78, "W": -3.32, "Y": -1.01, "V": -3.5}

# Fauchere Hydrophobicity scale
# Eur. J. Med. Chem. 18:369-375(1983).
# "Fauchere"
hphob_fauchere = {"A": 0.31, "R": -1.01, "N": -0.6, "D": -0.77, "C": 1.54,
      "Q": -0.22, "E": -0.64, "G": 0, "H": 0.13, "I": 1.8,
      "L": 1.7, "K": -0.99, "M": 1.23, "F": 1.79, "P": 0.72,
      "S": -0.04, "T": 0.26, "W": 2.25, "Y": 0.96, "V": 1.22}

# Goldsack & Chalifoux Free Energy of Mixing of the Hydrophobic Side Chains
# Journal of theoretical biology, 39(3), 645-651.
# "Goldsack"
hphob_goldsack = {"A": 0.75, "R": 0.75, "N": 0.69, "D": 0, "C": 1,
      "Q": 0.59, "E": 0, "G": 0, "H": 0, "I": 2.95,
      "L": 2.4, "K": 1.5, "M": 1.3, "F": 2.65, "P": 2.6,
      "S": 0, "T": 0.45, "W": 3, "Y": 2.85, "V": 1.7}

# Guy Hydrophobicity scale based on free energy of transfer (kcal/mole).
# Biophys J. 47:61-70(1985)
# "Guy"
hphob_guy = {"A": 0.1, "R": 1.91, "N": 0.48, "D": 0.78, "C": -1.42,
      "Q": 0.95, "E": 0.83, "G": 0.33, "H": -0.5, "I": -1.13,
      "L": -1.18, "K": 1.4, "M": -1.59, "F": -2.12, "P": 0.73,
      "S": 0.52, "T": 0.07, "W": -0.51, "Y": -0.21, "V": -1.27}

# Jones Hydrophobicity scale
# Journal of theoretical biology, 50(1), 167-183.
# "Jones"
hphob_jones = {"A": 0.87, "R": 0.85, "N": 0.09, "D": 0.66, "C": 1.52,
      "Q": 0, "E": 0.67, "G": 0.1, "H": 0.87, "I": 3.15,
      "L": 2.17, "K": 1.64, "M": 1.67, "F": 2.87, "P": 2.77,
      "S": 0.07, "T": 0.07, "W": 3.77, "Y": 2.67, "V": 1.87}

# Juretic Hydrophobicity scale
# Theoretical and computational chemistry, 5, 405-445.
# "Juretic"
hphob_juretic = {"A": 1.1, "R": -5.1, "N": -3.5, "D": -3.6, "C": 2.5,
      "Q": -3.68, "E": -3.2, "G": -0.64, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -4.11, "M": 1.9, "F": 2.8, "P": -1.9,
      "S": -0.5, "T": -0.7, "W": -0.46, "Y": -1.3, "V": 4.2}

# Kidera Hydrophobicity Factors
# Journal of Protein Chemistry, 4(1), 23-55.
# "Kidera"
hphob_kidera = {"A": -0.27, "R": 1.87, "N": 0.81, "D": 0.81, "C": -1.05,
      "Q": 1.1, "E": 1.17, "G": -0.16, "H": 0.28, "I": -0.77,
      "L": -1.1, "K": 1.7, "M": -0.73, "F": -1.43, "P": -0.75,
      "S": 0.42, "T": 0.63, "W": -1.57, "Y": -0.56, "V": -0.4}

# Miyazawa Hydrophobicity scale (contact energy derived from 3D data)
# Macromolecules 18:534-552(1985)
# "Miyazawa"
hphob_miyazawa = {"A": 5.33, "R": 4.18, "N": 3.71, "D": 3.59, "C": 7.93,
      "Q": 3.87, "E": 3.65, "G": 4.48, "H": 5.1, "I": 8.83,
      "L": 8.47, "K": 2.95, "M": 8.95, "F": 9.03, "P": 3.87,
      "S": 4.09, "T": 4.49, "W": 7.66, "Y": 5.89, "V": 7.63}

# Parker Hydrophilicity scale derived from HPLC peptide retention times
# Biochemistry 25:5425-5431(1986)
# "Parker"
hphob_parker = {"A": 2.1, "R": 4.2, "N": 7, "D": 10, "C": 1.4,
      "Q": 6, "E": 7.8, "G": 5.7, "H": 2.1, "I": -8,
      "L": -9.2, "K": 5.7, "M": -4.2, "F": -9.2, "P": 2.1,
      "S": 6.5, "T": 5.2, "W": -10, "Y": -1.9, "V": -3.7}

# Ponnuswamy Hydrophobic characteristics of folded proteins
# Progress in biophysics and molecular biology, 59(1), 57-103.
# "Ponnuswamy"
hphob_ponnuswamy = {"A": 0.85, "R": 0.2, "N": -0.48, "D": -1.1, "C": 2.1,
      "Q": -0.42, "E": -0.79, "G": 0, "H": 0.22, "I": 3.14,
      "L": 1.99, "K": -1.19, "M": 1.42, "F": 1.69, "P": -1.14,
      "S": -0.52, "T": -0.08, "W": 1.76, "Y": 1.37, "V": 2.53}

"""
Amino acid scale: Average surrounding hydrophobicity.
Author(s):
Manavalan P., Ponnuswamy P.K.
Reference: Nature 275:673-674(1978).
https://web.expasy.org/protscale/pscale/Hphob.Manavalan.html
"""
hphob_manavalan = {"A": 12.97, "R": 11.72, "N": 11.42, "D": 10.85, "C": 14.63,
      "Q": 11.76, "E": 11.89, "G": 12.43, "H": 12.16, "I": 15.67,
      "L": 14.9, "K": 11.36, "M": 14.39, "F": 14.0, "P": 11.37,
      "S": 11.23, "T": 11.69, "W": 13.93, "Y": 13.42, "V": 15.71}

# Rose Hydrophobicity scale
# Science 229:834-838(1985)
# "Rose"
hphob_rose = {"A": 0.74, "R": 0.64, "N": 0.63, "D": 0.62, "C": 0.91,
      "Q": 0.62, "E": 0.62, "G": 0.72, "H": 0.78, "I": 0.88,
      "L": 0.85, "K": 0.52, "M": 0.85, "F": 0.88, "P": 0.64,
      "S": 0.66, "T": 0.7, "W": 0.85, "Y": 0.76, "V": 0.86}

# Roseman Hydrophobicity scale
# J. Mol. Biol. 200:513-522(1988)
# "Roseman"
hphob_roseman = {"A": 0.39, "R": -3.95, "N": -1.91, "D": -3.81, "C": 0.25,
      "Q": -1.3, "E": -2.91, "G": 0, "H": -0.64, "I": 1.82,
      "L": 1.82, "K": -2.77, "M": 0.96, "F": 2.27, "P": 0.99,
      "S": -1.24, "T": -1, "W": 2.13, "Y": 1.47, "V": 1.3}

# Sweet Optimized Matchig Hydrophobicity (OMH)
# J. Mol. Biol. 171:479-488(1983).
# "Sweet
hphob_sweet = {"A": -0.4, "R": -0.59, "N": -0.92, "D": -1.31, "C": 0.17,
      "Q": -0.91, "E": -1.22, "G": -0.67, "H": -0.64, "I": 1.25,
      "L": 1.22, "K": -0.67, "M": 1.02, "F": 1.92, "P": -0.49,
      "S": -0.55, "T": -0.28, "W": 0.5, "Y": 1.67, "V": 0.91}

# Tanford Hydrophobicity scale
# J. Am. Chem. Soc. 84:4240-4274(1962)
# "Tanford"
hphob_tanford = {"A": 0.62, "R": -2.53, "N": -0.78, "D": -0.09, "C": 0.29,
      "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.4, "I": 1.38,
      "L": 1.53, "K": -1.5, "M": 0.64, "F": 1.19, "P": 0.12,
      "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.8}

# Wilson Hydrophobic constants derived from HPLC peptide retention times
# Biochem. J. 199:31-41(1981)
# "Wilson"
hphob_wilson = {"A": -0.3, "R": -1.1, "N": -0.2, "D": -1.4, "C": 6.3,
      "Q": -0.2, "E": 0, "G": 1.2, "H": -1.3, "I": 4.3,
      "L": 6.6, "K": -3.6, "M": 2.5, "F": 7.5, "P": 2.2,
      "S": -0.6, "T": -2.2, "W": 7.9, "Y": 7.1, "V": 5.9}

# Zimmerman Hydrophobicity scale
# Journal of theoretical biology, 21(2), 170-201.
# "Zimmerman"
hphob_zimmerman = {"A": 0.83, "R": 0.83, "N": 0.09, "D": 0.64, "C": 1.48,
      "Q": 0, "E": 0.65, "G": 0.1, "H": 1.1, "I": 3.07,
      "L": 2.52, "K": 1.6, "M": 1.4, "F": 2.75, "P": 2.7,
      "S": 0.14, "T": 0.54, "W": 0.31, "Y": 2.97, "V": 1.79}

"""
Amino acid scale: Proportion of residues 95% buried (in 12 proteins).
Author(s):
Chothia C.
Reference: J. Mol. Biol. 105:1-14(1976).
https://web.expasy.org/cgi-bin/protscale/protscale.pl
"""

hphob_chothia = {"A": 0.38, "R": 0.01, "N": 0.12, "D": 0.15, "C": 0.5,
      "Q": 0.07, "E": 0.18, "G": 0.36, "H": 0.17, "I": 0.6,
      "L": 0.45, "K": 0.03, "M": 0.4, "F": 0.5, "P": 0.18,
      "S": 0.22, "T": 0.23, "W": 0.27, "Y": 0.15, "V": 0.54}


""""
Amino acid scale: Free energy of transfer from inside to outside of a globular protein.
Author(s):
Janin J.
Reference: Nature 277:491-492(1979).
https://web.expasy.org/protscale/pscale/Hphob.Janin.html
"""

hphob_janin = {"A": 0.3, "R": -1.4, "N": -0.5, "D": -0.6, "C": 0.9,
      "Q": -0.7, "E": -0.7, "G": 0.3, "H": -0.1, "I": 0.7,
      "L": 0.5, "K": -1.8, "M": 0.4, "F": 0.5, "P": -0.3,
      "S": -0.1, "T": -0.2, "W": 0.3, "Y": -0.4, "V": 0.6}

"""
Amino acid scale: Hydration potential (kcal/mole) at 25�C.
Author(s):
Wolfenden R.V., Andersson L., Cullis P.M., Southgate C.C.F.
Reference: Biochemistry 20:849-855(1981).
https://web.expasy.org/protscale/pscale/Hphob.Wolfenden.html
"""

hphob_wolfenden = {"A": 1.940, "R": -19.920, "N": -9.680, "D": -10.950, "C": -1.240,
      "Q": -9.380, "E": -10.200, "G": 2.390, "H": -10.270, "I": 2.150,
      "L": 2.280, "K": -9.520, "M": -1.480, "F": -0.760, "P": 0.000,
      "S": -5.060, "T": -4.880, "W": -5.880, "Y": -6.110, "V": 1.990}


"""
mino acid scale: Hydrophobicity scale (pi-r).
Author(s):
Fauchere J.-L., Pliska V.E.
Reference: Eur. J. Med. Chem. 18:369-375(1983).
https://web.expasy.org/protscale/pscale/Hphob.Fauchere.html
"""

hphob_fauchere = {"A": 0.31, "R": -1.01, "N": -0.6, "D": -0.77, "C": 1.54,
      "Q": -0.22, "E": -0.64, "G": 0.0, "H": 0.13, "I": 1.8,
      "L": 1.7, "K": -0.99, "M": 1.23, "F": 1.79, "P": 0.72,
      "S": -0.04, "T": 0.26, "W": 2.25, "Y": 0.96, "V": 1.22}


"""
Amino acid scale: Antigenicity value X 10.
Author(s):
Welling G.W., Weijer W.J., Van der Zee R., Welling-Wester S.
Reference: FEBS Lett. 188:215-218(1985).
https://web.expasy.org/protscale/pscale/Hphob.Welling.html
"""

# Welling antigenicity scale
hphob_welling = {"A": 1.150, "R": 0.580, "N": -0.770, "D": 0.650, "C": -1.200,
      "Q": -0.110, "E": -0.710, "G": -1.840, "H": 3.120, "I": -2.920,
      "L": 0.750, "K": 2.060, "M": -3.850, "F": -1.410, "P": -0.530,
      "S": -0.260, "T": -0.450, "W": -1.140, "Y": 0.130, "V": -0.130}


hphob = hphob_kyte_doolittle

hydrophobicity_scales = {"KyteDoolitle": hphob_kyte_doolittle, "Aboderin": hphob_adoberin,
                "AbrahamLeo": hphob_abraham_leo, "Argos": hphob_agros,
                "BlackMould": hphob_black_mould, "BullBreese": hphob_bull_breese,
                "Casari": hphob_casari_sippl, "Cid": hphob_cid,
                "Cowan3.4": _hphob_cowan[3.4], "Cowan7.5": _hphob_cowan[7.5],
                "Eisenberg": hphob_eisenberg, "Engelman": hphob_engelman,
                "Fasman": hphob_fasman, "Fauchere": hphob_fauchere,
                "GoldSack": hphob_goldsack, "Guy": hphob_guy,
                "Jones": hphob_jones, "Juretic": hphob_juretic,
                "Kidera": hphob_kidera, "Miyazawa": hphob_miyazawa,
                "Parker": hphob_parker, "Ponnuswamy": hphob_ponnuswamy,
                "Rose": hphob_rose, "Roseman": hphob_roseman,
                "Sweet": hphob_sweet, "Tanford": hphob_tanford,
                "Wilson": hphob_wilson, "Zimmerman": hphob_zimmerman,
                "RaoArgos": hphob_rao_argos, "Manavalan": hphob_manavalan,
                  "Chothia": hphob_chothia, "Janin": hphob_janin,
                  "Wolfenden": hphob_wolfenden, "FaucherePiR": hphob_fauchere,
                  "Welling": hphob_welling, 'Default': hphob}


# Flexibility
# Normalized flexibility parameters (B-values), average
# Vihinen M., Torkkila E., Riikonen P. Proteins. 19(2):141-9(1994).
flexibility_vihinen = {"A": 0.984, "C": 0.906, "E": 1.094, "D": 1.068,
        "G": 1.031, "F": 0.915, "I": 0.927, "H": 0.950,
        "K": 1.102, "M": 0.952, "L": 0.935, "N": 1.048,
        "Q": 1.037, "P": 1.049, "S": 1.046, "R": 1.008,
        "T": 0.997, "W": 0.904, "V": 0.931, "Y": 0.929}

# Hydrophilicity
# 1 Hopp & Wood
# Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981).
hydrophilicity_hop_wood = {"A": -0.5, "R": 3.0, "N": 0.2, "D": 3.0, "C": -1.0,
      "Q": 0.2, "E": 3.0, "G": 0.0, "H": -0.5, "I": -1.8,
      "L": -1.8, "K": 3.0, "M": -1.3, "F": -2.5, "P": 0.0,
      "S": 0.3, "T": -0.4, "W": -3.4, "Y": -2.3, "V": -1.5}

# Surface accessibility
# Vergoten G & Theophanides T, Biomolecular Structure and Dynamics,
# pg.138 (1997).
# 1 Emini Surface fractional probability
surface_accessibility_vergoten = {"A": 0.815, "R": 1.475, "N": 1.296, "D": 1.283, "C": 0.394,
      "Q": 1.348, "E": 1.445, "G": 0.714, "H": 1.180, "I": 0.603,
      "L": 0.603, "K": 1.545, "M": 0.714, "F": 0.695, "P": 1.236,
      "S": 1.115, "T": 1.184, "W": 0.808, "Y": 1.089, "V": 0.606}

# 2 Janin Interior to surface transfer energy scale
surface_accessiblility_janin = {"A": 0.28, "R": -1.14, "N": -0.55, "D": -0.52, "C": 0.97,
      "Q": -0.69, "E": -1.01, "G": 0.43, "H": -0.31, "I": 0.60,
      "L": 0.60, "K": -1.62, "M": 0.43, "F": 0.46, "P": -0.42,
      "S": -0.19, "T": -0.32, "W": 0.29, "Y": -0.15, "V": 0.60}


# A two dimensional dictionary for calculating the instability index.
# Guruprasad K., Reddy B.V.B., Pandit M.W. Protein Engineering 4:155-161(1990).
# It is based on dipeptide values; therefore, the value for the dipeptide DG
# is DIWV['D']['G'].
DIWV = {"A": {"A": 1.0, "C": 44.94, "E": 1.0, "D": -7.49,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": -7.49,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 1.0, "P": 20.26, "S": 1.0, "R": 1.0,
              "T": 1.0, "W": 1.0, "V": 1.0, "Y": 1.0},
        "C": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 20.26,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 33.60,
              "K": 1.0, "M": 33.60, "L": 20.26, "N": 1.0,
              "Q": -6.54, "P": 20.26, "S": 1.0, "R": 1.0,
              "T": 33.60, "W": 24.68, "V": -6.54, "Y": 1.0},
        "E": {"A": 1.0, "C": 44.94, "E": 33.60, "D": 20.26,
              "G": 1.0, "F": 1.0, "I": 20.26, "H": -6.54,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 20.26, "P": 20.26, "S": 20.26, "R": 1.0,
              "T": 1.0, "W": -14.03, "V": 1.0, "Y": 1.0},
        "D": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": 1.0, "F": -6.54, "I": 1.0, "H": 1.0,
              "K": -7.49, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 1.0, "P": 1.0, "S": 20.26, "R": -6.54,
              "T": -14.03, "W": 1.0, "V": 1.0, "Y": 1.0},
        "G": {"A": -7.49, "C": 1.0, "E": -6.54, "D": 1.0,
              "G": 13.34, "F": 1.0, "I": -7.49, "H": 1.0,
              "K": -7.49, "M": 1.0, "L": 1.0, "N": -7.49,
              "Q": 1.0, "P": 1.0, "S": 1.0, "R": 1.0,
              "T": -7.49, "W": 13.34, "V": 1.0, "Y": -7.49},
        "F": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 13.34,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 1.0,
              "K": -14.03, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 1.0, "P": 20.26, "S": 1.0, "R": 1.0,
              "T": 1.0, "W": 1.0, "V": 1.0, "Y": 33.601},
        "I": {"A": 1.0, "C": 1.0, "E": 44.94, "D": 1.0,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 13.34,
              "K": -7.49, "M": 1.0, "L": 20.26, "N": 1.0,
              "Q": 1.0, "P": -1.88, "S": 1.0, "R": 1.0,
              "T": 1.0, "W": 1.0, "V": -7.49, "Y": 1.0},
        "H": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": -9.37, "F": -9.37, "I": 44.94, "H": 1.0,
              "K": 24.68, "M": 1.0, "L": 1.0, "N": 24.68,
              "Q": 1.0, "P": -1.88, "S": 1.0, "R": 1.0,
              "T": -6.54, "W": -1.88, "V": 1.0, "Y": 44.94},
        "K": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": -7.49, "F": 1.0, "I": -7.49, "H": 1.0,
              "K": 1.0, "M": 33.60, "L": -7.49, "N": 1.0,
              "Q": 24.64, "P": -6.54, "S": 1.0, "R": 33.60,
              "T": 1.0, "W": 1.0, "V": -7.49, "Y": 1.0},
        "M": {"A": 13.34, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 58.28,
              "K": 1.0, "M": -1.88, "L": 1.0, "N": 1.0,
              "Q": -6.54, "P": 44.94, "S": 44.94, "R": -6.54,
              "T": -1.88, "W": 1.0, "V": 1.0, "Y": 24.68},
        "L": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 1.0,
              "K": -7.49, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 33.60, "P": 20.26, "S": 1.0, "R": 20.26,
              "T": 1.0, "W": 24.68, "V": 1.0, "Y": 1.0},
        "N": {"A": 1.0, "C": -1.88, "E": 1.0, "D": 1.0,
              "G": -14.03, "F": -14.03, "I": 44.94, "H": 1.0,
              "K": 24.68, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": -6.54, "P": -1.88, "S": 1.0, "R": 1.0,
              "T": -7.49, "W": -9.37, "V": 1.0, "Y": 1.0},
        "Q": {"A": 1.0, "C": -6.54, "E": 20.26, "D": 20.26,
              "G": 1.0, "F": -6.54, "I": 1.0, "H": 1.0,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 20.26, "P": 20.26, "S": 44.94, "R": 1.0,
              "T": 1.0, "W": 1.0, "V": -6.54, "Y": -6.54},
        "P": {"A": 20.26, "C": -6.54, "E": 18.38, "D": -6.54,
              "G": 1.0, "F": 20.26, "I": 1.0, "H": 1.0,
              "K": 1.0, "M": -6.54, "L": 1.0, "N": 1.0,
              "Q": 20.26, "P": 20.26, "S": 20.26, "R": -6.54,
              "T": 1.0, "W": -1.88, "V": 20.26, "Y": 1.0},
        "S": {"A": 1.0, "C": 33.60, "E": 20.26, "D": 1.0,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 1.0,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 20.26, "P": 44.94, "S": 20.26, "R": 20.26,
              "T": 1.0, "W": 1.0, "V": 1.0, "Y": 1.0},
        "R": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": -7.49, "F": 1.0, "I": 1.0, "H": 20.26,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": 13.34,
              "Q": 20.26, "P": 20.26, "S": 44.94, "R": 58.28,
              "T": 1.0, "W": 58.28, "V": 1.0, "Y": -6.54},
        "T": {"A": 1.0, "C": 1.0, "E": 20.26, "D": 1.0,
              "G": -7.49, "F": 13.34, "I": 1.0, "H": 1.0,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": -14.03,
              "Q": -6.54, "P": 1.0, "S": 1.0, "R": 1.0,
              "T": 1.0, "W": -14.03, "V": 1.0, "Y": 1.0},
        "W": {"A": -14.03, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": -9.37, "F": 1.0, "I": 1.0, "H": 24.68,
              "K": 1.0, "M": 24.68, "L": 13.34, "N": 13.34,
              "Q": 1.0, "P": 1.0, "S": 1.0, "R": 1.0,
              "T": -14.03, "W": 1.0, "V": -7.49, "Y": 1.0},
        "V": {"A": 1.0, "C": 1.0, "E": 1.0, "D": -14.03,
              "G": -7.49, "F": 1.0, "I": 1.0, "H": 1.0,
              "K": -1.88, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 1.0, "P": 20.26, "S": 1.0, "R": 1.0,
              "T": -7.49, "W": 1.0, "V": 1.0, "Y": -6.54},
        "Y": {"A": 24.68, "C": 1.0, "E": -6.54, "D": 24.68,
              "G": -7.49, "F": 1.0, "I": 1.0, "H": 13.34,
              "K": 1.0, "M": 44.94, "L": 1.0, "N": 1.0,
              "Q": 1.0, "P": 13.34, "S": 1.0, "R": -15.91,
              "T": -7.49, "W": -9.37, "V": 1.0, "Y": 13.34},
        }

# Surface accessibility scales
surface_accessibility_scales = {"Emini": surface_accessibility_vergoten, "Janin": surface_accessiblility_janin}

# Hydrophilicity scales  
hydrophilicity_scales = {"HoppWood": hydrophilicity_hop_wood}

# Flexibility scales
flexibility_scales = {"Vihinen": flexibility_vihinen}


"""
https://www.peptideweb.com/images/pdf/pKa-and-pI-values-of-amino-acids.pdf
"""

pk_nterminal = { 
      "A": 9.69, "R": 9.04, "N": 8.80, "D": 9.82, "C": 10.78,
      "Q": 9.13, "E": 9.67, "G": 9.60, "H": 9.17, "I": 9.60, 
      "L": 9.60, "K": 8.95, "M": 9.21, "F": 9.13, "P": 10.60, 
      "S": 9.15, "T": 9.10, "W": 9.44, "Y": 9.11, "V": 9.62
}
pk_cterminal = {
      "A": 2.34, "R": 2.17, "N": 2.02, "D": 2.09, "C": 1.71,
      "Q": 2.19, "E": 2.17, "G": 2.34, "H": 1.82, "I": 2.36,
      "L": 2.36, "K": 2.18, "M": 2.28, "F": 1.83, "P": 1.99,
      "S": 2.21, "T": 2.09, "W": 2.43, "Y": 2.20, "V": 2.32
}

pk_sidechain = {
      "A": 0.0, "R": 12.48, "N": 0.0, "D": 3.86, "C": 8.33,
      "Q": 0.0, "E": 4.25, "G": 0.0, "H": 6.00, "I": 0.0,
      "L": 0.0, "K": 10.79, "M": 0.0, "F": 0.0, "P": 0.0,
      "S": 0.0, "T": 0.0, "W": 0.0, "Y": 10.07, "V": 0.0
}

charged_aas = ("K", "R", "H", "D", "E", "C", "Y")

"""
Amino acid scale: Conformational parameter for alpha helix.
Author(s):
Deleage G., Roux B.
Reference: Protein Engineering 1:289-294(1987).
https://web.expasy.org/protscale/pscale/alpha-helixRoux.html
"""

deleage_roux_alpha_helix = {
    "A": 1.489, "R": 1.224, "N": 0.772, "D": 0.924, "C": 0.966,
    "Q": 1.164, "E": 1.504, "G": 0.510, "H": 1.003, "I": 1.003,
    "L": 1.236, "K": 1.172, "M": 1.363, "F": 1.195, "P": 0.492,
    "S": 0.739, "T": 0.785, "W": 1.090, "Y": 0.787, "V": 0.990
}

"""
Amino acid scale: Conformational parameter for beta-sheet.
Author(s):
Deleage G., Roux B.
Reference: Protein Engineering 1:289-294(1987).
https://web.expasy.org/protscale/pscale/beta-sheetRoux.html
"""

deleage_roux_beta_sheet = {
    "A": 0.709, "R": 0.920, "N": 0.604, "D": 0.541, "C": 1.191,
    "Q": 0.840, "E": 0.567, "G": 0.657, "H": 0.863, "I": 1.799,
    "L": 1.261, "K": 0.721, "M": 1.210, "F": 1.393, "P": 0.354,
    "S": 0.928, "T": 1.221, "W": 1.306, "Y": 1.266, "V": 1.965
}

"""
Amino acid scale: Conformational parameter for beta-turn.
Author(s):
Deleage G., Roux B.
Reference: Protein Engineering 1:289-294(1987).
https://web.expasy.org/protscale/pscale/beta-turnRoux.html
"""

deleage_roux_beta_turn = {
    "A": 0.788, "R": 0.912, "N": 1.572, "D": 1.197, "C": 0.965,
    "Q": 0.997, "E": 1.149, "G": 1.860, "H": 0.970, "I": 0.240,
    "L": 0.670, "K": 1.302, "M": 0.436, "F": 0.624, "P": 1.415,
    "S": 1.316, "T": 0.739, "W": 0.546, "Y": 0.795, "V": 0.387
}

"""
Amino acid scale: Conformational parameter for coil.
Author(s):
Deleage G., Roux B.
Reference: Protein Engineering 1:289-294(1987).
https://web.expasy.org/protscale/pscale/CoilRoux.html
"""

deleage_roux_coil = {
    "A": 0.824, "R": 0.893, "N": 1.167, "D": 1.197, "C": 0.953,
    "Q": 0.947, "E": 0.761, "G": 1.251, "H": 1.068, "I": 0.886,
    "L": 0.810, "K": 0.897, "M": 0.810, "F": 0.797, "P": 1.540,
    "S": 1.130, "T": 1.148, "W": 0.941, "Y": 1.109, "V": 0.772
}


"""
Amino acid scale: Normalized frequency for alpha helix.
Author(s):
Levitt M.
Reference: Biochemistry 17:4277-4285(1978).
https://web.expasy.org/protscale/pscale/alpha-helixLevitt.html
"""

levitt_alpha_helix = {
    "A": 1.290, "R": 0.960, "N": 0.900, "D": 1.040, "C": 1.110,
    "Q": 1.270, "E": 1.440, "G": 0.560, "H": 1.220, "I": 0.970,
    "L": 1.300, "K": 1.230, "M": 1.470, "F": 1.070, "P": 0.520,
    "S": 0.820, "T": 0.820, "W": 0.990, "Y": 0.720, "V": 0.910
}

"""
Amino acid scale: Normalized frequency for beta-sheet.
Author(s):
Levitt M.
Reference: Biochemistry 17:4277-4285(1978).
https://web.expasy.org/protscale/pscale/beta-sheetLevitt.html
"""

levitt_beta_sheet = {
    "A": 0.900, "R": 0.990, "N": 0.760, "D": 0.720, "C": 0.740,
    "Q": 0.800, "E": 0.750, "G": 0.920, "H": 1.080, "I": 1.450,
    "L": 1.020, "K": 0.770, "M": 0.970, "F": 1.320, "P": 0.640,
    "S": 0.950, "T": 1.210, "W": 1.140, "Y": 1.250, "V": 1.490
}

"""
Amino acid scale: Normalized frequency for beta-turn.
Author(s):
Levitt M.
Reference: Biochemistry 17:4277-4285(1978).
https://web.expasy.org/protscale/pscale/beta-turnLevitt.html
"""

levitt_beta_turn = {
    "A": 0.770, "R": 0.880, "N": 1.280, "D": 1.410, "C": 0.810,
    "Q": 0.980, "E": 0.990, "G": 1.640, "H": 0.680, "I": 0.510,
    "L": 0.580, "K": 0.960, "M": 0.410, "F": 0.590, "P": 1.910,
    "S": 1.320, "T": 1.040, "W": 0.760, "Y": 1.050, "V": 0.470
}

"""
Amino acid scale: Conformational parameter for alpha helix (computed from 29 proteins).
Author(s):
Chou P.Y., Fasman G.D.
Reference: Adv. Enzym. 47:45-148(1978).
https://web.expasy.org/protscale/pscale/alpha-helixFasman.html
"""

chou_fasman_alpha_helix = {
    "A": 1.420, "R": 0.980, "N": 0.670, "D": 1.010, "C": 0.700,
    "Q": 1.110, "E": 1.510, "G": 0.570, "H": 1.000, "I": 1.080,
    "L": 1.210, "K": 1.160, "M": 1.450, "F": 1.130, "P": 0.570,
    "S": 0.770, "T": 0.830, "W": 1.080, "Y": 0.690, "V": 1.060
}

"""
Amino acid scale: Conformational parameter for beta-sheet (computed from 29 proteins).
Author(s):
Chou P.Y., Fasman G.D.
Reference: Adv. Enzym. 47:45-148(1978).
  https://web.expasy.org/protscale/pscale/beta-sheetFasman.html
"""

chou_fasman_beta_sheet = {
    "A": 0.830, "R": 0.930, "N": 0.890, "D": 0.540, "C": 1.190,
    "Q": 1.100, "E": 0.370, "G": 0.750, "H": 0.870, "I": 1.600,
    "L": 1.300, "K": 0.740, "M": 1.050, "F": 1.380, "P": 0.550,
    "S": 0.750, "T": 1.190, "W": 1.370, "Y": 1.470, "V": 1.700
}

"""
Amino acid scale: Conformational parameter for beta-turn (computed from 29 proteins).
Author(s):
Chou P.Y., Fasman G.D.
Reference: Adv. Enzym. 47:45-148(1978).
https://web.expasy.org/protscale/pscale/beta-turnFasman.html
"""

chou_fasman_beta_turn = {
    "A": 0.660, "R": 0.950, "N": 1.560, "D": 1.460, "C": 1.190,
    "Q": 0.980, "E": 0.740, "G": 1.560, "H": 0.950, "I": 0.470,
    "L": 0.590, "K": 1.010, "M": 0.600, "F": 0.600, "P": 1.520,
    "S": 1.430, "T": 0.960, "W": 0.960, "Y": 1.140, "V": 0.500
}

secondary_structure_scales = {
      "DeleageRouxAlphaHelix": deleage_roux_alpha_helix,
      "DeleageRouxBetaSheet": deleage_roux_beta_sheet,
      "DeleageRouxBetaTurn": deleage_roux_beta_turn,
      "DeleageRouxCoil": deleage_roux_coil,
      "LevittAlphaHelix": levitt_alpha_helix,
      "LevittBetaSheet": levitt_beta_sheet,
      "LevittBetaTurn": levitt_beta_turn,
      "ChouFasmanAlphaHelix": chou_fasman_alpha_helix,
      "ChouFasmanBetaSheet": chou_fasman_beta_sheet,
      "ChouFasmanBetaTurn": chou_fasman_beta_turn
      }

secondary_structure_scales_by_name = {
    'DeleageRoux' : {
        'alpha_helix': deleage_roux_alpha_helix,
        'beta_sheet': deleage_roux_beta_sheet,
        'beta_turn': deleage_roux_beta_turn,
        'coil': deleage_roux_coil
    },
    'Levitt' : {
        'alpha_helix': levitt_alpha_helix,
        'beta_sheet': levitt_beta_sheet,
        'beta_turn': levitt_beta_turn
    },
    'ChouFasman' : {
        'alpha_helix': chou_fasman_alpha_helix,
        'beta_sheet': chou_fasman_beta_sheet,
        'beta_turn': chou_fasman_beta_turn
    }
}


"""
Amino acid scale: Conformational preference for parallel beta strand.
Author(s):
Lifson S., Sander C.
Reference: Nature 282:109-111(1979).
https://web.expasy.org/protscale/pscale/Parallelbeta-strand.html
"""

beta_strand_parallel = {
    "A": 1.000, "R": 0.680, "N": 0.540, "D": 0.500, "C": 0.910,
    "Q": 0.280, "E": 0.590, "G": 0.790, "H": 0.380, "I": 2.600,
    "L": 1.420, "K": 0.590, "M": 1.490, "F": 1.300, "P": 0.350,
    "S": 0.700, "T": 0.590, "W": 0.890, "Y": 1.080, "V": 2.630
}

"""
Amino acid scale: Conformational preference for antiparallel beta strand.
Author(s):
Lifson S., Sander C.
Reference: Nature 282:109-111(1979).
https://web.expasy.org/protscale/pscale/Antiparallelbeta-strand.html
"""

beta_strand_antiparallel = {
    "A": 0.900, "R": 1.020, "N": 0.620, "D": 0.470, "C": 1.240,
    "Q": 1.180, "E": 0.620, "G": 0.560, "H": 1.120, "I": 1.540,
    "L": 1.260, "K": 0.740, "M": 1.090, "F": 1.230, "P": 0.420,
    "S": 0.870, "T": 1.300, "W": 1.750, "Y": 1.680, "V": 1.530
}

"""
Amino acid scale: Conformational preference for total beta strand (antiparallel+parallel).
Author(s):
Lifson S., Sander C.
Reference: Nature 282:109-111(1979).
https://web.expasy.org/protscale/pscale/Totalbeta-strand.html
"""

beta_strand_total = {
    "A": 0.920, "R": 0.930, "N": 0.600, "D": 0.480, "C": 1.160,
    "Q": 0.950, "E": 0.610, "G": 0.610, "H": 0.930, "I": 1.810,
    "L": 1.300, "K": 0.700, "M": 1.190, "F": 1.250, "P": 0.400,
    "S": 0.820, "T": 1.120, "W": 1.540, "Y": 1.530, "V": 1.810
}

"""
Amino acid scale: Atomic weight ratio of hetero elements in end group to C in side chain.
Author(s):
Grantham R.
Reference: Science 185:862-864(1974).
https://web.expasy.org/protscale/pscale/Ratioside.html
"""

ratioside = {
    "A": 0.000, "R": 0.650, "N": 1.330, "D": 1.380, "C": 2.750,
    "Q": 0.890, "E": 0.920, "G": 0.740, "H": 0.580, "I": 0.000,
    "L": 0.000, "K": 0.330, "M": 0.000, "F": 0.000, "P": 0.390,
    "S": 1.420, "T": 0.710, "W": 0.130, "Y": 0.200, "V": 0.000
}


"""
Amino acid scale: Polarity (p).
Author(s):
Grantham R.
Reference: Science 185:862-864(1974).
https://web.expasy.org/protscale/pscale/PolarityGrantham.html
"""

polarity_grantham = {
      "A": 8.100, "R": 10.500, "N": 11.600, "D": 13.000, "C": 5.500,
      "Q": 10.500, "E": 12.300, "G": 9.000, "H": 10.400, "I": 5.200,
      "L": 4.900, "K": 11.300, "M": 5.700, "F": 5.200, "P": 8.000,
      "S": 9.200, "T": 8.600, "W": 5.400, "Y": 6.200, "V": 5.900
      }

"""
Amino acid scale: Polarity.
Author(s):
Zimmerman J.M., Eliezer N., Simha R.
Reference: J. Theor. Biol. 21:170-201(1968).
https://web.expasy.org/protscale/pscale/PolarityZimmerman.html
"""

polarity_zimmerman = {
    "A": 0.000, "R": 52.000, "N": 3.380, "D": 49.700, "C": 1.480,
    "Q": 3.530, "E": 49.900, "G": 0.000, "H": 51.600, "I": 0.130,
    "L": 0.130, "K": 49.500, "M": 1.430, "F": 0.350, "P": 1.580,
    "S": 1.670, "T": 1.660, "W": 2.100, "Y": 1.610, "V": 0.130
}

polarity_scales = {
    "Grantham": polarity_grantham,
    "Zimmerman": polarity_zimmerman
}

"""
Amino acid scale: Relative mutability of amino acids (Ala=100).
Author(s):
Dayhoff M.O., Schwartz R.M., Orcutt B.C.
Reference: In "Atlas of Protein Sequence and Structure", Vol.5, Suppl.3 (1978).
https://web.expasy.org/protscale/pscale/Relativemutability.html
"""

mutability = {
      "A": 100.000, "R": 65.000, "N": 134.000, "D": 106.000, "C": 20.000,
      "Q": 93.000, "E": 102.000, "G": 49.000, "H": 66.000, "I": 96.000,
      "L": 40.000, "K": 56.000, "M": 94.000, "F": 41.000, "P": 56.000,
      "S": 120.000, "T": 97.000, "W": 18.000, "Y": 41.000, "V": 74.000
}

"""
Amino acid scale: Number of codon(s) coding for each amino acid in universal genetic code.
Author(s):
-
Reference: Most textbooks.
https://web.expasy.org/protscale/pscale/Numbercodons.html
"""

codons = {
    "A": 4.0, "R": 6.0, "N": 2.0, "D": 2.0, "C": 1.0,
    "Q": 2.0, "E": 2.0, "G": 4.0, "H": 2.0, "I": 3.0,
    "L": 6.0, "K": 2.0, "M": 1.0, "F": 2.0, "P": 4.0,
    "S": 6.0, "T": 4.0, "W": 1.0, "Y": 2.0, "V": 4.0
}

"""
Amino acid scale: Refractivity.
Author(s):
Jones. D.D.
Reference: J. Theor. Biol. 50:167-184(1975).
https://web.expasy.org/protscale/pscale/Refractivity.html
"""

refractivity = {
    "A": 4.340, "R": 26.660, "N": 13.280, "D": 12.000, "C": 35.770,
    "Q": 17.560, "E": 17.260, "G": 0.000, "H": 21.810, "I": 19.060,
    "L": 18.780, "K": 21.290, "M": 21.640, "F": 29.400, "P": 10.930,
    "S": 6.350, "T": 11.010, "W": 42.530, "Y": 31.530, "V": 13.920
}

"""
Amino acid scale: Bulkiness.
Author(s):
Zimmerman J.M., Eliezer N., Simha R.
Reference: J. Theor. Biol. 21:170-201(1968).

https://web.expasy.org/protscale/pscale/Bulkiness.html 
"""

bulkiness = {
    "A": 11.500, "R": 14.280, "N": 12.820, "D": 11.680, "C": 13.460,
    "Q": 14.450, "E": 13.570, "G": 3.400, "H": 13.690, "I": 21.400,
    "L": 21.400, "K": 15.710, "M": 16.250, "F": 19.800, "P": 17.430,
    "S": 9.470, "T": 15.770, "W": 21.670, "Y": 18.030, "V": 21.570
}

"""
Amino acid scale: Recognition factors.
Author(s):
Fraga S.
Reference: Can. J. Chem. 60:2606-2610(1982).
https://web.expasy.org/protscale/pscale/Recognitionfactors.html
"""

recognition_factors = {
    "A": 78.000, "R": 95.000, "N": 94.000, "D": 81.000, "C": 89.000,
    "Q": 87.000, "E": 78.000, "G": 84.000, "H": 84.000, "I": 88.000,
    "L": 85.000, "K": 87.000, "M": 80.000, "F": 81.000, "P": 91.000,
    "S": 107.000, "T": 93.000, "W": 104.000, "Y": 84.000, "V": 89.000
}

"""
Amino acid scale: Overall amino acid composition (%).
Author(s):
McCaldon P., Argos P.
Reference: Proteins: Structure, Function and Genetics 4:99-122(1988).
https://web.expasy.org/protscale/pscale/A.A.composition.html
"""

aa_composition_mccaldron = {
      "A": 8.300, "R": 5.700, "N": 4.400, "D": 5.300, "C": 1.700,
      "Q": 4.000, "E": 6.200, "G": 7.200, "H": 2.200, "I": 5.200,
      "L": 9.000, "K": 5.700, "M": 2.400, "F": 3.900, "P": 5.100,
      "S": 6.900, "T": 5.800, "W": 1.300, "Y": 3.200, "V": 6.600
      }

"""
Amino acid scale: Amino acid composition (%) in the UniProtKB/Swiss-Prot data bank.
Author(s):
Bairoch A.
Reference: Release notes for UniProtKB/Swiss-Prot release 2013_04 - April 2013.
https://web.expasy.org/protscale/pscale/A.A.Swiss-Prot.html 
"""

aa_composition_swissprot = {
    "A": 8.25, "R": 5.53, "N": 4.06, "D": 5.45, "C": 1.37,
    "Q": 3.93, "E": 6.75, "G": 7.07, "H": 2.27, "I": 5.96,
    "L": 9.66, "K": 5.84, "M": 2.42, "F": 3.86, "P": 4.70,
    "S": 6.56, "T": 5.34, "W": 1.08, "Y": 2.92, "V": 6.87
}

aa_composition_scales = {
    "McCaldron": aa_composition_mccaldron,
    "Swiss-Prot": aa_composition_swissprot
}

"""
Amino acid scale: Transmembrane tendency
Author(s):
Zhao, G., London E.
Reference: Protein Sci. 15:1987-2001(2006).
https://web.expasy.org/protscale/pscale/Transmembranetendency.html
"""

transmembrane_tendency = {
    "A": 0.380, "R": -2.570, "N": -1.620, "D": -3.270, "C": -0.300,
    "Q": -1.840, "E": -2.900, "G": -0.190, "H": -1.440, "I": 1.970,
    "L": 1.820, "K": -3.460, "M": 1.400, "F": 1.980, "P": -1.440,
    "S": -0.530, "T": -0.320, "W": 1.530, "Y": 0.490, "V": 1.460
}

"""
Amino acid scale: Molar fraction (%) of 3220 accessible residues.
Author(s):
Janin J.
Reference: Nature 277:491-492(1979).
https://web.expasy.org/protscale/pscale/accessibleresidues.html
"""

accessible_residues = {
    "A": 6.600, "R": 4.500, "N": 6.700, "D": 7.700, "C": 0.900,
    "Q": 5.200, "E": 5.700, "G": 6.700, "H": 2.500, "I": 2.800,
    "L": 4.800, "K": 10.300, "M": 1.000, "F": 2.400, "P": 4.800,
    "S": 9.400, "T": 7.000, "W": 1.400, "Y": 5.100, "V": 4.500
}

"""
Amino acid scale: Average area buried on transfer from standard state to folded protein.
Author(s):
Rose G.D., Geselowitz A.R., Lesser G.J., Lee R.H., Zehfus M.H.
Reference: Science 229:834-838(1985).
https://web.expasy.org/protscale/pscale/Averageburied.html
"""

average_buried_area = {
    "A": 86.600, "R": 162.200, "N": 103.300, "D": 97.800, "C": 132.300,
    "Q": 119.200, "E": 113.900, "G": 62.900, "H": 155.800, "I": 158.000,
    "L": 164.100, "K": 115.500, "M": 172.900, "F": 194.100, "P": 92.900,
    "S": 85.600, "T": 106.500, "W": 224.600, "Y": 177.700, "V": 141.000
}

"""
Amino acid scale: Molecular weight of each amino acid.
Author(s):
-
Reference: Most textbooks.
https://web.expasy.org/protscale/pscale/Molecularweight.html
"""

molecular_weights = {
    "A": 89.000, "R": 174.000, "N": 132.000, "D": 133.000, "C": 121.000,
    "Q": 146.000, "E": 147.000, "G": 75.000, "H": 155.000, "I": 131.000,
    "L": 131.000, "K": 146.000, "M": 149.000, "F": 165.000, "P": 115.000,
    "S": 105.000, "T": 119.000, "W": 204.000, "Y": 181.000, "V": 117.000
}

"""
Amino acid scale: Retention coefficient in HPLC, pH 2.1.
Author(s):
Meek J.L.
Reference: Proc. Natl. Acad. Sci. USA 77:1632-1636(1980).
https://web.expasy.org/protscale/pscale/HPLC2.1.html
"""

hplc_meek_2_1 = {
    "A": -0.100, "R": -4.500, "N": -1.600, "D": -2.800, "C": -2.200,
    "Q": -2.500, "E": -7.500, "G": -0.500, "H": 0.800, "I": 11.800,
    "L": 10.000, "K": -3.200, "M": 7.100, "F": 13.900, "P": 8.000,
    "S": -3.700, "T": 1.500, "W": 18.100, "Y": 8.200, "V": 3.300
}

"""
Amino acid scale: Retention coefficient in HFBA.
Author(s):
Browne C.A., Bennett H.P.J., Solomon S.
Reference: Anal. Biochem. 124:201-208(1982).
https://web.expasy.org/protscale/pscale/HPLCHFBA.html

"""

hplc_browne = {
    "A": 3.900, "R": 3.200, "N": -2.800, "D": -2.800, "C": -14.300,
    "Q": 1.800, "E": -7.500, "G": -2.300, "H": 2.000, "I": 11.000,
    "L": 15.000, "K": -2.500, "M": 4.100, "F": 14.700, "P": 5.600,
    "S": -3.500, "T": 1.100, "W": 17.800, "Y": 3.800, "V": 2.100
}

"""
Amino acid scale: Retention coefficient in HPLC, pH 7.4.
Author(s):
Meek J.L.
Reference: Proc. Natl. Acad. Sci. USA 77:1632-1636(1980).
https://web.expasy.org/protscale/pscale/HPLC7.4.html
"""

hplc_meek_7_4 = {
    "A": 0.500, "R": 0.800, "N": 0.800, "D": -8.200, "C": -6.800,
    "Q": -4.800, "E": -16.900, "G": 0.000, "H": -3.500, "I": 13.900,
    "L": 8.800, "K": 0.100, "M": 4.800, "F": 13.200, "P": 6.100,
    "S": 1.200, "T": 2.700, "W": 14.900, "Y": 6.100, "V": 2.700
}

"""
Amino acid scale: Retention coefficient in TFA.
Author(s):
Browne C.A., Bennett H.P.J., Solomon S.
Reference: Anal. Biochem. 124:201-208(1982).
https://web.expasy.org/protscale/pscale/HPLCTFA.html
"""

hplc_browne_tfa = {
    "A": 7.300, "R": -3.600, "N": -5.700, "D": -2.900, "C": -9.200,
    "Q": -0.300, "E": -7.100, "G": -1.200, "H": -2.100, "I": 6.600,
    "L": 20.000, "K": -3.700, "M": 5.600, "F": 19.200, "P": 5.100,
    "S": -4.100, "T": 0.800, "W": 16.300, "Y": 5.900, "V": 3.500
}

hplc_scales = {
    "Meek_2.1": hplc_meek_2_1,
    "Browne": hplc_browne,
    "Meek_7.4": hplc_meek_7_4,
    "Browne_TFA": hplc_browne_tfa
}


all_property_scales = {
    'hphob_kyte_doolittle': hphob_kyte_doolittle,
    "hphob_adoberin": hphob_adoberin,
    "hphob_abraham_leo": hphob_abraham_leo,
    "hphob_agros": hphob_agros,
    "hphob_rao_argos": hphob_rao_argos,
    "hphob_black_mould": hphob_black_mould,
    "hphob_bull_breese": hphob_bull_breese,
    "hphob_casari_sippl": hphob_casari_sippl,
    "hphob_cid": hphob_cid,
    "hphob_cowan_3_4": hphob_cowan_3_4,
    "hphob_cowan_7_5": hphob_cowan_7_5,
    "hphob_eisenberg": hphob_eisenberg,
    "hphob_engelman": hphob_engelman,
    "hphob_fasman": hphob_fasman,
    "hphob_fauchere": hphob_fauchere,
    "hphob_goldsack": hphob_goldsack,
    "hphob_guy": hphob_guy,
    "hphob_jones": hphob_jones,
    "hphob_juretic": hphob_juretic,
    "hphob_kidera": hphob_kidera,
    "hphob_miyazawa": hphob_miyazawa,
    "hphob_parker": hphob_parker,
    "hphob_ponnuswamy": hphob_ponnuswamy,
    "hphob_manavalan": hphob_manavalan,
    "hphob_rose": hphob_rose,
    "hphob_roseman": hphob_roseman,
    "hphob_sweet": hphob_sweet,
    "hphob_tanford": hphob_tanford,
    "hphob_wilson": hphob_wilson,
    "hphob_zimmerman": hphob_zimmerman,
    "hphob_chothia": hphob_chothia,
    "hphob_janin": hphob_janin,
    "hphob_wolfenden": hphob_wolfenden,
    "hphob_fauchere": hphob_fauchere,
    "hphob_welling": hphob_welling,
    "flexibility_vihinen": flexibility_vihinen,
    "hydrophilicity_hop_wood": hydrophilicity_hop_wood,
    "surface_accessibility_vergoten": surface_accessibility_vergoten,
    "surface_accessiblility_janin": surface_accessiblility_janin,
    "pk_nterminal": pk_nterminal,
    "pk_cterminal": pk_cterminal,
    "pk_sidechain": pk_sidechain,
    "deleage_roux_alpha_helix": deleage_roux_alpha_helix,
    "deleage_roux_beta_sheet": deleage_roux_beta_sheet,
    "deleage_roux_beta_turn": deleage_roux_beta_turn,
    "deleage_roux_coil": deleage_roux_coil,
    "levitt_alpha_helix": levitt_alpha_helix,
    "levitt_beta_sheet": levitt_beta_sheet,
    "levitt_beta_turn": levitt_beta_turn,
    "chou_fasman_alpha_helix": chou_fasman_alpha_helix,
    "chou_fasman_beta_sheet": chou_fasman_beta_sheet,
    "chou_fasman_beta_turn": chou_fasman_beta_turn,
    "beta_strand_parallel": beta_strand_parallel,
    "beta_strand_antiparallel": beta_strand_antiparallel,
    "beta_strand_total": beta_strand_total,
    "ratioside": ratioside,
    "polarity_grantham": polarity_grantham,
    "polarity_zimmerman": polarity_zimmerman,
    "mutability": mutability,
    "codons": codons,
    "refractivity": refractivity,
    "bulkiness": bulkiness,
    "recognition_factors": recognition_factors,
    "aa_composition_mccaldron": aa_composition_mccaldron,
    "aa_composition_swissprot": aa_composition_swissprot,
    "transmembrane_tendency": transmembrane_tendency,
    "accessible_residues": accessible_residues,
    "average_buried_area": average_buried_area,
    "molecular_weights": molecular_weights,
    "hplc_meek_2_1": hplc_meek_2_1,
    "hplc_browne": hplc_browne,
    "hplc_meek_7_4": hplc_meek_7_4,
    "hplc_browne_tfa": hplc_browne_tfa,
}
