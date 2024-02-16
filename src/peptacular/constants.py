PROTON_MASS = 1.00727646688
NEUTRON_MASS = 1.00866491597

MONOISOTOPIC_ATOMIC_MASSES = {
    'HYDROGEN': 1.007825035,
    'DEUTERIUM': 2.014101779,
    'LITHIUM': 7.016003,
    'BORON': 11.0093055,
    'CARBON': 12.0,
    'CARBON13': 13.00335483,
    'NITROGEN': 14.003074,
    'NITROGEN15': 15.00010897,
    'OXYGEN': 15.99491463,
    'OXYGEN18': 17.9991603,
    'FLUORINE': 18.99840322,
    'SODIUM': 22.9897677,
    'MAGNESIUM': 23.9850423,
    'ALUMINIUM': 26.9815386,
    'PHOSPHOROUS': 30.973762,
    'SULFUR': 31.9720707,
    'CHLORINE': 34.96885272,
    'POTASSIUM': 38.9637074,
    'CALCIUM': 39.9625906,
    'CHROMIUM': 51.9405098,
    'MANGANESE': 54.9380471,
    'IRON': 55.9349393,
    'NICKEL': 57.9353462,
    'COBALT': 58.9331976,
    'COPPER': 62.9295989,
    'ZINC': 63.9291448,
    'ARSENIC': 74.9215942,
    'BROMINE': 78.9183361,
    'SELENIUM': 79.9165196,
    'MOLYBDENUM': 97.9054073,
    'RUTHENIUM': 101.9043485,
    'PALLADIUM': 105.903478,
    'SILVER': 106.905092,
    'CADMIUM': 113.903357,
    'IODINE': 126.904473,
    'PLATINUM': 194.964766,
    'GOLD': 196.966543,
    'MERCURY': 201.970617,
}

AVERAGE_ATOMIC_MASSES = {
    'HYDROGEN': 1.00794,
    'DEUTERIUM': 2.014101779,
    'LITHIUM': 6.941,
    'BORON': 10.811,
    'CARBON': 12.0107,
    'CARBON13': 13.00335483,
    'NITROGEN': 14.0067,
    'NITROGEN15': 15.00010897,
    'OXYGEN': 15.9994,
    'OXYGEN18': 17.9991603,
    'FLUORINE': 18.9984032,
    'SODIUM': 22.98977,
    'MAGNESIUM': 24.305,
    'ALUMINIUM': 26.9815386,
    'PHOSPHOROUS': 30.973761,
    'SULFUR': 32.065,
    'CHLORINE': 35.453,
    'POTASSIUM': 39.0983,
    'CALCIUM': 40.078,
    'CHROMIUM': 51.9961,
    'MANGANESE': 54.938045,
    'IRON': 55.845,
    'NICKEL': 58.6934,
    'COBALT': 58.933195,
    'COPPER': 63.546,
    'ZINC': 65.409,
    'ARSENIC': 74.9215942,
    'BROMINE': 79.904,
    'SELENIUM': 78.96,
    'MOLYBDENUM': 95.94,
    'RUTHENIUM': 101.07,
    'PALLADIUM': 106.42,
    'SILVER': 107.8682,
    'CADMIUM': 112.411,
    'IODINE': 126.90447,
    'PLATINUM': 195.084,
    'GOLD': 196.96655,
    'MERCURY': 200.59,
}

AA_COMPOSITIONS = {
    "G": {"CARBON": 2, "HYDROGEN": 3, "NITROGEN": 1, "OXYGEN": 1},  # Glycine
    "A": {"CARBON": 3, "HYDROGEN": 5, "NITROGEN": 1, "OXYGEN": 1},  # Alanine
    "S": {"CARBON": 3, "HYDROGEN": 5, "NITROGEN": 1, "OXYGEN": 2},  # Serine
    "P": {"CARBON": 5, "HYDROGEN": 7, "NITROGEN": 1, "OXYGEN": 1},  # Proline
    "V": {"CARBON": 5, "HYDROGEN": 9, "NITROGEN": 1, "OXYGEN": 1},  # Valine
    "T": {"CARBON": 4, "HYDROGEN": 7, "NITROGEN": 1, "OXYGEN": 2},  # Threonine
    "C": {"CARBON": 3, "HYDROGEN": 5, "NITROGEN": 1, "OXYGEN": 1, "SULFUR": 1},  # Cysteine
    "I": {"CARBON": 6, "HYDROGEN": 11, "NITROGEN": 1, "OXYGEN": 1},  # Isoleucine
    "L": {"CARBON": 6, "HYDROGEN": 11, "NITROGEN": 1, "OXYGEN": 1},  # Leucine
    "N": {"CARBON": 4, "HYDROGEN": 6, "NITROGEN": 2, "OXYGEN": 2},  # Asparagine
    "D": {"CARBON": 4, "HYDROGEN": 5, "NITROGEN": 1, "OXYGEN": 3},  # Aspartic acid
    "Q": {"CARBON": 5, "HYDROGEN": 8, "NITROGEN": 2, "OXYGEN": 3},  # Glutamine
    "K": {"CARBON": 6, "HYDROGEN": 12, "NITROGEN": 2},  # Lysine
    "E": {"CARBON": 5, "HYDROGEN": 7, "NITROGEN": 1, "OXYGEN": 3},  # Glutamic acid
    "M": {"CARBON": 5, "HYDROGEN": 9, "NITROGEN": 1, "OXYGEN": 1, "SULFUR": 1},  # Methionine
    "H": {"CARBON": 6, "HYDROGEN": 7, "NITROGEN": 3, "OXYGEN": 1},  # Histidine
    "F": {"CARBON": 9, "HYDROGEN": 9, "NITROGEN": 1, "OXYGEN": 1},  # Phenylalanine
    "R": {"CARBON": 6, "HYDROGEN": 12, "NITROGEN": 4, "OXYGEN": 1},  # Arginine
    "Y": {"CARBON": 9, "HYDROGEN": 9, "NITROGEN": 1, "OXYGEN": 2},  # Tyrosine
    "W": {"CARBON": 11, "HYDROGEN": 10, "NITROGEN": 2},  # Tryptophan
    "U": {"CARBON": 3, "HYDROGEN": 5, "NITROGEN": 1, "OXYGEN": 1, "SELENIUM": 1},  # Selenocysteine
    "O": {"CARBON": 12, "HYDROGEN": 19, "NITROGEN": 3, "OXYGEN": 2},  # Pyrrolysine
}

MONOISOTOPIC_AA_MASSES, AVERAGE_AA_MASSES = {}, {}
for aa in AA_COMPOSITIONS:
    MONOISOTOPIC_AA_MASSES[aa] = sum([MONOISOTOPIC_ATOMIC_MASSES[k] * v for k, v in AA_COMPOSITIONS[aa].items()])
    AVERAGE_AA_MASSES[aa] = sum([AVERAGE_ATOMIC_MASSES[k] * v for k, v in AA_COMPOSITIONS[aa].items()])

AMINO_ACIDS = set(AA_COMPOSITIONS.keys())

PROTEASES = {'arg-c': 'R',
             'asp-n': '\\w(?=D)',
             'bnps-skatole': 'W',
             'caspase 1': '(?<=[FWYL]\\w[HAT])D(?=[^PEDQKR])',
             'caspase 2': '(?<=DVA)D(?=[^PEDQKR])',
             'caspase 3': '(?<=DMQ)D(?=[^PEDQKR])',
             'caspase 4': '(?<=LEV)D(?=[^PEDQKR])',
             'caspase 5': '(?<=[LW]EH)D',
             'caspase 6': '(?<=VE[HI])D(?=[^PEDQKR])',
             'caspase 7': '(?<=DEV)D(?=[^PEDQKR])',
             'caspase 8': '(?<=[IL]ET)D(?=[^PEDQKR])',
             'caspase 9': '(?<=LEH)D',
             'caspase 10': '(?<=IEA)D',
             'chymotrypsin high specificity': '([FY](?=[^P]))|(W(?=[^MP]))',
             'chymotrypsin low specificity': '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
             'chymotrypsin': '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
             'clostripain': 'R',
             'cnbr': 'M',
             'enterokinase': '(?<=[DE]{3})K',
             'factor xa': '(?<=[AFGILTVM][DE]G)R',
             'formic acid': 'D',
             'glutamyl endopeptidase': 'E',
             'glu-c': 'E',
             'granzyme b': '(?<=IEP)D',
             'hydroxylamine': 'N(?=G)',
             'iodosobenzoic acid': 'W',
             'lys-c': 'K',
             'lys-n': '\\w(?=K)',
             'ntcb': '\\w(?=C)',
             'pepsin ph1.3': '((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\\w[^P]))',
             'pepsin ph2.0': '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\\w[^P]))',
             'proline endopeptidase': '(?<=[HKR])P(?=[^P])',
             'proteinase k': '[AEFILTVWY]',
             'staphylococcal peptidase i': '(?<=[^E])E',
             'thermolysin': '[^DE](?=[AFILMV])',
             'thrombin': '((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
             'trypsin_full': '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
             'trypsin_exception': '((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))',
             'trypsin': '([KR](?=[^P]))',
             'trypsin/P': '([KR])',
             'non-specific': '()',
             'no-cleave': '_'}

MONOISOTOPIC_ION_ADJUSTMENTS = {
    'a': -MONOISOTOPIC_ATOMIC_MASSES['CARBON'] - MONOISOTOPIC_ATOMIC_MASSES['OXYGEN'],
    'b': 0,
    'c': MONOISOTOPIC_ATOMIC_MASSES['HYDROGEN'] * 3 + MONOISOTOPIC_ATOMIC_MASSES['NITROGEN'],
    'x': MONOISOTOPIC_ATOMIC_MASSES['CARBON'] + MONOISOTOPIC_ATOMIC_MASSES['OXYGEN'] * 2,
    'y': MONOISOTOPIC_ATOMIC_MASSES['HYDROGEN'] * 2 + MONOISOTOPIC_ATOMIC_MASSES['OXYGEN'],
    'z': MONOISOTOPIC_ATOMIC_MASSES['OXYGEN'] - MONOISOTOPIC_ATOMIC_MASSES['NITROGEN'] - MONOISOTOPIC_ATOMIC_MASSES[
        'HYDROGEN'],
    'I': - MONOISOTOPIC_ATOMIC_MASSES['CARBON'] - MONOISOTOPIC_ATOMIC_MASSES['OXYGEN'],
}

AVERAGE_ION_ADJUSTMENTS = {
    'a': AVERAGE_ATOMIC_MASSES['CARBON'] - AVERAGE_ATOMIC_MASSES['OXYGEN'],
    'b': 0,
    'c': AVERAGE_ATOMIC_MASSES['HYDROGEN'] * 3 + AVERAGE_ATOMIC_MASSES['NITROGEN'],
    'x': AVERAGE_ATOMIC_MASSES['CARBON'] + AVERAGE_ATOMIC_MASSES['OXYGEN'] * 2,
    'y': AVERAGE_ATOMIC_MASSES['HYDROGEN'] * 2 + AVERAGE_ATOMIC_MASSES['OXYGEN'],
    'z': AVERAGE_ATOMIC_MASSES['OXYGEN'] - AVERAGE_ATOMIC_MASSES['NITROGEN'] - AVERAGE_ATOMIC_MASSES['HYDROGEN'],
    'I': - AVERAGE_ATOMIC_MASSES['CARBON'] - AVERAGE_ATOMIC_MASSES['OXYGEN'],
}


FORWARD_IONS = {'a', 'b', 'c'}
BACKWARD_IONS = {'x', 'y', 'z'}
VALID_ION_TYPES = FORWARD_IONS.union(BACKWARD_IONS).union({'I'})

ISOTOPES = {"CARBON": [0.9893, 0.0107],
            "HYDROGEN": [0.999885, 0.000115],
            "OXYGEN": [0.99757, 0.00038, 0.00205],
            "NITROGEN": [0.99636, 0.00364],
            "SULFUR": [0.9499, 0.0075, 0.0425, 0.0001]}
