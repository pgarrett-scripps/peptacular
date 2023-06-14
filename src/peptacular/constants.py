amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

# Trypsin: https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#exceptions
TRYPTIC_COMPLEX_REGEXES = (
    [('([KR])([^P])', 1),
     ('([W])([K])([P])', 2),
     ('([M])([R])([P])', 2)],

    [('([CD])(K)([D])', 2),
     ('([C])([K])([HY])', 2),
     ('([C])([R])([K])', 2),
     ('([R])([R])([HR])', 2)])

TRYPTIC_SIMPLE_REGEXES = (
    [('([KR])([^P])', 1)],

    [])

THROMBIN_REGEXES = (
    [('([AFGILTVM])([AFGILTVWA])([P])([R])([^DE])([^DE])', 4),
     ('([G])([R])([G])', 2)],

    [])

THERMOLYSIN_REGEXES = (
    [('([^DE])([AFILMV])([^P])', 1)],

    [])

PROTEINASEK_REGEXES = (
    [('([AEFILTVWY])', 1)],

    [])

LYSC_REGEXES = ([('([K])', 1)],

                [])

LYSN_REGEXES = ([('([K])', 0)],

                [])

ENZYMES_OPTIONS = \
    {'Trypsin': TRYPTIC_SIMPLE_REGEXES,
     'Trypsin (expasy)': TRYPTIC_COMPLEX_REGEXES,
     'Thrombin': THROMBIN_REGEXES,
     'Thermolysin': THERMOLYSIN_REGEXES,
     'Protinase K': PROTEINASEK_REGEXES,
     'Lys C': LYSC_REGEXES,
     'Lys N': LYSN_REGEXES,
     }


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