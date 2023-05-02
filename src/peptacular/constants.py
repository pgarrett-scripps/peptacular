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
