from peptacular import constants
from peptacular.protein import digest_protein, combine_site_regexes

import unittest

"""
Results from: https://web.expasy.org/cgi-bin/peptide_cutter/peptidecutter.pl
"""

PROTEIN = 'MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKTRLDSSNLEKIYREVQLMKLLNHPHIIKLYQVMETKDMLYIVTE' \
          'FAKNGEMFDYLTSNGHLSENEARKKFWQILSAVEYCHDHHIVHRDLKTENLLLDGNMDIKLADFGFGNFYKSGEPLSTWCGSPPYAAPEVFEGKEYEGPQLDIW' \
          'SLGVVLYVLVCGSLPFDGPNLPTLRQRVLEGRFRIPFFMSQDCESLIRRMLVVDPARRITIAQIRQHRWMRAEPCLPGPACPAFSAHSYTSNLGDYDEQALGIM' \
          'QTLGVDRQRTVESLQNSSYNHFAAIYYLLLERLKEYRNAQCARPGPARQPRPRSSDLSGLEVPQEGLSTDPFRPALLCPQPQTLVQSVLQAEMDCELQSSLQWP' \
          'LFFPVDASCSGVFRPRPVSPSSLLDTAISEEARQGPGLEEEQDTQESLPSSTGRRHTLAEVSTRLSPLTAPCIVVSPSTTASPAEGTSSDSCLTFSASKSPAGL' \
          'SGTPATQGLLGACSPVRLASPFLGSQSATPVLQAQGGLGGAVLLPVSFQEGRRASDTSLTQGLKAFRQQLRKTTRTKGFLGLNKIKGLARQVCQVPASRASRGG' \
          'LSPFHAPAQSPGLHGGAAGSREGWSLLEEVLEQQRLLQLQHHPAAAPGCSQAPQPAPAPFVIAPCDGPGAAPLPSTLLTSGLPLLPPPLLQTGASPVASAAQLL' \
          'DTHLHIGTGPTALPAVPPPRLARLAPGCEPLGLLQGDCEMEDLMPCSLGTFVLVQ'


class TestProtein(unittest.TestCase):

    def test_trypsin_sites(self):
        cleavage_sites = combine_site_regexes(PROTEIN, constants.TRYPTIC_COMPLEX_REGEXES[0], constants.TRYPTIC_COMPLEX_REGEXES[1])

        sites = [23, 31, 35, 42, 45, 47, 50, 56, 60, 62, 70, 73, 79, 88, 96, 107, 127, 128, 129, 148, 151, 164, 175,
                 198, 233, 235, 240, 242, 256, 257, 265, 266, 273, 276, 279, 319, 321, 344, 346, 349, 360, 365, 449,
                 470, 480, 515, 537, 572, 573, 584, 587, 591, 592, 595, 597, 604, 606, 610, 619, 622, 645, 659, 748,
                 751]

        self.assertEqual(cleavage_sites, sites)

    def test_thermolysin_sites(self):
        cleavage_sites = combine_site_regexes(PROTEIN, constants.THERMOLYSIN_REGEXES[0], constants.THERMOLYSIN_REGEXES[1])
        sites = [1, 2, 3, 8, 11, 21, 23, 25, 32, 37, 38, 39, 40, 42, 43, 47, 52, 53, 54, 56, 57, 62, 67, 70, 76, 77,
                 79, 80, 85, 86, 88, 91, 92, 98, 100, 101, 105, 111, 114, 120, 129, 132, 133, 135, 136, 144, 145, 154,
                 155, 156, 160, 164, 165, 169, 172, 179, 189, 194, 204, 209, 211, 212, 213, 215, 216, 217, 223, 231,
                 235, 236, 240, 244, 245, 246, 253, 254, 257, 258, 259, 260, 263, 266, 268, 269, 271, 277, 279, 287,
                 290, 291, 293, 300, 307, 308, 310, 311, 314, 316, 322, 325, 333, 334, 335, 336, 339, 340, 341, 344,
                 350, 353, 358, 371, 378, 383, 386, 387, 388, 395, 396, 399, 400, 402, 412, 416, 417, 420, 427, 428,
                 433, 438, 439, 442, 443, 453, 473, 474, 480, 483, 488, 489, 490, 496, 499, 508, 510, 512, 517, 519,
                 524, 528, 529, 531, 535, 537, 538, 541, 542, 547, 550, 551, 553, 557, 560, 561, 562, 565, 567, 573,
                 578, 582, 584, 585, 589, 598, 599, 601, 604, 607, 608, 611, 616, 619, 624, 627, 631, 636, 640, 641,
                 649, 650, 654, 659, 660, 662, 667, 668, 683, 684, 685, 693, 700, 701, 707, 712, 713, 717, 720, 721,
                 723, 724, 726, 727, 731, 733, 739, 742, 748, 749, 751, 758, 760, 761, 775, 778, 779, 780, 781]
        self.assertEqual(cleavage_sites, sites)

    def test_thrombin_sites(self):
        cleavage_sites = combine_site_regexes(PROTEIN, constants.THROMBIN_REGEXES[0], constants.THROMBIN_REGEXES[1])
        sites = []
        self.assertEqual(cleavage_sites, sites)

    def test_proteinaseK_sites(self):
        cleavage_sites = combine_site_regexes(PROTEIN, constants.PROTEINASEK_REGEXES[0], constants.PROTEINASEK_REGEXES[1])
        sites = [2, 3, 6, 7, 9, 12, 22, 24, 26, 27, 29, 30, 32, 33, 38, 39, 40, 41, 43, 44, 48, 49, 51, 53, 54, 55, 57,
                 58, 61, 63, 68, 69, 71, 72, 74, 75, 77, 80, 81, 86, 87, 89, 90, 92, 94, 95, 99, 100, 101, 102, 103,
                 104, 105, 106, 110, 112, 114, 115, 116, 121, 123, 125, 126, 130, 131, 133, 134, 136, 137, 138, 139,
                 145, 146, 150, 152, 153, 155, 156, 157, 163, 165, 166, 168, 170, 173, 174, 178, 180, 182, 183, 189,
                 190, 191, 193, 194, 195, 196, 199, 200, 201, 205, 207, 208, 210, 212, 213, 214, 215, 216, 217, 218,
                 222, 224, 229, 231, 232, 236, 237, 238, 241, 243, 245, 246, 252, 254, 255, 259, 260, 261, 264, 267,
                 268, 269, 270, 272, 277, 280, 281, 284, 288, 291, 292, 294, 297, 298, 301, 304, 306, 308, 309, 311,
                 314, 315, 317, 322, 323, 324, 326, 331, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 345, 347,
                 348, 351, 354, 359, 369, 372, 373, 374, 377, 379, 381, 384, 387, 388, 389, 395, 396, 397, 400, 401,
                 403, 404, 408, 409, 413, 415, 417, 418, 419, 421, 423, 428, 429, 434, 439, 440, 442, 443, 444, 446,
                 447, 448, 454, 455, 456, 457, 460, 462, 464, 468, 473, 474, 475, 476, 477, 479, 481, 484, 485, 486,
                 489, 490, 491, 495, 496, 497, 500, 501, 503, 509, 510, 511, 513, 518, 520, 523, 525, 526, 529, 530,
                 532, 536, 538, 539, 542, 543, 548, 549, 551, 552, 554, 558, 561, 562, 563, 564, 566, 568, 570, 574,
                 577, 579, 580, 583, 585, 586, 590, 593, 594, 596, 599, 600, 602, 605, 608, 609, 612, 615, 617, 620,
                 625, 628, 630, 632, 637, 641, 642, 646, 648, 650, 651, 652, 653, 654, 655, 656, 660, 661, 663, 668,
                 669, 670, 676, 680, 682, 684, 685, 686, 687, 694, 695, 697, 700, 701, 702, 703, 706, 708, 709, 713,
                 714, 716, 718, 721, 722, 724, 725, 727, 728, 730, 732, 734, 736, 739, 740, 741, 743, 744, 749, 750,
                 752, 753, 757, 759, 761, 762, 767, 769, 771, 776, 778, 779, 780, 781, 782]

        self.assertEqual(cleavage_sites, sites)

    def test_lysC_sites(self):
        cleavage_sites = combine_site_regexes(PROTEIN, constants.LYSC_REGEXES[0], constants.LYSC_REGEXES[1])
        sites = [20, 35, 42, 50, 56, 60, 70, 79, 88, 96, 107, 128, 129, 151, 164, 175, 198, 346, 515, 584, 592, 597,
                 604, 606]

        self.assertEqual(cleavage_sites, sites)

    def test_lysN_sites(self):
        cleavage_sites = combine_site_regexes(PROTEIN, constants.LYSN_REGEXES[0], constants.LYSN_REGEXES[1])
        sites = [19, 34, 41, 49, 55, 59, 69, 78, 87, 95, 106, 127, 128, 150, 163, 174, 197, 345, 514, 583, 591, 596,
                 603, 605]

        self.assertEqual(cleavage_sites, sites)

    def test_digest_protein(self):

        peptides = set(digest_protein(protein_sequence='TIDERTIDEKTIDE',
                                      enzyme_regexes=constants.TRYPTIC_SIMPLE_REGEXES,
                                      missed_cleavages=2,
                                      min_len=0,
                                      max_len=100,
                                      non_enzymatic=False,
                                      semi_enzymatic=False))
        self.assertEqual(peptides,
                         {('TIDER', 0), ('TIDERTIDEK', 1), ('TIDERTIDEKTIDE', 2), ('TIDEK', 0), ('TIDEKTIDE', 1),
                          ('TIDE', 0)})

        peptides = set(digest_protein(protein_sequence='TIDERTIDEKTIDE',
                                      enzyme_regexes=constants.TRYPTIC_SIMPLE_REGEXES,
                                      missed_cleavages=1,
                                      min_len=0,
                                      max_len=100,
                                      non_enzymatic=False,
                                      semi_enzymatic=False))
        self.assertEqual(peptides, {('TIDER', 0), ('TIDERTIDEK', 1), ('TIDEK', 0), ('TIDEKTIDE', 1), ('TIDE', 0)})

        peptides = set(digest_protein(protein_sequence='KTIDERTIDEKTIDE',
                                      enzyme_regexes=constants.TRYPTIC_SIMPLE_REGEXES,
                                      missed_cleavages=1,
                                      min_len=0,
                                      max_len=100,
                                      non_enzymatic=False,
                                      semi_enzymatic=False))
        self.assertEqual(peptides,
                         {('K', 0), ('KTIDER', 1), ('TIDER', 0), ('TIDERTIDEK', 1), ('TIDEK', 0), ('TIDEKTIDE', 1),
                          ('TIDE', 0)})

        peptides = set(digest_protein(protein_sequence='TIDERTIDEKTIDEK',
                                      enzyme_regexes=constants.TRYPTIC_SIMPLE_REGEXES,
                                      missed_cleavages=1,
                                      min_len=0,
                                      max_len=100,
                                      non_enzymatic=False,
                                      semi_enzymatic=False))
        self.assertEqual(peptides, {('TIDER', 0), ('TIDERTIDEK', 1), ('TIDEK', 0), ('TIDEKTIDEK', 1), ('TIDEK', 0)})

        peptides = set(digest_protein(protein_sequence='TIDERTIDEKKTIDE',
                                      enzyme_regexes=constants.TRYPTIC_SIMPLE_REGEXES,
                                      missed_cleavages=1,
                                      min_len=0,
                                      max_len=100,
                                      non_enzymatic=False,
                                      semi_enzymatic=False))
        self.assertEqual(peptides,
                         {('TIDER', 0), ('TIDERTIDEK', 1), ('TIDEK', 0), ('TIDEKK', 1), ('K', 0), ('KTIDE', 1),
                          ('TIDE', 0)})

        peptides = set(digest_protein(protein_sequence='TIDERTIDEKTIDE',
                                      enzyme_regexes=constants.TRYPTIC_SIMPLE_REGEXES,
                                      missed_cleavages=0,
                                      min_len=0,
                                      max_len=100,
                                      non_enzymatic=False,
                                      semi_enzymatic=False))
        self.assertEqual(peptides, {('TIDER', 0), ('TIDEK', 0), ('TIDE', 0)})

        peptides = set(digest_protein(protein_sequence='TIDERTIDEKTIDE',
                                      enzyme_regexes=constants.TRYPTIC_SIMPLE_REGEXES,
                                      missed_cleavages=10,
                                      min_len=0,
                                      max_len=100,
                                      non_enzymatic=False,
                                      semi_enzymatic=False))
        self.assertEqual(peptides,
                         {('TIDER', 0), ('TIDERTIDEK', 1), ('TIDERTIDEKTIDE', 2), ('TIDEK', 0), ('TIDEKTIDE', 1),
                          ('TIDE', 0)})


if __name__ == "__main__":
    unittest.main()
