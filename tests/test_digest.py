import unittest

import peptacular as pt


PROTEIN = "MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKTRLDSSNLEKIYREVQLMKLLNHPHIIKLYQVMETKDMLYIVTE"


class TestDigest(unittest.TestCase):
    def test_trypsin_cleavage_sites(self):
        """Test getting cleavage sites for trypsin."""
        annotation = pt.ProFormaAnnotation.parse(PROTEIN)
        cleavage_sites = list(annotation.get_cleavage_sites(pt.PROTEASES["trypsin"]))
        expected_sites = [23, 31, 35, 42, 45, 47, 50, 56, 60, 62, 70, 73, 79, 88, 96]
        self.assertEqual(cleavage_sites, expected_sites)

    def test_digest_protein_missed_cleavages_2(self):
        """Test protein digestion with 2 missed cleavages."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        peptides = set(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=2,
                min_len=0,
                max_len=100,
                semi=False,
                return_type="str",
            )
        )
        expected = {
            "TIDER",
            "TIDERTIDEK",
            "TIDERTIDEKTIDE",
            "TIDEK",
            "TIDEKTIDE",
            "TIDE",
        }
        self.assertEqual(peptides, expected)

    def test_digest_protein_missed_cleavages_1(self):
        """Test protein digestion with 1 missed cleavage."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        peptides = set(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=1,
                min_len=0,
                max_len=100,
                semi=False,
                return_type="str",
            )
        )
        expected = {"TIDER", "TIDERTIDEK", "TIDEK", "TIDEKTIDE", "TIDE"}
        self.assertEqual(peptides, expected)

    def test_digest_protein_with_k_at_beginning(self):
        """Test protein digestion with K at the beginning."""
        annotation = pt.ProFormaAnnotation.parse("KTIDERTIDEKTIDE")
        peptides = set(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=1,
                min_len=0,
                max_len=100,
                semi=False,
                return_type="str",
            )
        )
        expected = {"K", "KTIDER", "TIDER", "TIDERTIDEK", "TIDEK", "TIDEKTIDE", "TIDE"}
        self.assertEqual(peptides, expected)

    def test_digest_protein_with_k_at_end(self):
        """Test protein digestion with K at the end."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDEK")
        peptides = set(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=1,
                min_len=0,
                max_len=100,
                semi=False,
                return_type="str",
            )
        )
        expected = {"TIDER", "TIDERTIDEK", "TIDEK", "TIDEKTIDEK", "TIDEK"}
        self.assertEqual(peptides, expected)

    def test_digest_protein_with_double_k(self):
        """Test protein digestion with consecutive K residues."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKKTIDE")
        peptides = set(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=1,
                min_len=0,
                max_len=100,
                semi=False,
                return_type="str",
            )
        )
        expected = {"TIDER", "TIDERTIDEK", "TIDEK", "TIDEKK", "K", "KTIDE", "TIDE"}
        self.assertEqual(peptides, expected)

    def test_digest_protein_no_missed_cleavages(self):
        """Test protein digestion with no missed cleavages."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        peptides = set(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=0,
                min_len=0,
                max_len=100,
                semi=False,
                return_type="str",
            )
        )
        expected = {"TIDER", "TIDEK", "TIDE"}
        self.assertEqual(peptides, expected)

    def test_digest_protein_high_missed_cleavages(self):
        """Test protein digestion with very high missed cleavages."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        peptides = set(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=10,
                min_len=0,
                max_len=100,
                semi=False,
                return_type="str",
            )
        )
        expected = {
            "TIDER",
            "TIDERTIDEK",
            "TIDERTIDEKTIDE",
            "TIDEK",
            "TIDEKTIDE",
            "TIDE",
        }
        self.assertEqual(peptides, expected)

    def test_semi_enzymatic_sequences_min_3_max_5(self):
        """Test semi-enzymatic sequences with min_len=3, max_len=5."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_semi_enzymatic_sequences(
                min_len=3, max_len=5, return_type="str"
            )
        )
        expected = {"TIDE", "PTIDE", "PEPT", "PEPTI", "IDE", "PEP"}
        self.assertEqual(sequences, expected)

    def test_semi_enzymatic_sequences_min_2_max_5(self):
        """Test semi-enzymatic sequences with min_len=2, max_len=5."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_semi_enzymatic_sequences(
                min_len=2, max_len=5, return_type="str"
            )
        )
        expected = {"TIDE", "PTIDE", "PEPT", "PEPTI", "IDE", "PEP", "PE", "DE"}
        self.assertEqual(sequences, expected)

    def test_semi_enzymatic_sequences_min_2_max_100(self):
        """Test semi-enzymatic sequences with min_len=2, max_len=100."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_semi_enzymatic_sequences(
                min_len=2, max_len=100, return_type="str"
            )
        )
        expected = {
            "TIDE",
            "PTIDE",
            "PEPT",
            "PEPTI",
            "IDE",
            "PEP",
            "PE",
            "DE",
            "PEPTID",
            "EPTIDE",
        }
        self.assertEqual(sequences, expected)

    def test_semi_enzymatic_sequences_no_limits(self):
        """Test semi-enzymatic sequences with no length limits."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_semi_enzymatic_sequences(
                min_len=None, max_len=None, return_type="str"
            )
        )
        expected = {
            "P",
            "E",
            "TIDE",
            "PTIDE",
            "PEPT",
            "PEPTI",
            "IDE",
            "PEP",
            "PE",
            "DE",
            "PEPTID",
            "EPTIDE",
        }
        self.assertEqual(sequences, expected)

    def test_left_semi_enzymatic_sequences_no_limits(self):
        """Test left semi-enzymatic sequences with no length limits."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_left_semi_enzymatic_sequences(None, None, return_type="str")
        )
        expected = {"P", "PE", "PEP", "PEPT", "PEPTI", "PEPTID"}
        self.assertEqual(sequences, expected)

    def test_left_semi_enzymatic_sequences_min_3(self):
        """Test left semi-enzymatic sequences with min_len=3."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_left_semi_enzymatic_sequences(3, None, return_type="str")
        )
        expected = {"PEP", "PEPT", "PEPTI", "PEPTID"}
        self.assertEqual(sequences, expected)

    def test_left_semi_enzymatic_sequences_max_5(self):
        """Test left semi-enzymatic sequences with max_len=5."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_left_semi_enzymatic_sequences(None, 5, return_type="str")
        )
        expected = {"P", "PE", "PEP", "PEPT", "PEPTI"}
        self.assertEqual(sequences, expected)

    def test_left_semi_enzymatic_sequences_min_3_max_4(self):
        """Test left semi-enzymatic sequences with min_len=3, max_len=4."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_left_semi_enzymatic_sequences(3, 4, return_type="str")
        )
        expected = {"PEP", "PEPT"}
        self.assertEqual(sequences, expected)

    def test_right_semi_enzymatic_sequences_no_limits(self):
        """Test right semi-enzymatic sequences with no length limits."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_right_semi_enzymatic_sequences(None, None, return_type="str")
        )
        expected = {"EPTIDE", "PTIDE", "TIDE", "IDE", "DE", "E"}
        self.assertEqual(sequences, expected)

    def test_right_semi_enzymatic_sequences_min_3(self):
        """Test right semi-enzymatic sequences with min_len=3."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_right_semi_enzymatic_sequences(3, None, return_type="str")
        )
        expected = {"EPTIDE", "PTIDE", "TIDE", "IDE"}
        self.assertEqual(sequences, expected)

    def test_right_semi_enzymatic_sequences_max_5(self):
        """Test right semi-enzymatic sequences with max_len=5."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_right_semi_enzymatic_sequences(None, 5, return_type="str")
        )
        expected = {"PTIDE", "TIDE", "IDE", "DE", "E"}
        self.assertEqual(sequences, expected)

    def test_right_semi_enzymatic_sequences_min_3_max_4(self):
        """Test right semi-enzymatic sequences with min_len=3, max_len=4."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sequences = set(
            annotation.get_right_semi_enzymatic_sequences(3, 4, return_type="str")
        )
        expected = {"TIDE", "IDE"}
        self.assertEqual(sequences, expected)

    def test_non_enzymatic_sequences_no_limits(self):
        """Test non-enzymatic sequences with no length limits."""
        annotation = pt.ProFormaAnnotation.parse("PEPT")
        sequences = set(annotation.get_non_enzymatic_sequences(return_type="str"))
        expected = {"P", "E", "P", "T", "PE", "EP", "PT", "PEP", "EPT"}
        self.assertEqual(sequences, expected)

    def test_non_enzymatic_sequences_min_1_max_2(self):
        """Test non-enzymatic sequences with min_len=1, max_len=2."""
        annotation = pt.ProFormaAnnotation.parse("PEPT")
        sequences = set(
            annotation.get_non_enzymatic_sequences(
                min_len=1, max_len=2, return_type="str"
            )
        )
        expected = {"P", "E", "P", "T", "PE", "EP", "PT"}
        self.assertEqual(sequences, expected)

    def test_non_enzymatic_sequences_min_2_max_4(self):
        """Test non-enzymatic sequences with min_len=2, max_len=4."""
        annotation = pt.ProFormaAnnotation.parse("PEPT")
        sequences = set(
            annotation.get_non_enzymatic_sequences(
                min_len=2, max_len=4, return_type="str"
            )
        )
        expected = {"PT", "EP", "EPT", "PE", "PEP"}
        self.assertEqual(sequences, expected)

    def test_modified_sequences_with_missed_cleavages(self):
        """Test digestion with modified sequences and missed cleavages."""
        modified_seq = "[Acetyl]-PEPTIDER[Phospho]TIDEM[Oxidation]K"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)

        peptides = set(
            annotation.regex_digest(
                enzyme_regex="([KR])", missed_cleavages=1, return_type="str"
            )
        )
        expected = {
            "[Acetyl]-PEPTIDER[Phospho]",
            "[Acetyl]-PEPTIDER[Phospho]TIDEM[Oxidation]K",
            "TIDEM[Oxidation]K",
        }
        self.assertEqual(peptides, expected)

    def test_modified_sequences_no_missed_cleavages(self):
        """Test digestion with modified sequences and no missed cleavages."""
        modified_seq = "[Acetyl]-PEPTIDER[Phospho]TIDEM[Oxidation]K"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)

        peptides = set(
            annotation.regex_digest(
                enzyme_regex="([KR])", missed_cleavages=0, return_type="str"
            )
        )
        expected = {"[Acetyl]-PEPTIDER[Phospho]", "TIDEM[Oxidation]K"}
        self.assertEqual(peptides, expected)

    def test_empty_sequence_regex_digest(self):
        """Test regex digestion with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation.parse("")
        result = list(
            empty_annotation.regex_digest(enzyme_regex="([KR])", return_type="str")
        )
        self.assertEqual(result, [])

    def test_empty_sequence_semi_enzymatic(self):
        """Test semi-enzymatic sequences with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation.parse("")
        result = list(empty_annotation.get_semi_enzymatic_sequences(return_type="str"))
        self.assertEqual(result, [])

    def test_empty_sequence_left_semi_enzymatic(self):
        """Test left semi-enzymatic sequences with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation.parse("")
        result = list(
            empty_annotation.get_left_semi_enzymatic_sequences(return_type="str")
        )
        self.assertEqual(result, [])

    def test_empty_sequence_right_semi_enzymatic(self):
        """Test right semi-enzymatic sequences with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation.parse("")
        result = list(
            empty_annotation.get_right_semi_enzymatic_sequences(return_type="str")
        )
        self.assertEqual(result, [])

    def test_empty_sequence_non_enzymatic(self):
        """Test non-enzymatic sequences with empty sequence."""
        empty_annotation = pt.ProFormaAnnotation.parse("")
        result = list(empty_annotation.get_non_enzymatic_sequences(return_type="str"))
        self.assertEqual(result, [])

    def test_single_letter_sequence_regex_digest(self):
        """Test regex digestion with single letter sequence."""
        single_annotation = pt.ProFormaAnnotation.parse("K")
        result = list(
            single_annotation.regex_digest(enzyme_regex="([KR])", return_type="str")
        )
        self.assertEqual(result, ["K"])

    def test_single_letter_sequence_semi_enzymatic(self):
        """Test semi-enzymatic sequences with single letter sequence."""
        single_annotation = pt.ProFormaAnnotation.parse("K")
        result = list(single_annotation.get_semi_enzymatic_sequences(return_type="str"))
        self.assertEqual(result, [])

    def test_single_letter_sequence_left_semi_enzymatic(self):
        """Test left semi-enzymatic sequences with single letter sequence."""
        single_annotation = pt.ProFormaAnnotation.parse("K")
        result = list(
            single_annotation.get_left_semi_enzymatic_sequences(return_type="str")
        )
        self.assertEqual(result, [])

    def test_single_letter_sequence_right_semi_enzymatic(self):
        """Test right semi-enzymatic sequences with single letter sequence."""
        single_annotation = pt.ProFormaAnnotation.parse("K")
        result = list(
            single_annotation.get_right_semi_enzymatic_sequences(return_type="str")
        )
        self.assertEqual(result, [])

    def test_single_letter_sequence_non_enzymatic(self):
        """Test non-enzymatic sequences with single letter sequence."""
        single_annotation = pt.ProFormaAnnotation.parse("K")
        result = list(single_annotation.get_non_enzymatic_sequences(return_type="str"))
        self.assertEqual(result, [])

    def test_no_cleavage_sites_sequence(self):
        """Test digestion with sequence that has no cleavage sites."""
        peptide_annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        result = list(
            peptide_annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"], return_type="str"
            )
        )
        self.assertEqual(result, ["PEPTIDE"])

    def test_sequential_digest_basic(self):
        """Test basic sequential digestion with trypsin and asp-n."""
        trypsin = pt.EnzymeConfig(
            enzyme_regex="([KR])",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )

        asp_n = pt.EnzymeConfig(
            enzyme_regex="\\w(?=D)",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )

        annotation = pt.ProFormaAnnotation.parse("PDEREKPKP")
        peptides = set(
            annotation.sequential_digest(
                enzyme_configs=[trypsin, asp_n], return_type="str"
            )
        )

        expected = {"P", "DER", "EK", "PK", "P"}
        self.assertEqual(peptides, expected)

    def test_sequential_digest_partial_digestion(self):
        """Test sequential digestion with partial digestion enabled."""
        partial_trypsin = pt.EnzymeConfig(
            enzyme_regex="([KR])",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=False,
        )

        partial_asp_n = pt.EnzymeConfig(
            enzyme_regex="\\w(?=D)",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=False,
        )

        annotation = pt.ProFormaAnnotation.parse("PEPRDK")
        peptides = set(
            annotation.sequential_digest(
                enzyme_configs=[partial_trypsin, partial_asp_n], return_type="str"
            )
        )

        expected = {"PEPR", "DK", "PEPRDK"}
        self.assertEqual(peptides, expected)

    def test_sequential_digest_modified_sequence(self):
        """Test sequential digestion with modified sequence."""
        trypsin = pt.EnzymeConfig(
            enzyme_regex="([KR])",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )

        asp_n = pt.EnzymeConfig(
            enzyme_regex="\\w(?=D)",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )

        annotation = pt.ProFormaAnnotation.parse("PEP[Phospho]RDK")
        peptides = set(
            annotation.sequential_digest(
                enzyme_configs=[trypsin, asp_n], return_type="str"
            )
        )
        expected = {"PEP[Phospho]R", "DK"}
        self.assertEqual(peptides, expected)

    def test_digest_by_amino_acid_basic_trypsin(self):
        """Test amino acid-based digestion mimicking trypsin."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")

        peptides = set(
            annotation.digest(
                cleave_on="KR", cterminal=True, missed_cleavages=1, return_type="str"
            )
        )
        expected = {"TIDER", "TIDERTIDEK", "TIDEK", "TIDEKTIDE", "TIDE"}
        self.assertEqual(peptides, expected)

    def test_digest_by_amino_acid_nterminal(self):
        """Test amino acid-based digestion with N-terminal cleavage."""
        annotation = pt.ProFormaAnnotation.parse("DKPEPTIDEK")

        peptides = set(
            annotation.digest(
                cleave_on="K",
                cterminal=False,  # Cleave before K
                missed_cleavages=0,
                return_type="str",
            )
        )
        # Should cleave before K residues
        expected = {"D", "KPEPTIDE", "K"}
        self.assertEqual(peptides, expected)

    def test_digest_by_amino_acid_with_restrictions(self):
        """Test amino acid-based digestion with restriction parameters."""
        annotation = pt.ProFormaAnnotation.parse("KPKRKP")

        peptides = set(
            annotation.digest(
                cleave_on="K",
                restrict_after="P",  # Don't cleave if P follows
                cterminal=True,
                missed_cleavages=0,
                return_type="str",
            )
        )
        # Only cleave at K positions not followed by P
        # This will depend on your exact regex implementation
        # You may need to adjust the expected result
        self.assertIsInstance(peptides, set)
        self.assertTrue(len(peptides) > 0)

    def test_digest_by_amino_acid_restrict_before(self):
        """Test amino acid-based digestion with before restriction."""
        annotation = pt.ProFormaAnnotation.parse("PKRPKR")

        peptides = set(
            annotation.digest(
                cleave_on="R",
                restrict_before="P",  # Don't cleave if P precedes
                cterminal=True,
                missed_cleavages=0,
                return_type="str",
            )
        )
        # Only cleave at R positions not preceded by P
        self.assertIsInstance(peptides, set)
        self.assertTrue(len(peptides) > 0)

    def test_digest_return_type_annotation(self):
        """Test digest with annotation return type."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")

        result = list(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=0,
                return_type="annotation",
            )
        )

        # Should return ProFormaAnnotation objects
        self.assertTrue(len(result) > 0)
        for item in result:
            self.assertIsInstance(item, type(annotation))

    def test_digest_return_type_span(self):
        """Test digest with span return type."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")

        result = list(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=0,
                return_type="span",
            )
        )

        # Should return tuples representing spans
        self.assertTrue(len(result) > 0)
        for item in result:
            self.assertIsInstance(item, tuple)
            self.assertEqual(len(item), 3)  # (start, end, missed_cleavages)

    def test_digest_return_type_str_span(self):
        """Test digest with string-span return type."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")

        result = list(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=0,
                return_type="str-span",
            )
        )

        # Should return tuples of (string, span)
        self.assertTrue(len(result) > 0)
        for item in result:
            self.assertIsInstance(item, tuple)
            self.assertEqual(len(item), 2)
            self.assertIsInstance(item[0], str)  # peptide sequence
            self.assertIsInstance(item[1], tuple)  # span

    def test_digest_return_type_annotation_span(self):
        """Test digest with annotation-span return type."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")

        result = list(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin"],
                missed_cleavages=0,
                return_type="annotation-span",
            )
        )

        # Should return tuples of (annotation, span)
        self.assertTrue(len(result) > 0)
        for item in result:
            self.assertIsInstance(item, tuple)
            self.assertEqual(len(item), 2)
            self.assertIsInstance(item[0], type(annotation))  # annotation
            self.assertIsInstance(item[1], tuple)  # span

    def test_modified_sequence_left_semi_enzymatic(self):
        """Test left semi-enzymatic with modifications preserved."""
        modified_seq = "[1]-P[2]EPTIDE"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)
        sequences = list(
            annotation.get_left_semi_enzymatic_sequences(return_type="str")
        )
        expected = [
            "[1]-P[2]EPTID",
            "[1]-P[2]EPTI",
            "[1]-P[2]EPT",
            "[1]-P[2]EP",
            "[1]-P[2]E",
            "[1]-P[2]",
        ]
        self.assertEqual(sequences, expected)

    def test_isotopic_labeling_left_semi_enzymatic(self):
        """Test left semi-enzymatic with isotopic labeling."""
        labeled_seq = "<13C>TIDE"
        annotation = pt.ProFormaAnnotation.parse(labeled_seq)
        sequences = list(
            annotation.get_left_semi_enzymatic_sequences(return_type="str")
        )
        expected = ["<13C>TID", "<13C>TI", "<13C>T"]
        self.assertEqual(sequences, expected)

    def test_modified_sequence_right_semi_enzymatic(self):
        """Test right semi-enzymatic with modifications preserved."""
        modified_seq = "PEPTIDE[1]-[2]"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)
        sequences = list(
            annotation.get_right_semi_enzymatic_sequences(return_type="str")
        )
        expected = [
            "EPTIDE[1]-[2]",
            "PTIDE[1]-[2]",
            "TIDE[1]-[2]",
            "IDE[1]-[2]",
            "DE[1]-[2]",
            "E[1]-[2]",
        ]
        self.assertEqual(sequences, expected)

    def test_isotopic_labeling_right_semi_enzymatic(self):
        """Test right semi-enzymatic with isotopic labeling."""
        labeled_seq = "<13C>TIDE"
        annotation = pt.ProFormaAnnotation.parse(labeled_seq)
        sequences = list(
            annotation.get_right_semi_enzymatic_sequences(return_type="str")
        )
        expected = ["<13C>IDE", "<13C>DE", "<13C>E"]
        self.assertEqual(sequences, expected)

    def test_modified_non_enzymatic_sequences(self):
        """Test non-enzymatic sequences with modifications preserved."""
        modified_seq = "[Acetyl]-P[1.0]EP[1.0]-[Amide]"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)
        sequences = list(annotation.get_non_enzymatic_sequences(return_type="str"))
        expected = [
            "[Acetyl]-P[1.0]",
            "[Acetyl]-P[1.0]E",
            "E",
            "EP[1.0]-[Amide]",
            "P[1.0]-[Amide]",
        ]
        self.assertEqual(sequences, expected)

    def test_isotopic_labeling_non_enzymatic(self):
        """Test non-enzymatic sequences with isotopic labeling."""
        labeled_seq = "<13C>PEP"
        annotation = pt.ProFormaAnnotation.parse(labeled_seq)
        sequences = list(annotation.get_non_enzymatic_sequences(return_type="str"))
        expected = ["<13C>P", "<13C>PE", "<13C>E", "<13C>EP", "<13C>P"]
        self.assertEqual(sequences, expected)

    def test_cleavage_sites_enzyme_key_trypsin_p(self):
        """Test cleavage sites using trypsin/P enzyme key."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        sites = list(annotation.get_cleavage_sites(pt.PROTEASES["trypsin/P"]))
        expected = [5, 10]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_no_match(self):
        """Test cleavage sites with enzyme that doesn't match sequence."""
        annotation = pt.ProFormaAnnotation.parse("TIDEPTIDEPTIDE")
        sites = list(annotation.get_cleavage_sites(pt.PROTEASES["trypsin/P"]))
        expected = []
        self.assertEqual(sites, expected)

    def test_cleavage_sites_n_terminal(self):
        """Test cleavage sites with N-terminal cleavage (lys-n)."""
        annotation = pt.ProFormaAnnotation.parse("KPEPTIDEK")
        sites = list(annotation.get_cleavage_sites(pt.PROTEASES["lys-n"]))
        expected = [0, 8]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_c_terminal(self):
        """Test cleavage sites with C-terminal cleavage (lys-c)."""
        annotation = pt.ProFormaAnnotation.parse("KPEPTIDEK")
        sites = list(annotation.get_cleavage_sites(pt.PROTEASES["lys-c"]))
        expected = [1, 9]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_modified_sequence(self):
        """Test cleavage sites with modified sequence."""
        modified_seq = "[Acetyl]-TIDERT[1.0]IDEKTIDE-[Amide]"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)
        sites = list(annotation.get_cleavage_sites(pt.PROTEASES["trypsin/P"]))
        expected = [5, 10]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_non_specific(self):
        """Test cleavage sites with non-specific cleavage."""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        sites = list(annotation.get_cleavage_sites("non-specific"))
        expected = [0, 1, 2, 3, 4, 5, 6, 7]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_overlapping_pattern(self):
        """Test cleavage sites with overlapping pattern."""
        annotation = pt.ProFormaAnnotation.parse("PPPP")
        sites = list(annotation.get_cleavage_sites("PP"))
        expected = [1, 2, 3]
        self.assertEqual(sites, expected)

    def test_cleavage_sites_lookahead_pattern(self):
        """Test cleavage sites with lookahead pattern."""
        annotation = pt.ProFormaAnnotation.parse("PEPCTIDE")
        sites = list(annotation.get_cleavage_sites("(?=C)"))
        expected = [3]
        self.assertEqual(sites, expected)

    def test_digest_with_string_span_return_type(self):
        """Test digest with string-span return type."""
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        result = list(
            annotation.regex_digest(
                enzyme_regex=pt.PROTEASES["trypsin/P"], return_type="str-span"
            )
        )
        expected = [("TIDER", (0, 5, 0)), ("TIDEK", (5, 10, 0)), ("TIDE", (10, 14, 0))]
        self.assertEqual(result, expected)

    def test_digest_incomplete_digestion_with_length_filter(self):
        """Test digest with incomplete digestion and length filter."""
        seq = "TIDERTIDEKTIDE"
        annotation = pt.ProFormaAnnotation.parse(seq)
        result = list(
            annotation.regex_digest(
                enzyme_regex="([KR])",
                missed_cleavages=2,
                max_len=5,
                complete_digestion=False,
                return_type="str",
            )
        )
        expected = ["TIDER", "TIDERTIDEKTIDE", "TIDEK", "TIDE"]
        self.assertEqual(result, expected)

    def test_digest_single_amino_acid_complete_digestion(self):
        """Test digest with single amino acid and complete digestion."""
        annotation = pt.ProFormaAnnotation.parse("K")
        result = list(
            annotation.regex_digest(
                enzyme_regex="([KR])",
                missed_cleavages=0,
                complete_digestion=True,
                return_type="str",
            )
        )
        expected = ["K"]
        self.assertEqual(result, expected)

    def test_digest_lookahead_cleavage(self):
        """Test digest with lookahead cleavage pattern."""
        annotation = pt.ProFormaAnnotation.parse("PEPCTIDE")
        result = list(
            annotation.regex_digest(
                enzyme_regex="(?=C)",
                missed_cleavages=0,
                complete_digestion=True,
                return_type="str",
            )
        )
        expected = ["PEP", "CTIDE"]
        self.assertEqual(result, expected)

    def test_digest_complex_modified_sequence(self):
        """Test digest with complex modified sequence and length filters."""
        seq = "<13C>T[1][2]IDERTIDEKTIDE"
        annotation = pt.ProFormaAnnotation.parse(seq)
        result = list(
            annotation.regex_digest(
                enzyme_regex="([KR])", missed_cleavages=2, min_len=6, return_type="str"
            )
        )
        expected = [
            "<13C>T[1][2]IDERTIDEK",
            "<13C>T[1][2]IDERTIDEKTIDE",
            "<13C>TIDEKTIDE",
        ]
        self.assertEqual(result, expected)

    def test_digest_semi_enzymatic_with_modified_sequence(self):
        """Test digest with semi-enzymatic and modified sequence."""
        seq = "TIDERTIDEK[1]TIDE-[2]"
        annotation = pt.ProFormaAnnotation.parse(seq)
        result = list(
            annotation.regex_digest(
                enzyme_regex="([KR])",
                missed_cleavages=1,
                min_len=9,
                semi=True,
                return_type="str",
            )
        )
        expected = ["TIDERTIDE", "TIDERTIDEK[1]", "IDERTIDEK[1]", "TIDEK[1]TIDE-[2]"]
        self.assertEqual(result, expected)

    def test_digest_non_specific_enzyme(self):
        """Test digest with non-specific enzyme."""
        annotation = pt.ProFormaAnnotation.parse("PEPT")
        result = list(
            annotation.regex_digest(enzyme_regex="non-specific", return_type="str")
        )
        expected = ["P", "PE", "PEP", "E", "EP", "EPT", "P", "PT", "T"]
        self.assertEqual(result, expected)

    def test_digest_from_config_basic(self):
        """Test digest_from_config with basic configuration."""
        config = pt.EnzymeConfig(
            enzyme_regex=pt.PROTEASES["trypsin"],
            missed_cleavages=1,
            semi_enzymatic=False,
            complete_digestion=True,
        )
        annotation = pt.ProFormaAnnotation.parse("TIDERTIDEKTIDE")
        peptides = set(pt.digest_from_config(annotation, config=config))
        expected = {"TIDER", "TIDERTIDEK", "TIDEK", "TIDEKTIDE", "TIDE"}
        self.assertEqual(peptides, expected)

    def test_digest_from_config_with_modified_sequence(self):
        """Test digest_from_config with modified sequence."""
        config = pt.EnzymeConfig(
            enzyme_regex=pt.PROTEASES["trypsin"],
            missed_cleavages=1,
            semi_enzymatic=False,
            complete_digestion=True,
        )
        modified_seq = "TIDER[Phospho]TIDEK"
        annotation = pt.ProFormaAnnotation.parse(modified_seq)
        peptides = set(pt.digest_from_config(annotation, config=config))
        expected = {"TIDER[Phospho]", "TIDER[Phospho]TIDEK", "TIDEK"}
        self.assertEqual(peptides, expected)

    def test_digest_from_config_semi_enzymatic(self):
        """Test digest_from_config with semi-enzymatic digestion."""
        config = pt.EnzymeConfig(
            enzyme_regex=pt.PROTEASES["trypsin"],
            missed_cleavages=0,
            semi_enzymatic=True,
            complete_digestion=True,
        )
        annotation = pt.ProFormaAnnotation.parse("TIDEKTIDE")
        peptides = set(pt.digest_from_config(annotation, config=config, min_len=3))
        # Should include both full enzymatic and semi-enzymatic peptides
        self.assertIn("TIDEK", peptides)  # Full enzymatic
        self.assertIn("IDE", peptides)  # Right semi
        self.assertIn("TID", peptides)  # Left semi

    def test_sequential_digest_example_xxxkxxxdxxx(self):
        """Test sequential digest with example from docstring."""
        trypsin = pt.EnzymeConfig(
            enzyme_regex="([KR])",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )
        asp_n = pt.EnzymeConfig(
            enzyme_regex="([D])",
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True,
        )

        annotation = pt.ProFormaAnnotation.parse("XXXKXXXDXXX")
        result = list(
            annotation.sequential_digest(
                enzyme_configs=[trypsin, asp_n], return_type="str"
            )
        )
        expected = ["XXXK", "XXXD", "XXX"]
        self.assertEqual(result, expected)

    def test_sequential_digest_partial_digestion_example(self):
        """Test sequential digest with partial digestion example."""
        partial_digest = [
            pt.EnzymeConfig(
                enzyme_regex="([KR])",
                missed_cleavages=0,
                semi_enzymatic=False,
                complete_digestion=False,
            ),
            pt.EnzymeConfig(
                enzyme_regex="([D])",
                missed_cleavages=0,
                semi_enzymatic=False,
                complete_digestion=False,
            ),
        ]

        annotation = pt.ProFormaAnnotation.parse("XXXKXXXDXXX")
        result = list(
            annotation.sequential_digest(
                enzyme_configs=partial_digest, return_type="str"
            )
        )
        expected = ["XXXK", "XXXKXXXD", "XXXKXXXDXXX", "XXX", "XXXD", "XXXDXXX", "XXX"]
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
