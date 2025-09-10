import unittest
from unittest.mock import patch
import warnings
import peptacular as pt
from peptacular.proforma.mod_builder import (
    get_mod_index_from_aa,
    get_mod_index_from_regex,
    get_mod_index,
    get_sites,
    ensure_single_static_mod,
    update_mod_list,
    apply_mods,
    build_mods,
)


class TestGetModIndexFromAA(unittest.TestCase):
    def test_single_match_beginning(self):
        """Test single match at beginning of peptide."""
        result = get_mod_index_from_aa("PEPTIDE", "P")
        self.assertEqual(result, {0, 2})

    def test_single_match_end(self):
        """Test single match at end of peptide."""
        result = get_mod_index_from_aa("PEPTIDE", "E")
        self.assertEqual(result, {1, 6})

    def test_single_match_middle(self):
        """Test single match in middle of peptide."""
        result = get_mod_index_from_aa("PEPTIDE", "T")
        self.assertEqual(result, {3})

    def test_multiple_amino_acids(self):
        """Test multiple amino acids in mod pattern."""
        result = get_mod_index_from_aa("PEPTIDE", "PE")
        self.assertEqual(result, {0, 1, 2, 6})

    def test_no_match(self):
        """Test when no matches are found."""
        result = get_mod_index_from_aa("PEPTIDE", "X")
        self.assertEqual(result, set())

    def test_empty_peptide(self):
        """Test with empty peptide sequence."""
        result = get_mod_index_from_aa("", "P")
        self.assertEqual(result, set())

    def test_empty_mod_aa(self):
        """Test with empty modification amino acids."""
        result = get_mod_index_from_aa("PEPTIDE", "")
        self.assertEqual(result, set())

    def test_all_positions_match(self):
        """Test when all positions match."""
        result = get_mod_index_from_aa("PPPP", "P")
        self.assertEqual(result, {0, 1, 2, 3})

    def test_single_character_peptide(self):
        """Test with single character peptide."""
        result = get_mod_index_from_aa("P", "P")
        self.assertEqual(result, {0})


class TestGetModIndexFromRegexImplemented(unittest.TestCase):
    def test_single_character_regex_match(self):
        """Test regex matching single character."""
        result = get_mod_index_from_regex("PEPTIDE", r"P")
        self.assertEqual(result, {0})  # Should match first P only

    def test_regex_match_at_end(self):
        """Test regex matching at end of sequence."""
        result = get_mod_index_from_regex("PEPTIDE", r"E$")
        self.assertEqual(result, {6})  # Should match E at end

    def test_regex_match_at_beginning(self):
        """Test regex matching at beginning of sequence."""
        result = get_mod_index_from_regex("PEPTIDE", r"^P")
        self.assertEqual(result, {0})  # Should match P at start

    def test_regex_no_match(self):
        """Test regex with no matches."""
        result = get_mod_index_from_regex("PEPTIDE", r"X")
        self.assertEqual(result, set())

    def test_regex_multiple_matches(self):
        """Test regex that could match multiple positions."""
        result = get_mod_index_from_regex("PEPTIDE", r"E")
        self.assertEqual(result, {1})  # Should match first E only

    def test_regex_character_class(self):
        """Test regex with character class."""
        result = get_mod_index_from_regex("PEPTIDE", r"[PT]")
        self.assertEqual(result, {0})  # Should match first P

    def test_regex_lookahead_pattern(self):
        """Test regex with lookahead pattern."""
        result = get_mod_index_from_regex("PEPTIDE", r"P(?=E)")
        self.assertEqual(result, {0})  # P followed by E

    def test_regex_lookbehind_pattern(self):
        """Test regex with lookbehind pattern."""
        result = get_mod_index_from_regex("PEPTIDE", r"(?<=P)E")
        self.assertEqual(result, {1})  # E preceded by P

    def test_regex_empty_sequence(self):
        """Test regex with empty sequence."""
        result = get_mod_index_from_regex("", r"P")
        self.assertEqual(result, set())

    def test_regex_case_sensitive(self):
        """Test regex case sensitivity."""
        result = get_mod_index_from_regex("peptide", r"P")
        self.assertEqual(result, set())  # No match due to case

    def test_regex_case_insensitive_flag(self):
        """Test regex with case insensitive flag."""
        result = get_mod_index_from_regex("peptide", r"(?i)P")
        self.assertEqual(result, {0})  # Should match with flag

    def test_regex_word_boundary(self):
        """Test regex with word boundaries."""
        # This tests regex functionality but may not be practically useful for peptides
        result = get_mod_index_from_regex("PEPTIDE", r"\bP")
        # Result depends on regex implementation of word boundaries with amino acids


class TestGetSites(unittest.TestCase):
    def setUp(self):
        self.phospho_mod = pt.Mod("Phospho", 1)
        self.oxidation_mod = pt.Mod("Oxidation", 1)

    def test_empty_mods(self):
        """Test with empty modifications dictionary."""
        result = get_sites("PEPTIDE", {}, is_regex=False)
        self.assertEqual(result, {})

    def test_single_mod_single_aa(self):
        """Test single modification on single amino acid."""
        mods = {"P": [self.phospho_mod]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {0: [self.phospho_mod], 2: [self.phospho_mod]}
        self.assertEqual(result, expected)

    def test_single_mod_multiple_aa(self):
        """Test single modification on multiple amino acids."""
        mods = {"PE": [self.oxidation_mod]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {
            0: [self.oxidation_mod],
            1: [self.oxidation_mod],
            2: [self.oxidation_mod],
            6: [self.oxidation_mod],
        }
        self.assertEqual(result, expected)

    def test_multiple_mods_same_aa(self):
        """Test multiple modifications on same amino acid."""
        mods = {"P": [self.phospho_mod, self.oxidation_mod]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {
            0: [self.phospho_mod, self.oxidation_mod],
            2: [self.phospho_mod, self.oxidation_mod],
        }
        self.assertEqual(result, expected)

    def test_multiple_mod_types(self):
        """Test multiple different modification types."""
        mods = {"P": [self.phospho_mod], "E": [self.oxidation_mod]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {
            0: [self.phospho_mod],
            1: [self.oxidation_mod],
            2: [self.phospho_mod],
            6: [self.oxidation_mod],
        }
        self.assertEqual(result, expected)

    def test_numeric_modifications(self):
        """Test with numeric modification values."""
        mods = {"P": [123.456]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {0: [123.456], 2: [123.456]}
        self.assertEqual(result, expected)

    def test_string_modifications(self):
        """Test with string modification values."""
        mods = {"P": ["CustomMod"]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {0: ["CustomMod"], 2: ["CustomMod"]}
        self.assertEqual(result, expected)

    def test_mixed_modification_types(self):
        """Test with mixed modification value types."""
        mods = {"P": [self.phospho_mod, 123.456, "StringMod"]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        expected = {
            0: [self.phospho_mod, 123.456, "StringMod"],
            2: [self.phospho_mod, 123.456, "StringMod"],
        }
        self.assertEqual(result, expected)

    def test_no_matching_sites(self):
        """Test when no sites match the modification pattern."""
        mods = {"X": [self.phospho_mod]}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        self.assertEqual(result, {})

    def test_empty_mod_list(self):
        """Test with empty modification list."""
        mods = {"P": []}
        result = get_sites("PEPTIDE", mods, is_regex=False)
        self.assertEqual(result, {0: [], 2: []})


class TestEnsureSingleStaticMod(unittest.TestCase):
    def setUp(self):
        self.phospho_mod = pt.Mod("Phospho", 1)
        self.oxidation_mod = pt.Mod("Oxidation", 1)

    def test_single_mod_per_site_no_warning(self):
        """Test no warning when single mod per site."""
        mods = {0: [self.phospho_mod], 1: [self.oxidation_mod]}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ensure_single_static_mod(mods)
            self.assertEqual(len(w), 0)

    def test_multiple_mods_per_site_warning(self):
        """Test warning when multiple mods per site."""
        mods = {0: [self.phospho_mod, self.oxidation_mod]}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ensure_single_static_mod(mods)
            self.assertEqual(len(w), 1)
            self.assertIn("Multiple static modifications", str(w[0].message))

    def test_empty_mods_no_warning(self):
        """Test no warning with empty modifications."""
        mods = {}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ensure_single_static_mod(mods)
            self.assertEqual(len(w), 0)


class TestUpdateModList(unittest.TestCase):
    def setUp(self):
        self.phospho_mod = pt.Mod("Phospho", 1)
        self.oxidation_mod = pt.Mod("Oxidation", 1)

    def test_update_strategy_clears_existing(self):
        """Test update strategy clears existing modifications."""
        mod_list = pt.ModList([self.phospho_mod])
        result = update_mod_list(mod_list, "update", [self.oxidation_mod])
        self.assertEqual(list(result), [self.oxidation_mod])

    def test_merge_strategy_extends_existing(self):
        """Test merge strategy extends existing modifications."""
        mod_list = pt.ModList([self.phospho_mod])
        result = update_mod_list(mod_list, "merge", [self.oxidation_mod])
        self.assertEqual(list(result), [self.phospho_mod, self.oxidation_mod])

    def test_error_strategy_with_empty_list(self):
        """Test error strategy with empty modification list."""
        mod_list = pt.ModList()
        result = update_mod_list(mod_list, "error", [self.phospho_mod])
        self.assertEqual(list(result), [self.phospho_mod])

    def test_error_strategy_with_existing_mods_raises(self):
        """Test error strategy raises when modifications already exist."""
        mod_list = pt.ModList([self.phospho_mod])
        with self.assertRaises(ValueError):
            update_mod_list(mod_list, "error", [self.oxidation_mod])

    def test_unknown_strategy_raises(self):
        """Test unknown strategy raises ValueError."""
        mod_list = pt.ModList()
        with self.assertRaises(ValueError):
            update_mod_list(mod_list, "unknown", [self.phospho_mod])


class TestApplyMods(unittest.TestCase):
    def setUp(self):
        self.annotation = pt.ProFormaAnnotation(sequence="PEPTIDE")
        self.phospho_mod = pt.Mod("Phospho", 1)
        self.oxidation_mod = pt.Mod("Oxidation", 1)
        self.acetyl_mod = pt.Mod("Acetyl", 1)

    def test_apply_internal_mods_only(self):
        """Test applying only internal modifications."""
        internal_mods = {"P": [self.phospho_mod]}
        result = apply_mods(self.annotation, internal=internal_mods, inplace=False)
        self.assertEqual(result.serialize(), "P[Phospho]EP[Phospho]TIDE")

    def test_apply_nterm_mods_only(self):
        """Test applying only N-terminal modifications."""
        nterm_mods = {"P": [self.acetyl_mod]}
        result = apply_mods(self.annotation, nterm=nterm_mods, inplace=False)
        self.assertEqual(result.serialize(), "[Acetyl]-PEPTIDE")

    def test_apply_cterm_mods_only(self):
        """Test applying only C-terminal modifications."""
        cterm_mods = {"E": [self.oxidation_mod]}
        result = apply_mods(self.annotation, cterm=cterm_mods, inplace=False)
        self.assertEqual(result.serialize(), "PEPTIDE-[Oxidation]")

    def test_apply_all_mod_types(self):
        """Test applying all modification types together."""
        nterm_mods = {"P": [self.acetyl_mod]}
        cterm_mods = {"E": [self.oxidation_mod]}
        internal_mods = {"T": [self.phospho_mod]}

        result = apply_mods(
            self.annotation,
            nterm=nterm_mods,
            cterm=cterm_mods,
            internal=internal_mods,
            inplace=False,
        )
        self.assertEqual(result.serialize(), "[Acetyl]-PEPT[Phospho]IDE-[Oxidation]")

    def test_inplace_true_returns_same_object(self):
        """Test inplace=True returns the same annotation object."""
        internal_mods = {"P": [self.phospho_mod]}
        result = apply_mods(self.annotation, internal=internal_mods, inplace=True)
        self.assertIs(result, self.annotation)

    def test_inplace_false_returns_copy(self):
        """Test inplace=False returns a different annotation object."""
        internal_mods = {"P": [self.phospho_mod]}
        result = apply_mods(self.annotation, internal=internal_mods, inplace=False)
        self.assertIsNot(result, self.annotation)
        self.assertEqual(self.annotation.serialize(), "PEPTIDE")  # Original unchanged

    def test_merge_strategy_update(self):
        """Test merge strategy 'update' clears existing mods."""
        # Pre-populate annotation with modifications
        self.annotation.get_internal_mod_dict()[0] = pt.ModList([self.oxidation_mod])

        internal_mods = {"P": [self.phospho_mod]}
        result = apply_mods(
            self.annotation, internal=internal_mods, inplace=True, merge_strat="update"
        )

        # Should only have phospho mod, oxidation should be cleared
        self.assertEqual(result.serialize(), "P[Phospho]EP[Phospho]TIDE")

    def test_merge_strategy_merge(self):
        """Test merge strategy 'merge' combines with existing mods."""
        # Pre-populate annotation with modifications
        self.annotation.get_internal_mod_dict()[0] = pt.ModList([self.oxidation_mod])

        internal_mods = {"P": [self.phospho_mod]}
        result = apply_mods(
            self.annotation, internal=internal_mods, inplace=True, merge_strat="merge"
        )

        # Should have both oxidation and phospho mods
        mod_list = result.get_internal_mod_dict()[0]
        self.assertIn(self.oxidation_mod, mod_list)
        self.assertIn(self.phospho_mod, mod_list)

    def test_none_modifications(self):
        """Test with None modification parameters."""
        result = apply_mods(
            self.annotation, nterm=None, cterm=None, internal=None, inplace=False
        )
        self.assertEqual(result.serialize(), "PEPTIDE")

    def test_empty_modification_dicts(self):
        """Test with empty modification dictionaries."""
        result = apply_mods(
            self.annotation, nterm={}, cterm={}, internal={}, inplace=False
        )
        self.assertEqual(result.serialize(), "PEPTIDE")


class TestBuildMods(unittest.TestCase):
    def setUp(self):
        self.annotation = pt.ProFormaAnnotation(sequence="PEPTIDE")
        self.phospho_mod = pt.Mod("Phospho", 1)
        self.oxidation_mod = pt.Mod("Oxidation", 1)
        self.acetyl_mod = pt.Mod("Acetyl", 1)

    def test_no_modifications_single_result(self):
        """Test build_mods with no modifications returns single result."""
        results = list(build_mods(self.annotation, max_variable_mods=0))
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "PEPTIDE")

    def test_static_modifications_only(self):
        """Test build_mods with only static modifications."""
        results = list(
            build_mods(
                self.annotation,
                internal_static={"P": [self.phospho_mod]},
                max_variable_mods=0,
            )
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "P[Phospho]EP[Phospho]TIDE")

    def test_variable_modifications_only(self):
        """Test build_mods with only variable modifications."""
        results = list(
            build_mods(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
            )
        )
        # Should generate: no mods, P@0, P@2
        self.assertEqual(len(results), 3)
        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)
        self.assertIn("P[Phospho]EPTIDE", serialized)
        self.assertIn("PEP[Phospho]TIDE", serialized)

    def test_static_and_variable_modifications(self):
        """Test build_mods with both static and variable modifications."""
        results = list(
            build_mods(
                self.annotation,
                internal_static={"T": [self.phospho_mod]},
                internal_variable={"P": [self.oxidation_mod]},
                max_variable_mods=1,
            )
        )
        # Should generate: static only, static + P@0, static + P@2
        self.assertEqual(len(results), 3)
        serialized = [r.serialize() for r in results]
        self.assertIn("PEPT[Phospho]IDE", serialized)
        self.assertIn("P[Oxidation]EPT[Phospho]IDE", serialized)
        self.assertIn("PEP[Oxidation]T[Phospho]IDE", serialized)

    def test_nterm_static_modifications(self):
        """Test N-terminal static modifications."""
        results = list(
            build_mods(
                self.annotation,
                nterm_static={"P": [self.acetyl_mod]},
                max_variable_mods=0,
            )
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "[Acetyl]-PEPTIDE")

    def test_cterm_static_modifications(self):
        """Test C-terminal static modifications."""
        results = list(
            build_mods(
                self.annotation,
                cterm_static={"E": [self.oxidation_mod]},
                max_variable_mods=0,
            )
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "PEPTIDE-[Oxidation]")

    def test_nterm_variable_modifications(self):
        """Test N-terminal variable modifications."""
        results = list(
            build_mods(
                self.annotation,
                nterm_variable={"P": [self.acetyl_mod]},
                max_variable_mods=1,
            )
        )
        self.assertEqual(len(results), 2)
        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)
        self.assertIn("[Acetyl]-PEPTIDE", serialized)

    def test_cterm_variable_modifications(self):
        """Test C-terminal variable modifications."""
        results = list(
            build_mods(
                self.annotation,
                cterm_variable={"E": [self.oxidation_mod]},
                max_variable_mods=1,
            )
        )
        self.assertEqual(len(results), 2)
        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)
        self.assertIn("PEPTIDE-[Oxidation]", serialized)

    def test_labile_static_modifications(self):
        """Test labile static modifications."""
        results = list(
            build_mods(
                self.annotation,
                labile_static={"P": [self.phospho_mod]},
                max_variable_mods=0,
            )
        )
        self.assertEqual(len(results), 1)
        # Labile mods should appear in serialization
        serialized = results[0].serialize()
        self.assertIn("{Phospho}^2PEPTIDE", serialized)

    def test_labile_variable_modifications(self):
        """Test labile variable modifications."""
        results = list(
            build_mods(
                self.annotation,
                labile_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
            )
        )
        # Should generate combinations with labile mods
        self.assertEqual(len(results), 3)
        serialized = [r.serialize() for r in results]

        serialized_counter = {s: serialized.count(s) for s in serialized}
        self.assertEqual(serialized_counter.get("PEPTIDE", 0), 1)
        self.assertEqual(serialized_counter.get("{Phospho}PEPTIDE", 0), 2)

    def test_max_variable_mods_limit(self):
        """Test max_variable_mods parameter limits combinations."""
        # With max_variable_mods=1, should get 3 results for 2 possible sites
        results_1 = list(
            build_mods(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
            )
        )
        self.assertEqual(len(results_1), 3)  # no mods, P@0, P@2

        # With max_variable_mods=2, should get 4 results
        results_2 = list(
            build_mods(
                self.annotation,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=2,
            )
        )
        serialized = [r.serialize() for r in results_2]
        self.assertEqual(len(results_2), 4)  # no mods, P@0, P@2, P@0+P@2
        self.assertTrue("PEPTIDE" in serialized)
        self.assertTrue("P[Phospho]EPTIDE" in serialized)
        self.assertTrue("PEP[Phospho]TIDE" in serialized)
        self.assertTrue("P[Phospho]EP[Phospho]TIDE" in serialized)

    def test_site_conflict_prevention(self):
        """Test that site conflicts are prevented."""
        # Two different mods for the same amino acid should not conflict
        # since they're applied to the same sites
        results = list(
            build_mods(
                self.annotation,
                internal_variable={"P": [self.phospho_mod, self.oxidation_mod]},
                max_variable_mods=1,
            )
        )
        # Should generate: no mods, phospho@0, phospho@2, oxidation@0, oxidation@2
        self.assertEqual(len(results), 5)
        serializxed = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serializxed)
        self.assertIn("P[Phospho]EPTIDE", serializxed)
        self.assertIn("PEP[Phospho]TIDE", serializxed)
        self.assertIn("P[Oxidation]EPTIDE", serializxed)
        self.assertIn("PEP[Oxidation]TIDE", serializxed)

    def test_multiple_variable_mod_types(self):
        """Test with multiple types of variable modifications."""
        results = list(
            build_mods(
                self.annotation,
                nterm_variable={"P": [self.acetyl_mod]},
                internal_variable={"T": [self.phospho_mod]},
                cterm_variable={"E": [self.oxidation_mod]},
                max_variable_mods=2,
            )
        )
        # Should generate many combinations
        self.assertEqual(len(results), 7)
        serialized = [r.serialize() for r in results]
        self.assertIn("PEPTIDE", serialized)
        self.assertIn("[Acetyl]-PEPTIDE", serialized)
        self.assertIn("PEPT[Phospho]IDE", serialized)
        self.assertIn("PEPTIDE-[Oxidation]", serialized)
        self.assertIn("[Acetyl]-PEPT[Phospho]IDE", serialized)
        self.assertIn("[Acetyl]-PEPTIDE-[Oxidation]", serialized)
        self.assertIn("PEPT[Phospho]IDE-[Oxidation]", serialized)

    def test_use_regex_parameter(self):
        """Test use_regex parameter is passed through."""
        # This mainly tests that the parameter doesn't cause errors
        results = list(
            build_mods(
                self.annotation,
                internal_static={"P": [self.phospho_mod]},
                max_variable_mods=0,
                use_regex=False,
            )
        )
        self.assertEqual(len(results), 1)

    def test_empty_sequence_annotation(self):
        """Test with empty sequence annotation."""
        empty_annotation = pt.ProFormaAnnotation(sequence="")
        results = list(build_mods(empty_annotation, max_variable_mods=0))
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].serialize(), "")

    def test_single_amino_acid_sequence(self):
        """Test with single amino acid sequence."""
        single_aa = pt.ProFormaAnnotation(sequence="P")
        results = list(
            build_mods(
                single_aa,
                internal_variable={"P": [self.phospho_mod]},
                max_variable_mods=1,
            )
        )
        self.assertEqual(len(results), 2)
        serialized = [r.serialize() for r in results]
        self.assertIn("P", serialized)
        self.assertIn("P[Phospho]", serialized)


if __name__ == "__main__":
    unittest.main()
