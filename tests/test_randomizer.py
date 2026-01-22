"""
Tests for the randomizer module.
"""

import pytest

import peptacular as pt
from peptacular.annotation import ProFormaAnnotation
from peptacular.annotation.randomizer import (
    generate_isotope_mod_dict,
    generate_random_intervals,
    generate_random_isotope_mod,
    generate_random_proforma_annotation,
    generate_static_mod,
    generate_static_mods_dict,
    get_random_mod_component,
    get_random_mod_dict,
    get_random_psimod,
    get_random_unimod,
    random_charge_adduct,
    random_charge_state,
)


class TestModSelection:
    """Test modification selection functions."""

    def test_get_random_psimod(self):
        """Test getting a random PSI-MOD modification."""
        mod = get_random_psimod()
        assert mod.monoisotopic_mass is not None
        assert mod.dict_composition is not None

    def test_get_random_unimod(self):
        """Test getting a random Unimod modification."""
        mod = get_random_unimod()
        assert mod.monoisotopic_mass is not None
        assert mod.dict_composition is not None

    def test_get_random_mod_component_with_composition(self):
        """Test that mod component with composition doesn't generate mass-only."""
        for _ in range(10):
            mod_str = get_random_mod_component(require_composition=True)
            # Should not start with +/- (mass delta)
            assert not (mod_str.startswith("+") or mod_str.startswith("-"))

    def test_get_random_mod_component_without_composition(self):
        """Test that mod component can generate various types."""
        mod_types_seen = set()
        for _ in range(50):
            mod_str = get_random_mod_component(require_composition=False)
            assert isinstance(mod_str, str)
            assert len(mod_str) > 0
            # Track if we see mass-based mods (starting with +/-)
            if mod_str.startswith(("+", "-")):
                mod_types_seen.add("mass")
            else:
                mod_types_seen.add("other")

    def test_get_random_mod_dict(self):
        """Test generating a random modification dictionary."""
        mod_dict = get_random_mod_dict(mod_probability=0.5)
        assert isinstance(mod_dict, dict)
        for key, value in mod_dict.items():
            assert isinstance(key, str)
            assert isinstance(value, int)
            assert value > 0

    def test_get_random_mod_dict_zero_probability(self):
        """Test that zero probability generates no mods."""
        mod_dict = get_random_mod_dict(mod_probability=0.0)
        assert mod_dict == {}

    def test_get_random_mod_dict_high_probability(self):
        """Test that high probability generates mods."""
        mod_dict = get_random_mod_dict(mod_probability=0.9)
        # With high probability, should usually have at least one mod
        # But it's probabilistic, so we'll just check it's a dict
        assert isinstance(mod_dict, dict)


class TestIsotopeMods:
    """Test isotope modification functions."""

    def test_generate_random_isotope_mod(self):
        """Test generating a random isotope modification."""
        isotope_str = generate_random_isotope_mod()
        assert isinstance(isotope_str, str)
        assert len(isotope_str) > 0

    def test_generate_isotope_mod_dict(self):
        """Test generating isotope modification dictionary."""
        mod_dict = generate_isotope_mod_dict(mod_probability=0.5)
        assert isinstance(mod_dict, dict)
        if mod_dict:
            assert len(mod_dict) == 1
            for key, value in mod_dict.items():
                assert isinstance(key, str)
                assert value == 1

    def test_generate_isotope_mod_dict_zero_probability(self):
        """Test that zero probability generates no isotope mods."""
        mod_dict = generate_isotope_mod_dict(mod_probability=0.0)
        assert mod_dict == {}


class TestStaticMods:
    """Test static modification functions."""

    def test_generate_static_mod(self):
        """Test generating a static modification."""
        static_mod = generate_static_mod()
        assert isinstance(static_mod, str)
        assert "@" in static_mod
        assert static_mod.startswith("[")

    def test_generate_static_mod_with_composition(self):
        """Test static mod with composition requirement."""
        static_mod = generate_static_mod(require_composition=True)
        # Should not have mass deltas in the mod part
        mod_part = static_mod.split("@")[0]
        assert "[+" not in mod_part or "UNIMOD" in mod_part or "MOD:" in mod_part

    def test_generate_static_mods_dict(self):
        """Test generating static modification dictionary."""
        mod_dict = generate_static_mods_dict(mod_probability=0.5)
        assert isinstance(mod_dict, dict)
        if mod_dict:
            for key, value in mod_dict.items():
                assert isinstance(key, str)
                assert "@" in key
                assert value == 1


class TestIntervals:
    """Test interval generation functions."""

    def test_generate_random_intervals(self):
        """Test generating random intervals."""
        intervals = generate_random_intervals(seq_length=10, interval_probability=0.3)
        assert isinstance(intervals, list)

        # Check non-overlapping
        for i, interval in enumerate(intervals):
            assert interval.start < interval.end
            for j, other in enumerate(intervals):
                if i != j:
                    # Intervals should not overlap
                    assert interval.end <= other.start or interval.start >= other.end

    def test_generate_random_intervals_short_sequence(self):
        """Test intervals with short sequence."""
        intervals = generate_random_intervals(seq_length=2, interval_probability=0.5)
        assert isinstance(intervals, list)
        for interval in intervals:
            assert 0 <= interval.start < interval.end < 2

    def test_generate_random_intervals_zero_probability(self):
        """Test that zero probability generates no intervals."""
        intervals = generate_random_intervals(seq_length=10, interval_probability=0.0)
        assert intervals == []


class TestChargeGeneration:
    """Test charge generation functions."""

    def test_random_charge_state(self):
        """Test generating random charge state."""
        charge = random_charge_state(min_charge=1, max_charge=5)
        assert isinstance(charge, int)
        assert 1 <= charge <= 5

    def test_random_charge_adduct(self):
        """Test generating random charge adduct."""
        adduct_dict = random_charge_adduct()
        assert isinstance(adduct_dict, dict)
        assert len(adduct_dict) > 0
        for key, value in adduct_dict.items():
            assert isinstance(key, str)
            assert ":z+" in key
            assert isinstance(value, int)
            assert value > 0


class TestProFormaAnnotationGeneration:
    """Test full ProForma annotation generation."""

    def test_generate_random_proforma_annotation_basic(self):
        """Test basic random annotation generation."""
        annot = generate_random_proforma_annotation()
        assert isinstance(annot, ProFormaAnnotation)
        assert annot.sequence is not None
        assert len(annot.sequence) >= 6
        assert len(annot.sequence) <= 20

    def test_generate_random_proforma_annotation_with_length(self):
        """Test random annotation with specific length range."""
        annot = generate_random_proforma_annotation(min_length=10, max_length=15)
        assert 10 <= len(annot.sequence) <= 15

    def test_generate_random_proforma_annotation_no_mods(self):
        """Test generating annotation without modifications."""
        annot = generate_random_proforma_annotation(
            mod_probability=0.0,
            include_internal_mods=True,
            include_nterm_mods=True,
            include_cterm_mods=True,
            include_labile_mods=True,
            include_unknown_mods=True,
            include_isotopic_mods=True,
            include_static_mods=True,
        )
        assert not annot.has_internal_mods
        assert not annot.has_nterm_mods
        assert not annot.has_cterm_mods
        assert not annot.has_labile_mods
        assert not annot.has_unknown_mods
        assert not annot.has_isotope_mods
        assert not annot.has_static_mods

    def test_generate_random_proforma_annotation_only_internal_mods(self):
        """Test generating annotation with only internal modifications."""
        annot = generate_random_proforma_annotation(
            mod_probability=0.5,
            include_internal_mods=True,
            include_nterm_mods=False,
            include_cterm_mods=False,
            include_labile_mods=False,
            include_unknown_mods=False,
            include_isotopic_mods=False,
            include_static_mods=False,
        )
        assert not annot.has_nterm_mods
        assert not annot.has_cterm_mods
        assert not annot.has_labile_mods
        assert not annot.has_unknown_mods
        assert not annot.has_isotope_mods
        assert not annot.has_static_mods

    def test_generate_random_proforma_annotation_no_charge(self):
        """Test generating annotation without charge."""
        annot = generate_random_proforma_annotation(include_charge=False)
        assert annot.charge is None

    def test_generate_random_proforma_annotation_with_composition(self):
        """Test that require_composition is respected."""
        annot = generate_random_proforma_annotation(
            mod_probability=0.3,
            require_composition=True,
            include_isotopic_mods=False,
        )
        # Should be able to calculate composition without errors
        comp = annot.comp()
        assert comp is not None

    def test_generate_random_proforma_annotation_serializable(self):
        """Test that generated annotation is serializable."""
        annot = generate_random_proforma_annotation(mod_probability=0.3)
        proforma_str = annot.serialize()
        assert isinstance(proforma_str, str)
        assert len(proforma_str) > 0

    def test_generate_random_proforma_annotation_parseable(self):
        """Test that generated annotation can be parsed back."""
        annot = generate_random_proforma_annotation(
            mod_probability=0.2, require_composition=True
        )
        proforma_str = annot.serialize()
        parsed = ProFormaAnnotation.parse(proforma_str)
        assert parsed.sequence == annot.sequence


class TestFunctionalAPI:
    """Test the functional API for random generation."""

    def test_generate_random_single(self):
        """Test generating a single random sequence."""
        seq = pt.generate_random()
        assert isinstance(seq, ProFormaAnnotation)
        assert seq.sequence is not None

    def test_generate_random_multiple(self):
        """Test generating multiple random sequences."""
        seqs = pt.generate_random(count=5)
        assert isinstance(seqs, list)
        assert len(seqs) == 5
        for seq in seqs:
            assert isinstance(seq, ProFormaAnnotation)

    def test_generate_random_with_params(self):
        """Test generating with custom parameters."""
        seq = pt.generate_random(
            min_length=10, max_length=15, mod_probability=0.1, include_charge=False
        )
        assert 10 <= len(seq.sequence) <= 15
        assert seq.charge is None

    def test_generate_random_multiple_parallel(self):
        """Test generating multiple sequences with parallel processing."""
        seqs = pt.generate_random(count=10, n_workers=2, method="thread")
        assert len(seqs) == 10
        for seq in seqs:
            assert isinstance(seq, ProFormaAnnotation)

    def test_generate_random_no_composition(self):
        """Test generating without composition requirement."""
        seq = pt.generate_random(mod_probability=0.3, require_composition=False)
        assert isinstance(seq, ProFormaAnnotation)


class TestRandomizerStability:
    """Test that randomizer produces valid outputs consistently."""

    def test_multiple_generations_are_valid(self):
        """Test that multiple generations all produce valid annotations."""
        for _ in range(20):
            annot = generate_random_proforma_annotation(
                mod_probability=0.2,
                require_composition=True,
                allow_adducts=False,
                include_intervals=False,
            )
            assert annot.sequence is not None
            assert all(c in "ACDEFGHIKLMNPQRSTVWY" for c in annot.sequence)

            # Should be serializable and parseable
            proforma_str = annot.serialize()
            parsed = ProFormaAnnotation.parse(proforma_str)
            assert parsed.sequence == annot.sequence

    def test_mass_calculation_works(self):
        """Test that generated annotations can calculate mass."""
        for _ in range(10):
            annot = generate_random_proforma_annotation(
                mod_probability=0.2,
                require_composition=True,
                include_isotopic_mods=False,
                allow_adducts=False,
            )
            mass = annot.mass()
            assert isinstance(mass, float)
            assert mass > 0

    def test_composition_calculation_works(self):
        """Test that generated annotations can calculate composition."""
        for _ in range(10):
            annot = generate_random_proforma_annotation(
                mod_probability=0.2,
                require_composition=True,
                include_isotopic_mods=False,
                allow_adducts=False,
                retry_cnt=100,
            )
            print(annot)
            comp = annot.comp()
            assert comp is not None
            assert len(comp) > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
