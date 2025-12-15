import pytest
import peptacular as pt


class TestApplyMods:
    def test_basic_mass(self):
        seq = pt.ProFormaAnnotation.parse("PEPTIDE")
        mass = seq.base_mass(monoisotopic=True)
        assert mass == pytest.approx(781.349)  # type: ignore

    def test_basic_average_mass(self):
        seq = pt.ProFormaAnnotation.parse("PEPTIDE")
        mass = seq.base_mass(monoisotopic=False)
        assert mass == pytest.approx(781.808)  # type: ignore

    def test_mass_with_internal_mod(self):
        seq = pt.ProFormaAnnotation.parse("PEP[+10]TIDE")
        mass = seq.base_mass(monoisotopic=True)
        assert mass == pytest.approx(791.349)  # type: ignore

    def test_mass_with_nterm_mod(self):
        seq = pt.ProFormaAnnotation.parse("[+10]-PEPTIDE")
        mass = seq.base_mass(monoisotopic=True)
        assert mass == pytest.approx(791.349)  # type: ignore

    def test_mass_with_cterm_mod(self):
        seq = pt.ProFormaAnnotation.parse("PEPTIDE-[+10]")
        mass = seq.base_mass(monoisotopic=True)
        assert mass == pytest.approx(791.349)  # type: ignore

    def test_mass_with_multiple_mods(self):
        seq = pt.ProFormaAnnotation.parse("[+5]-PEP[+10]TIDE-[-15]")
        mass = seq.base_mass(monoisotopic=True)
        assert mass == pytest.approx(781.349)  # type: ignore

    def test_mass_with_labile_mod(self):
        seq = pt.ProFormaAnnotation.parse("{+10}PEPTIDE")
        mass = seq.base_mass(monoisotopic=True)
        assert mass == pytest.approx(791.349)  # type: ignore

    def test_mass_with_unknown_mod(self):
        seq = pt.ProFormaAnnotation.parse("[+10]?PEPTIDE")
        mass = seq.base_mass(monoisotopic=True)
        assert mass == pytest.approx(791.349)  # type: ignore

    def test_mass_with_interval_mod(self):
        seq = pt.ProFormaAnnotation.parse("PEP(TID)[+10]E")
        mass = seq.base_mass(monoisotopic=True)
        assert mass == pytest.approx(791.349)  # type: ignore
