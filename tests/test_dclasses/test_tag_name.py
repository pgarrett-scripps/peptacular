"""
Tests for parsing named modification tags.
"""

import peptacular as pt


class TestTagName:
    """Tests for parsing named modifications"""

    def test_simple_name(self):
        """Test parsing a simple modification name"""
        result = pt.ModificationTags.from_string("Oxidation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Oxidation"
        assert result.cv is None

    def test_various_names(self):
        """Test parsing various modification names"""
        names = ["Phospho", "Acetyl", "Amidated", "Methylation"]
        for name in names:
            result = pt.ModificationTags.from_string(name).tags[0]
            assert isinstance(result, pt.TagName)
            assert result.name == name

    def test_case_sensitive_name(self):
        """Test that modification names are case-sensitive"""
        result1 = pt.ModificationTags.from_string("Phospho").tags[0]
        assert isinstance(result1, pt.TagName)
        assert result1.name == "Phospho"

        result2 = pt.ModificationTags.from_string("phospho").tags[0]
        assert isinstance(result2, pt.TagName)
        assert result2.name == "phospho"

    def test_unimod_shorthand_with_name(self):
        """Test parsing UNIMOD shorthand U:Oxidation"""
        result = pt.ModificationTags.from_string("U:Oxidation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Oxidation"
        assert result.cv == pt.CV.UNIMOD

    def test_unimod_shorthand_lowercase(self):
        """Test parsing UNIMOD shorthand u:Phospho (case-insensitive)"""
        result = pt.ModificationTags.from_string("u:Phospho").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Phospho"
        assert result.cv == pt.CV.UNIMOD

    def test_psi_mod_shorthand_with_name(self):
        """Test parsing PSI-MOD shorthand M:Phospho"""
        result = pt.ModificationTags.from_string("M:Phospho").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Phospho"
        assert result.cv == pt.CV.PSI_MOD

    def test_psi_mod_shorthand_lowercase(self):
        """Test parsing PSI-MOD shorthand m:Acetyl (case-insensitive)"""
        result = pt.ModificationTags.from_string("m:Acetyl").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Acetyl"
        assert result.cv == pt.CV.PSI_MOD

    def test_resid_shorthand_with_name(self):
        """Test parsing RESID shorthand R:Methylation"""
        result = pt.ModificationTags.from_string("R:Methylation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Methylation"
        assert result.cv == pt.CV.RESID

    def test_resid_shorthand_lowercase(self):
        """Test parsing RESID shorthand r:Biotinylation (case-insensitive)"""
        result = pt.ModificationTags.from_string("r:Biotinylation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Biotinylation"
        assert result.cv == pt.CV.RESID

    def test_gnome_shorthand_with_name(self):
        """Test parsing GNO shorthand G:Glycan"""
        result = pt.ModificationTags.from_string("G:Glycosylation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Glycosylation"
        assert result.cv == pt.CV.GNOME

    def test_gnome_shorthand_lowercase(self):
        """Test parsing GNO shorthand g:Mannosylation (case-insensitive)"""
        result = pt.ModificationTags.from_string("g:Mannosylation").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "Mannosylation"
        assert result.cv == pt.CV.GNOME

    def test_xlmod_shorthand_with_name(self):
        """Test parsing XLMOD shorthand X:Crosslink"""
        result = pt.ModificationTags.from_string("X:DSS").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "DSS"
        assert result.cv == pt.CV.XL_MOD

    def test_xlmod_shorthand_lowercase(self):
        """Test parsing XLMOD shorthand x:BS3 (case-insensitive)"""
        result = pt.ModificationTags.from_string("x:BS3").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "BS3"
        assert result.cv == pt.CV.XL_MOD

    # test serialization
    def test_serialization(self):
        """Test serialization and deserialization of TagName"""
        original = pt.TagName(name="Phospho", cv=pt.CV.UNIMOD)
        serialized = original.serialize()
        assert serialized == "Phospho"

    # test mass and composition retrieval
    def test_mass(self):
        """Test retrieving mass for a known modification"""
        tag = pt.TagName(name="Oxidation", cv=pt.CV.UNIMOD)
        mass = tag.get_mass(monoisotopic=True)
        assert mass is not None
        assert abs(mass - 15.9949) < 0.0001  # Approximate monoisotopic mass of Oxidation

    def test_composition(self):
        """Test retrieving composition for a known modification"""
        tag = pt.TagName(name="Phospho", cv=pt.CV.UNIMOD)
        composition = tag.get_composition()
        assert composition is not None
