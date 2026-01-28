"""
Tests for parsing CV accession modification tags.
"""

import peptacular as pt


class TestTagAccession:
    """Tests for parsing CV accession modifications"""

    def test_unimod_accession(self):
        """Test parsing UNIMOD accession"""
        result = pt.ModificationTags.from_string("UNIMOD:35").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "35"
        assert result.cv == pt.CV.UNIMOD

    def test_unimod_case_insensitive(self):
        """Test that UNIMOD is case-insensitive"""
        result = pt.ModificationTags.from_string("unimod:35").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.cv == pt.CV.UNIMOD

    def test_psi_mod_accession(self):
        """Test parsing PSI-MOD accession (use MOD:)"""
        result = pt.ModificationTags.from_string("MOD:00719").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "00719"
        assert result.cv == pt.CV.PSI_MOD

    def test_resid_accession(self):
        """Test parsing RESID accession"""
        result = pt.ModificationTags.from_string("RESID:AA0037").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "AA0037"
        assert result.cv == pt.CV.RESID

    def test_gnome_accession(self):
        """Test parsing GNO accession (use GNO:)"""
        result = pt.ModificationTags.from_string("GNO:G12345").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "G12345"
        assert result.cv == pt.CV.GNOME

    def test_xlmod_accession(self):
        """Test parsing XLMOD accession (full form only)"""
        result = pt.ModificationTags.from_string("XLMOD:02001").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "02001"
        assert result.cv == pt.CV.XL_MOD

    def test_xlmod_case_insensitive(self):
        """Test that XLMOD is case-insensitive"""
        result = pt.ModificationTags.from_string("xlmod:02001").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "02001"
        assert result.cv == pt.CV.XL_MOD
