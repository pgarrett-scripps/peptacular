"""
Tests for parsing CV accession modification tags.
"""

import peptacular as pt


class TestTagAccession:
    """Tests for parsing CV accession modifications"""

    def test_unimod_accession(self):
        """Test parsing UNIMOD accession"""
        result = pt.parse_modification_tag("UNIMOD:35")
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "35"
        assert result.cv == pt.CV.UNIMOD

    def test_unimod_case_insensitive(self):
        """Test that UNIMOD is case-insensitive"""
        result = pt.parse_modification_tag("unimod:35")
        assert isinstance(result, pt.TagAccession)
        assert result.cv == pt.CV.UNIMOD

    def test_psi_mod_accession(self):
        """Test parsing PSI-MOD accession (use MOD:)"""
        result = pt.parse_modification_tag("MOD:00719")
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "00719"
        assert result.cv == pt.CV.PSI_MOD

    def test_resid_accession(self):
        """Test parsing RESID accession"""
        result = pt.parse_modification_tag("RESID:AA0037")
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "AA0037"
        assert result.cv == pt.CV.RESID

    def test_gnome_accession(self):
        """Test parsing GNO accession (use GNO:)"""
        result = pt.parse_modification_tag("GNO:G12345")
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "G12345"
        assert result.cv == pt.CV.GNOME

    def test_xlmod_accession(self):
        """Test parsing XLMOD accession (full form only)"""
        result = pt.parse_modification_tag("XLMOD:02001")
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "02001"
        assert result.cv == pt.CV.XL_MOD

    def test_xlmod_case_insensitive(self):
        """Test that XLMOD is case-insensitive"""
        result = pt.parse_modification_tag("xlmod:02001")
        assert isinstance(result, pt.TagAccession)
        assert result.accession == "02001"
        assert result.cv == pt.CV.XL_MOD
