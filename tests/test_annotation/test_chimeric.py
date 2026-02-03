import peptacular as pt


class TestAnnotationConstruction:
    """Test creating ProFormaAnnotation objects directly"""

    def test_parse_chimeric(self):
        """Test creating an empty annotation"""
        annots = pt.parse_chimeric("PEPTIDE+PEPTIDE")
        assert len(annots) == 2

        serialzied_chimeric = pt.serialize_chimeric(annots)
        assert serialzied_chimeric == "PEPTIDE+PEPTIDE"

    def test_serialize_chimeric(self):
        """Test creating an empty annotation"""
        sequence = "PEPTIDE+PEPTIDE"
        annots = pt.parse_chimeric([sequence, sequence, sequence, sequence])
        assert len(annots) == 4
        for annot in annots:
            assert len(annot) == 2

        serialized_chimerics = pt.serialize_chimeric(annots)
        for serialized in serialized_chimerics:
            assert serialized == sequence
