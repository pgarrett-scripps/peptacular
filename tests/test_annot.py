import pytest

import peptacular as pt


class TestAnnotationConstruction:
    """Test creating ProFormaAnnotation objects directly"""

    def test_labile_before_global_annotation(self):
        """Test creating an empty annotation"""
        with pytest.raises(ValueError):
            _ = pt.ProFormaAnnotation.parse("{Oxidation}<blabla@M>PEPTIDE")

    def test_complex_annotation(self):
        """Test creating a complex annotation directly"""
        sequence = "(>>>name3)<[Acetyl]@C><[Oxidation]@M>(>>name2)(>name1){Glycan}[Phospho?][Methyl]-PE[Methyl]PTI[Phospho]DE-[Amide][+57]/[Na:z+1^3]"
        annot = pt.ProFormaAnnotation.parse(sequence)
        assert annot.serialize() == sequence
