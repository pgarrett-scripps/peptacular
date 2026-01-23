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

    def test_internal_charged_annotation(self):
        """Test creating a complex annotation directly"""
        sequence = "PEP[Formula:Na:z+1]TIDE/2"
        annot = pt.ProFormaAnnotation.parse(sequence)
        assert annot.serialize() == sequence
        assert annot.charge_state == 2
        assert annot.frag().charge_state == 3

    def test_labile_mod_dropped_from_annotation1(self):
        """Test creating a complex annotation directly"""
        sequence = "{+100}PEPTIDE"
        base_mass = pt.mass("PEPTIDE", ion_type="y")
        annot = pt.ProFormaAnnotation.parse(sequence)
        assert annot.serialize() == sequence
        assert annot.frag(ion_type="y").mass == base_mass

    def test_labile_mod_dropped_from_annotation2(self):
        """Test creating a complex annotation directly"""
        sequence = "{+100}PEPTIDE"
        base_mass = pt.mass("PEPTIDE", ion_type="p")
        annot = pt.ProFormaAnnotation.parse(sequence)
        assert annot.serialize() == sequence
        assert annot.frag(ion_type="p").mass == base_mass + 100.0

    def test_labile_mod_dropped_from_annotation3(self):
        """Test creating a complex annotation directly"""
        sequence = "{+100}PEPTIDE"
        base_mass = pt.mass("PEPTIDE", ion_type="n")
        annot = pt.ProFormaAnnotation.parse(sequence)
        assert annot.serialize() == sequence
        assert annot.frag(ion_type="n").mass == base_mass + 100.0

    def test_cant_fragment_with_unknown_mod(self):
        """Test creating a complex annotation directly"""
        sequence = "[Oxidation]?PEMAT"
        annot = pt.ProFormaAnnotation.parse(sequence)
        with pytest.raises(ValueError):
            annot.fragment()

    def test_cant_fragment_range(self):
        """Test creating a complex annotation directly"""
        sequence = "PRT(ESFRMS)[+19.0523]ISK"
        annot = pt.ProFormaAnnotation.parse(sequence)
        with pytest.raises(ValueError):
            annot.fragment()

    def test_cant_composition_with_iso_and_delta_mass(self):
        """Test creating a complex annotation directly"""
        sequence = "<13C>PEPTI[+100]ISK"
        annot = pt.ProFormaAnnotation.parse(sequence)
        with pytest.raises(ValueError):
            annot.comp()

    def test_can_mass_with_iso_and_delta_mass(self):
        """Test creating a complex annotation directly"""
        sequence = "<13C>PEPTI[+100]DE"
        annot = pt.ProFormaAnnotation.parse(sequence)
        annot.mass()
        m: int | float = pt.mass("<13C>PEPTIDE") + 100.0
        assert annot.mass() == m
