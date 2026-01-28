import peptacular as pt


class TestComp:
    def test_comp(self):
        seq = pt.ProFormaAnnotation.parse("PEPTIDE")
        comp = seq.comp()
        comp = {str(e): c for e, c in comp.items()}
        assert comp == {"C": 34, "H": 53, "N": 7, "O": 15}

    def test_comp_with_isotope(self):
        seq = pt.ProFormaAnnotation.parse("<15N>PEPTIDE")
        comp = seq.comp()
        comp = {str(e): c for e, c in comp.items()}
        assert comp == {"C": 34, "H": 53, "15N": 7, "O": 15}

    def test_comp_with_multiple_isotopes(self):
        seq = pt.ProFormaAnnotation.parse("<15N><13C>PEPTIDE")
        comp = seq.comp()
        comp = {str(e): c for e, c in comp.items()}
        assert comp == {"13C": 34, "H": 53, "15N": 7, "O": 15}

    def test_comp_with_mods(self):
        seq = pt.ProFormaAnnotation.parse("PEPTIDE[Phospho]")
        comp = seq.comp()
        comp = {str(e): c for e, c in comp.items()}
        assert comp == {"C": 34, "H": 54, "N": 7, "O": 18, "P": 1}

    def test_comp_with_mods_and_isotopes(self):
        seq = pt.ProFormaAnnotation.parse("<15N>PEPTIDE[Phospho]")
        comp = seq.comp()
        comp = {str(e): c for e, c in comp.items()}
        assert comp == {"C": 34, "H": 54, "15N": 7, "O": 18, "P": 1}
