import pytest
import peptacular as pt


class TestProFormaAnnotationCopy:
    """Test suite for pt.ProFormaAnnotation copy methods."""

    @pytest.fixture
    def simple_annotation(self):
        """Create a simple annotation for testing."""
        return pt.ProFormaAnnotation.parse("PEPTIDE")

    @pytest.fixture
    def complex_annotation(self):
        """Create a complex annotation with various modifications."""
        return pt.ProFormaAnnotation.parse(
            "[Acetyl]-PEP[Phospho]T[+79.966]IDE[Oxidation]-[Amidation]/3"
        )

    @pytest.fixture
    def annotation_with_internals(self):
        """Create annotation with internal modifications."""
        annot = pt.ProFormaAnnotation(sequence="PEPTIDE")
        annot.append_internal_mod_at_index(0, 79.966)
        annot.append_internal_mod_at_index(3, 15.995)
        return annot

    def test_deep_copy_creates_independent_object(
        self, complex_annotation: pt.ProFormaAnnotation
    ):
        """Test that deep copy creates completely independent object."""
        original = complex_annotation
        copy = original.copy()

        # Objects should be equal but not the same
        assert copy == original
        assert copy is not original

        # Modify copy
        copy.append_nterm_mod(42.0)

        # Original should be unchanged
        assert copy != original


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
