import pytest
import peptacular as pt
from peptacular.proforma.annotation import ProFormaAnnotation


class TestProFormaAnnotationCopy:
    """Test suite for ProFormaAnnotation copy methods."""
    
    @pytest.fixture
    def simple_annotation(self):
        """Create a simple annotation for testing."""
        return ProFormaAnnotation.parse("PEPTIDE")
    
    @pytest.fixture
    def complex_annotation(self):
        """Create a complex annotation with various modifications."""
        return ProFormaAnnotation.parse(
            "[Acetyl]-PEP[Phospho]T[+79.966]IDE[Oxidation]-[Amidation]/3"
        )
    
    @pytest.fixture
    def annotation_with_internals(self):
        """Create annotation with internal modifications."""
        annot = ProFormaAnnotation(sequence="PEPTIDE")
        annot.append_internal_mod_at_index(0, 79.966)
        annot.append_internal_mod_at_index(3, 15.995)
        return annot
    
    def test_deep_copy_creates_independent_object(self, complex_annotation):
        """Test that deep copy creates completely independent object."""
        original = complex_annotation
        copy = original.copy(deep=True)
        
        # Objects should be equal but not the same
        assert copy == original
        assert copy is not original
        
        # Modify copy
        copy.append_nterm_mod(42.0)
        
        # Original should be unchanged
        assert copy != original
        assert len(original._nterm_mod_list) < len(copy._nterm_mod_list)
    
    def test_shallow_copy_creates_independent_object(self, complex_annotation):
        """Test that shallow copy creates independent object."""
        original = complex_annotation
        copy = original.copy(deep=False)
        
        # Objects should be equal but not the same
        assert copy == original
        assert copy is not original
        
        # Modify copy
        copy.append_nterm_mod(42.0)
        
        # Original should be unchanged
        assert copy != original
        assert len(original._nterm_mod_list) < len(copy._nterm_mod_list)
    
    def test_shallow_copy_containers_are_independent(self, annotation_with_internals):
        """Test that shallow copy creates new ModList/ModDict containers."""
        original = annotation_with_internals
        copy = original.copy(deep=False)
        
        # Should initially be equal
        assert copy == original
        assert copy is not original
        
        # Sequences should be the same reference (strings are immutable)
        assert copy.sequence is original.sequence
        
        # ModLists should be DIFFERENT objects (shallow copied)
        assert copy.get_nterm_mod_list() is not original.get_nterm_mod_list()
        assert copy.get_internal_mod_dict() is not original.get_internal_mod_dict()
        
        # But since Mod is immutable, the Mod objects themselves can be shared
        # This is safe and faster!
    
    def test_copy_allows_independent_nterm_modification(self, simple_annotation):
        """Test that copy allows independent N-term modification."""
        original = simple_annotation
        original.append_nterm_mod(42.0)
        
        copy_deep = original.copy(deep=True)
        copy_shallow = original.copy(deep=False)
        
        # Modify copies
        copy_deep.append_nterm_mod(57.021)
        copy_shallow.append_nterm_mod(79.966)
        
        # Original should only have one mod
        assert len(original._nterm_mod_list) == 1
        assert 42.0 in original._nterm_mod_list
        assert 57.021 not in original._nterm_mod_list
        assert 79.966 not in original._nterm_mod_list
        
        # Copies should have different mods
        assert len(copy_deep._nterm_mod_list) == 2
        assert len(copy_shallow._nterm_mod_list) == 2
        assert 57.021 in copy_deep._nterm_mod_list
        assert 79.966 in copy_shallow._nterm_mod_list
    
    def test_copy_allows_independent_cterm_modification(self, simple_annotation):
        """Test that copy allows independent C-term modification."""
        original = simple_annotation
        original.append_cterm_mod(17.027)
        
        copy_deep = original.copy(deep=True)
        copy_shallow = original.copy(deep=False)
        
        # Modify copies
        copy_deep.append_cterm_mod(15.995)
        copy_shallow.append_cterm_mod(79.966)
        
        # Original should only have one mod
        assert len(original._cterm_mod_list) == 1
        assert 17.027 in original._cterm_mod_list
        
        # Copies should have different mods
        assert len(copy_deep._cterm_mod_list) == 2
        assert len(copy_shallow._cterm_mod_list) == 2
    
    def test_copy_allows_independent_internal_modification(self, annotation_with_internals):
        """Test that copy allows independent internal modification."""
        original = annotation_with_internals
        
        copy_deep = original.copy(deep=True)
        copy_shallow = original.copy(deep=False)
        
        # Modify copies
        copy_deep.append_internal_mod_at_index(5, 42.0)
        copy_shallow.append_internal_mod_at_index(6, 57.021)
        
        # Original should not have mods at index 5 or 6
        assert 5 not in original._internal_mod_dict
        assert 6 not in original._internal_mod_dict
        
        # Copies should have different mods
        assert 5 in copy_deep._internal_mod_dict
        assert 6 in copy_shallow._internal_mod_dict
        assert 5 not in copy_shallow._internal_mod_dict
        assert 6 not in copy_deep._internal_mod_dict
    
    def test_copy_allows_independent_labile_modification(self, simple_annotation):
        """Test that copy allows independent labile modification."""
        original = simple_annotation
        original.append_labile_mod(98.0)
        
        copy_deep = original.copy(deep=True)
        copy_shallow = original.copy(deep=False)
        
        # Modify copies
        copy_deep.append_labile_mod(80.0)
        copy_shallow.append_labile_mod(79.966)
        
        # Original should only have one mod
        assert len(original._labile_mod_list) == 1
        assert 98.0 in original._labile_mod_list
        
        # Copies should have different mods
        assert len(copy_deep._labile_mod_list) == 2
        assert len(copy_shallow._labile_mod_list) == 2
    
    def test_copy_preserves_charge(self, simple_annotation):
        """Test that copy preserves charge state."""
        original = simple_annotation
        original.charge = 2
        
        copy_deep = original.copy(deep=True)
        copy_shallow = original.copy(deep=False)
        
        assert copy_deep.charge == 2
        assert copy_shallow.charge == 2
        
        # Modify copies
        copy_deep.charge = 3
        copy_shallow.charge = 4
        
        # Original should be unchanged
        assert original.charge == 2
    
    def test_copy_empty_annotation(self):
        """Test copying an empty annotation."""
        original = ProFormaAnnotation()
        
        deep = original.copy(deep=True)
        shallow = original.copy(deep=False)
        
        assert deep == original
        assert shallow == original
        assert deep is not original
        assert shallow is not original
        assert len(deep.sequence) == 0
        assert len(shallow.sequence) == 0
    
    def test_copy_preserves_all_mod_types(self, complex_annotation):
        """Test that all modification types are preserved in copy."""
        original = complex_annotation
        
        deep = original.copy(deep=True)
        shallow = original.copy(deep=False)
        
        for copy in [deep, shallow]:
            # Check all mod types are preserved
            assert copy.has_nterm_mods == original.has_nterm_mods
            assert copy.has_cterm_mods == original.has_cterm_mods
            assert copy.has_internal_mods == original.has_internal_mods
            assert copy.has_charge == original.has_charge
            
            if original.has_nterm_mods:
                assert copy.nterm_mods == original.nterm_mods
            if original.has_cterm_mods:
                assert copy.cterm_mods == original.cterm_mods
    
    def test_multiple_copies_independence(self, annotation_with_internals):
        """Test creating multiple copies maintains independence."""
        original = annotation_with_internals
        
        copy1 = original.copy(deep=True)
        copy2 = original.copy(deep=True)
        copy3 = original.copy(deep=False)
        copy4 = original.copy(deep=False)
        
        # Modify each copy differently
        copy1.append_nterm_mod(42.0)
        copy2.append_cterm_mod(17.027)
        copy3.append_labile_mod(98.0)
        copy4.append_static_mod(57.021)
        
        # All should be different from each other
        assert copy1 != copy2
        assert copy1 != copy3
        assert copy1 != copy4
        assert copy2 != copy3
        
        # Original should be unchanged
        assert not original.has_nterm_mods
        assert not original.has_cterm_mods
        assert not original.has_labile_mods
        assert not original.has_static_mods
    
    def test_copy_with_intervals(self):
        """Test copying annotations with interval modifications."""
        original = ProFormaAnnotation(sequence="PEPTIDE")
        original.append_interval((0, 3, [79.966], 0))
        
        deep = original.copy(deep=True)
        shallow = original.copy(deep=False)
        
        for copy in [deep, shallow]:
            assert copy.has_intervals
            assert copy.intervals == original.intervals
            
            # Modify copy
            copy.append_interval((3, 7, [15.995], 0))
            
            # Original should only have one interval
            assert len(original._interval_list) == 1
            assert len(copy._interval_list) == 2
    
    def test_copy_preserves_sequence(self, complex_annotation):
        """Test that sequence is preserved exactly in copy."""
        original = complex_annotation
        
        deep = original.copy(deep=True)
        shallow = original.copy(deep=False)
        
        assert deep.sequence == original.sequence
        assert shallow.sequence == original.sequence
        assert len(deep) == len(original)
        assert len(shallow) == len(original)
    
    def test_hash_consistency_after_copy(self, complex_annotation):
        """Test that hash is consistent between original and copies."""
        original = complex_annotation
        
        deep = original.copy(deep=True)
        shallow = original.copy(deep=False)
        
        # Same content should have same hash
        assert hash(deep) == hash(original)
        assert hash(shallow) == hash(original)
        
        # Can use in sets
        annotation_set = {original, deep, shallow}
        assert len(annotation_set) == 1  # Should be treated as same
    
    def test_copy_with_static_mods(self):
        """Test copying with static modifications."""
        original = ProFormaAnnotation(
            sequence="PEPTIDE",
            static_mods=[57.021, 79.966]
        )
        
        deep = original.copy(deep=True)
        shallow = original.copy(deep=False)
        
        for copy in [deep, shallow]:
            assert copy.has_static_mods
            assert copy.static_mods == original.static_mods
            
            # Modify copy
            copy.append_static_mod(15.995)
            
            # Original should have 2 static mods
            assert len(original._static_mod_list) == 2
            assert len(copy._static_mod_list) == 3
    
    def test_shallow_copy_shares_immutable_mod_objects(self, simple_annotation):
        """Test that shallow copy safely shares immutable Mod objects."""
        original = simple_annotation
        original.append_nterm_mod(42.0)
        
        shallow = original.copy(deep=False)
        
        # The ModLists should be different objects
        assert shallow._nterm_mod_list is not original._nterm_mod_list
        
        # But the Mod objects inside can be the same (they're immutable, so safe!)
        # This is an optimization and doesn't affect correctness
        # We just verify that modifications to one list don't affect the other
        shallow.append_nterm_mod(57.021)
        
        assert len(original._nterm_mod_list) == 1
        assert len(shallow._nterm_mod_list) == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])