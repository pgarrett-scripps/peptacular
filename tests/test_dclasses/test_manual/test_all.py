import peptacular as pt


class TestAllComps:
    def test_modification_tag(self):
        """Test that all components can be created without error"""
        # Test ModificationTag
        mod_tag = pt.ModificationTags.from_string("Phospho")
        assert isinstance(mod_tag, pt.ModificationTags)
