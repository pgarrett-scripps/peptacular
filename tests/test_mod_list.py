import unittest

from peptacular.proforma_dataclasses import Mod
from peptacular.proforma.mod_list import ModList


class TestModList(unittest.TestCase):
    """Test cases for the ModList class."""

    def setUp(self):
        """Set up test fixtures."""
        self.mod1 = Mod("Oxidation", 1)
        self.mod2 = Mod(15.994915, 1)
        self.mod3 = Mod("Phosphorylation", 2)
        
    def test_init_empty(self):
        """Test initialization with no modifications."""
        mod_list = ModList()
        self.assertEqual(len(mod_list), 0)
        
    def test_init_with_mods(self):
        """Test initialization with modifications."""
        mod_list = ModList([self.mod1, self.mod2])
        self.assertEqual(len(mod_list), 2)
        self.assertEqual(mod_list[0], self.mod1)
        self.assertEqual(mod_list[1], self.mod2)

    def test_init_with_convertible_types(self):
        """Test initialization with convertible types."""
        mod_list = ModList(["Oxidation", 15.994915, 42])
        self.assertEqual(len(mod_list), 3)
        self.assertIsInstance(mod_list[0], Mod)
        self.assertIsInstance(mod_list[1], Mod)
        self.assertIsInstance(mod_list[2], Mod)

    def test_append_convertible_types(self):
        """Test appending convertible types."""
        mod_list = ModList()
        mod_list.append("Oxidation")
        mod_list.append(15.994915)
        
        self.assertEqual(len(mod_list), 2)
        self.assertIsInstance(mod_list[0], Mod)
        self.assertIsInstance(mod_list[1], Mod)

    def test_extend_convertible_types(self):
        """Test extending with convertible types."""
        mod_list = ModList([self.mod1])
        mod_list.extend(["Oxidation", 15.994915])
        self.assertEqual(len(mod_list), 3)
        self.assertIsInstance(mod_list[1], Mod)
        self.assertIsInstance(mod_list[2], Mod)

    def test_setitem_convertible_type(self):
        """Test setting with convertible types."""
        mod_list = ModList([self.mod1])
        mod_list[0] = "Phosphorylation"
        self.assertIsInstance(mod_list[0], Mod)
        self.assertEqual(mod_list[0].val, "Phosphorylation")

    def test_setitem_slice(self):
        """Test setting a slice with convertible types."""
        mod_list = ModList([self.mod1, self.mod2, self.mod3])
        mod_list[1:3] = ["Acetylation", 42.010565]
        self.assertEqual(len(mod_list), 3)
        self.assertIsInstance(mod_list[1], Mod)
        self.assertIsInstance(mod_list[2], Mod)

    def test_contains_with_conversion(self):
        """Test membership with automatic conversion."""
        mod_list = ModList([Mod("Oxidation", 1), self.mod2])
        self.assertTrue("Oxidation" in mod_list)
        self.assertTrue(15.994915 in mod_list)
        self.assertFalse("NonExistent" in mod_list)

    def test_contains_invalid_conversion(self):
        """Test membership with non-convertible types."""
        mod_list = ModList([self.mod1])
        self.assertFalse(None in mod_list)
        self.assertFalse([] in mod_list)

    def test_remove_with_conversion(self):
        """Test removing with conversion."""
        mod_list = ModList([Mod("Oxidation", 1), self.mod2])
        mod_list.remove("Oxidation")
        self.assertEqual(len(mod_list), 1)
        self.assertEqual(mod_list[0], self.mod2)

    def test_count_with_conversion(self):
        """Test counting with conversion."""
        mod_list = ModList([Mod("Oxidation", 1), self.mod2, Mod("Oxidation", 1)])
        count = mod_list.count("Oxidation")
        self.assertEqual(count, 2)

    def test_index_with_conversion(self):
        """Test finding index with conversion."""
        oxidation_mod = Mod("Oxidation", 1)
        mod_list = ModList([oxidation_mod, self.mod2])
        
        index = mod_list.index("Oxidation")
        self.assertEqual(index, 0)

    def test_add_operations(self):
        """Test addition operations."""
        mod_list1 = ModList([self.mod1])
        mod_list2 = ModList([self.mod2])
        
        # Test + operator
        result = mod_list1 + mod_list2
        self.assertIsInstance(result, ModList)
        self.assertEqual(len(result), 2)
        
        # Test += operator
        original_id = id(mod_list1)
        mod_list1 += [self.mod3]
        self.assertEqual(id(mod_list1), original_id)  # Should be same object
        self.assertEqual(len(mod_list1), 2)

    def test_add_with_convertible_types(self):
        """Test addition with convertible types."""
        mod_list = ModList([self.mod1])
        result = mod_list + ["Oxidation", 15.994915]
        
        self.assertIsInstance(result, ModList)
        self.assertEqual(len(result), 3)
        self.assertIsInstance(result[1], Mod)
        self.assertIsInstance(result[2], Mod)

    def test_add_invalid_type(self):
        """Test adding with invalid type raises TypeError."""
        mod_list = ModList([self.mod1])
        with self.assertRaises(TypeError):
            mod_list + "invalid_string"

    def test_validation_error(self):
        """Test that invalid items during construction raise ValueError."""
        with self.assertRaises(ValueError):
            # This should fail in convert_to_mod if we pass something invalid
            ModList([None])  # None should not be convertible

    def test_basic_list_operations(self):
        """Test essential inherited list operations."""
        mod_list = ModList([self.mod1, self.mod2, self.mod3])
        
        # Test indexing and slicing
        self.assertEqual(mod_list[0], self.mod1)
        self.assertEqual(mod_list[-1], self.mod3)
        
        subset = mod_list[1:3]
        self.assertEqual(len(subset), 2)
        
        # Test iteration
        mods_from_iteration = list(mod_list)
        self.assertEqual(len(mods_from_iteration), 3)
        
        # Test pop and clear
        popped = mod_list.pop()
        self.assertEqual(popped, self.mod3)
        self.assertEqual(len(mod_list), 2)
        
        mod_list.clear()
        self.assertEqual(len(mod_list), 0)

    def test_repr(self):
        """Test string representation."""
        mod_list = ModList([self.mod1, self.mod2])
        repr_str = repr(mod_list)
        self.assertTrue(repr_str.startswith("ModList(["))
        self.assertTrue(repr_str.endswith("])"))
        
        # Test empty
        empty_list = ModList()
        self.assertEqual(repr(empty_list), "ModList([])")


if __name__ == "__main__":
    unittest.main()
