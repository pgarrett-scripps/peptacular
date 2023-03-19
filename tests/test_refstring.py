import unittest

from peptacular.refseq import RefString


class TestRefString(unittest.TestCase):
    def setUp(self):
        self.rs = RefString("Hello World")

    def test_init(self):
        self.assertIsInstance(self.rs._value, memoryview)

        with self.assertRaises(ValueError):
            RefString(123)

    def test_repr(self):
        self.assertEqual(repr(self.rs), "RefString(Hello World)")

    def test_str(self):
        self.assertEqual(str(self.rs), "Hello World")
        self.assertEqual(self.rs, "Hello World")
        self.assertEqual(self.rs, RefString("Hello World"))
        self.assertEqual(self.rs[1:-2], RefString("ello Wor"))
        self.assertEqual(self.rs[1:-2], "ello Wor")
        self.assertEqual(self.rs[1:9], "Hello World"[1:9])

    def test_getitem(self):
        self.assertIsInstance(self.rs[0], RefString)
        self.assertEqual(self.rs[0], "H")
        self.assertIsInstance(self.rs[:5], RefString)
        self.assertEqual(self.rs[:5], "Hello")

    def test_len(self):
        self.assertEqual(len(self.rs), 11)

    def test_eq(self):
        self.assertEqual(self.rs, "Hello World")
        self.assertNotEqual(self.rs, "Goodbye World")

    def test_startswith(self):
        self.assertTrue(self.rs.startswith("Hello"))
        self.assertFalse(self.rs.startswith("Goodbye"))

    def test_endswith(self):
        self.assertTrue(self.rs.endswith("World"))
        self.assertFalse(self.rs.endswith("Universe"))

    def test_find(self):
        self.assertEqual(self.rs.find("World"), 6)
        self.assertEqual(self.rs.find("Universe"), -1)

    def test_rfind(self):
        self.assertEqual(self.rs.rfind("o"), 7)
        self.assertEqual(self.rs.rfind("Universe"), -1)

    def test_count(self):
        self.assertEqual(self.rs.count("o"), 2)
        self.assertEqual(self.rs.count("Universe"), 0)

    def test_split(self):
        self.assertEqual(self.rs.split(), ["Hello", "World"])
        self.assertEqual(self.rs.split("o"), ["Hell", " W", "rld"])

    def test_replace(self):
        self.assertEqual(self.rs.replace("Hello", "Goodbye"), "Goodbye World")
        self.assertEqual(self.rs.replace("o", "z", 1), "Hellz World")

    def test_join(self):
        s = '-'
        rs = RefString(s)
        words = ["apple", "banana", "cherry"]
        result = rs.join(words)
        expected = s.join(words)
        self.assertEqual(result, expected)
        self.assertEqual(result, RefString(expected))

    def test_lower(self):
        s = 'Hello World'
        rs = RefString(s)
        result = rs.lower()
        expected = s.lower()
        self.assertEqual(result, expected)
        self.assertEqual(result, RefString(expected))

    def test_upper(self):
        s = 'Hello World'
        rs = RefString(s)
        result = rs.upper()
        expected = s.upper()
        self.assertEqual(result, expected)
        self.assertEqual(result, RefString(expected))

    def test_strip(self):
        s = "   Hello World   "
        rs = RefString(s)
        result = rs.strip()
        expected = s.strip()
        self.assertEqual(result, expected)
        self.assertEqual(result, RefString(expected))

    def test_lstrip(self):
        s = "   Hello World   "
        rs = RefString(s)
        result = rs.lstrip()
        expected = s.lstrip()
        self.assertEqual(result, expected)
        self.assertEqual(result, RefString(expected))


    def test_rstrip(self):
        s = "   Hello World   "
        rs = RefString(s)
        result = rs.rstrip()
        expected = s.rstrip()
        self.assertEqual(result, expected)
        self.assertEqual(result, RefString(expected))

    def test_ljust(self):
        s = "Hello"
        rs = RefString(s)
        result = rs.ljust(10, "-")
        expected = s.ljust(10, "-")
        self.assertEqual(result, expected)
        self.assertEqual(result, RefString(expected))

    def test_rjust(self):
        s = "Hello"
        rs = RefString(s)
        result = rs.rjust(10, "-")
        expected = s.rjust(10, "-")
        self.assertEqual(result, expected)
        self.assertEqual(result, RefString(expected))

    def test_center(self):
        s = "Hello"
        rs = RefString(s)
        result = rs.center(10, "-")
        expected = s.center(10, '-')
        self.assertEqual(result, expected)
        self.assertEqual(result, RefString(expected))

    def test_hash(self):
        self.assertEqual(hash(RefString('Hello')), hash(RefString('Hello')))
        self.assertEqual(hash(RefString('Hello')), hash('Hello'))
        self.assertNotEqual(hash(RefString('World')), hash(RefString('Hello')))
        self.assertNotEqual(hash(RefString('World')), hash('Hello'))


if __name__ == '__main__':
    unittest.main()
