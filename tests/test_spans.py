import unittest

from peptacular.spans import build_non_enzymatic_spans, build_left_semi_spans, build_right_semi_spans, \
    build_enzymatic_spans, build_semi_spans

from peptacular import span_to_sequence


class TestSpans(unittest.TestCase):

    def test_get_non_enzymatic_spans(self):
        span = (0, 4, 1)
        expected_output = [(0, 1, 0), (0, 2, 0), (0, 3, 0), (1, 2, 0), (1, 3, 0), (1, 4, 0), (2, 3, 0),
                           (2, 4, 0), (3, 4, 0)]
        self.assertEqual(list(build_non_enzymatic_spans(span)), expected_output)

        span = (0, 3, 2)
        min_len = 1
        max_len = 2
        expected_output = [(0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 3, 0), (2, 3, 0)]
        self.assertEqual(list(build_non_enzymatic_spans(span, min_len, max_len)), expected_output)

    def test_get_left_spans(self):
        span = (0, 4, 1)
        expected_output = [(0, 3, 1), (0, 2, 1), (0, 1, 1)]
        self.assertEqual(list(build_left_semi_spans(span)), expected_output)

        span = (0, 3, 2)
        min_len = 1
        max_len = 2
        expected_output = [(0, 2, 2), (0, 1, 2)]
        self.assertEqual(list(build_left_semi_spans(span, min_len, max_len)), expected_output)

    def test_get_right_spans(self):
        span = (0, 4, 1)
        expected_output = [(1, 4, 1), (2, 4, 1), (3, 4, 1)]
        self.assertEqual(list(build_right_semi_spans(span)), expected_output)

        span = (0, 3, 2)
        min_len = 1
        max_len = 2
        expected_output = [(1, 3, 2), (2, 3, 2)]
        self.assertEqual(list(build_right_semi_spans(span, min_len, max_len)), expected_output)

    def testspan_to_sequence(self):
        sequence = "ABCDEFGH"
        span = (0, 4, 1)
        expected_output = "ABCD"
        self.assertEqual(span_to_sequence(sequence, span), expected_output)

        sequence = "IJKLMNOP"
        span = (2, 6, 1)
        expected_output = "KLMN"
        self.assertEqual(span_to_sequence(sequence, span), expected_output)

    def test_get_non_enzymatic_spans_no_sub_spans(self):
        span = (0, 0, 1)
        expected_output = []
        self.assertEqual(list(build_non_enzymatic_spans(span)), expected_output)

    def test_get_left_spans_span_of_one(self):
        span = (0, 0, 1)
        expected_output = []
        self.assertEqual(list(build_left_semi_spans(span)), expected_output)

    def test_get_right_spans_span_of_one(self):
        span = (0, 0, 1)
        expected_output = []
        self.assertEqual(list(build_right_semi_spans(span)), expected_output)

    def testspan_to_sequence_empty_sequence(self):
        sequence = ""
        span = (0, 0, 1)
        expected_output = ""
        self.assertEqual(span_to_sequence(sequence, span), expected_output)

    def testspan_to_sequence_span_beyond_sequence_length(self):
        sequence = "ABC"
        span = (0, 10, 1)
        expected_output = "ABC"
        self.assertEqual(span_to_sequence(sequence, span), expected_output)

    def testspan_to_sequence_negative_span(self):
        sequence = "ABCDEFGH"
        span = (-4, -1, 1)
        self.assertEqual(span_to_sequence(sequence, span), 'EFG')

    def testspan_to_sequence_none_span(self):
        sequence = "ABCDEFGH"
        span = (-4, None, 1)
        self.assertEqual(span_to_sequence(sequence, span), 'EFGH')

    def test_get_non_enzymatic_spans_non_integer_input(self):
        span = (0.5, 4.5, 1)
        with self.assertRaises(TypeError):
            list(build_non_enzymatic_spans(span))

    def test_get_left_spans_non_integer_input(self):
        span = (0.5, 4.5, 1)
        with self.assertRaises(TypeError):
            list(build_left_semi_spans(span))

    def test_get_right_spans_non_integer_input(self):
        span = (0.5, 4.5, 1)
        with self.assertRaises(TypeError):
            list(build_right_semi_spans(span))

    def testspan_to_sequence_non_integer_input(self):
        sequence = "ABCDEFGH"
        span = (0.5, 4.5, 1)
        with self.assertRaises(TypeError):
            span_to_sequence(sequence, span)

    def test_get_enzymatic_spans(self):
        end_index = 10
        enzyme_sites = [2, 5, 8]
        missed_cleavages = 1
        min_len = 2
        max_len = 6

        expected_spans = [(0, 2, 0), (0, 5, 1), (2, 5, 0), (2, 8, 1), (5, 8, 0), (5, 10, 1), (8, 10, 0)]
        spans = list(build_enzymatic_spans(end_index, enzyme_sites, missed_cleavages, min_len, max_len))
        self.assertEqual(spans, expected_spans)

    def test_get_enzymatic_spans2(self):
        end_index = 10
        enzyme_sites = [2, 5, 8]
        missed_cleavages = 1
        min_len = 2
        max_len = 5

        expected_spans = [(0, 2, 0), (0, 5, 1), (2, 5, 0), (5, 8, 0), (5, 10, 1), (8, 10, 0)]
        spans = list(build_enzymatic_spans(end_index, enzyme_sites, missed_cleavages, min_len, max_len))
        self.assertEqual(spans, expected_spans)

    def test_get_enzymatic_spans_oob(self):
        end_index = 10
        enzyme_sites = [2, 5, 8]
        missed_cleavages = 3
        min_len = None
        max_len = None

        expected_spans = [(0, 2, 0), (0, 5, 1), (0, 8, 2), (0, 10, 3), (2, 5, 0), (2, 8, 1), (2, 10, 2), (5, 8, 0),
                          (5, 10, 1), (8, 10, 0)]
        spans = list(build_enzymatic_spans(end_index, enzyme_sites, missed_cleavages, min_len, max_len))
        self.assertEqual(spans, expected_spans)

    def test_get_semi_spans(self):
        spans = [(0, 5, 0), (0, 8, 1)]
        min_len = 2
        max_len = 5

        expected_semi_spans = [(0, 4, 0), (0, 3, 0), (0, 2, 0), (1, 5, 0), (2, 5, 0), (3, 5, 0), (3, 8, 1),
                               (4, 8, 1), (5, 8, 1), (6, 8, 1)]
        semi_spans = list(build_semi_spans(spans, min_len, max_len))
        self.assertEqual(semi_spans, expected_semi_spans)

    def test_get_semi_spans2(self):
        spans = [(10, 11, 0), (10, 15, 1), (10, 18, 2)]
        min_len = 2
        max_len = 5

        expected_semi_spans = [(10, 14, 1), (10, 13, 1), (10, 12, 1), (11, 15, 1), (12, 15, 1), (13, 15, 1),
                               (13, 18, 2),
                               (14, 18, 2), (15, 18, 2), (16, 18, 2)]
        semi_spans = list(build_semi_spans(spans, min_len, max_len))
        self.assertEqual(semi_spans, expected_semi_spans)

    def test_get_semi_spans3(self):
        spans = [(5, 15, 2), (10, 15, 1), (14, 15, 0)]
        min_len = 2
        max_len = 5

        expected_semi_spans = [(5, 10, 2), (5, 9, 2), (5, 8, 2), (5, 7, 2), (10, 14, 1), (10, 13, 1), (10, 12, 1),
                               (11, 15, 1), (12, 15, 1), (13, 15, 1)]
        semi_spans = list(build_semi_spans(spans, min_len, max_len))
        self.assertEqual(semi_spans, expected_semi_spans)


if __name__ == "__main__":
    unittest.main()
