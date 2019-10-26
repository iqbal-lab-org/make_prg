import unittest
from hypothesis import given
from hypothesis.strategies import text
from make_prg.utils import *


class TestRemoveGaps(unittest.TestCase):
    def test_empty_string_returns_empty(self):
        seq = ""

        actual = remove_gaps(seq)
        expected = ""

        self.assertEqual(actual, expected)

    def test_string_with_no_gaps_returns_original(self):
        seq = "ACGT"

        actual = remove_gaps(seq)
        expected = seq

        self.assertEqual(actual, expected)

    def test_string_with_one_gaps_returns_original_without_gap(self):
        seq = "ACGT-"

        actual = remove_gaps(seq)
        expected = "ACGT"

        self.assertEqual(actual, expected)

    def test_string_with_many_gaps_returns_original_without_gap(self):
        seq = "A-CGT--"

        actual = remove_gaps(seq)
        expected = "ACGT"

        self.assertEqual(actual, expected)

    @given(text())
    def test_all_input_space_doesnt_break(self, seq):
        actual = remove_gaps(seq)

        self.assertFalse("-" in actual)
