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


class TestIntegerEncoder(unittest.TestCase):
    def test_DNAtoInt_empty_string_raises_assert_error(self):
        encoder = IntegerEncoder()
        char = ""

        with self.assertRaises(ConversionError) as context:
            encoder._DNA_to_int(char)

        self.assertTrue("Char '' is not in" in str(context.exception))

    def test_DNAtoInt_char_not_valid_raises_assert_error(self):
        encoder = IntegerEncoder()
        char = "F"

        with self.assertRaises(ConversionError) as context:
            encoder._DNA_to_int(char)

        self.assertTrue("Char 'F' is not in" in str(context.exception))

    def test_DNAtoInt_char_valid_returns_int(self):
        encoder = IntegerEncoder()
        char = "A"

        actual = encoder._DNA_to_int(char)
        expected = 1

        self.assertEqual(actual, expected)

    def test_DNAtoInt_char_valid_but_lowercase_returns_int(self):
        encoder = IntegerEncoder()
        char = "a"

        actual = encoder._DNA_to_int(char)
        expected = 1

        self.assertEqual(actual, expected)
