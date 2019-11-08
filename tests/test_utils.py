import unittest
from unittest.mock import patch, Mock, call

from hypothesis import given
from hypothesis.strategies import text, characters, integers

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
    def test_dnaToInt_empty_string_raises_assert_error(self):
        encoder = IntegerEncoder()
        char = ""

        with self.assertRaises(ConversionError) as context:
            encoder._dna_to_int(char)

        self.assertTrue("Char '' is not in" in str(context.exception))

    @given(text(alphabet=characters(blacklist_characters="ACGTacgt")))
    def test_dnaToInt_char_not_valid_raises_assert_error(self, char):
        encoder = IntegerEncoder()

        with self.assertRaises(ConversionError) as context:
            encoder._dna_to_int(char)

        self.assertTrue("Char '{}' is not in".format(char) in str(context.exception))

    def test_dnaToInt_char_valid_returns_int(self):
        encoder = IntegerEncoder()
        char = "A"

        actual = encoder._dna_to_int(char)
        expected = 1

        self.assertEqual(actual, expected)

    def test_dnaToInt_char_valid_but_lowercase_returns_int(self):
        encoder = IntegerEncoder()
        char = "a"

        actual = encoder._dna_to_int(char)
        expected = 1

        self.assertEqual(actual, expected)

    def test_encode_unit_empty_string_raises_error(self):
        encoder = IntegerEncoder()
        unit = ""

        with self.assertRaises(EncodeError) as err:
            encoder._encode_unit(unit)

        self.assertTrue("Cannot encode an empty string")

    @patch.object(IntegerEncoder, "_dna_to_int", side_effect=[1, 2, 3, 4])
    def test_encode_unit_dna_returns_list_of_ints_between_1_and_4(
        self, mock_method: Mock
    ):
        encoder = IntegerEncoder()
        unit = "ACGT"

        actual = encoder._encode_unit(unit)
        expected = [1, 2, 3, 4]

        self.assertEqual(actual, expected)

        calls = [call(c) for c in "ACGT"]
        mock_method.assert_has_calls(calls)

        self.assertEqual(mock_method.call_count, 4)

    @given(integers())
    def test_encode_unit_integer_string_returns_list_with_just_that_int(self, unit):
        encoder = IntegerEncoder()

        actual = encoder._encode_unit(str(unit))
        expected = [unit]

        self.assertEqual(actual, expected)

