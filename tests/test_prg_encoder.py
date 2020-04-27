import unittest
from unittest.mock import patch, Mock, call
from io import BytesIO

from hypothesis import given
from hypothesis.strategies import text, characters, integers, from_regex

from make_prg.prg_encoder import PrgEncoder, ConversionError, EncodeError


class TestPrgEncoder(unittest.TestCase):
    def test_dnaToInt_empty_string_raises_assert_error(self):
        encoder = PrgEncoder()
        char = ""

        with self.assertRaises(ConversionError) as context:
            encoder._dna_to_int(char)

        self.assertTrue("Char '' is not in" in str(context.exception))

    @given(characters(blacklist_characters="ACGTacgt", max_codepoint=10000))
    def test_dnaToInt_char_not_valid_raises_assert_error(self, char):
        encoder = PrgEncoder()

        with self.assertRaises(ConversionError):
            encoder._dna_to_int(char)

    def test_dnaToInt_char_valid_returns_int(self):
        encoder = PrgEncoder()
        char = "A"

        actual = encoder._dna_to_int(char)
        expected = 1

        self.assertEqual(actual, expected)

    def test_dnaToInt_char_valid_but_lowercase_returns_int(self):
        encoder = PrgEncoder()
        char = "a"

        actual = encoder._dna_to_int(char)
        expected = 1

        self.assertEqual(actual, expected)

    def test_dnaToInt_char_valid_non_default_encoding(self):
        encoder = PrgEncoder(encoding={"A": 7})
        char = "a"

        actual = encoder._dna_to_int(char)
        expected = 7

        self.assertEqual(actual, expected)

    def test_encode_unit_empty_string_raises_error(self):
        encoder = PrgEncoder()
        unit = ""

        with self.assertRaises(EncodeError):
            encoder._encode_unit(unit)

    @patch.object(PrgEncoder, "_dna_to_int", side_effect=[1, 2, 3, 4])
    def test_encode_unit_dna_returns_list_of_ints_between_1_and_4(
        self, mock_method: Mock
    ):
        encoder = PrgEncoder()
        unit = "ACGT"

        actual = encoder._encode_unit(unit)
        expected = [1, 2, 3, 4]

        self.assertEqual(actual, expected)

        calls = [call(c) for c in "ACGT"]
        mock_method.assert_has_calls(calls)

        self.assertEqual(mock_method.call_count, 4)

    @given(integers(min_value=0))
    def test_encode_unit_integer_string_returns_list_with_just_that_int(self, unit):
        encoder = PrgEncoder()

        actual = encoder._encode_unit(str(unit))
        expected = [unit]

        self.assertEqual(actual, expected)

    def test_encode_unit_prg_with_invalid_chars_raises_error(self):
        encoder = PrgEncoder()
        unit = "foo"

        with self.assertRaises(EncodeError) as context:
            encoder._encode_unit(unit)

        self.assertTrue(
            "Unit {} contains invalid characters".format(unit) in str(context.exception)
        )

    def test_encode_empty_string_returns_empty(self):
        encoder = PrgEncoder()
        prg = ""

        actual = encoder.encode(prg)
        expected = []

        self.assertEqual(actual, expected)

    def test_encode_prg_one_site_deletion(self):
        encoder = PrgEncoder()
        prg = " 5  6 C 5 "

        actual = encoder.encode(prg)
        expected = [5, 6, 2, 5]

        self.assertEqual(actual, expected)

    def test_encode_prg_one_site_deletion2(self):
        encoder = PrgEncoder()
        prg = " 5 C 6  5 "

        actual = encoder.encode(prg)
        expected = [5, 2, 6, 5]

        self.assertEqual(actual, expected)

    def test_encode_prg_with_one_site_and_one_alt_no_spaces_at_ends(self):
        encoder = PrgEncoder()
        prg = "5  6 C 5"

        actual = encoder.encode(prg)
        expected = [5, 6, 2, 5]

        self.assertEqual(actual, expected)

    def test_encode_prg_with_one_site_and_two_alts(self):
        encoder = PrgEncoder()
        prg = "5 A 6 C 5"

        actual = encoder.encode(prg)
        expected = [5, 1, 6, 2, 5]

        self.assertEqual(actual, expected)

    def test_encode_prg_with_one_site_and_two_alts_longer_than_one_base(self):
        encoder = PrgEncoder()
        prg = "5 GA 6 CT 5"

        actual = encoder.encode(prg)
        expected = [5, 3, 1, 6, 2, 4, 5]

        self.assertEqual(actual, expected)

    def test_encode_prg_with_long_site_numbers_and_two_alts_longer_than_one_base(self):
        encoder = PrgEncoder()
        prg = "55 GA 63 Ct 55"

        actual = encoder.encode(prg)
        expected = [55, 3, 1, 63, 2, 4, 55]

        self.assertEqual(actual, expected)

    def test_encode_prg_with_space_at_start_and_only_letters(self):
        encoder = PrgEncoder()
        prg = " a "

        actual = encoder.encode(prg)
        expected = [1]

        self.assertEqual(actual, expected)

    @given(from_regex(r" ?[0-9]* [ACGTacgt]* [0-9]* ", fullmatch=True))
    def test_encode_permutations_of_valid_input(self, prg):
        encoder = PrgEncoder()

        encoder.encode(prg)

        self.assertTrue(True)

    def test_write_encoding_to_empty_encoding_writes_nothing(self):
        encoding = []
        ostream = BytesIO()
        PrgEncoder.write(encoding, ostream)
        ostream.seek(0)

        self.assertEqual(ostream.read(), b"")

    def test_write_encoding_to_encoding_with_one_int(self):
        encoding = [1]
        write_to = BytesIO()
        PrgEncoder.write(encoding, write_to)
        write_to.seek(0)

        actual = write_to.read()
        expected = int(1).to_bytes(4, "little")

        self.assertEqual(actual, expected)

    def test_write_encoding_to_encoding_with_two_ints(self):
        encoding = [1, 4]
        write_to = BytesIO()
        PrgEncoder.write(encoding, write_to)
        write_to.seek(0)

        actual = write_to.read()
        expected = int(1).to_bytes(4, "little") + int(4).to_bytes(4, "little")

        self.assertEqual(actual, expected)
