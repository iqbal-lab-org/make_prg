from unittest import TestCase
from unittest.mock import patch, Mock, call, MagicMock
from io import BytesIO
import random

from hypothesis import given
from hypothesis.strategies import characters, integers, from_regex

from make_prg.prg_encoder import PrgEncoder, ConversionError, EncodeError, to_bytes


class TestDnaToInt(TestCase):
    def test_empty_string_raises_assert_error(self):
        encoder = PrgEncoder()
        char = ""

        with self.assertRaises(ConversionError) as context:
            encoder._dna_to_int(char)

        self.assertTrue("Char '' is not in" in str(context.exception))

    @given(characters(blacklist_characters="ACGTacgt", max_codepoint=10000))
    def test_char_not_valid_raises_assert_error(self, char):
        encoder = PrgEncoder()

        with self.assertRaises(ConversionError):
            encoder._dna_to_int(char)

    def test_default_encoding_int(self):
        encoder = PrgEncoder()
        uppercase = encoder._dna_to_int("A")
        expected = 1
        self.assertEqual(uppercase, expected)

        lowercase = encoder._dna_to_int("a")
        self.assertEqual(lowercase, expected)

    def test_custom_encoding(self):
        encoder = PrgEncoder(encoding={"A": 7})
        char = "a"

        actual = encoder._dna_to_int(char)
        expected = 7

        self.assertEqual(actual, expected)


class TestInvalidPrgs(TestCase):
    def test_empty_string_fails(self):
        encoder = PrgEncoder()
        with self.assertRaises(EncodeError):
            encoder._encode_unit("")

    def test_repeated_odd_marker_fails(self):
        prg = "5 A 6 C 5 AT 5 T 6 G 5"

        with self.assertRaises(ValueError):
            PrgEncoder().encode(prg)


class TestEncodeUnit(TestCase):
    @patch.object(PrgEncoder, "_dna_to_int", side_effect=[1, 2, 3, 4])
    def test_dna_returns_list_of_ints_between_1_and_4(self, mock_method: Mock):
        encoder = PrgEncoder()
        actual = encoder._encode_unit("ACGT")
        expected = [1, 2, 3, 4]
        self.assertEqual(actual, expected)

        self.assertEqual(mock_method.call_args_list, [call(c) for c in "ACGT"])

    @given(integers(min_value=0))
    def test_single_numeric_chars_converted_to_ints(self, integer):
        encoder = PrgEncoder()
        actual = encoder._encode_unit(str(integer))
        expected = [integer]

        self.assertEqual(actual, expected)

    def test_invalid_string_fails(self):
        encoder = PrgEncoder()
        with self.assertRaises(EncodeError):
            encoder._encode_unit("foo")

    def test_encode_empty_string_returns_empty(self):
        encoder = PrgEncoder()
        actual = encoder.encode("")

        self.assertEqual(actual, [])

    def test_encode_prg_with_one_snp(self):
        encoder = PrgEncoder()
        prg = "5 A 6 C 5"

        actual = encoder.encode(prg)
        expected = [5, 1, 6, 2, 6]

        self.assertEqual(actual, expected)

    def test_encode_prg_one_site_deletion(self):
        encoder = PrgEncoder()
        prg = " 5  6 C 5 "

        actual = encoder.encode(prg)
        expected = [5, 6, 2, 6]

        self.assertEqual(actual, expected)

    def test_encode_prg_one_site_deletion2(self):
        encoder = PrgEncoder()
        prg = " 5 C 6  5 "

        actual = encoder.encode(prg)
        expected = [5, 2, 6, 6]

        self.assertEqual(actual, expected)

    def test_encode_prg_multi_base_alleles(self):
        encoder = PrgEncoder()
        prg = "5 GA 6 CT 5"

        actual = encoder.encode(prg)
        expected = [5, 3, 1, 6, 2, 4, 6]

        self.assertEqual(actual, expected)

    def test_encode_prg_multiple_alleles(self):
        encoder = PrgEncoder()
        prg = "5 GA 6 CT 6 TA 5"

        actual = encoder.encode(prg)
        expected = [5, 3, 1, 6, 2, 4, 6, 4, 1, 6]

        self.assertEqual(actual, expected)

    def test_encode_prg_nested_variation(self):
        encoder = PrgEncoder()
        prg = "5 A 7 C 8 T 8 A 7 6 CT 6 TA 5"

        actual = encoder.encode(prg)
        expected = [5, 1, 7, 2, 8, 4, 8, 1, 8, 6, 2, 4, 6, 4, 1, 6]

        self.assertEqual(actual, expected)

    def test_encode_prg_nonlinear_markers(self):
        encoder = PrgEncoder()
        prg = "55 GA 63 Ct 55"

        actual = encoder.encode(prg)
        expected = [55, 3, 1, 63, 2, 4, 56]

        self.assertEqual(actual, expected)

    def test_encode_prg_spacing_no_variants(self):
        encoder = PrgEncoder()
        prg = " a "

        actual = encoder.encode(prg)
        expected = [1]

        self.assertEqual(actual, expected)

    @given(from_regex(r" ?[0-9]* [ACGTacgt]* [0-9]* ", fullmatch=True))
    def test_permutations_of_valid_input_passes(self, prg):
        encoder = PrgEncoder()
        encoder.encode(prg)


class TestWritePrgInts(TestCase):
    def test_empty_encoding_writes_nothing(self):
        encoding = []
        ostream = BytesIO()
        PrgEncoder.write(encoding, ostream)
        ostream.seek(0)

        self.assertEqual(ostream.read(), b"")

    def test_write_single_int(self):
        prg_ints = [1]
        write_to = BytesIO()
        PrgEncoder.write(prg_ints, write_to)
        write_to.seek(0)

        actual = write_to.read()
        expected = to_bytes(1)
        self.assertEqual(actual, expected)

    def test_write_two_ints(self):
        encoding = [1, 4]
        write_to = BytesIO()
        PrgEncoder.write(encoding, write_to)
        write_to.seek(0)

        actual = write_to.read()
        expected = to_bytes(1) + to_bytes(4)

        self.assertEqual(actual, expected)

    def test_write_multiple_ints_function_calls(self):
        num_elems = 100
        prg_ints = [random.randint for _ in range(num_elems)]

        ostream = BytesIO()
        ostream.write = MagicMock(spec=True)

        # Check call to write triggers the production of bytes
        with patch(
            "make_prg.prg_encoder.to_bytes", return_value=b""
        ) as mocked_to_bytes:
            PrgEncoder.write(prg_ints, ostream)
            expected_bytes_calls = [call(i) for i in prg_ints]
            self.assertEqual(expected_bytes_calls, mocked_to_bytes.call_args_list)

            # Check ostream write operation only called once
            ostream.write.assert_called_once()
