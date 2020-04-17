import unittest

from hypothesis import given
from hypothesis.strategies import text
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_prg.seq_utils import remove_gaps, get_interval_seqs


class TestGetIntervals(unittest.TestCase):
    def test_ambiguous_bases_one_seq(self):
        alignment = AlignIO.MultipleSeqAlignment([SeqRecord(Seq("RWAAT"))])
        result = get_interval_seqs(alignment)
        expected = {"GAAAT", "AAAAT", "GTAAT", "ATAAT"}
        self.assertEqual(set(result), expected)

    def test_ambiguous_bases_one_seq_with_repeated_base(self):
        alignment = AlignIO.MultipleSeqAlignment([SeqRecord(Seq("RRAAT"))])
        result = get_interval_seqs(alignment)
        expected = {"GAAAT", "AAAAT", "GGAAT", "AGAAT"}
        self.assertEqual(set(result), expected)


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
