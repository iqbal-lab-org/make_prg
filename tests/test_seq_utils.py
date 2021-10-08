from unittest import TestCase

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tests.from_msa import make_alignment
from make_prg.seq_utils import (
    get_alignment_seqs,
    count,
    get_expanded_sequences,
    has_empty_sequence,
)


class TestSeqIteration(TestCase):
    input_seqs = ["TTAGGTTT", "TTA--TTT", "GGA-TTTT"]

    def test_get_seqs_from_alignment(self):
        msa = make_alignment(self.input_seqs)
        result = list(get_alignment_seqs(msa))
        self.assertEqual(result, self.input_seqs)

    def test_count_alignment_seqs(self):
        msa = make_alignment(self.input_seqs)
        result = count(get_alignment_seqs(msa))
        self.assertEqual(result, 3)


class TestGetIntervals(TestCase):
    def test_ambiguous_bases_one_seq(self):
        alignment = AlignIO.MultipleSeqAlignment([SeqRecord(Seq("RWAAT"))])
        result = get_expanded_sequences(alignment)
        expected = {"GAAAT", "AAAAT", "GTAAT", "ATAAT"}
        self.assertEqual(set(result), expected)

    def test_ambiguous_bases_one_seq_with_repeated_base(self):
        alignment = AlignIO.MultipleSeqAlignment([SeqRecord(Seq("RRAAT"))])
        result = get_expanded_sequences(alignment)
        expected = {"GAAAT", "AAAAT", "GGAAT", "AGAAT"}
        self.assertEqual(set(result), expected)

    def test_first_sequence_in_is_first_sequence_out(self):
        alignment = make_alignment(["TTTT", "AAAA", "CC-C"])
        result = get_expanded_sequences(alignment)
        expected = ["TTTT", "AAAA", "CCC"]
        self.assertEqual(expected, result)


class TestEmptySeq(TestCase):
    def test_sub_alignment_with_empty_sequence(self):
        msa = make_alignment(["TTAGGTTT", "TTA--TTT", "GGA-TTTT"])
        self.assertTrue(has_empty_sequence(msa, [3, 4]))
