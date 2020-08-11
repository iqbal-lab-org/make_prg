import unittest

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_prg.utils import get_interval_seqs


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

    def test_first_sequence_in_is_first_sequence_out(self):
        alignment = AlignIO.MultipleSeqAlignment(
            [SeqRecord(Seq("TTTT")), SeqRecord(Seq("AAAA")), SeqRecord(Seq("CC-C")),]
        )
        result = get_interval_seqs(alignment)
        expected = ["TTTT", "AAAA", "CCC"]
        self.assertEqual(expected, result)
