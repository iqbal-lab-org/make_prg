import os
import random
from unittest import TestCase, mock, skip

from Bio.AlignIO import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_prg.make_prg_from_msa import AlignedSeq
from make_prg.exceptions import ClusteringError
from make_prg.seq_utils import standard_bases

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, "tests", "data", "make_prg_from_msa")


@mock.patch.object(AlignedSeq, "check_nonmatch_intervals")
@mock.patch.object(AlignedSeq, "get_prg")
@mock.patch.object(AlignedSeq, "get_consensus")
class TestIntervalPartitioning(TestCase):
    def test_all_non_match(self, get_consensus, _, __):
        get_consensus.return_value = "******"
        tester = AlignedSeq("_", alignment="_", min_match_length=3)
        match, non_match = tester.interval_partition()
        self.assertEqual(match, [])
        self.assertEqual(non_match, [[0, 5]])

    def test_all_match(self, get_consensus, _, __):
        get_consensus.return_value = "ATATAAA"
        tester = AlignedSeq("_", alignment="_", min_match_length=3)
        match, non_match = tester.interval_partition()
        self.assertEqual(match, [[0, 6]])
        self.assertEqual(non_match, [])

    def test_short_match_counted_as_non_match(self, get_consensus, _, __):
        get_consensus.return_value = "AT***"
        tester = AlignedSeq("_", alignment="_", min_match_length=3)
        match, non_match = tester.interval_partition()
        self.assertEqual(match, [])
        self.assertEqual(non_match, [[0, 4]])

    def test_match_non_match_match(self, get_consensus, _, __):
        get_consensus.return_value = "ATT**AAAC"
        tester = AlignedSeq("_", alignment="_", min_match_length=3)
        match, non_match = tester.interval_partition()
        self.assertEqual(match, [[0, 2], [5, 8]])
        self.assertEqual(non_match, [[3, 4]])

    def test_end_in_non_match(self, get_consensus, _, __):
        get_consensus.return_value = "**ATT**AAA*C"
        tester = AlignedSeq("_", alignment="_", min_match_length=3)
        match, non_match = tester.interval_partition()
        self.assertEqual(match, [[2, 4], [7, 9]])
        self.assertEqual(non_match, [[0, 1], [5, 6], [10, 11]])


class TestKmeansClusters(TestCase):
    def test_one_seq_returns_single_id(self):
        alignment = MultipleSeqAlignment([SeqRecord(Seq("AAAT"), id="s1")])
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        self.assertEqual(result, [["s1"]])

    def test_two_seqs_one_below_min_match_len_separate_clusters(self):
        alignment = MultipleSeqAlignment(
            [SeqRecord(Seq("AATTTAT"), id="s1"), SeqRecord(Seq("AA---AT"), id="s2")]
        )
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 5], alignment, 5)
        self.assertEqual(result, [["s1"], ["s2"]])

    def test_two_identical_seqs_returns_two_ids_clustered(self):
        alignment = MultipleSeqAlignment(
            [SeqRecord(Seq("AAAT"), id="s1"), SeqRecord(Seq("AAAT"), id="s2"),]
        )
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        self.assertEqual(result, [["s1", "s2"]])

    def test_sequences_in_short_interval_separate_clusters(self):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAAT"), id="s1"),
                SeqRecord(Seq("AATT"), id="s2"),
                SeqRecord(Seq("AAGT"), id="s3"),
            ]
        )
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 3], alignment, 5)
        # Each sequence is below min_match_len (5), so goes into own cluster
        self.assertEqual([["s1"], ["s2"], ["s3"]], result)

    @skip(
        "This fails, probably because kmean clustering should never run with this input"
    )
    def test_ambiguous_sequences_in_short_interval_separate_clusters(self):
        alignment = MultipleSeqAlignment(
            [SeqRecord(Seq("ARAT"), id="s1"), SeqRecord(Seq("WAAT"), id="s2"),]
        )
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 3], alignment, 5)
        self.assertEqual([["s1"], ["s2"]], result)

    def test_two_identical_sequences_clustered_together(self):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAAT"), id="s1"),
                SeqRecord(Seq("AAAT"), id="s2"),
                SeqRecord(Seq("C-CC"), id="s3"),
            ]
        )
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        self.assertEqual([["s1", "s2"], ["s3"]], result)

    def test_all_sequences_below_min_match_len(self):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AA---AT"), id="s1"),
                SeqRecord(Seq("AA---TT"), id="s2"),
                SeqRecord(Seq("CA--CAT"), id="s3"),
            ]
        )
        result = AlignedSeq.kmeans_cluster_seqs_in_interval(
            [0, len(alignment[0])], alignment, 6
        )
        self.assertEqual([["s1"], ["s2"], ["s3"]], result)


class TestKMeansOrdering(TestCase):
    def test_first_sequence_placed_in_first_cluster(self):
        """
        Runs kmeans clustering on randomly generated multiple sequence alignments
        """
        seq_len = 20
        num_seqs = 20
        bases = list(standard_bases)
        # Function has different behaviour at below and above seq_len
        for seq_len in [seq_len - 1, seq_len + 1]:
            with self.subTest(min_match_len=seq_len):
                for _ in range(20):  # Run on a number of random alignments
                    records = []
                    for i in range(num_seqs):
                        rand_seq = "".join(
                            [random.choice(bases) for _ in range(seq_len)]
                        )
                        records.append(SeqRecord(Seq(rand_seq), id=f"s{i}"))
                    alignment = MultipleSeqAlignment(records)
                    result = AlignedSeq.kmeans_cluster_seqs_in_interval(
                        [0, seq_len - 1], alignment, 1
                    )
                    self.assertTrue(result[0][0] == "s0")

    def test_one_long_one_short_sequence_separate_and_ordered_clusters(self):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AATTAATTATATAATAAC"), id="s1"),
                SeqRecord(Seq("A--------------AAT"), id="s2"),
            ]
        )
        order_1 = AlignedSeq.kmeans_cluster_seqs_in_interval(
            [0, len(alignment[0])], alignment, 5
        )
        self.assertEqual(order_1, [["s1"], ["s2"]])

        order_2 = AlignedSeq.kmeans_cluster_seqs_in_interval(
            [0, len(alignment[0])], alignment[::-1], 5
        )
        self.assertEqual(order_2, [["s2"], ["s1"]])


def msas_equal(al1: MultipleSeqAlignment, al2: MultipleSeqAlignment):
    if len(al1) != len(al2):
        return False
    for i in range(len(al1)):
        if al1[i].seq != al2[i].seq:
            return False
        if al1[i].id != al2[i].id:
            return False
    return True


class TestSubAlignments(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAAT"), id="s1"),
                SeqRecord(Seq("C--C"), id="s2"),
                SeqRecord(Seq("AATT"), id="s3"),
                SeqRecord(Seq("GNGG"), id="s4"),
            ]
        )

    def test_get_subalignment_sequence_order_maintained2(self):
        result = AlignedSeq.get_sub_alignment_by_list_id(["s1", "s3"], self.alignment)
        expected = MultipleSeqAlignment([self.alignment[0], self.alignment[2]])
        self.assertTrue(msas_equal(expected, result))

    def test_get_subalignment_sequence_order_maintained(self):
        """
        Sequences given rearranged are still output in input order
        """
        result = AlignedSeq.get_sub_alignment_by_list_id(["s3", "s1"], self.alignment)
        expected = MultipleSeqAlignment([self.alignment[0], self.alignment[2]])
        self.assertTrue(msas_equal(expected, result))

    def test_get_subalignment_with_interval(self):
        result = AlignedSeq.get_sub_alignment_by_list_id(
            ["s2", "s3"], self.alignment, [0, 2]
        )
        expected = MultipleSeqAlignment(
            [SeqRecord(Seq("C--"), id="s2"), SeqRecord(Seq("AAT"), id="s3"),]
        )
        self.assertTrue(msas_equal(expected, result))


class TestMakePrgFromMsaFile_IntegrationTests(TestCase):
    def test_answers(self):
        infile = os.path.join(data_dir, "match.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "ACGTGTTTTGTAACTGTGCCACACTCTCGAGACTGCATATGTGTC")

        infile = os.path.join(data_dir, "nonmatch.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, " 5 AAACGTGGTT 6 CCCCCCCCCC 5 ")

        infile = os.path.join(data_dir, "match.nonmatch.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 TGGTT 6 CCCCC 5 ")

        infile = os.path.join(data_dir, "nonmatch.match.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, " 5 AAACGT 6 CCCCCC 5 GGTT")

        infile = os.path.join(data_dir, "match.nonmatch.match.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = os.path.join(data_dir, "shortmatch.nonmatch.match.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, " 5 AAACGT 6 ATTTTC 5 GGTT")

        infile = os.path.join(data_dir, "match.nonmatch.shortmatch.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAAC 5 GTGGTT 6 CCCCCT 5 ")

        infile = os.path.join(data_dir, "match.staggereddash.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACGTGGTT")

        infile = os.path.join(data_dir, "contains_n.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = os.path.join(data_dir, "contains_RYKMSW.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = os.path.join(data_dir, "contains_n_and_RYKMSW.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = os.path.join(data_dir, "contains_n_and_RYKMSW_no_variants.fa")
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACGTGGTT")

        # with pytest.raises(Exception):
        #    aseq = AlignedSeq("test/fails.fa")
