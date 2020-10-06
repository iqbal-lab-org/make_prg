from pathlib import Path
import random
from unittest import TestCase, skip

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_prg.from_msa.aligned_seq import AlignedSeq
from make_prg.seq_utils import standard_bases
from tests.from_msa import make_alignment, MSA

this_dir = Path(__file__).resolve().parent
data_dir = this_dir.parent / "data" / "make_prg_from_msa"


class TestConsensusString(TestCase):
    def test_all_match(self):
        alignment = make_alignment(["AATTA", "AATTA"])
        result = AlignedSeq.get_consensus(alignment)
        self.assertEqual(result, "AATTA")

    def test_mixed_match_nonmatch(self):
        alignment = make_alignment(["AAGTA", "CATTA"])
        result = AlignedSeq.get_consensus(alignment)
        self.assertEqual(result, "*A*TA")

    def test_indel_nonmatch(self):
        alignment = make_alignment(["AAAA", "A--A"])
        result = AlignedSeq.get_consensus(alignment)
        self.assertEqual(result, "A**A")

    def test_IUPACAmbiguous_nonmatch(self):
        alignment = make_alignment(["RYA", "RTA"])
        result = AlignedSeq.get_consensus(alignment)
        self.assertEqual(result, "**A")

    def test_N_special_treatment(self):
        """
        i)A and N at pos 2 are different, but still consensus
        ii)N and N at pos 0 are same, but not consensus"""
        alignment = make_alignment(["NTN", "NTA"])
        result = AlignedSeq.get_consensus(alignment)
        self.assertEqual(result, "*TA")

    def test_all_gap_nonmatch(self):
        alignment = make_alignment(["A--A", "A--A"])
        result = AlignedSeq.get_consensus(alignment)
        self.assertEqual(result, "A**A")


class TestKmeansClusters(TestCase):
    def test_one_seq_returns_single_id(self):
        alignment = MSA([SeqRecord(Seq("AAAT"), id="s1")])
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        self.assertEqual(result, [["s1"]])

    def test_two_seqs_one_below_min_match_len_separate_clusters(self):
        alignment = MSA(
            [SeqRecord(Seq("AATTTAT"), id="s1"), SeqRecord(Seq("AA---AT"), id="s2")]
        )
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 5], alignment, 5)
        self.assertEqual(result, [["s1"], ["s2"]])

    def test_two_identical_seqs_returns_two_ids_clustered(self):
        alignment = MSA(
            [SeqRecord(Seq("AAAT"), id="s1"), SeqRecord(Seq("AAAT"), id="s2"),]
        )
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        self.assertEqual(result, [["s1", "s2"]])

    def test_sequences_in_short_interval_separate_clusters(self):
        alignment = MSA(
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
        alignment = MSA(
            [SeqRecord(Seq("ARAT"), id="s1"), SeqRecord(Seq("WAAT"), id="s2"),]
        )
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 3], alignment, 5)
        self.assertEqual([["s1"], ["s2"]], result)

    def test_two_identical_sequences_clustered_together(self):
        alignment = MSA(
            [
                SeqRecord(Seq("AAAT"), id="s1"),
                SeqRecord(Seq("AAAT"), id="s2"),
                SeqRecord(Seq("C-CC"), id="s3"),
            ]
        )
        result = AlignedSeq.kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        self.assertEqual([["s1", "s2"], ["s3"]], result)

    def test_all_sequences_below_min_match_len(self):
        alignment = MSA(
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
                    alignment = MSA(records)
                    result = AlignedSeq.kmeans_cluster_seqs_in_interval(
                        [0, seq_len - 1], alignment, 1
                    )
                    self.assertTrue(result[0][0] == "s0")

    def test_one_long_one_short_sequence_separate_and_ordered_clusters(self):
        alignment = MSA(
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


def msas_equal(al1: MSA, al2: MSA):
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
        cls.alignment = MSA(
            [
                SeqRecord(Seq("AAAT"), id="s1"),
                SeqRecord(Seq("C--C"), id="s2"),
                SeqRecord(Seq("AATT"), id="s3"),
                SeqRecord(Seq("GNGG"), id="s4"),
            ]
        )

    def test_get_subalignment_sequence_order_maintained2(self):
        result = AlignedSeq.get_sub_alignment_by_list_id(["s1", "s3"], self.alignment)
        expected = MSA([self.alignment[0], self.alignment[2]])
        self.assertTrue(msas_equal(expected, result))

    def test_get_subalignment_sequence_order_maintained(self):
        """
        Sequences given rearranged are still output in input order
        """
        result = AlignedSeq.get_sub_alignment_by_list_id(["s3", "s1"], self.alignment)
        expected = MSA([self.alignment[0], self.alignment[2]])
        self.assertTrue(msas_equal(expected, result))

    def test_get_subalignment_with_interval(self):
        result = AlignedSeq.get_sub_alignment_by_list_id(
            ["s2", "s3"], self.alignment, [0, 2]
        )
        expected = MSA(
            [SeqRecord(Seq("C--"), id="s2"), SeqRecord(Seq("AAT"), id="s3"),]
        )
        self.assertTrue(msas_equal(expected, result))


class TestMakePrgFromMsaFile_IntegrationTests(TestCase):
    def test_answers(self):
        infile = data_dir / "match.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "ACGTGTTTTGTAACTGTGCCACACTCTCGAGACTGCATATGTGTC")

        infile = data_dir / "nonmatch.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, " 5 AAACGTGGTT 6 CCCCCCCCCC 5 ")

        infile = data_dir / "match.nonmatch.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 TGGTT 6 CCCCC 5 ")

        infile = data_dir / "nonmatch.match.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, " 5 AAACGT 6 CCCCCC 5 GGTT")

        infile = data_dir / "match.nonmatch.match.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = data_dir / "shortmatch.nonmatch.match.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, " 5 AAACGT 6 ATTTTC 5 GGTT")

        infile = data_dir / "match.nonmatch.shortmatch.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAAC 5 GTGGTT 6 CCCCCT 5 ")

        infile = data_dir / "match.staggereddash.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACGTGGTT")

        infile = data_dir / "contains_n.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = data_dir / "contains_RYKMSW.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = data_dir / "contains_n_and_RYKMSW.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = data_dir / "contains_n_and_RYKMSW_no_variants.fa"
        aseq = AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACGTGGTT")

        # with pytest.raises(Exception):
        #    aseq = AlignedSeq("test/fails.fa")
