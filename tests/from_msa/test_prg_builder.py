from pathlib import Path
from unittest import TestCase

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_prg.from_msa.prg_builder import PrgBuilder
from make_prg.from_msa.interval_partition import IntervalType, Interval
from tests.from_msa import make_alignment, MSA

this_dir = Path(__file__).resolve().parent
data_dir = this_dir.parent / "data" / "make_prg_from_msa"


class TestConsensusString(TestCase):
    def test_all_match(self):
        alignment = make_alignment(["AATTA", "AATTA"])
        result = PrgBuilder.get_consensus(alignment)
        self.assertEqual(result, "AATTA")

    def test_mixed_match_nonmatch(self):
        alignment = make_alignment(["AAGTA", "CATTA"])
        result = PrgBuilder.get_consensus(alignment)
        self.assertEqual(result, "*A*TA")

    def test_indel_nonmatch(self):
        alignment = make_alignment(["AAAA", "A--A"])
        result = PrgBuilder.get_consensus(alignment)
        self.assertEqual(result, "A**A")

    def test_IUPACAmbiguous_nonmatch(self):
        alignment = make_alignment(["RYA", "RTA"])
        result = PrgBuilder.get_consensus(alignment)
        self.assertEqual(result, "**A")

    def test_N_special_treatment(self):
        """
        i)A and N at pos 2 are different, but still consensus
        ii)N and N at pos 0 are same, but not consensus"""
        alignment = make_alignment(["NTN", "NTA"])
        result = PrgBuilder.get_consensus(alignment)
        self.assertEqual(result, "*TA")

    def test_all_gap_nonmatch(self):
        alignment = make_alignment(["A--A", "A--A"])
        result = PrgBuilder.get_consensus(alignment)
        self.assertEqual(result, "A**A")


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
        cls.alignment = make_alignment(
            ["AAAT", "C--C", "AATT", "GNGG"], ["s1", "s2", "s3", "s4"]
        )

    def test_GivenOrderedIds_SubalignmentInSequenceOrder(self):
        result = PrgBuilder.get_sub_alignment_by_list_id(["s1", "s3"], self.alignment)
        expected = MSA([self.alignment[0], self.alignment[2]])
        self.assertTrue(msas_equal(expected, result))

    def test_GivenUnorderedIds_SubalignmentStillInSequenceOrder(self):
        """
        Sequences given rearranged are still output in input order
        """
        result = PrgBuilder.get_sub_alignment_by_list_id(["s3", "s1"], self.alignment)
        expected = MSA([self.alignment[0], self.alignment[2]])
        self.assertTrue(msas_equal(expected, result))

    def test_get_subalignment_with_interval(self):
        result = PrgBuilder.get_sub_alignment_by_list_id(
            ["s2", "s3"], self.alignment, [0, 2]
        )
        expected = make_alignment(["C--", "AAT"], ["s2", "s3"])
        self.assertTrue(msas_equal(expected, result))


class TestSkipClustering(TestCase):
    """
    Test the conditions under which clustering is not performed.
    """

    def setUp(self):
        """
        Set of parameters whereby clustering is to be performed.
        We'll modify each of them in turn
        """
        self.aligned_seqs = ["ATTTTTTA", "A--TTTTA", "ATTTCTTA"]
        self.tested_params = {
            "interval": Interval(IntervalType.Match, 0, 7),
            "max_nesting": 2,
            "nesting_level": 1,
            "min_match_length": 2,
            "alignment": make_alignment(self.aligned_seqs),
        }

    def test_original_params_no_skip_clustering(self):
        self.assertFalse(PrgBuilder.skip_clustering(**self.tested_params))

    def test_max_nesting_reached_skip_clustering(self):
        self.tested_params["nesting_level"] = 2
        self.assertTrue(PrgBuilder.skip_clustering(**self.tested_params))

    def test_small_interval_skip_clustering(self):
        self.tested_params["interval"].stop = 1
        self.assertTrue(PrgBuilder.skip_clustering(**self.tested_params))

    def test_too_few_seqs_skip_clustering(self):
        self.tested_params["alignment"] = self.tested_params["alignment"][0:1]
        self.assertTrue(PrgBuilder.skip_clustering(**self.tested_params))

    def test_ambiguous_alignment_skip_clustering(self):
        """
        `added_seq` below is an equally valid alignment as "A--TTTTA" to the sequence
        "ATTAATTA"
        If we have such ambiguous alignments (defined as more than one gapped alignment
        corresponding to the same ungapped sequence), we choose not to cluster the
        alignment, as it can create ambiguous graphs (whereby different paths spell same sequence)
        """
        added_seq = "ATTTT--A"
        self.tested_params["alignment"] = make_alignment(
            self.aligned_seqs + [added_seq]
        )
        self.assertTrue(PrgBuilder.skip_clustering(**self.tested_params))


class Test_Integration_FullBuilds(TestCase):
    def test_answers_non_nested(self):
        infile = data_dir / "match.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, "ACGTGTTTTGTAACTGTGCCACACTCTCGAGACTGCATATGTGTC")

        infile = data_dir / "nonmatch.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, " 5 AAACGTGGTT 6 CCCCCCCCCC 5 ")

        infile = data_dir / "match.nonmatch.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, "AAACG 5 TGGTT 6 CCCCC 5 ")

        infile = data_dir / "nonmatch.match.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, " 5 AAACGT 6 CCCCCC 5 GGTT")

        infile = data_dir / "match.nonmatch.match.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = data_dir / "shortmatch.nonmatch.match.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, " 5 AAACGT 6 ATTTTC 5 GGTT")

        infile = data_dir / "match.nonmatch.shortmatch.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, "AAAC 5 GTGGTT 6 CCCCCT 5 ")

        infile = data_dir / "match.staggereddash.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, "AAACGTGGTT")

        infile = data_dir / "contains_n.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = data_dir / "contains_RYKMSW.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = data_dir / "contains_n_and_RYKMSW.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")

        infile = data_dir / "contains_n_and_RYKMSW_no_variants.fa"
        aseq = PrgBuilder(infile)
        self.assertEqual(aseq.prg, "AAACGTGGTT")

        # with pytest.raises(Exception):
        #    aseq = AlignedSeq("test/fails.fa")

    def test_nested_snp_backgrounds(self):
        infile = data_dir / "nested_snps_seq_backgrounds.fa"
        aseq = PrgBuilder(infile, min_match_length=3)
        self.assertEqual(
            aseq.prg, " 5 AAAA 7 T 8 C 7 AAAAAA 6 CCCC 9 T 10 G 9 CCCCCC 5 "
        )

    def test_nested_snps_under_del(self):
        infile = data_dir / "nested_snps_deletion.fa"
        aseq = PrgBuilder(infile, min_match_length=1)
        self.assertEqual(aseq.prg, "A 5 AA 7 C 8 T 7 AAAA 9 T 10 G 9 AA 6 A 5 AA")
