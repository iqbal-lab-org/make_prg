from pathlib import Path
from unittest import TestCase

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_prg.from_msa.prg_builder import PrgBuilder
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
        cls.alignment = MSA(
            [
                SeqRecord(Seq("AAAT"), id="s1"),
                SeqRecord(Seq("C--C"), id="s2"),
                SeqRecord(Seq("AATT"), id="s3"),
                SeqRecord(Seq("GNGG"), id="s4"),
            ]
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
        expected = MSA(
            [SeqRecord(Seq("C--"), id="s2"), SeqRecord(Seq("AAT"), id="s3"),]
        )
        self.assertTrue(msas_equal(expected, result))


class TestMakePrgFromMsaFile_IntegrationTests(TestCase):
    def test_answers(self):
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
