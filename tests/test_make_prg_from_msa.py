import os
from unittest import TestCase, mock

from make_prg.make_prg_from_msa import AlignedSeq

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
