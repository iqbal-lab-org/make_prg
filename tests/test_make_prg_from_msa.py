import os
import unittest

from make_prg import make_prg_from_msa

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, "tests", "data", "make_prg_from_msa")


class TestMakePrgFromMsa(unittest.TestCase):
    def test_answers(self):
        infile = os.path.join(data_dir, "match.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, "ACGTGTTTTGTAACTGTGCCACACTCTCGAGACTGCATATGTGTC")

        infile = os.path.join(data_dir, "nonmatch.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, " 5 AAACGTGGTT 6 CCCCCCCCCC 5 ")

        infile = os.path.join(data_dir, "match.nonmatch.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 TGGTT 6 CCCCC 5 ")

        infile = os.path.join(data_dir, "nonmatch.match.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, " 5 AAACGT 6 CCCCCC 5 GGTT")

        infile = os.path.join(data_dir, "match.nonmatch.match.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 C 6 T 5 GGTT")

        infile = os.path.join(data_dir, "shortmatch.nonmatch.match.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, " 5 AAACGT 6 ATTTTC 5 GGTT")

        infile = os.path.join(data_dir, "match.nonmatch.shortmatch.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAAC 5 GTGGTT 6 CCCCCT 5 ")

        infile = os.path.join(data_dir, "match.staggereddash.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACGTGGTT")

        infile = os.path.join(data_dir, "contains_n.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 C 6 T 5 GGTT")

        infile = os.path.join(data_dir, "contains_RYKMSW.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 C 6 T 5 GGTT")

        infile = os.path.join(data_dir, "contains_n_and_RYKMSW.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACG 5 C 6 T 5 GGTT")

        infile = os.path.join(data_dir, "contains_n_and_RYKMSW_no_variants.fa")
        aseq = make_prg_from_msa.AlignedSeq(infile)
        self.assertEqual(aseq.prg, "AAACGTGGTT")

        # with pytest.raises(Exception):
        #    aseq = make_prg_from_msa.AlignedSeq("test/fails.fa")
