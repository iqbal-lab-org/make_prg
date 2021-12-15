from unittest import TestCase
from make_prg.prg_builder import PrgBuilder
from make_prg.utils.recursive_tree_drawer import RecursiveTreeDrawer
import filecmp
from pathlib import Path

workdir = Path("tests/data/utils/recursive_tree_drawer")


class Test_RecursiveTreeDrawer(TestCase):
    """
    TODO: it might not be needed to test RecursiveTreeDrawer, the recursive tree draw is used solely for debug purposes.
    TODO: if tests are needed, see the TODOs below.
    TODO: this should not be a big bang test, but several small unit tests.
    TODO: refactor this by adding several smaller unit tests.
    TODO: this is a regression test, need to check if the truth recursive tree drawing is correct.
    """
    def test___big_bang(self):
        prg_builder = PrgBuilder(locus_name="sample", msa_file=workdir / "sample_msa.fa", alignment_format="fasta",
                   max_nesting=5, min_match_length=7)
        recursive_tree_drawer = RecursiveTreeDrawer(prg_builder.root)
        recursive_tree_drawer.output_graph(workdir / "sample.recursion_tree.png")
        self.assertTrue(filecmp.cmp(workdir / "sample.recursion_tree.png",
                                    workdir / "sample.recursion_tree.truth.png"))
