import filecmp
from pathlib import Path
from unittest import TestCase

from make_prg.utils.gfa import GFA_Output
from tests.test_helpers import sample_prg

workdir = Path("tests/data/utils/gfa")


class Test_GFA_Output(TestCase):
    """
    TODO: this should not be a big bang test, but several small unit tests.
    TODO: refactor this by adding several smaller unit tests.
    TODO: this is a regression test, need to check if the truth GFA is correct.
    """

    def test___big_bang(self):
        GFA_Output.write_gfa(str(workdir / "big_bang"), sample_prg)
        self.assertTrue(
            filecmp.cmp(workdir / "big_bang.gfa", workdir / "big_bang.truth.gfa")
        )
