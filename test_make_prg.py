#!/usr/bin/env python3

import pytest
from make_prg_from_msa import AlignedSeq


def test_answers():
    aseq = AlignedSeq("test/match.fa")
    assert aseq.prg == "ACGTGTTTTGTAACTGTGCCACACTCTCGAGACTGCATATGTGTC"

    aseq = AlignedSeq("test/nonmatch.fa")
    assert aseq.prg == " 5 AAACGTGGTT 6 CCCCCCCCCC 5 "

    aseq = AlignedSeq("test/match.nonmatch.fa")
    assert aseq.prg == "AAACG 5 TGGTT 6 CCCCC 5 "

    aseq = AlignedSeq("test/nonmatch.match.fa")
    assert aseq.prg == " 5 AAACGT 6 CCCCCC 5 GGTT"

    aseq = AlignedSeq("test/match.nonmatch.match.fa")
    assert aseq.prg == "AAACG 5 T 6 C 5 GGTT"

    aseq = AlignedSeq("test/shortmatch.nonmatch.match.fa")
    assert aseq.prg == " 5 AAACGT 6 ATTTTC 5 GGTT"

    aseq = AlignedSeq("test/match.nonmatch.shortmatch.fa")
    assert aseq.prg == "AAAC 5 GTGGTT 6 CCCCCT 5 "

    aseq = AlignedSeq("test/match.staggereddash.fa")
    assert aseq.prg == "AAACGTGGTT"

    aseq = AlignedSeq("test/contains_n.fa")
    assert aseq.prg == "AAACG 5 T 6 C 5 GGTT"

    aseq = AlignedSeq("test/contains_RYKMSW.fa")
    assert aseq.prg == "AAACG 5 T 6 C 5 GGTT"

    aseq = AlignedSeq("test/contains_n_and_RYKMSW.fa")
    assert aseq.prg == "AAACG 5 T 6 C 5 GGTT"

    aseq = AlignedSeq("test/contains_n_and_RYKMSW_no_variants.fa")
    assert aseq.prg == "AAACGTGGTT"

    with pytest.raises(Exception):
        aseq = AlignedSeq("test/fails.fa")
    print("Done")