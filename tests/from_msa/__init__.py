from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_prg.from_msa import MSA


def make_alignment(seqs: List[str], ids: List[str] = None) -> MSA:
    seq_lengths = set(map(len, seqs))
    assert (
        len(seq_lengths) == 1
    ), "Sequences are not the same length, does not represent an alignment"
    if ids is None:
        seqrecords = [SeqRecord(Seq(seq), id=f"s{i}") for i, seq in enumerate(seqs)]
    else:
        seqrecords = [SeqRecord(Seq(seq), id=ID) for seq, ID in zip(seqs, ids)]
    return MSA(seqrecords)
