from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_prg import MSA


def make_alignment(seqs: List[str], ids: List[str] = None) -> MSA:
    if ids is None:
        seqrecords = [SeqRecord(Seq(seq), id=f"s_{i}") for i, seq in enumerate(seqs)]
    else:
        seqrecords = [SeqRecord(Seq(seq), id=ID) for seq, ID in zip(seqs, ids)]
    return MSA(seqrecords)
