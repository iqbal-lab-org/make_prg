import logging
from typing import Generator, Sequence
import itertools

from Bio import AlignIO


def remove_duplicates(seqs: Sequence) -> Generator:
    seen = set()
    for x in seqs:
        if x in seen:
            continue
        seen.add(x)
        yield x


def remove_gaps(sequence: str) -> str:
    return sequence.replace("-", "")


iupac = {
    "R": "GA",
    "Y": "TC",
    "K": "GT",
    "M": "AC",
    "S": "GC",
    "W": "AT",
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
}
allowed_bases = set(iupac.keys())


def get_interval_seqs(interval_alignment: AlignIO.MultipleSeqAlignment):
    """Replace - with nothing, remove seqs containing N or other non-allowed letters
    and duplicate sequences containing RYKMSW, replacing with AGCT alternatives """
    gapless_seqs = [str(record.seq.ungap("-")) for record in interval_alignment]

    callback_seqs = []
    expanded_set = set()
    for seq in remove_duplicates(gapless_seqs):
        if len(expanded_set) == 0:
            callback_seqs.append(seq)
        if not set(seq).issubset(allowed_bases):
            continue
        alternatives = [iupac[base] for base in seq]
        for expanded_seq in itertools.product(*alternatives):
            expanded_set.add("".join(expanded_seq))

    if len(expanded_set) == 0:
        logging.warning(
            "WARNING: Every sequence must have contained an N in this slice - redo sequence curation because this is nonsense"
        )
        logging.warning(f'Sequences were: {" ".join(callback_seqs)}')
        logging.warning(
            "Using these sequences anyway, and should be ignored downstream"
        )
        return sorted(callback_seqs)
    return sorted(list(expanded_set))
