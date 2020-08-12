import logging
from typing import Generator, Sequence
import itertools

from Bio import AlignIO

NONMATCH = "*"


def is_non_match(letter: str):
    return letter == NONMATCH


def remove_duplicates(seqs: Sequence) -> Generator:
    seen = set()
    for x in seqs:
        if x in seen:
            continue
        seen.add(x)
        yield x


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
standard_bases = {"A", "C", "G", "T"}
ambiguous_bases = allowed_bases.difference(standard_bases)


def get_interval_seqs(interval_alignment: AlignIO.MultipleSeqAlignment):
    """
    Replace - with nothing, remove seqs containing N or other non-allowed letters
    and duplicate sequences containing RYKMSW, replacing with AGCT alternatives

    The sequences are deliberately returned in the order they are received
    """
    gapless_seqs = [str(record.seq.ungap("-")) for record in interval_alignment]

    callback_seqs, expanded_seqs = [], []
    expanded_set = set()
    for seq in remove_duplicates(gapless_seqs):
        if len(expanded_set) == 0:
            callback_seqs.append(seq)
        if not set(seq).issubset(allowed_bases):
            continue
        alternatives = [iupac[base] for base in seq]
        for tuple_product in itertools.product(*alternatives):
            expanded_str = "".join(tuple_product)
            if expanded_str not in expanded_set:
                expanded_set.add(expanded_str)
                expanded_seqs.append(expanded_str)

    if len(expanded_set) == 0:
        logging.warning(
            "WARNING: Every sequence must have contained an N in this slice - redo sequence curation"
        )
        logging.warning(f'Sequences were: {" ".join(callback_seqs)}')
        logging.warning(
            "Using these sequences anyway, and should be ignored downstream"
        )
        return callback_seqs
    return expanded_seqs
