import logging
from typing import Generator, Sequence
import itertools

from Bio import AlignIO

from make_prg import MSA


class PartitioningError(Exception):
    pass


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
            "WARNING: Every sequence must have contained an N in this slice - redo sequence curation because this is nonsense"
        )
        logging.warning(f'Sequences were: {" ".join(callback_seqs)}')
        logging.warning(
            "Using these sequences anyway, and should be ignored downstream"
        )
        return callback_seqs
    return expanded_seqs


def enforce_multisequence_nonmatch_intervals(
    match_intervals, non_match_intervals, alignment: MSA
):
    """
    Goes through non-match intervals and makes sure there is more than one sequence there, else makes it a match
    interval.
    Example reasons for such a conversion to occur:
        - 'N' in a sequence causes it to be filtered out, and left with a single useable sequence
        - '-' in sequences causes them to appear different, but they are the same
    """
    for i in reversed(range(len(non_match_intervals))):
        interval = non_match_intervals[i]
        interval_alignment = alignment[:, interval[0] : interval[1] + 1]
        interval_seqs = get_interval_seqs(interval_alignment)
        if len(interval_seqs) < 2:
            match_intervals.append(non_match_intervals[i])
            non_match_intervals.pop(i)
    match_intervals.sort()


def enforce_alignment_interval_bijection(
    match_intervals, non_match_intervals, alignment_length: int
):
    """
    Check each position in an alignment is in one, and one only, (match or non_match) interval
    """
    for i in range(alignment_length):
        count_match = 0
        for interval in match_intervals:
            if interval[0] <= i <= interval[1]:
                count_match += 1
        count_non_match = 0
        for interval in non_match_intervals:
            if interval[0] <= i <= interval[1]:
                count_non_match += 1

        if count_match > 1 or count_non_match > 1:
            raise PartitioningError(
                f"Failed interval partitioning: position {i}"
                " appears in more than one interval"
            )
        if not count_match ^ count_non_match:  # test fails if they are the same integer
            msg = ["neither", "nor"] if count_match == 0 else ["both", "and"]
            raise PartitioningError(
                "Failed interval partitioning: alignment position %d"
                "classified as %s match %s non-match " % (i, *msg)
            )
