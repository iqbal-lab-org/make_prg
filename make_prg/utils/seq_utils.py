import copy
import itertools
from typing import Generator, List, Tuple

import numpy as np
from Bio import pairwise2
from Bio.Seq import Seq

from make_prg import MSA

NONMATCH = "*"
GAP = "-"
Sequence = str
Sequences = List[str]


def is_non_match(letter: str):
    return letter == NONMATCH


def is_gap(letter: str):
    return letter == GAP


def remove_duplicates(seqs: Sequences) -> Generator:
    seen = set()
    for x in seqs:
        if x in seen:
            continue
        seen.add(x)
        yield x


def has_empty_sequence(alignment: MSA, interval: Tuple[int, int]) -> bool:
    sub_alignment = alignment[:, interval[0] : interval[1] + 1]
    for record in sub_alignment:
        if all(map(is_gap, record.seq)):
            return True
    return False


def ungap(seq: str) -> str:
    return seq.replace(GAP, "")


def count(iterable) -> int:
    return sum(1 for _ in iterable)


def get_alignment_seqs(alignment: MSA) -> Generator:
    for record in alignment:
        yield str(record.seq)


def get_number_of_unique_ungapped_sequences(sub_alignment: MSA) -> int:
    alignment_seqs = get_alignment_seqs(sub_alignment)
    ungapped_sequences = map(ungap, alignment_seqs)
    deduplicated_ungapped_sequences = remove_duplicates(ungapped_sequences)
    number_of_unique_nongapped_sequences = count(deduplicated_ungapped_sequences)
    return number_of_unique_nongapped_sequences


def get_number_of_unique_gapped_sequences(sub_alignment: MSA) -> int:
    alignment_seqs = get_alignment_seqs(sub_alignment)
    deduplicated_gapped_sequences = remove_duplicates(alignment_seqs)
    number_of_unique_gapped_sequences = count(deduplicated_gapped_sequences)
    return number_of_unique_gapped_sequences


class SequenceCurationError(Exception):
    pass


class SequenceExpander:
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
    expandable_bases = set(iupac.keys())
    allowed_bases = expandable_bases.union({"N"})
    standard_bases = {"A", "C", "G", "T"}
    ambiguous_bases = expandable_bases.difference(standard_bases)

    @classmethod
    def check_if_there_is_sequence_with_disallowed_bases(cls, sequences: List[str]):
        for sequence in sequences:
            if not set(sequence).issubset(cls.allowed_bases):
                raise SequenceCurationError(
                    "A slice of a sequence has a disallowed base.\n"
                    f"Allowed bases: {cls.allowed_bases}.\n"
                    f"Sequence: {sequence}\n"
                    f"Redo sequence curation.\n"
                )

    @classmethod
    def check_all_sequences_are_composed_of_ACGT(cls, sequences: List[str]):
        for sequence in sequences:
            sequence_is_composed_of_ACGT_only = set(sequence).issubset(
                cls.standard_bases
            )
            assert (
                sequence_is_composed_of_ACGT_only
            ), f"Sequence ({sequence}) should be composed of ACTG only."

    @classmethod
    def get_expanded_sequences(cls, sequences: List[str]) -> Sequences:
        """
        Expand sequences in the given list of sequences, following the translation
        table in SequenceExpander.iupac.
        It does the following steps:
            1. Check that we don't have disallowed bases;
            2. Remove sequences with N;
            3. Duplicate sequences containing RYKMSW, replacing with AGCT alternatives;
        Note 1: The sequences are deliberately returned in the order they are received.
        Note 2: Returned sequences are composed of ACGT only
        """

        cls.check_if_there_is_sequence_with_disallowed_bases(sequences)

        expanded_seqs = []
        expanded_set = set()

        for seq in remove_duplicates(sequences):
            if "N" in seq:
                continue

            alternatives = [cls.iupac[base] for base in seq]
            for tuple_product in itertools.product(*alternatives):
                expanded_str = "".join(tuple_product)
                if expanded_str not in expanded_set:
                    expanded_set.add(expanded_str)
                    expanded_seqs.append(expanded_str)

        all_sequences_contained_N = len(expanded_seqs) == 0
        if all_sequences_contained_N:
            raise SequenceCurationError(
                "All sequences in this slice contained N. Redo sequence curation.\n"
                f"Sequences: {sequences}"
            )

        cls.check_all_sequences_are_composed_of_ACGT(expanded_seqs)
        return expanded_seqs

    @classmethod
    def get_expanded_sequences_from_MSA(cls, alignment: MSA) -> Sequences:
        gapless_seqs = list(map(ungap, get_alignment_seqs(alignment)))
        return cls.get_expanded_sequences(gapless_seqs)


def align(
    seq1: str,
    seq2: str,
    match_score: float = 2,
    mismatch_score: float = -0.9,
    gap_open_score: float = -1.1,
    gap_extend_score: float = -1,
) -> Tuple[str, str]:

    seq1_is_empty = len(seq1) == 0
    seq2_is_empty = len(seq2) == 0

    if seq1_is_empty:
        return GAP * len(seq2), seq2
    elif seq2_is_empty:
        return seq1, GAP * len(seq1)
    else:
        # TODO: replace alignment library from BioPython to e.g. edlib due to
        #  performance reasons?
        alignment = pairwise2.align.globalms(
            seq1,
            seq2,
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score,
            one_alignment_only=True,
        )
        alignment = alignment[0]
        return alignment.seqA, alignment.seqB


def remove_columns_full_of_gaps_from_MSA(alignment: MSA) -> MSA:
    """
    Note: This code is long and a bit convoluted because it is optimised (it was too
    slow if done in the most intuitive
    way).
    """
    alignment_as_array = np.array([list(rec) for rec in alignment], str, order="F")
    gapless_sequences = [[] for _ in range(len(alignment))]
    for column_index in range(alignment.get_alignment_length()):
        column_bases = alignment_as_array[:, column_index]
        column_bases_deduplicated = list(set(column_bases))
        just_gaps = column_bases_deduplicated == [GAP]
        if not just_gaps:
            for gapless_sequence, base in zip(gapless_sequences, column_bases):
                gapless_sequence.append(base)

    gapless_records = []
    for gapless_sequence, previous_record in zip(gapless_sequences, alignment):
        new_record = copy.deepcopy(previous_record)
        new_record.seq = Seq("".join(gapless_sequence))
        gapless_records.append(new_record)

    gapless_alignment = MSA(gapless_records)
    return gapless_alignment


def get_consensus_from_MSA(alignment: MSA) -> str:
    """Produces a 'consensus string' from an MSA: at each position of the
    MSA, the string has a base if all aligned sequences agree, and a "*" if not.
    IUPAC ambiguous bases result in non-consensus and are later expanded in the prg.
    N results in consensus at that position unless they are all N."""
    consensus_string_as_list = []
    for i in range(alignment.get_alignment_length()):
        column = set([record.seq[i] for record in alignment])
        column = column.difference({"N"})
        column_is_ambiguous = (
            len(SequenceExpander.ambiguous_bases.intersection(column)) > 0
            or len(column) != 1
        )
        column_is_all_gaps = column == {"-"}
        not_a_match = column_is_ambiguous or column_is_all_gaps
        if not_a_match:
            consensus_string_as_list.append(NONMATCH)
        else:
            consensus_string_as_list.append(column.pop())
    consensus_string = "".join(consensus_string_as_list)
    return consensus_string
