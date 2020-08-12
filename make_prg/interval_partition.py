"""
Code responsible for converting a consensus string into a set of disjoint
match/non_match intervals.
"""
from enum import Enum, auto
from typing import List, Tuple

from make_prg import MSA

from make_prg.seq_utils import get_interval_seqs


class IntervalType(Enum):
    Match = auto()
    NonMatch = auto()


def get_type(letter: str) -> IntervalType:
    if letter == "*":
        return IntervalType.NonMatch
    else:
        return IntervalType.Match


def is_type(letter: str, interval_type: IntervalType) -> bool:
    if get_type(letter) is interval_type:
        return True
    else:
        return False


class Interval:
    def __init__(self, it_type: IntervalType, start: int, stop: int = None):
        self.type = it_type
        self.start = start
        self.stop = stop if stop is not None else start

    def modify_by(self, left_delta: int, right_delta: int):
        self.start += left_delta
        self.stop += right_delta

    def __len__(self) -> int:
        return self.stop - self.start + 1

    def __lt__(self, other: "Interval") -> bool:
        return self.start < other.start


Intervals = List[Interval]


class IntervalPartitioner:
    """Produces a list of intervals in which we have
    consensus sequence longer than min_match_length, and
    a list of the non-match intervals left."""

    def __init__(self, consensus_string: str, min_match_length: int, alignment: MSA):
        self._match_intervals: Intervals = list()
        self._non_match_intervals: Intervals = list()
        self.mml = min_match_length

        cur_interval = self._new_interval(consensus_string[0], 0)

        for i, letter in enumerate(consensus_string[1:], start=1):
            if is_type(letter, cur_interval.type):
                cur_interval.modify_by(0, 1)  # simple interval extension
            else:
                self._add_interval(cur_interval, alignment)
                cur_interval = self._new_interval(letter, i)
        self._add_interval(cur_interval, alignment)

    def get_intervals(self) -> Tuple[Intervals, Intervals, Intervals]:
        return (
            sorted(self._match_intervals),
            sorted(self._non_match_intervals),
            sorted(self._match_intervals + self._non_match_intervals),
        )

    def _new_interval(self, letter: str, start_pos: int) -> Interval:
        return Interval(get_type(letter), start_pos)

    def _append(self, interval: Interval):
        if interval.type is IntervalType.Match:
            self._match_intervals.append(interval)
        else:
            self._non_match_intervals.append(interval)

    def _pop(self, it_type: IntervalType) -> Interval:
        if it_type is IntervalType.Match:
            return self._match_intervals.pop()
        else:
            return self._non_match_intervals.pop()

    def _add_interval(self, interval: Interval, alignment: MSA):
        if interval.type is IntervalType.Match:
            if len(interval) < self.mml:
                try:
                    last_non_match = self._pop(IntervalType.NonMatch)
                    last_non_match.modify_by(0, len(interval))
                except IndexError:
                    last_non_match = Interval(
                        interval.start, interval.stop, IntervalType.NonMatch
                    )
                self.append(last_non_match)
                return
        else:
            pass
        self._append(interval)


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


class PartitioningError(Exception):
    pass


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
