"""
Code responsible for converting a consensus string into a set of disjoint
match/non_match intervals.
"""
from enum import Enum, auto
from typing import List, Tuple, Optional

from make_prg.from_msa import MSA

from make_prg.seq_utils import get_expanded_sequences, is_non_match, has_empty_sequence


class PartitioningError(Exception):
    pass


class IntervalType(Enum):
    Match = auto()
    NonMatch = auto()

    @classmethod
    def from_char(cls, letter: str) -> "IntervalType":
        if letter == "*":
            return IntervalType.NonMatch
        else:
            return IntervalType.Match


def is_type(letter: str, interval_type: IntervalType) -> bool:
    if IntervalType.from_char(letter) is interval_type:
        return True
    else:
        return False


class Interval:
    """Stores a closed interval [a,b]"""

    def __init__(self, it_type: IntervalType, start: int, stop: int = None):
        self.type = it_type
        self.start = start
        if stop is not None:
            assert stop >= start
        self.stop = stop if stop is not None else start

    def modify_by(self, left_delta: int, right_delta: int):
        self.start += left_delta
        self.stop += right_delta

    def contains(self, position: int):
        return self.start <= position <= self.stop

    def __len__(self) -> int:
        return self.stop - self.start + 1

    def __lt__(self, other: "Interval") -> bool:
        return self.start < other.start

    def __eq__(self, other: "Interval") -> bool:
        return (
            self.start == other.start
            and self.stop == other.stop
            and self.type is other.type
        )

    def __repr__(self):
        return f"[{self.start}, {self.stop}]"


Intervals = List[Interval]


class IntervalPartitioner:
    """Produces a list of intervals in which we have
    consensus sequence longer than min_match_length, and
    a list of the non-match intervals left."""

    def __init__(self, consensus_string: str, min_match_length: int, alignment: MSA):
        self._match_intervals: Intervals = list()
        self._non_match_intervals: Intervals = list()
        self.mml = min_match_length

        if len(consensus_string) < self.mml:
            # In this case, a match of less than the min_match_length gets counted
            # as a match (usually, it counts as a non_match)
            it_type = IntervalType.Match
            if any(map(is_non_match, consensus_string)):
                it_type = IntervalType.NonMatch
            self._append(Interval(it_type, 0, len(consensus_string) - 1))

        else:
            cur_interval = self._new_interval(consensus_string[0], 0)

            for i, letter in enumerate(consensus_string[1:], start=1):
                if is_type(letter, cur_interval.type):
                    cur_interval.modify_by(0, 1)  # simple interval extension
                else:
                    new_interval = self._add_interval(cur_interval, alignment)
                    if new_interval is None:
                        cur_interval = self._new_interval(letter, i)
                    else:
                        cur_interval = new_interval
            self._add_interval(cur_interval, alignment, end=True)

        self.enforce_multisequence_nonmatch_intervals(
            self._match_intervals, self._non_match_intervals, alignment
        )
        self.enforce_alignment_interval_bijection(
            self._match_intervals,
            self._non_match_intervals,
            alignment.get_alignment_length(),
        )

    def get_intervals(self) -> Tuple[Intervals, Intervals, Intervals]:
        return (
            sorted(self._match_intervals),
            sorted(self._non_match_intervals),
            sorted(self._match_intervals + self._non_match_intervals),
        )

    def _new_interval(self, letter: str, start_pos: int) -> Interval:
        return Interval(IntervalType.from_char(letter), start_pos)

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

    def _add_interval(
        self, interval: Interval, alignment: MSA, end: bool = False
    ) -> Optional[Interval]:
        """
        i)If we are given a match interval < min_match_length, we return an extended non_match interval
        ii)If we are given a non_match interval containing 1+ empty sequence, we pad it with
           previous match_interval, if any, to avoid empty alleles in resulting prg.
        """
        if interval.type is IntervalType.Match:
            # The +1 is because we also extend the non_match interval
            if len(interval) < self.mml:
                try:
                    last_non_match = self._pop(IntervalType.NonMatch)
                    last_non_match.modify_by(0, len(interval) + 1)
                except IndexError:
                    last_non_match = Interval(
                        IntervalType.NonMatch, interval.start, interval.stop + 1
                    )
                if end:  # If this is final call, go to append the interval
                    last_non_match.modify_by(0, -1)
                    self._append(last_non_match)
                return last_non_match
        else:
            if len(self._match_intervals) > 0 and has_empty_sequence(
                alignment, (interval.start, interval.stop)
            ):
                # Pad interval with sequence to avoid empty alleles
                len_match = len(self._match_intervals[-1])
                if len_match - 1 < self.mml:
                    # Case: match is now too small, converted to non_match
                    self._match_intervals.pop()
                    interval.modify_by(-1 * len_match, 0)
                    if len(self._non_match_intervals) > 0:
                        # Case: merge previous non_match with this non_match
                        self._non_match_intervals[-1].modify_by(0, len(interval))
                        return None
                else:
                    self._match_intervals[-1].modify_by(0, -1)
                    interval.modify_by(-1, 0)

        self._append(interval)
        return None

    @classmethod
    def enforce_multisequence_nonmatch_intervals(
        cls, match_intervals: Intervals, non_match_intervals: Intervals, alignment: MSA
    ) -> None:
        """
        Goes through non-match intervals and makes sure there is more than one sequence there, else makes it a match
        interval.
        Modifies the intervals in-place.
        Example reasons for such a conversion to occur:
            - 'N' in a sequence causes it to be filtered out, and left with a single useable sequence
            - '-' in sequences causes them to appear different, but they are the same
        """
        if len(alignment) == 0:  # For testing convenience
            return
        for i in reversed(range(len(non_match_intervals))):
            interval = non_match_intervals[i]
            interval_alignment = alignment[:, interval.start : interval.stop + 1]
            interval_seqs = get_expanded_sequences(interval_alignment)
            if len(interval_seqs) < 2:
                changed_interval = non_match_intervals[i]
                match_intervals.append(
                    Interval(
                        IntervalType.Match,
                        changed_interval.start,
                        changed_interval.stop,
                    )
                )
                non_match_intervals.pop(i)

    @classmethod
    def enforce_alignment_interval_bijection(
        cls,
        match_intervals: Intervals,
        non_match_intervals: Intervals,
        alignment_length: int,
    ):
        """
        Check each position in an alignment is in one, and one only, (match or non_match) interval
        """
        for i in range(alignment_length):
            count_match = 0
            for interval in match_intervals:
                if interval.contains(i):
                    count_match += 1
            count_non_match = 0
            for interval in non_match_intervals:
                if interval.contains(i):
                    count_non_match += 1

            if count_match > 1 or count_non_match > 1:
                raise PartitioningError(
                    f"Failed interval partitioning: position {i}"
                    " appears in more than one interval"
                )
            if (
                not count_match ^ count_non_match
            ):  # test fails if they are the same integer
                msg = ["neither", "nor"] if count_match == 0 else ["both", "and"]
                raise PartitioningError(
                    "Failed interval partitioning: alignment position %d"
                    "classified as %s match %s non-match " % (i, msg[0], msg[1])
                )
