from unittest import TestCase
from typing import List

from tests import make_alignment, MSA
from make_prg.interval_partition import (
    IntervalType,
    Interval,
    IntervalPartitioner,
    PartitioningError,
)

Lists = List[List[int]]

Match = IntervalType.Match
NonMatch = IntervalType.NonMatch


def make_typed_intervals(lists: Lists, it_type: IntervalType):
    result = list()
    for elem in lists:
        result.append(Interval(it_type, elem[0], elem[1]))
    return result


def make_intervals(match_lists: Lists, non_match_lists: Lists):
    return (
        make_typed_intervals(match_lists, Match),
        make_typed_intervals(non_match_lists, NonMatch),
    )


class TestIntervalConsistency(TestCase):
    def test_nonmatch_interval_switching_indels(self):
        """Because the sequences are the same, despite different alignment"""
        alignment = make_alignment(["A---A", "A-A--"])
        match_intervals, non_match_intervals = make_intervals([], [[0, 5]])
        IntervalPartitioner.enforce_multisequence_nonmatch_intervals(
            match_intervals, non_match_intervals, alignment
        )
        self.assertEqual(match_intervals, make_typed_intervals([[0, 5]], Match))
        self.assertEqual(non_match_intervals, [])

    def test_nonmatch_interval_switching_Ns(self):
        """'N's make sequences get removed"""
        alignment = make_alignment(["ANAAA", "ATAAT"])
        match_intervals, non_match_intervals = make_intervals([], [[0, 5]])
        IntervalPartitioner.enforce_multisequence_nonmatch_intervals(
            match_intervals, non_match_intervals, alignment
        )
        self.assertEqual(match_intervals, make_typed_intervals([[0, 5]], Match))
        self.assertEqual(non_match_intervals, [])

    def test_position_in_several_intervals_fails(self):
        match_intervals = make_typed_intervals([[0, 1], [1, 2]], Match)
        with self.assertRaises(PartitioningError):
            IntervalPartitioner.enforce_alignment_interval_bijection(
                match_intervals, [], 3
            )

    def test_position_in_no_interval_fails(self):
        match_intervals = make_typed_intervals([[0, 1]], Match)
        with self.assertRaises(PartitioningError):
            IntervalPartitioner.enforce_alignment_interval_bijection(
                match_intervals, [], 3
            )

    def test_position_in_match_and_nonmatch_intervals_fails(self):
        match_intervals, nmatch_intervals = make_intervals([[0, 2]], [[2, 3]])
        with self.assertRaises(PartitioningError):
            IntervalPartitioner.enforce_alignment_interval_bijection(
                match_intervals, nmatch_intervals, 4
            )

    def test_bijection_respected_passes(self):
        match_intervals, nmatch_intervals = make_intervals([[0, 2], [5, 10]], [[3, 4]])
        IntervalPartitioner.enforce_alignment_interval_bijection(
            match_intervals, nmatch_intervals, 11
        )


class TestIntervalPartitioning(TestCase):
    def test_all_non_match(self):
        tester = IntervalPartitioner("******", min_match_length=3, alignment=MSA([]))
        match, non_match, _ = tester.get_intervals()
        self.assertEqual(match, [])
        self.assertEqual(non_match, make_typed_intervals([[0, 5]], NonMatch))

    def test_all_match(self):
        tester = IntervalPartitioner("ATATAAA", min_match_length=3, alignment=MSA([]))
        match, non_match, _ = tester.get_intervals()
        self.assertEqual(match, make_typed_intervals([[0, 6]], Match))
        self.assertEqual(non_match, [])

    def test_short_match_counted_as_non_match(self):
        tester = IntervalPartitioner("AT***", min_match_length=3, alignment=MSA([]))
        match, non_match, _ = tester.get_intervals()
        self.assertEqual(match, [])
        self.assertEqual(non_match, make_typed_intervals([[0, 4]], NonMatch))

    def test_match_non_match_match(self):
        tester = IntervalPartitioner("ATT**AAAC", min_match_length=3, alignment=MSA([]))
        match, non_match, _ = tester.get_intervals()
        self.assertEqual(match, make_typed_intervals([[0, 2], [5, 8]], Match))
        self.assertEqual(non_match, make_typed_intervals([[3, 4]], NonMatch))

    def test_end_in_non_match(self):
        tester = IntervalPartitioner(
            "**ATT**AAA*C", min_match_length=3, alignment=MSA([])
        )
        match, non_match, _ = tester.get_intervals()
        self.assertEqual(match, make_typed_intervals([[2, 4], [7, 9]], Match))
        self.assertEqual(
            non_match, make_typed_intervals([[0, 1], [5, 6], [10, 11]], NonMatch)
        )

    def test_consensus_smaller_than_min_match_len(self):
        """
        Usually, a match smaller than min_match_length counts as non-match,
        but if the whole string is smaller than min_match_length, counts as match.
        """
        tester1 = IntervalPartitioner("TTATT", min_match_length=7, alignment=MSA([]))
        match, non_match, _ = tester1.get_intervals()
        self.assertEqual(match, make_typed_intervals([[0, 4]], Match))
        self.assertEqual(non_match, [])

        tester2 = IntervalPartitioner("T*ATT", min_match_length=7, alignment=MSA([]))
        match, non_match, _ = tester2.get_intervals()
        self.assertEqual(match, [])
        self.assertEqual(non_match, make_typed_intervals([[0, 4]], NonMatch))
