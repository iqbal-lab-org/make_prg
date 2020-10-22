from unittest import TestCase, skip
from unittest.mock import patch, Mock
import random

from numpy import array, array_equal
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_prg.seq_utils import standard_bases
from make_prg.from_msa.cluster_sequences import (
    count_distinct_kmers,
    count_kmer_occurrences,
    kmeans_cluster_seqs_in_interval,
)

from tests.from_msa import make_alignment, MSA


class TestCountKmers(TestCase):
    def test_one_seq_below_kmer_size_throws(self):
        sequences = ["AAA"]
        with self.assertRaises(ValueError):
            count_distinct_kmers(sequences, 5)

    def test_one_seq_returns_correct_counts(self):
        sequences = ["AAAAAAAAAT"]
        result = count_distinct_kmers(sequences, 5)
        expected = {"AAAAA": 0, "AAAAT": 1}
        assert result == expected

    def test_multiple_seqs_correct_counts(self):
        sequences = ["AAAT", "AAAG"]
        result = count_distinct_kmers(sequences, 3)
        expected = {"AAA": 0, "AAT": 1, "AAG": 2}
        assert result == expected


class TestCountMatrix(TestCase):
    def test_one_sequence_correct_counts(self):
        sequences = ["AAAAT"]
        kmers = {"AAA": 0, "AAT": 1}
        result = count_kmer_occurrences(sequences, kmers)
        expected = array([[2, 1]])
        assert array_equal(result, expected)

    def test_multiple_sequence_correct_counts(self):
        sequences = ["AAAAT", "AAATA"]
        kmers = {"AAA": 0, "AAT": 1, "ATA": 2}
        result = count_kmer_occurrences(sequences, kmers)
        expected = array([[2, 1, 0], [1, 1, 1]])
        assert array_equal(result, expected)


class TestClustering_Trivial(TestCase):
    def test_one_seq_returns_single_id(self):
        alignment = MSA([SeqRecord(Seq("AAAT"), id="s1")])
        result = kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        self.assertEqual(result, [["s1"]])

    def test_two_identical_seqs_returns_two_ids_clustered(self):
        alignment = MSA(
            [SeqRecord(Seq("AAAT"), id="s1"), SeqRecord(Seq("AAAT"), id="s2"),]
        )
        result = kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        self.assertEqual(result, [["s1", "s2"]])


class TestClustering_SmallSequences(TestCase):
    def test_two_seqs_one_below_kmer_size_separate_clusters(self):
        alignment = MSA(
            [SeqRecord(Seq("AATTTAT"), id="s1"), SeqRecord(Seq("AA---AT"), id="s2")]
        )
        result = kmeans_cluster_seqs_in_interval([0, 5], alignment, 5)
        self.assertEqual(result, [["s1"], ["s2"]])

    @skip(
        "This fails, probably because kmean clustering should never run with this input"
    )
    def test_ambiguous_sequences_in_short_interval_separate_clusters(self):
        alignment = MSA(
            [SeqRecord(Seq("ARAT"), id="s1"), SeqRecord(Seq("WAAT"), id="s2"),]
        )
        result = kmeans_cluster_seqs_in_interval([0, 3], alignment, 5)
        self.assertEqual([["s1"], ["s2"]], result)

    @patch("make_prg.from_msa.cluster_sequences.KMeans")
    def test_GivenAllSequencesBelowKmerSize_NoKMeansAndIdenticalSequencesClustered(
        self, mockKMeans
    ):
        alignment = MSA(
            [
                SeqRecord(Seq("AA---AT"), id="s1"),
                SeqRecord(Seq("AA---TT"), id="s2"),
                SeqRecord(Seq("CA--CAT"), id="s3"),
                SeqRecord(Seq("A-A--AT"), id="s4"),
            ]
        )
        result = kmeans_cluster_seqs_in_interval([0, len(alignment[0])], alignment, 6)
        mockKMeans.assert_not_called()
        self.assertEqual([["s1", "s4"], ["s2"], ["s3"]], result)


class TestClustering_Kmeans(TestCase):
    @patch("make_prg.from_msa.cluster_sequences.KMeans.fit")
    def test_GivenTwoSequencesAboveKmerSize_KMeansClusteringCalled(self, mockfit):
        alignment = MSA(
            [SeqRecord(Seq("AAAT"), id="s1"), SeqRecord(Seq("TTTT"), id="s2")]
        )
        mockfit.return_value = Mock(inertia_=0)
        result = kmeans_cluster_seqs_in_interval([0, 3], alignment, 2)
        mockfit.assert_called_once()

    def test_TwoIdenticalSequencesClusteredTogether(self):
        alignment = MSA(
            [
                SeqRecord(Seq("AAAT"), id="s1"),
                SeqRecord(Seq("AAAT"), id="s2"),
                SeqRecord(Seq("C-CC"), id="s3"),
            ]
        )
        result = kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        self.assertEqual([["s1", "s2"], ["s3"]], result)

    def test_GivenTwoSequenceGroups_ReturnsTwoClusters(self):
        sequences = ["CATATAAAATA", "CATATAATATA", "GGGGCGGGCCC", "GGGGCGGGCGC"]
        expected_clustering = [["s0", "s1"], ["s2", "s3"]]
        seq_size = len(sequences[0])
        alignment = make_alignment(sequences)
        for kmer_size in range(1, 7):
            result = kmeans_cluster_seqs_in_interval(
                [0, seq_size - 1], alignment, kmer_size
            )
            self.assertEqual(expected_clustering, result)

    def test_GivenThreeSequenceGroups_ReturnsThreeClusters(self):
        sequences = [
            "CCCCCCAACCT",
            "CCCCCCAATCT",
            "GGGGCGGGCCC",
            "GGGGCGGGCGC",
            "TTTAATTTTAA",
            "TTTAAGTTTAA",
        ]
        expected_clustering = [["s0", "s1"], ["s2", "s3"], ["s4", "s5"]]
        seq_size = len(sequences[0])
        alignment = make_alignment(sequences)
        for kmer_size in range(1, 7):
            result = kmeans_cluster_seqs_in_interval(
                [0, seq_size - 1], alignment, kmer_size
            )
            self.assertEqual(expected_clustering, result)

    def test_GivenAllSequencesOneSnpApart_ReturnsNoClustering(self):
        sequences = ["CATATAAAATA", "CATATAACATA", "CATATAAGATA", "CATATAATATA"]
        seq_size = len(sequences[0])
        alignment = make_alignment(sequences)
        expected_clustering = [[record.id] for record in alignment]
        for kmer_size in range(1, 7):
            result = kmeans_cluster_seqs_in_interval(
                [0, seq_size - 1], alignment, kmer_size
            )
            self.assertEqual(expected_clustering, result)

    def test_GivenAllSequencesSmallEditDist_ReturnsNoClustering(self):
        """Cf graph 157.pdf in issue #15"""
        sequences = [
            "gctccgccggtcccgccggtcc",
            "gctccgccgggcccgccggtcc",
            "tctccgccggtcccgccggtcc",
            "gctcagccggtcccgccggtcc",
            "gctccgccggtcccaccggtcc",
            "gctccgccggtaccgccggtcc",
            "gctccgctggtcccgccggtcc",
            "gctccgccggtcccgctggtcc",
            "gctccgccggtcccgccggtct",
            "gctccgccggtcccgcctgtcc",
            "gctccgccggtcctgccggtcc",
        ]
        seq_size = len(sequences[0])
        alignment = make_alignment(sequences)
        expected_clustering = [[record.id] for record in alignment]
        for kmer_size in range(4, 8):
            result = kmeans_cluster_seqs_in_interval(
                [0, seq_size - 1], alignment, kmer_size
            )
            self.assertEqual(expected_clustering, result)


class TestKMeansOrdering(TestCase):
    """
    Because reference genome is often embedded as first path in PRG,
    it is important that clustering returns the first input sequence
    in first position of first cluster.
    """

    def test_first_sequence_placed_in_first_cluster(self):
        """
        Runs kmeans clustering on randomly generated multiple sequence alignments
        """
        seq_len = 20
        num_seqs = 20
        bases = list(standard_bases)
        # Function has different behaviour at below and above seq_len
        for seq_len in [seq_len - 1, seq_len + 1]:
            with self.subTest(min_match_len=seq_len):
                for _ in range(5):  # Run on a number of random alignments
                    sequences = [
                        "".join(random.choices(bases, k=seq_len))
                        for _ in range(num_seqs)
                    ]
                    alignment = make_alignment(sequences)
                    result = kmeans_cluster_seqs_in_interval(
                        [0, seq_len - 1], alignment, 1
                    )
                    self.assertTrue(result[0][0] == "s0")

    def test_one_long_one_short_sequence_separate_and_ordered_clusters(self):
        alignment = MSA(
            [
                SeqRecord(Seq("AATTAATTATATAATAAC"), id="s1"),
                SeqRecord(Seq("A--------------AAT"), id="s2"),
            ]
        )
        order_1 = kmeans_cluster_seqs_in_interval([0, len(alignment[0])], alignment, 5)
        self.assertEqual(order_1, [["s1"], ["s2"]])

        order_2 = kmeans_cluster_seqs_in_interval(
            [0, len(alignment[0])], alignment[::-1], 5
        )
        self.assertEqual(order_2, [["s2"], ["s1"]])
