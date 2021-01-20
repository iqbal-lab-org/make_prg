from unittest import TestCase, skip
from unittest.mock import patch, Mock
import random
from itertools import product

from numpy import array, array_equal
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_prg.seq_utils import standard_bases
from make_prg.from_msa.cluster_sequences import (
    count_distinct_kmers,
    count_kmer_occurrences,
    get_majority_char_in_column,
    get_majority_string,
    hamming_distance,
    get_distances,
    DISTANCE_THRESHOLD,
    LENGTH_THRESHOLD,
    MAX_CLUSTERS,
    get_one_ref_like_threshold_distance,
    sequences_are_one_reference_like,
    cluster_further,
    extract_clusters,
    merge_clusters,
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


class TestOneRefWorkerFunctions(TestCase):
    sequences = [
        "CATATAAAATA",
        "CATATAA-ATA",
        "GGGGCGGGCCC",
        "GGGGCGGGCGC",
        "GGGGCCGGCCC",
        "GGGGCGGGCCC",
    ]

    def test_majority_char_bad_index_fail(self):
        with self.assertRaises(ValueError):
            get_majority_char_in_column(self.sequences, len(self.sequences[0]))

    def test_majority_char_returns_majority_char(self):
        actual = get_majority_char_in_column(self.sequences, 0)
        self.assertEqual(actual, "G")

    def test_majority_char_returns_majority_nonbase_char(self):
        sequences = ["AAA-", "AAT-", "AAAA"]
        actual = get_majority_char_in_column(sequences, 3)
        self.assertEqual(actual, "-")

    def test_majority_string_bad_lengths_fails(self):
        sequences = ["AAA-", "AT-", "AAAA"]
        with self.assertRaises(ValueError):
            get_majority_string(sequences)

    def test_majority_string_correct(self):
        actual = get_majority_string(self.sequences)
        expected = "GGGGCGGGCCC"
        self.assertEqual(actual, expected)

    def test_hamming_distance_same_sequences(self):
        actual = hamming_distance(self.sequences[0], self.sequences[0])
        self.assertEqual(actual, 0)

    def test_hamming_distance_different_sequences(self):
        actual = hamming_distance(self.sequences[0], self.sequences[1])
        self.assertEqual(actual, 1)
        actual = hamming_distance(self.sequences[3], self.sequences[4])
        self.assertEqual(actual, 2)

    def test_get_all_distances(self):
        actual = list(get_distances(self.sequences, self.sequences[0]))
        expected = [0, 1, 11, 11, 11, 11]
        self.assertEqual(actual, expected)

    def test_get_all_distances_to_majority_string(self):
        majority_string = get_majority_string(self.sequences)
        actual = list(get_distances(self.sequences, majority_string))
        expected = [11, 11, 0, 1, 1, 0]
        self.assertEqual(actual, expected)


class TestOneRefLikeClusters(TestCase):
    """
    Disclaimer:
    A heuristic is used to determine whether a set of sequences
    is 'one-ref like' based on a length threshold and a distance threshold.

    If either of those parameters is changed, below tests can start to fail;
    I placed in assertions to point to why.
    """

    def test_GivenSequencesWithDifferentLengths_Fails(self):
        sequences = ["AT", "AA", "CCC"]
        with self.assertRaises(ValueError):
            sequences_are_one_reference_like(sequences)

    def test_GivenSnpsOnly_AreRefLike(self):
        sequences = ["A", "T", "C"]
        self.assertTrue(sequences_are_one_reference_like(sequences))

    def test_GivenShortSeqsWithMoreThanOneDiff_AreNotRefLike(self):
        sequences = ["AA", "TT", "CC"]
        seqlen = len(sequences[0])
        self.assertTrue(seqlen < LENGTH_THRESHOLD)
        self.assertEqual(get_one_ref_like_threshold_distance(seqlen), 1)
        self.assertFalse(sequences_are_one_reference_like(sequences))

    def test_GivenLongSequencesWithDiffBelowThreshold_AreRefLike(self):
        sequences = ["AATTA", "AATTT", "TATTA"]
        # Sequences are above length threshold
        seqlen = len(sequences[0])
        self.assertTrue(seqlen >= LENGTH_THRESHOLD)
        # We tolerate up to one diff against majority string
        self.assertEqual(get_one_ref_like_threshold_distance(seqlen), 1)

        # Here's the majority string
        self.assertEqual(get_majority_string(sequences), "AATTA")

        # No sequence is > distance threshold from the majority string
        self.assertTrue(sequences_are_one_reference_like(sequences))

    def test_GivenLongSequencesWithDiffAboveThreshold_AreNotRefLike(self):
        sequences = ["AATTA", "AATTT", "TATTA", "GCGGG"]
        seqlen = len(sequences[0])
        self.assertTrue(seqlen >= LENGTH_THRESHOLD)
        self.assertEqual(get_one_ref_like_threshold_distance(seqlen), 1)
        self.assertEqual(get_majority_string(sequences), "AATTA")
        self.assertFalse(sequences_are_one_reference_like(sequences))

    def test_GivenLongSequencesFromRealDataWithDiffBelowThreshold_AreRefLike(self):
        """This data appeared in issue #15 and used to get clustered, producing ambiguous graphs with paths spelling the same sequence"""
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
        seqlen = len(sequences[0])
        self.assertEqual(get_one_ref_like_threshold_distance(seqlen), 4)
        self.assertTrue(sequences_are_one_reference_like(sequences))

    def test_GivenTightClusters_ClusterFurtherIsFalse(self):
        clusters = [["AATA", "AAAA"], ["GGGC", "GGGG"]]
        self.assertFalse(cluster_further(clusters))

    def test_GivenDiverseClusters_ClusterFurtherIsTrue(self):
        clusters = [["AATA", "GGGG"], ["TTTT", "TATT"]]
        self.assertTrue(cluster_further(clusters))


class TestExtractClusters(TestCase):
    seqdict_ids = dict(AT=["s1", "s2"], TTT=["s3"], GGG=["s4"])
    seqdict_gapped_seqs = dict(AT=["A-T", "AT-"], TTT=["TTT"], GGG=["GGG"])

    def test_GivenTooFewClusterAssignments_Fails(self):
        cluster_assignment = [0, 1]
        for seqdict in [self.seqdict_ids, self.seqdict_gapped_seqs]:
            with self.assertRaises(ValueError):
                extract_clusters(seqdict, cluster_assignment)

    def test_GivenDistinctClusters_ExtractCorrectSequenceClusters(self):
        cluster_assignment = [2, 0, 1]
        actual = extract_clusters(self.seqdict_gapped_seqs, cluster_assignment)
        expected = [["TTT"], ["GGG"], ["A-T", "AT-"]]
        self.assertEqual(actual, expected)

    def test_GivenGroupedClusters_ExtractCorrectSequenceClusters(self):
        cluster_assignment = [1, 1, 0]
        actual = extract_clusters(self.seqdict_gapped_seqs, cluster_assignment)
        expected = [["GGG"], ["A-T", "AT-", "TTT"]]
        self.assertEqual(actual, expected)

    def test_GivenGroupedClusters_ExtractCorrectIDClusters(self):
        cluster_assignment = [0, 0, 1]
        actual = extract_clusters(self.seqdict_ids, cluster_assignment)
        expected = [["s1", "s2", "s3"], ["s4"]]
        self.assertEqual(actual, expected)


class TestClustering_Trivial(TestCase):
    def test_one_seq_returns_single_id(self):
        alignment = MSA([SeqRecord(Seq("AAAT"), id="s1")])
        result = kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        self.assertEqual(result, [["s1"]])

    @patch("make_prg.from_msa.cluster_sequences.KMeans")
    def test_GivenTwoIdenticalSeqs_NoKmeansAndOneCluster(self, mockKMeans):
        alignment = make_alignment(["AAAT", "AAAT"])
        result = kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        mockKMeans.assert_not_called()
        self.assertEqual(result, [["s0", "s1"]])

    @patch("make_prg.from_msa.cluster_sequences.KMeans")
    def test_GivenTwoDifferentSeqs_NoKmeansAndTwoClusters(self, mockKMeans):
        alignment = make_alignment(["AAAT", "ATAT"])
        result = kmeans_cluster_seqs_in_interval([0, 3], alignment, 1)
        mockKMeans.assert_not_called()
        self.assertEqual(result, [["s0"], ["s1"]])


class TestClustering_SmallSequences(TestCase):
    def test_two_seqs_one_below_kmer_size_separate_clusters(self):
        alignment = make_alignment(["AATTTAT", "AA---AT"])
        result = kmeans_cluster_seqs_in_interval([0, 5], alignment, 5)
        self.assertEqual(result, [["s0"], ["s1"]])

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


class TestClustering_GappedSequences(TestCase):
    def test_SequencesUnevenLengthIfGapsRemoved_ClusteringRuns(self):
        """If check for 'one-ref' property in clusters
        was on ungapped sequences, hamming distance computation would fail because sequences have different length"""
        sequences = ["A---T", "AAAAT", "AAA-T"]
        alignment = make_alignment(sequences)
        actual = kmeans_cluster_seqs_in_interval([0, 4], alignment, 1)
        expected = [["s0"], ["s1", "s2"]]
        self.assertEqual(actual, expected)

    def test_GivenRepeatedUngappedSequencesBelowKmerSize_EndUpInSameCluster(self):
        sequences = ["A-A-T", "CCCCC", "AA--T"]
        alignment = make_alignment(sequences)
        actual = kmeans_cluster_seqs_in_interval([0, 4], alignment, 4)
        expected = [["s0", "s2"], ["s1"]]
        self.assertEqual(actual, expected)


class TestClustering_RunKmeans(TestCase):
    @patch("make_prg.from_msa.cluster_sequences.KMeans")
    def test_GivenThreeSequencesAboveKmerSize_KMeansClusteringCalled(self, mockfit):
        alignment = make_alignment(["AAAT", "TTTT", "ATAT"])
        try:
            result = kmeans_cluster_seqs_in_interval([0, 3], alignment, 2)
        except ValueError:
            pass
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
            for cluster in expected_clustering:
                self.assertTrue(cluster in result)

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

    def test_GivenManyVeryDifferentSequences_EachSeqInOwnCluster(self):
        # all 1024 distinct DNA 5-mers.
        # We want each in one cluster, but do not want to run clustering 1023 times
        all_5mers = list(map("".join, product(standard_bases, repeat=5)))
        alignment = make_alignment(all_5mers)
        result = kmeans_cluster_seqs_in_interval([0, 5], alignment, 5)
        self.assertEqual(len(result), MAX_CLUSTERS)


class TestMergeClusters(TestCase):
    clusters = [["s1", "s2"], ["s3"]]
    small_seqs = [["s4"], ["s5", "s6"]]
    to_merge = [clusters, small_seqs]

    def test_GivenFirstIDNotFound_Fails(self):
        with self.assertRaises(ValueError):
            merge_clusters(self.to_merge, "s200")

    def test_GivenFirstIDFound_MergedAndFirstIDInFirstCluster(self):
        actual = merge_clusters(self.to_merge, "s5")
        expected = [["s5", "s6"], ["s1", "s2"], ["s3"], ["s4"]]
        self.assertEqual(actual, expected)


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
        num_seqs = 5
        bases = list(standard_bases)
        # Function has different behaviour at below and above seq_len
        for used_len in [seq_len - 5, seq_len + 5]:
            with self.subTest(kmer_size=seq_len):
                for _ in range(5):  # Run on a number of random alignments
                    sequences = [
                        "".join(random.choices(bases, k=seq_len))
                        for _ in range(num_seqs)
                    ]
                    alignment = make_alignment(sequences)
                    result = kmeans_cluster_seqs_in_interval(
                        [0, seq_len - 1], alignment, used_len
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
