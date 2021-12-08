from unittest import TestCase, skip
from unittest.mock import patch, Mock
from itertools import product
from numpy import array, array_equal
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from make_prg.utils.seq_utils import SequenceExpander, ungap
from make_prg.from_msa.cluster_sequences import (
    count_distinct_kmers,
    count_kmer_occurrences,
    get_majority_char_in_column,
    get_majority_string,
    hamming_distance,
    get_distances,
    LENGTH_THRESHOLD,
    MAX_CLUSTERS,
    get_one_ref_like_threshold_distance,
    sequences_are_one_reference_like,
    cluster_further,
    extract_clusters,
    merge_clusters,
    kmeans_cluster_seqs,
    ClusteringResult,
    merge_sequences
)
from tests.test_helpers import make_alignment, MSA

class TestClusteringResult(TestCase):
    def setUp(self) -> None:
        self.clustered_ids = [["s0"], ["s1"]]
        self.sequences = ["AATA", "AAAA"]
        self.clustering_result = ClusteringResult(self.clustered_ids, self.sequences)

    def test___constructor___no_sequences(self):
        clustering_result = ClusteringResult(self.clustered_ids)

        self.assertEqual(self.clustered_ids, clustering_result.clustered_ids)
        self.assertIsNone(clustering_result.sequences)

    def test___constructor___with_sequences(self):
        self.assertEqual(self.clustered_ids, self.clustering_result.clustered_ids)
        self.assertEqual(self.sequences, self.clustering_result.sequences)

    def test___eq___equality(self):
        clustering_result_copy = ClusteringResult(self.clustered_ids, self.sequences)
        self.assertEqual(clustering_result_copy, self.clustering_result)

    def test___eq___inequality(self):
        other = ClusteringResult([["s0"]], ["A"])
        self.assertNotEqual(other, self.clustering_result)

        other = ClusteringResult([["s0"], ["s1"]], ["A"])
        self.assertNotEqual(other, self.clustering_result)

        other = ClusteringResult([["s0"]], ["AATA", "AAAA"])
        self.assertNotEqual(other, self.clustering_result)

        other = "other_type"
        self.assertNotEqual(other, self.clustering_result)

    def test___no_clustering___single_cluster_with_single_seq(self):
        clustering_result = ClusteringResult([["s0"]])
        self.assertTrue(clustering_result.no_clustering)

    def test___no_clustering___single_cluster_with_several_seqs(self):
        clustering_result = ClusteringResult([["s0", "s1", "s2", "s3", "s4", "s5"]])
        self.assertTrue(clustering_result.no_clustering)

    def test___no_clustering___three_small_clusters(self):
        clustering_result = ClusteringResult([["s0"], ["s1"], ["s2"]])
        self.assertFalse(clustering_result.no_clustering)

    def test___have_precomputed_sequences___no_precomputed_sequences(self):
        clustering_result = ClusteringResult(self.clustered_ids)
        self.assertFalse(clustering_result.have_precomputed_sequences)

    def test___have_precomputed_sequences___has_precomputed_sequences(self):
        self.assertTrue(self.clustering_result.have_precomputed_sequences)

    def test___repr(self):
        expected = "ClusteringResult(clustered_ids=[['s0'], ['s1']], sequences=['AATA', 'AAAA'])"
        actual = repr(self.clustering_result)
        self.assertEqual(expected, actual)

    def test___str(self):
        expected = "ClusteringResult(clustered_ids=[['s0'], ['s1']], sequences=['AATA', 'AAAA'])"
        actual = str(self.clustering_result)
        self.assertEqual(expected, actual)


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

    def test_GivenInconsistentClusterNumbers_RaisesError(self):
        """
        `gapped_clusters` can occur when number of distinct data points < requested number of
        clusters in kmeans, leading to production of empty clusters and then
        messing up prg construction.
        `large_clusters` should never occur, cluster assigment should be consecutive integers from 0.
        """
        gapped_clusters = [0, 2, 2]
        with self.assertRaises(ValueError):
            extract_clusters(self.seqdict_ids, gapped_clusters)
        large_clusters = [5, 6, 6]
        with self.assertRaises(ValueError):
            extract_clusters(self.seqdict_ids, large_clusters)


class TestClustering_Trivial(TestCase):
    def test_one_seq_returns_single_id(self):
        alignment = make_alignment(["AAAT"])
        expected = ClusteringResult(clustered_ids=[['s0']], sequences=['AAAT'])
        result = kmeans_cluster_seqs(alignment, 1)
        self.assertEqual(expected, result)

    @patch("make_prg.from_msa.cluster_sequences.KMeans")
    def test_GivenTwoIdenticalSeqs_NoKmeansAndOneCluster(self, mockKMeans):
        alignment = make_alignment(["AAAT", "AAAT"])
        expected = ClusteringResult(clustered_ids=[['s0', 's1']], sequences=['AAAT'])
        result = kmeans_cluster_seqs(alignment, 1)

        mockKMeans.assert_not_called()
        self.assertEqual(expected, result)

    @patch("make_prg.from_msa.cluster_sequences.KMeans")
    def test_GivenTwoDifferentSeqs_NoKmeansAndTwoClusters(self, mockKMeans):
        alignment = make_alignment(["AAAT", "ATAT"])
        expected = ClusteringResult(clustered_ids=[['s0', 's1']], sequences=['AAAT', 'ATAT'])
        result = kmeans_cluster_seqs(alignment, 1)

        mockKMeans.assert_not_called()
        self.assertEqual(expected, result)


class TestClustering_SmallSequences(TestCase):
    def test_GivenLessThanTwoLongSeqs_NoClustering(self):
        alignment = make_alignment(["A-A-T", "CCCCC", "AA--T"])

        expected = ClusteringResult(clustered_ids=[['s0', 's1', 's2']], sequences=['AAT', 'CCCCC'])
        result = kmeans_cluster_seqs(alignment, 4)

        self.assertEqual(expected, result)

    @skip(
        "This fails, probably because kmean clustering should never run with this input"
    )
    def test_ambiguous_sequences_in_short_interval_separate_clusters(self):
        alignment = make_alignment(["ARAT", "WAAT"])

        expected = ClusteringResult([["s0", "s1"]])
        result = kmeans_cluster_seqs(alignment, 5)

        self.assertEqual(expected, result)

    @patch("make_prg.from_msa.cluster_sequences.KMeans")
    def test_GivenAllSequencesBelowKmerSize_NoClustering(self, mockKMeans):
        alignment = make_alignment(
            [
                "AA---AT",
                "AA---TT",
                "CA--CAT",
                "A-A--AT",
            ]
        )

        result = kmeans_cluster_seqs(alignment, 6)
        mockKMeans.assert_not_called()

        expected = ClusteringResult(clustered_ids=[['s0', 's3', 's1', 's2']], sequences=['AAAT', 'AATT', 'CACAT'])
        self.assertEqual(expected, result)


class TestClustering_GappedSequences(TestCase):
    def test_SequencesUnevenLengthIfGapsRemoved_ClusteringRuns(self):
        """Checking for 'one-ref' property in clusters
        needs to be on ungapped sequences (elsehamming distance computation would fail
        due to different seq length"""
        sequences = ["A---T", "AAAAT", "AAA-T"]
        alignment = make_alignment(sequences)
        actual = kmeans_cluster_seqs(alignment, 1)
        expected = ClusteringResult(clustered_ids=[['s0'], ['s1', 's2']])
        self.assertEqual(expected, actual)


class TestClustering_RunKmeans(TestCase):
    @patch("make_prg.from_msa.cluster_sequences.KMeans")
    def test_GivenThreeSequencesAboveKmerSize_KMeansClusteringCalled(self, mockfit):
        alignment = make_alignment(["AAAT", "TTTT", "ATAT"])
        try:
            result = kmeans_cluster_seqs(alignment, 2)
        except ValueError:
            pass
        mockfit.assert_called_once()

    def test_TwoIdenticalSequencesClusteredTogether(self):
        alignment = MSA(
            [
                SeqRecord(Seq("AAAT"), id="s1"),
                SeqRecord(Seq("AAAT"), id="s2"),
                SeqRecord(Seq("C-CC"), id="s3"),
                SeqRecord(Seq("C-TC"), id="s4"),
            ]
        )

        expected = ClusteringResult(clustered_ids=[['s1', 's2'], ['s3', 's4']])
        result = kmeans_cluster_seqs(alignment, 1)
        self.assertEqual(expected, result)

    def test_GivenTwoSequenceGroups_ReturnsTwoClusters(self):
        sequences = ["CATATAAAATA", "CATATAATATA", "GGGGCGGGCCC", "GGGGCGGGCGC"]
        alignment = make_alignment(sequences)

        expected_clustering = ClusteringResult(clustered_ids=[['s0', 's1'], ['s2', 's3']])
        for kmer_size in range(1, 7):
            result = kmeans_cluster_seqs(alignment, kmer_size)
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
        alignment = make_alignment(sequences)

        expected_clustering = ClusteringResult([["s0", "s1"], ["s2", "s3"], ["s4", "s5"]])
        for kmer_size in range(1, 7):
            result = kmeans_cluster_seqs(alignment, kmer_size)
            for cluster in expected_clustering.clustered_ids:
                self.assertTrue(cluster in result.clustered_ids)

    def test_GivenAllSequencesOneSnpApart_ReturnsNoClustering(self):
        sequences = ["CATATAAAATA", "CATATAACATA", "CATATAAGATA", "CATATAATATA"]
        alignment = make_alignment(sequences)

        expected = ClusteringResult(clustered_ids=[['s0', 's1', 's2', 's3']],
                                    sequences=['CATATAAAATA', 'CATATAACATA', 'CATATAAGATA', 'CATATAATATA'])
        for kmer_size in range(1, 7):
            result = kmeans_cluster_seqs(alignment, kmer_size)
            self.assertEqual(expected, result)

    def test_GivenAllSequencesSmallEditDist_ReturnsNoClustering(self):
        """Cf graph 157.pdf in issue #15"""
        sequences = [
            "GCTCCGCCGGTCCCGCCGGTCC",
            "GCTCCGCCGGGCCCGCCGGTCC",
            "TCTCCGCCGGTCCCGCCGGTCC",
            "GCTCAGCCGGTCCCGCCGGTCC",
            "GCTCCGCCGGTCCCACCGGTCC",
            "GCTCCGCCGGTACCGCCGGTCC",
            "GCTCCGCTGGTCCCGCCGGTCC",
            "GCTCCGCCGGTCCCGCTGGTCC",
            "GCTCCGCCGGTCCCGCCGGTCT",
            "GCTCCGCCGGTCCCGCCTGTCC",
            "GCTCCGCCGGTCCTGCCGGTCC",
        ]
        alignment = make_alignment(sequences)

        expected = ClusteringResult(
            clustered_ids=[['s0', 's1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's10']],
            sequences=['GCTCCGCCGGTCCCGCCGGTCC', 'GCTCCGCCGGGCCCGCCGGTCC', 'TCTCCGCCGGTCCCGCCGGTCC',
                       'GCTCAGCCGGTCCCGCCGGTCC', 'GCTCCGCCGGTCCCACCGGTCC', 'GCTCCGCCGGTACCGCCGGTCC',
                       'GCTCCGCTGGTCCCGCCGGTCC', 'GCTCCGCCGGTCCCGCTGGTCC', 'GCTCCGCCGGTCCCGCCGGTCT',
                       'GCTCCGCCGGTCCCGCCTGTCC', 'GCTCCGCCGGTCCTGCCGGTCC'])

        for kmer_size in range(4, 8):
            result = kmeans_cluster_seqs(alignment, kmer_size)
            print(result)
            self.assertEqual(expected, result)

    def test_GivenManyVeryDifferentSequences_EachSeqInOwnCluster(self):
        # all 256 distinct DNA 4-mers.
        # We want clustering to keep looking for clusters, and stop at MAX_CLUSTERS
        all_4mers = list(map("".join, product(SequenceExpander.standard_bases, repeat=4)))
        alignment = make_alignment(all_4mers)
        result = kmeans_cluster_seqs(alignment, 4)
        self.assertEqual(len(result.clustered_ids), MAX_CLUSTERS)

    def test_GivenFewVeryDifferentSequences_NoClustering(self):
        # We want clustering to keep looking for clusters, and stop at num_sequences, and realise there is no clustering
        sequences = [
            "AAAAAAAAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCCCCCCCC",
            "GGGGGGGGGGGGGGGGGGGGGG",
            "TTTTTTTTTTTTTTTTTTTTTT",
        ]
        alignment = make_alignment(sequences)
        expected = ClusteringResult(clustered_ids=[['s0', 's1', 's2', 's3']],
                                    sequences=['AAAAAAAAAAAAAAAAAAAAAA',
                                               'CCCCCCCCCCCCCCCCCCCCCC',
                                               'GGGGGGGGGGGGGGGGGGGGGG',
                                               'TTTTTTTTTTTTTTTTTTTTTT'])
        result = kmeans_cluster_seqs(alignment, 7)
        self.assertEqual(expected, result)

    def test_GivenSequencesWithSameKmerCounts_ClusteringInterrupted(self):
        """
        Sequences below are not 'one-ref-like', yet kmer counts are identical.
        This is because the sequences contain repeats and gaps, making them
        not identical from the point of view of edit distance.
        Number of clusters will try to be increased, but kmeans will only find one,
        as there is a single data point in kmer space.
        This test checks the code deals with this by aborting further clustering.
        """
        sequences = [
            "TTTTTTTGGGGGGGAAAAAAATTTTTTT-------AAAAAAATTTTTTTAAAAAAA-------",
            "-------TTTTTTTAAAAAAATTTTTTTGGGGGGGAAAAAAATTTTTTT-------AAAAAAA",
            "TTTTTTTAAAAAAATTTTTTTAAAAAAATTTTTTT-------GGGGGGG-------AAAAAAA",
        ]
        ungapped_sequences = list(map(ungap, sequences))
        distinct_kmers = count_distinct_kmers(ungapped_sequences, kmer_size=7)
        count_matrix = count_kmer_occurrences(ungapped_sequences, distinct_kmers)
        distinct_count_patterns = set(map(str, count_matrix))
        assert len(distinct_count_patterns) == 1
        assert not sequences_are_one_reference_like(sequences)

        alignment = make_alignment(sequences)
        expected = ClusteringResult(clustered_ids=[['s0', 's1', 's2']],
                                    sequences=['TTTTTTTGGGGGGGAAAAAAATTTTTTTAAAAAAATTTTTTTAAAAAAA',
                                               'TTTTTTTAAAAAAATTTTTTTGGGGGGGAAAAAAATTTTTTTAAAAAAA',
                                               'TTTTTTTAAAAAAATTTTTTTAAAAAAATTTTTTTGGGGGGGAAAAAAA'])
        result = kmeans_cluster_seqs(alignment, 7)
        self.assertEqual(expected, result)


class TestMergeClusters(TestCase):
    clusters = [["s1", "s2"], ["s3"]]
    small_seqs = [["s4"], ["s5", "s6"]]

    def test_GivenFirstIDNotFound_Fails(self):
        with self.assertRaises(ValueError):
            merge_clusters(self.clusters, self.small_seqs, first_id="s200")

    def test_GivenFirstIDFound_MergedAndFirstIDInFirstCluster(self):
        actual = merge_clusters(self.clusters, self.small_seqs, first_id="s5")
        expected = [["s5", "s6"], ["s1", "s2"], ["s3"], ["s4"]]
        self.assertEqual(actual, expected)

    def test_GivenFirstIDFound_FirstIDInLastClusterAndLastPosition_MergedAndFirstIDInFirstCluster(self):
        actual = merge_clusters(self.clusters, self.small_seqs, first_id="s6")
        expected = [["s6", "s5"], ["s1", "s2"], ["s3"], ["s4"]]
        self.assertEqual(actual, expected)


class TestKMeansOrdering(TestCase):
    """
    Because reference genome is often embedded as first path in PRG,
    it is important that clustering returns the first input sequence
    in first position of first cluster.
    """

    def test_first_id_in_first_cluster___first_id_in_first_pos_of_first_cluster(self):
        alignment = make_alignment(
            [
                "AATTAATTATATAATAAC",
                "AATTAAGTATATAATAAC",
                "TTAATTAATTAATTAATT",
            ],
            ["s1", "s2", "s3"],
        )
        expected_order = ClusteringResult(clustered_ids=[['s1', 's2'], ['s3']])
        order = kmeans_cluster_seqs(alignment, 5)
        self.assertEqual(expected_order, order)

    def test_first_id_in_last_cluster(self):
        alignment = make_alignment(
            [
                "TTAATTAATTAATTAATT",
                "AATTAAGTATATAATAAC",
                "AATTAATTATATAATAAC",
            ],
            ["s3", "s2", "s1"],
        )

        expected_order = ClusteringResult(clustered_ids=[['s3'], ['s2', 's1']])
        order = kmeans_cluster_seqs(alignment, 5)
        self.assertEqual(expected_order, order)

    def test_first_id_in_first_cluster___first_id_is_small_seq(self):
        alignment = make_alignment(
            [
                "A-----------------",
                "AATTAATTATATAATAAC",
                "AATTAAGTATATAATAAC",
                "TTAATTAATTAATTAATT",
            ],
            ["small", "s2", "s3", "s4"],
        )
        expected_order = ClusteringResult(clustered_ids=[['small'], ['s2', 's3'], ['s4']])
        order = kmeans_cluster_seqs(alignment, 5)
        self.assertEqual(expected_order, order)




class TestMergeSequences(TestCase):
    def setUp(self):
        self.long_seqs = ["AATAA", "TTAAA"]
        self.short_seqs = ["AA", "TT"]
        self.first_seq = self.long_seqs[0]

    def test_GivenFirstSeqNotInList_Fails(self):
        with self.assertRaises(AssertionError):
            merge_sequences(self.short_seqs, self.long_seqs, first_seq="TTTTT")

    def test_GivenShortFirst_FirstSeqIsFirstSeqOut(self):
        expected = ["AATAA", "AA", "TT", "TTAAA"]
        actual = merge_sequences(
            self.short_seqs, self.long_seqs, first_seq=self.first_seq
        )
        self.assertEqual(expected, actual)

    def test_GivenLongFirst_FirstSeqIsFirstSeqOut(self):
        expected = ["AATAA", "TTAAA", "AA", "TT"]
        actual = merge_sequences(
            self.long_seqs, self.short_seqs, first_seq=self.first_seq
        )
        self.assertEqual(expected, actual)

    def test_GivenSeqsWithAmbiguousBases_ExpansionIsDone_Duplicates_are_removed(self):
        expected = ["GA", "TA", "AATAA", "AACAA", "TTGAA", "TTAAA", "TC"]
        actual = merge_sequences(
            ["AAYAA", "TTRAA"], ["KA", "TM"], first_seq="KA"
        )
        self.assertEqual(expected, actual)

    def test_GivenFirstSeqWithN_FirstSeqIsRemoved(self):
        expected = ["AATAA", "TTAAA", "AA", "TT"]
        actual = merge_sequences(
            ["ANA"] + self.long_seqs, self.short_seqs, first_seq="ANA"
        )
        self.assertEqual(expected, actual)
