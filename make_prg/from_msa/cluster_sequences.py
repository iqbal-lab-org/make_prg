import logging
from collections import defaultdict
from typing import List, Dict, Iterator, Union
from itertools import starmap, repeat, chain
from collections import Counter
import time

import numpy as np
from sklearn.cluster import KMeans

from make_prg.from_msa import MSA

IDs = List[str]
SeqToIDs = Dict[str, IDs]
ClusteredIDs = List[IDs]
Sequence = str
Sequences = List[str]
ClusteredSeqs = List[Sequences]
KmerIDs = Dict[str, int]

DISTANCE_THRESHOLD: float = 0.2
LENGTH_THRESHOLD: int = 5
MAX_CLUSTERS: int = 10


def count_distinct_kmers(seqs: Sequences, kmer_size: int) -> KmerIDs:
    for seq in seqs:
        if len(seq) < kmer_size:
            raise ValueError(f"Input sequence {seq} has length < kmer size {kmer_size}")
    result = dict()
    kmer_id = 0
    for seq in seqs:
        for start_pos in range(len(seq) - kmer_size + 1):
            kmer = seq[start_pos : start_pos + kmer_size]
            if kmer not in result:
                result[kmer] = kmer_id
                kmer_id += 1
    return result


def count_kmer_occurrences(seqs: Sequences, kmers: KmerIDs):
    """
    Computes a count matrix of kmers in :param seqs
    :param kmers: all possible kmers assumed to have been counted in seqs
    :return The counts of all kmers in each input sequence
    """
    num_kmers = len(kmers)
    kmer_size = len(next(iter(kmers)))
    result = np.zeros(shape=(len(seqs), num_kmers))
    for j, seq in enumerate(seqs):
        counts = result[j]
        for i in range(len(seq) - kmer_size + 1):
            kmer = seq[i : i + kmer_size]
            counts[kmers[kmer]] += 1
        result[j] = counts
    return result


def get_majority_char_in_column(sequences: Sequences, col_idx: int) -> str:
    max_idx = len(sequences[0]) - 1
    if not (0 <= col_idx <= max_idx):
        raise ValueError(f"Column index {col_idx} not in range(0,{max_idx})")

    column = [seq[col_idx] for seq in sequences]
    return Counter(column).most_common(1)[0][0]


def get_majority_string(sequences: Sequences) -> str:
    if len(sequences) == 1:
        return sequences[0]
    seqlen = len(sequences[0])
    if not all([len(seq) == seqlen for seq in sequences]):
        raise ValueError("Not all sequences have the same length")

    arguments = zip(repeat(sequences), range(seqlen))
    return "".join(starmap(get_majority_char_in_column, arguments))


def hamming_distance(seq1: Sequence, seq2: Sequence) -> int:
    seqlen = len(seq1)
    result = 0
    for i in range(seqlen):
        if seq1[i] != seq2[i]:
            result += 1
    return result


def get_distances(sequences: Sequences, consensus_string: Sequence) -> Iterator[int]:
    arguments = zip(sequences, repeat(consensus_string))
    return starmap(hamming_distance, arguments)


def get_one_ref_like_threshold_distance(seqlen: int) -> int:
    if seqlen < LENGTH_THRESHOLD:
        return 1
    else:
        return int(DISTANCE_THRESHOLD * seqlen)


def sequences_are_one_reference_like(sequences: Sequences) -> bool:
    majority_string = get_majority_string(sequences)
    distance_iterator = get_distances(sequences, majority_string)
    min_thresh = get_one_ref_like_threshold_distance(len(sequences[0]))
    return all(map(lambda dist: dist <= min_thresh, distance_iterator))


def cluster_further(clusters: ClusteredSeqs) -> bool:
    not_ref_like = lambda sequences: not sequences_are_one_reference_like(sequences)
    return any(map(not_ref_like, clusters))


def extract_clusters(
    seqdict: SeqToIDs, cluster_assignment: List[int], extract_IDs: bool = False
) -> Union[ClusteredSeqs, ClusteredIDs]:
    if extract_IDs:
        value_pool = list(seqdict.values())
    else:
        value_pool = [[seq] for seq in seqdict.keys()]
    num_elems = len(cluster_assignment)
    if num_elems != len(value_pool):
        raise ValueError(
            f"Mismatch between number of sequences/ID lists and number of cluster assignments"
        )
    num_clusters = max(cluster_assignment) + 1
    result = [list() for _ in range(num_clusters)]
    for cluster_num, clustered_elem in zip(cluster_assignment, value_pool):
        result[cluster_num].extend(clustered_elem)
    return result


def merge_clusters(clusters: List[ClusteredIDs], first_id: str) -> ClusteredIDs:
    merged_clusters = list()
    first_id_cluster = []
    for cluster in chain.from_iterable(clusters):
        if first_id in cluster:
            first_id_cluster = cluster
        else:
            merged_clusters.append(cluster)
    if len(first_id_cluster) == 0:
        raise ValueError(f"Could not find {first_id} in any cluster")
    return [first_id_cluster] + merged_clusters


def kmeans_cluster_seqs_in_interval(
    interval: List[int], alignment: MSA, kmer_size: int,
) -> ClusteredIDs:
    """Divide sequences in interval into subgroups of similar sequences."""
    interval_alignment = alignment[:, interval[0] : interval[1] + 1]

    logging.debug("Get kmeans partition of interval [%d, %d]", interval[0], interval[1])

    # Find unique sequences for clustering, but keep each sequence's IDs
    seq_to_ids: SeqToIDs = defaultdict(list)
    small_seq_to_ids: SeqToIDs = defaultdict(list)

    for record in interval_alignment:
        seq = str(record.seq.ungap("-"))
        if len(seq) >= kmer_size:
            seq_to_ids[seq].append(record.id)
        else:
            small_seq_to_ids[seq].append(record.id)

    num_clusters = 1
    num_sequences = len(seq_to_ids)
    if num_sequences > 2:
        distinct_sequences = list(seq_to_ids)
        distinct_kmers = count_distinct_kmers(distinct_sequences, kmer_size)
        count_matrix = count_kmer_occurrences(distinct_sequences, distinct_kmers)
        cluster_assignment = [0 for _ in range(len(seq_to_ids))]
        seqclustering: ClusteredSeqs = extract_clusters(seq_to_ids, cluster_assignment)

        while cluster_further(seqclustering):
            num_clusters += 1
            if num_clusters >= MAX_CLUSTERS:
                num_clusters = num_sequences
            if num_clusters == num_sequences:
                break
            start = time.time()
            kmeans = KMeans(n_clusters=num_clusters, random_state=2).fit(count_matrix)
            cluster_assignment = list(kmeans.predict(count_matrix))
            seqclustering = extract_clusters(seq_to_ids, cluster_assignment)

    if num_clusters == 1 or num_clusters == num_sequences:
        cluster_assignment = list(range(num_sequences))
    id_clustering = []
    if num_sequences > 0:
        id_clustering: ClusteredIDs = extract_clusters(
            seq_to_ids, cluster_assignment, extract_IDs=True
        )

    first_id = interval_alignment[0].id
    result = merge_clusters([id_clustering, small_seq_to_ids.values()], first_id)

    assert len(interval_alignment) == sum(
        [len(i) for i in result]
    ), "Each input sequence should be in a cluster"
    return result
