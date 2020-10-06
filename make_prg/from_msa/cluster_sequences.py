import logging
from collections import defaultdict
from typing import List, Dict

import numpy as np
from sklearn.cluster import KMeans

from make_prg.from_msa import MSA

IDs = List[str]
Clustering = List[IDs]
Sequence = str
Sequences = List[str]
KmerIDs = Dict[str, int]


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


def kmeans_cluster_seqs_in_interval(
    interval: List[int], alignment: MSA, kmer_size: int,
) -> Clustering:
    """Divide sequences in interval into subgroups of similar sequences."""
    interval_alignment = alignment[:, interval[0] : interval[1] + 1]

    logging.debug("Get kmeans partition of interval [%d, %d]", interval[0], interval[1])

    # Find unique sequences for clustering, but keep each sequence's IDs
    seq_to_ids = defaultdict(list)
    small_seq_to_ids = defaultdict(list)

    for record in interval_alignment:
        seq = str(record.seq.ungap("-"))
        if len(seq) >= kmer_size:
            seq_to_ids[seq].append(record.id)
        else:
            small_seq_to_ids[seq].append(record.id)

    clustered_ids = []
    if len(seq_to_ids) < 2:
        clustered_ids.extend(seq_to_ids.values())

    else:
        distinct_sequences = list(seq_to_ids)
        distinct_kmers = count_distinct_kmers(distinct_sequences, kmer_size)
        count_matrix = count_kmer_occurrences(distinct_sequences, distinct_kmers)

        number_of_clusters = 1
        kmeans = KMeans(n_clusters=number_of_clusters, random_state=2).fit(count_matrix)
        pre_cluster_inertia = kmeans.inertia_

        cluster_inertia = pre_cluster_inertia
        while (
            cluster_inertia > 0
            and cluster_inertia > pre_cluster_inertia / 2
            and number_of_clusters < len(seq_to_ids)
        ):
            number_of_clusters += 1
            kmeans = KMeans(n_clusters=number_of_clusters, random_state=2).fit(
                count_matrix
            )
            cluster_inertia = kmeans.inertia_

        # convert cluster numbers to sequence record IDs
        if pre_cluster_inertia > 0:
            cluster_ids = list(kmeans.predict(count_matrix))
            for i in range(max(cluster_ids) + 1):
                clustered_ids.append([])
            all_ids = list(seq_to_ids.values())
            for i, cluster_id in enumerate(cluster_ids):
                clustered_ids[cluster_id].extend(all_ids[i])
        else:
            clustered_ids = list(seq_to_ids.values())

    first_id = interval_alignment[0].id
    id_lists = [[]]  # Reserve space for first seq id
    for cluster in clustered_ids:
        if first_id in set(cluster):
            id_lists[0] = cluster
        else:
            id_lists.append(cluster)

    for ids in small_seq_to_ids.values():
        if first_id in set(ids):
            id_lists[0] = ids
        else:
            id_lists.append(ids)

    assert len(interval_alignment) == sum(
        [len(i) for i in id_lists]
    ), "I seem to have lost (or gained?) some sequences in the process of clustering"
    return id_lists
