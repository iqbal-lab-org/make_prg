import logging
from collections import defaultdict
from typing import List
from itertools import chain

import numpy as np
from Bio.AlignIO import MultipleSeqAlignment
from sklearn.cluster import KMeans

from make_prg.io_utils import load_alignment_file
from make_prg.seq_utils import (
    remove_duplicates,
    remove_gaps,
    get_interval_seqs,
    ambiguous_bases,
)


class AlignedSeq(object):
    """
    Object based on a set of aligned sequences. Note min_match_length must be strictly greater than max_nesting + 1.
    """

    def __init__(
        self,
        msa_file,
        alignment_format="fasta",
        max_nesting=2,
        nesting_level=1,
        min_match_length=3,
        site=5,
        alignment=None,
        interval=None,
        prg_file=None,
    ):
        self.msa_file = msa_file
        self.alignment_format = alignment_format
        self.max_nesting = max_nesting
        self.nesting_level = nesting_level
        self.min_match_length = min_match_length
        self.site = site
        self.alignment = alignment
        if self.alignment is None:
            self.alignment = load_alignment_file(msa_file, alignment_format)

        self.interval = interval
        self.consensus = self.get_consensus()
        self.length = len(self.consensus)
        (self.match_intervals, self.non_match_intervals) = self.interval_partition()
        self.check_nonmatch_intervals()
        self.all_intervals = self.match_intervals + self.non_match_intervals
        logging.info("Non match intervals: %s", self.non_match_intervals)
        self.all_intervals.sort()

        # properties for stats
        self.subAlignedSeqs = {}

        # make prg
        self.delim_char = " "
        self.prg = ""
        if prg_file:
            logging.info(
                "Reading from a PRG file which already exists. To regenerate, delete it."
            )
            with open(prg_file, "r") as f:
                self.prg = f.read()
        else:
            self.prg = self.get_prg()

    def get_consensus(self):
        """Given a set of aligment records from AlignIO, creates
        a consensus string.
        IUPAC ambiguous bases result in non-consensus and are later expanded in the prg.
        N results in consensus at that position unless they are all N."""
        first_string = str(self.alignment[0].seq)
        consensus_string = ""
        for i, letter in enumerate(first_string):
            consensus = True
            for record in self.alignment:
                if letter == "N" or record.seq[i] == "N":
                    if letter == "N" and record.seq[i] != "N":
                        letter = record.seq[i]
                    continue
                if letter != record.seq[i] or record.seq[i] in ambiguous_bases:
                    consensus = False
                    break
            if consensus and letter != "N":
                consensus_string += letter
            else:
                consensus_string += "*"
        assert len(first_string) == len(consensus_string)
        return consensus_string

    def interval_partition(self):
        """Return a list of intervals in which we have
        consensus sequence longer than min_match_length, and
        a list of the non-match intervals left."""
        match_intervals = []
        non_match_intervals = []
        match_count, match_start, non_match_start = 0, 0, 0

        logging.debug("consensus: %s" % self.consensus)
        for i in range(self.length):
            letter = self.consensus[i]
            if letter != "*":
                # In a match region.
                if match_count == 0:
                    match_start = i
                match_count += 1
            elif match_count > 0:
                # Have reached a non-match. Check if previous match string is long enough to add to match_regions
                match_string = remove_gaps(
                    self.consensus[match_start : match_start + match_count]
                )
                match_len = len(match_string)
                logging.debug("have match string %s" % match_string)

                if match_len >= self.min_match_length:
                    if non_match_start < match_start:
                        non_match_intervals.append([non_match_start, match_start - 1])
                        logging.debug(
                            f"add non-match interval [{non_match_start},{match_start - 1}]"
                        )
                    end = match_start + match_count - 1
                    match_intervals.append([match_start, end])
                    logging.debug(f"add match interval [{match_start},{end}]")
                    non_match_start = i
                match_count = 0
                match_start = non_match_start

        end = self.length - 1
        if self.length < self.min_match_length:
            # Special case: a short sequence can still get classified as a match interval
            added_interval = "match" if "*" in self.consensus else "non_match"
            if added_interval == "match":
                match_intervals.append([0, end])
            else:
                non_match_intervals.append([0, end])
            logging.debug(f"add whole short {added_interval} interval [0,{end}]")
            match_count = 0
            non_match_start = end + 1

        # At end add last intervals
        if match_count > 0:
            if match_count >= self.min_match_length:
                match_intervals.append([match_start, end])
                logging.debug(f"add final match interval [{match_start},{end}]")
                if non_match_start < match_start:
                    end = match_start - 1
        if match_count != self.length and non_match_start <= end:
            non_match_intervals.append([non_match_start, end])
            logging.debug(f"add non-match interval [{non_match_start},{end}]")

        # check all stretches of consensus are in an interval, and intervals don't overlap
        for i in range(self.length):
            count_match = 0
            for interval in match_intervals:
                if interval[0] <= i <= interval[1]:
                    count_match += 1
            count_non_match = 0
            for interval in non_match_intervals:
                if interval[0] <= i <= interval[1]:
                    count_non_match += 1

            assert count_match | count_non_match, (
                "Failed to correctly identify match intervals: position %d "
                "appeared in both/neither match and non-match intervals" % i
            )
            assert count_match + count_non_match == 1, (
                "Failed to correctly identify match intervals: position "
                "%d appeared in %d intervals" % (i, count_match + count_non_match)
            )

        return match_intervals, non_match_intervals

    def check_nonmatch_intervals(self):
        """
        Goes through non-match intervals and makes sure there is more than one sequence there, else makes it a match
        interval.
        Example reasons for such a conversion to occur:
            - 'N' in a sequence causes it to be filtered out, and left with a single useable sequence
            - '-' in sequences causes them to appear different, but they are the same
        """
        for i in reversed(range(len(self.non_match_intervals))):
            interval = self.non_match_intervals[i]
            interval_alignment = self.alignment[:, interval[0] : interval[1] + 1]
            interval_seqs = get_interval_seqs(interval_alignment)
            if len(interval_seqs) < 2:
                self.match_intervals.append(self.non_match_intervals[i])
                self.non_match_intervals.pop(i)
        self.match_intervals.sort()

    @classmethod
    def kmeans_cluster_seqs_in_interval(
        self,
        interval: List[int],
        alignment: MultipleSeqAlignment,
        min_match_length: int,
    ) -> List[List[str]]:
        """Divide sequences in interval into subgroups of similar sequences."""
        interval_alignment = alignment[:, interval[0] : interval[1] + 1]

        if interval[1] - interval[0] <= min_match_length:
            logging.info("Small variation site in interval %s \n", interval)
            interval_seqs = get_interval_seqs(interval_alignment)
            id_lists = [
                [
                    record.id
                    for record in alignment
                    if record.seq[interval[0] : interval[1] + 1].ungap("-") == seq
                ]
                for seq in interval_seqs
            ]
            return id_lists

        logging.debug(
            "Get kmeans partition of interval [%d, %d]", interval[0], interval[1]
        )

        seq_to_ids = defaultdict(list)
        small_seq_to_ids = defaultdict(list)

        for record in interval_alignment:
            seq = str(record.seq.ungap("-"))
            if len(seq) >= min_match_length:
                seq_to_ids[seq].append(record.id)
            else:
                small_seq_to_ids[seq].append(record.id)
        logging.debug(
            f"Add classes for {len(seq_to_ids)} long "
            f"and {len(small_seq_to_ids)} small sequences"
        )

        # The clustering is performed on unique sequences
        clustered_ids = []
        if len(seq_to_ids) <= 1:
            logging.info("<= 1 sequence >= min_match_len, no clustering to perform.")
            if len(seq_to_ids) == 1:
                clustered_ids = [list(seq_to_ids.values())[0]]

        else:
            # first transform sequences into kmer occurrence vectors using a dict
            logging.debug("First transform sequences into kmer occurrence vectors")

            # collect all kmers
            kmer_dict = {}
            n = 0
            for seq in seq_to_ids:
                for i in range(len(seq) - min_match_length + 1):
                    kmer = seq[i : i + min_match_length]
                    if kmer not in kmer_dict:
                        kmer_dict[kmer] = n
                        n += 1
            logging.debug(f"Found {n} kmers")

            # count all kmers
            seq_kmer_counts = np.zeros(shape=(len(seq_to_ids), n))
            for j, seq in enumerate(seq_to_ids):
                counts = np.zeros(n)
                for i in range(len(seq) - min_match_length + 1):
                    kmer = seq[i : i + min_match_length]
                    counts[kmer_dict[kmer]] += 1
                seq_kmer_counts[j] = counts

            # cluster sequences using kmeans
            logging.debug("Now cluster:")
            kmeans = KMeans(n_clusters=1, random_state=2).fit(seq_kmer_counts)
            pre_cluster_inertia = kmeans.inertia_

            cluster_inertia = pre_cluster_inertia
            number_of_clusters = 1
            logging.debug(f"initial inertia: {cluster_inertia}")
            while (
                cluster_inertia > 0
                and cluster_inertia > pre_cluster_inertia / 2
                and number_of_clusters < len(seq_to_ids)
            ):
                number_of_clusters += 1
                kmeans = KMeans(n_clusters=number_of_clusters, random_state=2).fit(
                    seq_kmer_counts
                )
                cluster_inertia = kmeans.inertia_
                logging.debug(
                    "number of clusters: %d, inertia: %f",
                    number_of_clusters,
                    cluster_inertia,
                )

            # convert cluster numbers to sequence record IDs
            logging.debug("Extract equivalence classes from this partition")
            if pre_cluster_inertia > 0:
                cluster_ids = list(kmeans.predict(seq_kmer_counts))
                for i in range(max(cluster_ids) + 1):
                    clustered_ids.append([])
                all_ids = list(seq_to_ids.values())
                for i, cluster_id in enumerate(cluster_ids):
                    clustered_ids[cluster_id].extend(all_ids[i])
            else:
                logging.debug("pre_cluster_inertia is 0! No clustering.")
                clustered_ids = list(seq_to_ids.values())

        logging.debug("Merge id lists for the partitions")
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

    @classmethod
    def get_sub_alignment_by_list_id(
        self, id_list: List[str], alignment: MultipleSeqAlignment, interval=None
    ):
        list_records = [record for record in alignment if record.id in id_list]
        sub_alignment = MultipleSeqAlignment(list_records)
        if interval:
            sub_alignment = sub_alignment[:, interval[0] : interval[1] + 1]
        return sub_alignment

    def get_prg(self):
        prg = ""
        # last_char = None
        # skip_char = False

        for interval in self.all_intervals:
            if interval in self.match_intervals:
                # all seqs are not necessarily exactly the same: some can have 'N'
                # thus still process all of them, to get the one with no 'N'.
                sub_alignment = self.alignment[:, interval[0] : interval[1] + 1]
                seqs = get_interval_seqs(sub_alignment)
                assert 0 < len(seqs) <= 1, "Got >1 filtered sequences in match interval"
                seq = seqs[0]
                prg += seq

            else:
                # Define variant site number and increment for next available
                site_num = self.site
                self.site += 2
                variant_prgs = []

                # Define the variant seqs to add
                if (self.nesting_level == self.max_nesting) or (
                    interval[1] - interval[0] <= self.min_match_length
                ):
                    # Have reached max nesting level, just add all variants in interval.
                    logging.debug(
                        "Have reached max nesting level or have a small variant site, so add all variant "
                        "sequences in interval."
                    )
                    sub_alignment = self.alignment[:, interval[0] : interval[1] + 1]
                    logging.debug(
                        "Variant seqs found: %s"
                        % list(
                            remove_duplicates(
                                [str(record.seq) for record in sub_alignment]
                            )
                        )
                    )
                    variant_prgs = get_interval_seqs(sub_alignment)
                    logging.debug("Which is equivalent to: %s" % variant_prgs)
                else:
                    # divide sequences into subgroups and define prg for each subgroup.
                    logging.debug(
                        "Divide sequences into subgroups and define prg for each subgroup."
                    )
                    recur = True
                    id_lists = self.kmeans_cluster_seqs_in_interval(
                        interval, self.alignment, self.min_match_length
                    )
                    list_sub_alignments = [
                        self.get_sub_alignment_by_list_id(
                            id_list, self.alignment, interval
                        )
                        for id_list in id_lists
                    ]
                    num_clusters = len(id_lists)

                    if len(list_sub_alignments) == self.num_seqs:
                        logging.debug(
                            "Clustering did not group any sequences together, each seq is a cluster"
                        )
                        recur = False
                    elif interval[0] not in self.subAlignedSeqs:
                        self.subAlignedSeqs[interval[0]] = []
                        logging.debug(
                            "subAlignedSeqs now has keys: %s",
                            list(self.subAlignedSeqs.keys()),
                        )
                    else:
                        logging.debug(
                            "subAlignedSeqs already had key %d in keys: %s. This shouldn't happen.",
                            interval[0],
                            list(self.subAlignedSeqs.keys()),
                        )

                    while len(list_sub_alignments) > 0:
                        sub_alignment = list_sub_alignments.pop(0)
                        sub_aligned_seq = AlignedSeq(
                            msa_file=self.msa_file,
                            alignment_format=self.alignment_format,
                            max_nesting=self.max_nesting,
                            nesting_level=self.nesting_level + 1,
                            min_match_length=self.min_match_length,
                            site=self.site,
                            alignment=sub_alignment,
                            interval=interval,
                        )
                        variant_prgs.append(sub_aligned_seq.prg)
                        self.site = sub_aligned_seq.site

                        if recur:
                            # logging.debug("None not in snp_scores - try to add sub__aligned_seq to list in
                            # dictionary")
                            self.subAlignedSeqs[interval[0]].append(sub_aligned_seq)
                            # logging.debug("Length of subAlignedSeqs[%d] is %d", interval[0],
                            # len(self.subAlignedSeqs[interval[0]]))
                    assert num_clusters == len(variant_prgs), (
                        "I don't seem to have a sub-prg sequence for all parts of the partition - there are %d "
                        "classes in partition, and %d variant seqs"
                        % (num_clusters, len(variant_prgs))
                    )
                assert len(variant_prgs) > 1, "Only have one variant seq"

                assert len(variant_prgs) == len(
                    list(remove_duplicates(variant_prgs))
                ), "have repeat variant seqs"

                # Add the variant seqs to the prg
                prg += "%s%d%s" % (
                    self.delim_char,
                    site_num,
                    self.delim_char,
                )  # considered making it so start of prg was not delim_char,
                # but that would defeat the point if it
                while len(variant_prgs) > 1:
                    prg += variant_prgs.pop(0)
                    prg += "%s%d%s" % (self.delim_char, site_num + 1, self.delim_char)
                prg += variant_prgs.pop()
                prg += "%s%d%s" % (self.delim_char, site_num, self.delim_char)

        return prg

    @property
    def max_nesting_level_reached(self):
        max_nesting = []
        if self.subAlignedSeqs == {}:
            logging.debug(
                "self.subAlignedSeqs == {} at nesting level %d for interval %s",
                self.nesting_level,
                self.interval,
            )
            max_nesting.append(self.nesting_level)
        else:
            logging.debug(
                "self.subAlignedSeqs.keys(): %s", list(self.subAlignedSeqs.keys())
            )
            logging.debug(
                "self.subAlignedSeqs[self.subAlignedSeqs.keys()[0]]: %s",
                self.subAlignedSeqs[list(self.subAlignedSeqs.keys())[0]],
            )
            for interval_start in list(self.subAlignedSeqs.keys()):
                logging.debug("interval start: %d", interval_start)
                for subaseq in self.subAlignedSeqs[interval_start]:
                    logging.debug(
                        "type of subAlignedSeqs object in list: %s", type(subaseq)
                    )
                    recur = subaseq.max_nesting_level_reached
                    logging.debug(
                        "recur max level nesting returned: %d, which has type %s",
                        recur,
                        type(recur),
                    )
                    max_nesting.append(recur)
        m = max(max_nesting)
        logging.debug("found the max of %s is %d", max_nesting, m)
        return m

    @property
    def prop_in_match_intervals(self):
        length_match_intervals = 0
        for interval in self.match_intervals:
            length_match_intervals += interval[1] - interval[0] + 1
        return length_match_intervals / float(self.length)

    @property
    def num_seqs(self):
        return len(self.alignment)
