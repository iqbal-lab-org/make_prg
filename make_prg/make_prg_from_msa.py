import gzip
import logging
import os

import numpy as np
from Bio import AlignIO
from Bio.AlignIO import MultipleSeqAlignment
from sklearn.cluster import KMeans

from make_prg.utils import remove_duplicates, remove_gaps


def get_interval_seqs(interval_alignment):
    """Replace - with nothing, remove seqs containing N or other non-allowed letters
    and duplicate sequences containing RYKMSW, replacing with AGCT alternatives """
    allowed_bases = {"A", "C", "G", "T", "R", "Y", "K", "M", "S", "W"}
    iupac = {
        "R": ["G", "A"],
        "Y": ["T", "C"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "S": ["G", "C"],
        "W": ["A", "T"],
    }
    seqs = []
    gapless_seqs = [
        remove_gaps(str(record.seq)).upper() for record in interval_alignment
    ]
    unique_seqs = remove_duplicates(gapless_seqs)

    for seq in unique_seqs:
        if allowed_bases.issuperset(seq):
            new_seqs = [seq]
            for letter in iupac.keys():
                letter_seqs = []
                for t in new_seqs:
                    if letter in t:
                        letter_seqs.append(t.replace(letter, iupac[letter][0]))
                        letter_seqs.append(t.replace(letter, iupac[letter][1]))
                    else:
                        letter_seqs.append(t)
                new_seqs = letter_seqs
            seqs.extend(new_seqs)
    ret_list = list(set(seqs))
    if len(ret_list) == 0:
        logging.warning(
            "WARNING: Every sequence must have contained an N in this slice - redo sequence curation because this is nonsense"
        )
        logging.warning(
            "Sequences were",
            " ".join(
                list(
                    remove_duplicates(
                        [
                            remove_gaps(str(record.seq)).upper()
                            for record in interval_alignment
                        ]
                    )
                )
            ),
        )
        logging.warning(
            "Using these sequences anyway, and should be ignored downstream"
        )
        seqs = list(
            remove_duplicates(
                [remove_gaps(str(record.seq)).upper() for record in interval_alignment]
            )
        )
    return sorted(list(set(seqs)))


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
        staggered_start=False,
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
        self.staggered_start = staggered_start
        self.alignment = alignment
        if not self.alignment:
            logging.info("Read from MSA file %s", self.msa_file)
            if ".gz" in self.msa_file:
                logging.debug("MSA is gzipped")
                handle = gzip.open(self.msa_file, "rt")
                self.alignment = AlignIO.read(handle, self.alignment_format)
                handle.close()
            else:
                self.alignment = AlignIO.read(self.msa_file, self.alignment_format)
        self.interval = interval
        self.num_seqs = len(self.alignment)
        self.consensus = self.get_consensus()
        self.length = len(self.consensus)
        (self.match_intervals, self.non_match_intervals) = self.get_match_intervals
        self.check_nonmatch_intervals()
        self.all_intervals = self.match_intervals + self.non_match_intervals
        logging.info("Non match intervals: %s", self.non_match_intervals)
        self.all_intervals.sort()
        if self.nesting_level == 1:
            self.length_match_intervals = 0
            for interval in self.match_intervals:
                self.length_match_intervals += interval[1] - interval[0] + 1
            self.prop_in_match_intervals = self.length_match_intervals / float(
                self.length
            )

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
        self.kmer_dict = {}

    def get_consensus(self):
        """Given a set of aligment records from AlignIO, creates
        a consensus string.
        Lower and upper case are equivalent
        Non AGCT symbols RYKMSW result in non-consensus and are substituted in graph
        N results in consensus at that position unless they are all N."""
        missing = ["N"]
        if self.staggered_start:
            missing = ["N", "-"]
        first_string = str(self.alignment[0].seq)
        consensus_string = ""
        for i, letter in enumerate(first_string):
            consensus = True
            if letter.upper() not in missing:
                missing = ["N"]
            for record in self.alignment:
                if (record.seq[i].upper() not in missing and letter.upper() not in missing) and (
                    record.seq[i].upper() != letter.upper()
                    or record.seq[i].upper() in ["R", "Y", "K", "M", "S", "W"]
                ):
                    consensus = False
                    break
                if letter.upper() in missing and record.seq[i].upper() not in missing:
                    letter = record.seq[i].upper()
            if consensus and letter.upper() not in missing:
                consensus_string += letter
            else:
                consensus_string += "*"
        assert len(first_string) == len(consensus_string)
        return consensus_string

    @property
    def get_match_intervals(self):
        """Return a list of intervals in which we have
        consensus sequence longer than min_match_length, and
        a list of the non-match intervals left."""
        match_intervals = []
        non_match_intervals = []
        match_count = 0
        match_start = 0
        non_match_start = 0

        logging.debug("consensus: %s" % self.consensus)
        if len(remove_gaps(self.consensus)) < self.min_match_length:
            # It makes no sense to classify a fully consensus sequence as
            # a non-match just because it is too short.
            if "*" in self.consensus:
                interval_alignment = self.alignment[:, 0 : self.length]
                interval_seqs = get_interval_seqs(interval_alignment)
                if len(interval_seqs) > 1:
                    logging.debug(
                        "add short non-match whole interval [%d,%d]"
                        % (0, self.length - 1)
                    )
                    non_match_intervals.append([0, self.length - 1])
                else:
                    logging.debug(
                        "add short match whole interval [%d,%d]" % (0, self.length - 1)
                    )
                    match_intervals.append([0, self.length - 1])
            else:
                match_intervals.append([0, self.length - 1])
                logging.debug(
                    "add short match whole interval [%d,%d]" % (0, self.length - 1)
                )
        else:
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
                        # if the non_match sequences in the interval are really the same, add a match interval
                        interval_alignment = self.alignment[
                            :, non_match_start : match_start + 1
                        ]
                        interval_seqs = get_interval_seqs(interval_alignment)
                        if non_match_start < match_start and len(interval_seqs) > 1:
                            non_match_intervals.append(
                                [non_match_start, match_start - 1]
                            )
                            logging.debug(
                                "add non-match interval as have alts [%d,%d]"
                                % (non_match_start, match_start - 1)
                            )
                        elif non_match_start < match_start:
                            match_intervals.append([non_match_start, match_start - 1])
                            logging.debug(
                                "add match interval as only one seq [%d,%d]"
                                % (non_match_start, match_start - 1)
                            )
                        match_intervals.append(
                            [match_start, match_start + match_count - 1]
                        )
                        logging.debug(
                            "add match interval to complete step [%d,%d]"
                            % (match_start, match_start + match_count - 1)
                        )
                        non_match_start = i
                    match_count = 0
                    match_start = non_match_start

            # At end add last intervals
            match_string = remove_gaps(
                self.consensus[match_start : match_start + match_count]
            )
            match_len = len(match_string)
            logging.debug("at end have match string %s" % match_string)
            if 0 < match_len < self.min_match_length:
                logging.debug(
                    "have short match region at end, so include it in non-match-region before - "
                    "match count was %d" % match_count
                )
                match_count = 0
                match_start = non_match_start
                logging.debug("match count is now %d" % match_count)

            if match_count > 0:
                interval_alignment = self.alignment[
                    :, non_match_start : match_start + 1
                ]
            else:
                interval_alignment = self.alignment[:, non_match_start : self.length]
            interval_seqs = get_interval_seqs(interval_alignment)
            if len(interval_seqs) == 1:
                match_intervals.append([non_match_start, self.length - 1])
                logging.debug(
                    "add match interval at end as only one seq [%d,%d]"
                    % (non_match_start, self.length - 1)
                )
            elif len(interval_seqs) > 1 and non_match_start < match_start:
                non_match_intervals.append([non_match_start, match_start - 1])
                logging.debug(
                    "add non-match interval at end as have alts [%d,%d]"
                    % (non_match_start, match_start - 1)
                )
                match_intervals.append([match_start, self.length - 1])
                logging.debug(
                    "add match interval at end [%d,%d]" % (match_start, self.length - 1)
                )
            else:
                non_match_intervals.append([non_match_start, self.length - 1])
                logging.debug(
                    "add only non-match interval at end as have alts [%d,%d]"
                    % (non_match_start, self.length - 1)
                )

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
        """Goes through non-match intervals and makes sure there is more than one sequence there, else makes it a match
        interval."""
        for i in reversed(range(len(self.non_match_intervals))):
            interval = self.non_match_intervals[i]
            interval_alignment = self.alignment[:, interval[0] : interval[1] + 1]
            interval_seqs = get_interval_seqs(interval_alignment)
            if len(interval_seqs) < 2:
                self.match_intervals.append(self.non_match_intervals[i])
                self.non_match_intervals.pop(i)
        self.match_intervals.sort()

    def kmeans_cluster_seqs_in_interval(
        self, interval
    ):  # , kmer_size=self.min_match_length):
        """Divide sequences in interval into subgroups of similar
           sequences. Return a list of lists of ids."""
        if interval[1] - interval[0] <= self.min_match_length:
            logging.info("Small variation site in interval %s \n", interval)
            logging.debug(
                "interval[1] - interval[0] <= self.min_match_length: %d <= %d",
                interval[1] - interval[0],
                self.min_match_length,
            )
            interval_alignment = self.alignment[:, interval[0] : interval[1] + 1]
            interval_seqs = get_interval_seqs(interval_alignment)
            assert len(interval_seqs) == len(
                list(remove_duplicates(interval_seqs))
            ), "should not have duplicate alternative allele sequences"
            return_id_lists = [
                [
                    record.id
                    for record in self.alignment
                    if remove_gaps(str(record.seq[interval[0] : interval[1] + 1]))
                    == seq
                ]
                for seq in interval_seqs
            ]
        else:
            logging.debug(
                "Get kmeans partition of interval [%d, %d]", interval[0], interval[1]
            )
            interval_alignment = self.alignment[:, interval[0] : interval[1] + 1]
            seq_dict_keys = get_interval_seqs(interval_alignment)

            interval_seq_dict = {}
            small_interval_seq_dict = {}

            for seq in seq_dict_keys:
                if len(seq) >= self.min_match_length:
                    interval_seq_dict[seq] = []
                else:
                    small_interval_seq_dict[seq] = []

            for record in interval_alignment:
                seq = remove_gaps(str(record.seq))
                if seq in list(interval_seq_dict.keys()):
                    interval_seq_dict[seq].append(record.id)
                elif seq in list(small_interval_seq_dict.keys()):
                    small_interval_seq_dict[seq].append(record.id)

            assert len(seq_dict_keys) == len(
                list(remove_duplicates(seq_dict_keys))
            ), "error, have duplicate dictionary keys"
            assert (
                len(
                    [
                        key
                        for key in list(interval_seq_dict.keys())
                        if key in list(small_interval_seq_dict.keys())
                    ]
                )
                == 0
            ), "error, should have no overlap of keys"
            assert (
                len(
                    [
                        key
                        for key in list(small_interval_seq_dict.keys())
                        if key in list(interval_seq_dict.keys())
                    ]
                )
                == 0
            ), "error, should have no overlap of keys"

            logging.debug(
                "Add classes corresponding to %d small sequences"
                % len(list(small_interval_seq_dict.keys()))
            )

            logging.debug(
                "Now add classes corresponding to %d longer sequences"
                % len(list(interval_seq_dict.keys()))
            )
            interval_seqs = list(interval_seq_dict.keys())
            big_return_id_lists = []
            if len(interval_seqs) > 1:
                # first transform sequences into kmer occurance vectors using a dict
                logging.debug("First transform sequences into kmer occurance vectors")

                # make dict based on number of kmers in all sequences
                self.kmer_dict = {}
                n = 0
                for j, seq in enumerate(interval_seqs):
                    for i in range(len(seq) - self.min_match_length + 1):
                        if seq not in list(self.kmer_dict.keys()):
                            self.kmer_dict[seq[i : i + self.min_match_length]] = n
                            n += 1
                logging.debug("These vectors have length %d" % n)

                # transform to vectors using dict
                seq_kmer_counts = np.zeros(shape=(len(interval_seqs), n))
                for j, seq in enumerate(interval_seqs):
                    counts = np.zeros(n)
                    for i in range(len(seq) - self.min_match_length + 1):
                        counts[self.kmer_dict[seq[i : i + self.min_match_length]]] += 1
                    seq_kmer_counts[j] = counts

                # cluster sequences using kmeans
                logging.debug("Now cluster:")
                kmeans = KMeans(n_clusters=1, random_state=2).fit(seq_kmer_counts)
                pre_cluster_inertia = kmeans.inertia_

                if pre_cluster_inertia == 0:
                    logging.debug("pre_cluster_intertia is 0!")
                    for key in list(interval_seq_dict.keys()):
                        logging.debug(
                            "seq: %s, num_seqs with this seq: %d",
                            key,
                            len(interval_seq_dict[key]),
                        )

                cluster_inertia = pre_cluster_inertia
                number_of_clusters = 1
                logging.debug(
                    "number of clusters: %d, inertia: %f",
                    number_of_clusters,
                    cluster_inertia,
                )
                while (
                    cluster_inertia > 0
                    and cluster_inertia > pre_cluster_inertia / 2
                    and number_of_clusters <= len(interval_seqs)
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

                # now extract the equivalence class details from this partition and return
                logging.debug("Extract equivalence classes from this partition")
                if pre_cluster_inertia > 0:
                    equiv_class_ids = list(kmeans.predict(seq_kmer_counts))
                    for i in range(max(equiv_class_ids) + 1):
                        big_return_id_lists.append([])
                    for i, val in enumerate(equiv_class_ids):
                        big_return_id_lists[val].extend(
                            interval_seq_dict[interval_seqs[i]]
                        )
                else:
                    logging.debug("default to not clustering")
                    big_return_id_lists = [
                        interval_seq_dict[key] for key in interval_seq_dict.keys()
                    ]
            elif len(interval_seqs) == 1:
                big_return_id_lists = [interval_seq_dict[interval_seqs[0]]]

            # now merge big and small return_id_lists so as to maintain the order of seqs before
            logging.debug("Merge return id lists for the partitions")
            return_id_lists = []
            added_ids = []
            big_keys = list(interval_seq_dict.keys())
            small_keys = list(small_interval_seq_dict.keys())
            for seq in seq_dict_keys:
                if seq in small_keys:
                    logging.debug(
                        "add (small) return ids: %s" % small_interval_seq_dict[seq]
                    )
                    return_id_lists.append(small_interval_seq_dict[seq])
                elif seq in big_keys:
                    not_added = [
                        nid for nid in interval_seq_dict[seq] if nid not in added_ids
                    ]
                    if len(not_added) == len(interval_seq_dict[seq]):
                        logging.debug(
                            "want to add (big) return ids: %s" % interval_seq_dict[seq]
                        )
                        for i in range(len(big_return_id_lists)):
                            if interval_seq_dict[seq][0] in big_return_id_lists[i]:
                                logging.debug(
                                    "add (big) return ids %d: %s"
                                    % (i, big_return_id_lists[i])
                                )
                                return_id_lists.append(big_return_id_lists[i])
                                added_ids.extend(return_id_lists[-1])
                                break
                    else:
                        assert (
                            len(not_added) == 0
                        ), "Equivalent sequences should be in same part of partition and are not"
                else:
                    logging.warning(
                        "Key %s doesn't seem to be in either big keys or small keys"
                    )
        assert (
            len(return_id_lists) > 1
        ), "should have some alternate alleles, not only one sequence, this is a non-match interval"
        return return_id_lists

    def get_sub_alignment_by_list_id(self, list_of_id, interval=None):
        list_records = [record for record in self.alignment if record.id in list_of_id]
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
                # WLOG can take first sequence as all same in this interval
                sub_alignment = self.alignment[:, interval[0] : interval[1] + 1]
                seqs = get_interval_seqs(sub_alignment)
                assert len(seqs) > 0
                seq = seqs[0]
                prg += seq

            else:
                # Define variant site number and increment for next available
                site_num = self.site
                self.site += 2
                variant_seqs = []

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
                    variant_seqs = get_interval_seqs(sub_alignment)
                    logging.debug("Which is equivalent to: %s" % variant_seqs)
                else:
                    # divide sequences into subgroups and define prg for each subgroup.
                    logging.debug(
                        "Divide sequences into subgroups and define prg for each subgroup."
                    )
                    recur = True
                    list_list_id = self.kmeans_cluster_seqs_in_interval(interval)
                    list_sub_alignments = [
                        self.get_sub_alignment_by_list_id(list_id, interval)
                        for list_id in list_list_id
                    ]
                    num_classes_in_partition = len(list_list_id)

                    if len(list_sub_alignments) == self.num_seqs:
                        logging.debug(
                            "Partition does not group any sequences together, all seqs get unique class in partition"
                        )
                        recur = False
                    elif interval[0] not in list(self.subAlignedSeqs.keys()):
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
                        sub__aligned_seq = AlignedSeq(
                            msa_file=self.msa_file,
                            alignment_format=self.alignment_format,
                            max_nesting=self.max_nesting,
                            nesting_level=self.nesting_level + 1,
                            min_match_length=self.min_match_length,
                            site=self.site,
                            alignment=sub_alignment,
                            interval=interval,
                        )
                        variant_seqs.append(sub__aligned_seq.prg)
                        self.site = sub__aligned_seq.site

                        if recur:
                            # logging.debug("None not in snp_scores - try to add sub__aligned_seq to list in
                            # dictionary")
                            self.subAlignedSeqs[interval[0]].append(sub__aligned_seq)
                            # logging.debug("Length of subAlignedSeqs[%d] is %d", interval[0],
                            # len(self.subAlignedSeqs[interval[0]]))
                    assert num_classes_in_partition == len(variant_seqs), (
                        "I don't seem to have a sub-prg sequence for all parts of the partition - there are %d "
                        "classes in partition, and %d variant seqs"
                        % (num_classes_in_partition, len(variant_seqs))
                    )
                assert len(variant_seqs) > 1, "Only have one variant seq"

                if len(variant_seqs) != len(list(remove_duplicates(variant_seqs))):
                    print("variant_seqs: ")
                    for s in variant_seqs:
                        print(s)
                        print(", ")

                assert len(variant_seqs) == len(
                    list(remove_duplicates(variant_seqs))
                ), "have repeat variant seqs"

                # Add the variant seqs to the prg
                prg += "%s%d%s" % (
                    self.delim_char,
                    site_num,
                    self.delim_char,
                )  # considered making it so start of prg was not delim_char,
                # but that would defeat the point if it
                while len(variant_seqs) > 1:
                    prg += variant_seqs.pop(0)
                    prg += "%s%d%s" % (self.delim_char, site_num + 1, self.delim_char)
                prg += variant_seqs.pop()
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
