import logging
from typing import List

from make_prg.from_msa import MSA
from make_prg.io_utils import load_alignment_file
from make_prg.from_msa.cluster_sequences import kmeans_cluster_seqs_in_interval
from make_prg.seq_utils import (
    ambiguous_bases,
    remove_duplicates,
    get_interval_seqs,
    NONMATCH,
)
from make_prg.from_msa.interval_partition import IntervalPartitioner


class AlignedSeq(object):
    """
    Object based on a set of aligned sequences.
    Note min_match_length must be strictly greater than max_nesting + 1.
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
        self.alignment: MSA = alignment
        if self.alignment is None:
            self.alignment = load_alignment_file(msa_file, alignment_format)

        self.interval = interval
        self.consensus = self.get_consensus(self.alignment)
        self.length = len(self.consensus)
        (
            self.match_intervals,
            self.non_match_intervals,
            self.all_intervals,
        ) = IntervalPartitioner(
            self.consensus, self.min_match_length, self.alignment
        ).get_intervals()
        logging.info(
            "match intervals: %s; non_match intervals: %s",
            self.match_intervals,
            self.non_match_intervals,
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

    @classmethod
    def get_consensus(cls, alignment: MSA):
        """ Produces a 'consensus string' from an MSA: at each position of the
        MSA, the string has a base if all aligned sequences agree, and a "*" if not.
        IUPAC ambiguous bases result in non-consensus and are later expanded in the prg.
        N results in consensus at that position unless they are all N."""
        consensus_string = ""
        for i in range(alignment.get_alignment_length()):
            column = set([record.seq[i] for record in alignment])
            column = column.difference({"N"})
            if (
                len(ambiguous_bases.intersection(column)) > 0
                or len(column) != 1
                or column == {"-"}
            ):
                consensus_string += NONMATCH
            else:
                consensus_string += column.pop()

        return consensus_string

    @classmethod
    def get_sub_alignment_by_list_id(
        self, id_list: List[str], alignment: MSA, interval=None
    ):
        list_records = [record for record in alignment if record.id in id_list]
        sub_alignment = MSA(list_records)
        if interval:
            sub_alignment = sub_alignment[:, interval[0] : interval[1] + 1]
        return sub_alignment

    def get_prg(self):
        prg = ""

        for interval in self.all_intervals:
            if interval in self.match_intervals:
                # all seqs are not necessarily exactly the same: some can have 'N'
                # thus still process all of them, to get the one with no 'N'.
                sub_alignment = self.alignment[:, interval.start : interval.stop + 1]
                seqs = get_interval_seqs(sub_alignment)
                assert len(seqs) == 1, "Got >1 filtered sequences in match interval"
                seq = seqs[0]
                prg += seq

            else:
                # Define variant site number and increment for next available
                site_num = self.site
                self.site += 2
                variant_prgs = []

                # Define the variant seqs to add
                if (self.nesting_level == self.max_nesting) or (
                    interval.stop - interval.start <= self.min_match_length
                ):
                    logging.debug(
                        "Have reached max nesting level or have a small variant site, so add all variant "
                        "sequences in interval."
                    )
                    sub_alignment = self.alignment[
                        :, interval.start : interval.stop + 1
                    ]
                    variant_prgs = get_interval_seqs(sub_alignment)
                    logging.debug(f"Variant seqs found: {variant_prgs}")
                else:
                    logging.debug(
                        "Divide sequences into subgroups and define prg for each subgroup."
                    )
                    recur = True
                    id_lists = kmeans_cluster_seqs_in_interval(
                        [interval.start, interval.stop],
                        self.alignment,
                        self.min_match_length,
                    )
                    list_sub_alignments = [
                        self.get_sub_alignment_by_list_id(
                            id_list, self.alignment, [interval.start, interval.stop]
                        )
                        for id_list in id_lists
                    ]
                    num_clusters = len(id_lists)

                    if len(list_sub_alignments) == self.num_seqs:
                        logging.debug(
                            "Clustering did not group any sequences together, each seq is a cluster"
                        )
                        recur = False
                    elif interval.start not in self.subAlignedSeqs:
                        self.subAlignedSeqs[interval.start] = []
                        logging.debug(
                            "subAlignedSeqs now has keys: %s",
                            list(self.subAlignedSeqs.keys()),
                        )
                    else:
                        logging.debug(
                            "subAlignedSeqs already had key %d in keys: %s. This shouldn't happen.",
                            interval.start,
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
                            self.subAlignedSeqs[interval.start].append(sub_aligned_seq)
                    assert num_clusters == len(variant_prgs), (
                        "I don't seem to have a sub-prg sequence for all parts of the partition - there are %d "
                        "classes in partition, and %d variant seqs"
                        % (num_clusters, len(variant_prgs))
                    )
                assert len(variant_prgs) > 1, "Only have one variant seq"

                assert len(variant_prgs) == len(
                    list(remove_duplicates(variant_prgs))
                ), "have repeat variant seqs"

                # Add the variant seqs to the prg.
                prg += f"{self.delim_char}{site_num}{self.delim_char}"
                while len(variant_prgs) > 1:
                    prg += variant_prgs.pop(0)
                    prg += f"{self.delim_char}{site_num + 1}{self.delim_char}"
                prg += variant_prgs.pop()
                prg += f"{self.delim_char}{site_num}{self.delim_char}"

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
            length_match_intervals += interval.stop - interval.start + 1
        return length_match_intervals / float(self.length)

    @property
    def num_seqs(self):
        return len(self.alignment)
