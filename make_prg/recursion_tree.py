"""
This module contains the classes to explicitly represent the make_prg process recursion tree
"""

from typing import List, Set, Optional, Tuple
from loguru import logger
from make_prg.from_msa import MSA
from make_prg.from_msa.cluster_sequences import kmeans_cluster_seqs, ClusteringResult
from make_prg.utils.seq_utils import (
    SequenceExpander,
    remove_columns_full_of_gaps_from_MSA,
    get_consensus_from_MSA,
    get_number_of_unique_ungapped_sequences,
    get_number_of_unique_gapped_sequences
)
from make_prg.from_msa.interval_partition import IntervalPartitioner, Interval, Intervals
from make_prg.update.denovo_variants import UpdateData
from make_prg.update.MLPath import MLPathError
from abc import ABC, abstractmethod

SubMSAs = List[MSA]


class RecursiveTreeNode(ABC):
    def __init__(self, nesting_level: int, alignment: MSA,
                 parent: Optional["RecursiveTreeNode"], prg_builder: "PrgBuilder",
                 children_subalignments: SubMSAs):
        """
        Builds a new recursive tree node.
        Notes:
            1. To build a node you need to know in advance how it is going to cluster (horizontally or vertically).
                This is given in the children_subalignments param. NodeFactory.build() takes care of this (see below).
        """
        self.nesting_level: int = nesting_level
        self.alignment: MSA = remove_columns_full_of_gaps_from_MSA(alignment)
        self.parent: Optional["RecursiveTreeNode"] = parent
        self.prg_builder = prg_builder
        self.id: int = self.prg_builder.get_next_node_id()

        # generate recursion tree
        self._children: List["RecursiveTreeNode"] = self._get_children(children_subalignments)

        self.log_that_node_was_created()

    @property
    def children(self):
        return self._children

    def _get_children(self, children_subalignments: SubMSAs) -> List["RecursiveTreeNode"]:
        return [NodeFactory.build(alignment, self.prg_builder, self) for alignment in children_subalignments]

    @abstractmethod
    def preorder_traversal_to_build_prg(self, prg_as_list: List[str], delim_char: str = " "):
        pass

    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def is_root(self) -> bool:
        return self.parent is None

    def replace_child(self, old_child: "RecursiveTreeNode", new_child: "RecursiveTreeNode"):
        old_child_is_one_of_the_children = old_child in self.children
        assert old_child_is_one_of_the_children, f"Failure to replace a child, {old_child} does not exist"

        old_child_index = self.children.index(old_child)
        self.children[old_child_index] = new_child

    # Note: trivial method, untested
    def log_that_node_was_created(self):
        logger.trace("Created node:\n" + str(self))

    def __repr__(self):
        return f"{self.__class__.__name__}:\n" \
               f"Id = {self.id}\n" \
               f"Nesting level = {self.nesting_level}\n" \
               f"Parent = {'None' if self.parent is None else f'Id = {self.parent.id}'}\n" \
               f"Children = [{', '.join(f'Id = {child.id}' for child in self.children)}]\n" \
               f"Alignment:\n{format(self.alignment, 'fasta')}"

    def __str__(self):
        return repr(self)


class MultiIntervalNode(RecursiveTreeNode):
    """
    Represent nodes that are clustered vertically, i.e. in intervals
    """
    def __init__(self, nesting_level: int, alignment: MSA, parent: Optional["RecursiveTreeNode"],
                 prg_builder: "PrgBuilder", interval_subalignments: SubMSAs):
        super().__init__(nesting_level, alignment, parent, prg_builder, interval_subalignments)
        assert not self.is_leaf(), "MultiIntervalNodes should never be leaves"

    def preorder_traversal_to_build_prg(self, prg_as_list: List[str], delim_char: str = " "):
        for child in self.children:
            child.preorder_traversal_to_build_prg(prg_as_list, delim_char)


class MultiClusterNode(RecursiveTreeNode):
    """
    Represent nodes that are clustered horizontally, i.e. in sequence clusters
    """
    def __init__(self, nesting_level: int, alignment: MSA, parent: Optional["RecursiveTreeNode"],
                 prg_builder: "PrgBuilder", cluster_subalignments: SubMSAs):
        super().__init__(nesting_level, alignment, parent, prg_builder, cluster_subalignments)
        assert not self.is_leaf(), "MultiClusterNodes should never be leaves"

    def preorder_traversal_to_build_prg(self, prg_as_list: List[str], delim_char: str = " "):
        site_num = self.prg_builder.get_next_site_num()
        prg_as_list.extend(f"{delim_char}{site_num}{delim_char}")

        for child_index, child in enumerate(self.children):
            site_num_to_separate_alleles = (
                (site_num + 1) if (child_index < len(self.children) - 1) else site_num
            )
            child.preorder_traversal_to_build_prg(prg_as_list, delim_char)
            prg_as_list.extend(
                f"{delim_char}{site_num_to_separate_alleles}{delim_char}"
            )


class UpdateError(Exception):
    pass


class LeafNode(RecursiveTreeNode):
    """
    Represent nodes that are never clustered in either direction.
    These nodes are the only ones that can get updated and indexed.
    """
    def __init__(self, nesting_level: int, alignment: MSA, parent: Optional["RecursiveTreeNode"],
                 prg_builder: "PrgBuilder"):
        super().__init__(nesting_level, alignment, parent, prg_builder, [])
        assert self.is_leaf(), f"Leaf node ({self}) is not a leaf"
        self.new_sequences: Set[str] = set()
        self.indexed_PRG_intervals: Set[Tuple[int, int]] = set()

    def preorder_traversal_to_build_prg(self, prg_as_list: List[str], delim_char: str = " ", do_indexing=True):
        # Note: here we could have several intervals, not anymore
        expanded_sequences = SequenceExpander.get_expanded_sequences_from_MSA(self.alignment)

        single_seq = len(expanded_sequences) == 1
        if single_seq:
            start_index = len(prg_as_list)
            prg_as_list.extend(expanded_sequences[0])
            end_index = len(prg_as_list)
            if do_indexing:
                self.prg_builder.update_PRG_index(start_index, end_index, node=self)
        else:
            # Add the variant seqs to the prg
            site_num = self.prg_builder.get_next_site_num()
            prg_as_list.extend(f"{delim_char}{site_num}{delim_char}")
            for seq_index, seq in enumerate(expanded_sequences):
                site_num_for_this_seq = (
                    (site_num + 1) if (seq_index < len(expanded_sequences) - 1) else site_num
                )
                start_index = len(prg_as_list)
                prg_as_list.extend(seq)
                end_index = len(prg_as_list)

                prg_as_list.extend(
                    f"{delim_char}{site_num_for_this_seq}{delim_char}"
                )

                if do_indexing:
                    self.prg_builder.update_PRG_index(start_index, end_index, node=self)

    ##################################################################################
    # update methods
    def add_data_to_batch_update(self, update_data: UpdateData):
        update_data_PRG_interval = update_data.ml_path_node_key
        update_data_PRG_interval_is_indexed = update_data_PRG_interval in self.indexed_PRG_intervals
        if not update_data_PRG_interval_is_indexed:
            raise UpdateError(f"PRG interval {update_data_PRG_interval} not found in indexed "
                              f"PRG intervals for node: {self.indexed_PRG_intervals}")

        seq_to_add_as_list = []
        for PRG_interval in sorted(self.indexed_PRG_intervals):
            ML_sequence_should_be_replaced = PRG_interval == update_data_PRG_interval
            if ML_sequence_should_be_replaced:
                seq_to_add_as_list.append(update_data.new_node_sequence)
            else:
                try:
                    ml_path_node = update_data.ml_path.get_node_given_interval_in_PRG_space(PRG_interval)
                    seq_to_add_as_list.append(ml_path_node.sequence)
                except MLPathError:
                    pass
        seq_to_add = "".join(seq_to_add_as_list)

        there_has_been_padding = seq_to_add != update_data.new_node_sequence
        if there_has_been_padding:
            logger.trace(f"Sequence {update_data.new_node_sequence} padded to {seq_to_add} on "
                         f"add_data_to_batch_update() for {self.prg_builder.locus_name}")

        self.new_sequences.add(seq_to_add)

    def add_indexed_PRG_interval(self, interval: Tuple[int, int]):
        self.indexed_PRG_intervals.add(interval)

    def batch_update(self):
        no_update_to_be_done = len(self.new_sequences) == 0
        if no_update_to_be_done:
            return
        self._update_leaf()

    def _update_leaf(self):
        logger.trace(f"Updating subMSA for {self.prg_builder.locus_name}")
        logger.trace(f"Node: {str(self)}")
        logger.trace(f"Sequences added to update: {self.new_sequences}")

        an_aligner_was_given = self.prg_builder.aligner is not None
        assert an_aligner_was_given, "Cannot make updates without a Multiple Sequence Aligner."

        updated_alignment = self.prg_builder.aligner.get_updated_alignment(
            current_alignment=self.alignment,
            new_sequences=self.new_sequences
        )

        # create a new, updated node, and add it to the tree
        updated_child = NodeFactory.build(updated_alignment, self.prg_builder, self.parent)

        # in some cases, we might have just a single node, which is both the root and a leaf
        need_to_replace_the_root = self.is_root()
        if need_to_replace_the_root:
            self.prg_builder.replace_root(updated_child)
        else:
            self.parent.replace_child(self, updated_child)

        # whenever we do an update in a node, we will mess up the prg builder index for the sites to the right of the
        # updated site in the PRG string. Just to be safe, we clear the prg index, as it will be regenerated again
        # when rebuilding the PRG
        self.prg_builder.clear_PRG_index()

    def clear_PRG_interval_index(self):
        self.indexed_PRG_intervals.clear()
    ##################################################################################


class NodeFactory:
    """
    Class responsible to build the different RecursiveTreeNodes depending on its alignment and other parameters.
    """
    @staticmethod
    def build(alignment: MSA, prg_builder: "PrgBuilder", parent_node: Optional[RecursiveTreeNode] = None) -> RecursiveTreeNode:
        """
        Method that builds the correct node given the alignment and other parameters.
        Note that the root has a None parent. For the root also, we need to first try to build a MultiIntervalNode and
        then a MulticlusterNode,as we first want to split its intervals, and then cluster.
        """
        parent_node_is_leaf = isinstance(parent_node, LeafNode)
        assert not parent_node_is_leaf, "Error on building a Recursive Tree node: parent node should never be a leaf"

        min_match_length = prg_builder.min_match_length
        nesting_level = NodeFactory._get_nesting_level(parent_node)

        # if parent is multi interval, children should not be
        skip_building_multi_interval_node = isinstance(parent_node, MultiIntervalNode)
        if not skip_building_multi_interval_node:
            all_intervals, match_intervals = NodeFactory._get_all_and_match_intervals(alignment, min_match_length)

            has_a_single_interval = NodeFactory._infer_if_has_single_interval(all_intervals=all_intervals, match_intervals=match_intervals,
                                                                nesting_level=nesting_level,
                                                                max_nesting=prg_builder.max_nesting,
                                                                alignment=alignment,
                                                                min_match_length=min_match_length)

            if not has_a_single_interval:
                interval_subalignments = NodeFactory._partition_alignment_into_interval_subalignments(alignment,
                                                                                                      all_intervals)
                return MultiIntervalNode(nesting_level, alignment, parent_node, prg_builder, interval_subalignments)

        # if parent is multi cluster, children should not be
        skip_building_multi_cluster_node = isinstance(parent_node, MultiClusterNode)
        if not skip_building_multi_cluster_node:
            clustering_result = kmeans_cluster_seqs(alignment, min_match_length)
            cluster_further = NodeFactory._infer_if_we_should_cluster_further(alignment, clustering_result)

            if cluster_further:
                cluster_subalignments = NodeFactory._get_subalignments_by_clustering(alignment, clustering_result)
                return MultiClusterNode(nesting_level, alignment, parent_node, prg_builder,
                                        cluster_subalignments)

        # here the node is a leaf
        return LeafNode(nesting_level, alignment, parent_node, prg_builder)

    #####################################################################################################
    #  helper methods
    @staticmethod
    def _alignment_has_issues(alignment: MSA) -> bool:
        num_unique_nongapped_seqs = get_number_of_unique_ungapped_sequences(alignment)
        too_few_unique_sequences = num_unique_nongapped_seqs <= 2
        if too_few_unique_sequences:
            return True

        num_unique_gapped_seqs = get_number_of_unique_gapped_sequences(alignment)

        # this is an assert as if it happens, something very wrong with the two previous functions is happening
        assert num_unique_nongapped_seqs <= num_unique_gapped_seqs

        alignment_has_ambiguity = num_unique_nongapped_seqs < num_unique_gapped_seqs
        if alignment_has_ambiguity:
            # TODO: fix alignment by deduplicating nongapped seqs and realigning them
            # TODO: the annoying thing is that this deduplicated nongapped alignment will differ from the input alignment
            # TODO: but I think it can be an ok solution
            return True

        return False

    @staticmethod
    def _get_nesting_level(parent_node: Optional[RecursiveTreeNode]) -> int:
        is_root = parent_node is None
        if is_root:
            nesting_level = 0
        else:
            need_to_go_down_one_level = isinstance(parent_node, MultiClusterNode)
            if need_to_go_down_one_level:
                nesting_level = parent_node.nesting_level + 1
            else:
                nesting_level = parent_node.nesting_level
        return nesting_level
    #####################################################################################################

    #####################################################################################################
    #  interval methods
    @staticmethod
    def _infer_if_has_single_interval(all_intervals: Intervals, match_intervals: Intervals,
                                      nesting_level: int, max_nesting: int,
                                      alignment: MSA, min_match_length: int) -> bool:
        single_match_interval = (len(all_intervals) == 1) and (all_intervals[0] in match_intervals)
        if single_match_interval:
            return True

        max_nesting_level_reached = nesting_level == max_nesting
        if max_nesting_level_reached:
            return True

        small_variant_site = alignment.get_alignment_length() < min_match_length
        if small_variant_site:
            return True

        return False

    @staticmethod
    def _get_all_and_match_intervals(alignment: MSA, min_match_length: int) -> Tuple[Intervals, Intervals]:
        consensus = get_consensus_from_MSA(alignment)
        interval_partitioner = IntervalPartitioner(consensus, min_match_length, alignment)
        (
            match_intervals,
            non_match_intervals,
            all_intervals,
        ) = interval_partitioner.get_intervals()
        return all_intervals, match_intervals

    @staticmethod
    def _partition_alignment_into_interval_subalignments(alignment: MSA, all_intervals: Intervals) -> List[MSA]:
        return [alignment[:, interval.start: interval.stop + 1] for interval in all_intervals]
    #####################################################################################################

    #####################################################################################################
    #  clustering methods
    @staticmethod
    def _infer_if_we_should_cluster_further(alignment: MSA, clustering_result: ClusteringResult) -> bool:
        if clustering_result.no_clustering:
            return False

        alignment_has_issues = NodeFactory._alignment_has_issues(alignment)
        if alignment_has_issues:
            return False

        return True

    @staticmethod
    def _get_subalignments_by_clustering(alignment: MSA, clustering_result: ClusteringResult) -> SubMSAs:
        list_sub_alignments = [
            NodeFactory._get_sub_alignment_by_list_id(alignment, clustered_id) for clustered_id in clustering_result.clustered_ids
        ]
        return list_sub_alignments

    @staticmethod
    def _get_sub_alignment_by_list_id(alignment: MSA, id_list: List[str]) -> MSA:
        list_records = [record for record in alignment if record.id in id_list]
        sub_alignment = MSA(list_records)
        return sub_alignment
    #####################################################################################################
