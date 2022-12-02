"""
This module contains classes to explicitly represent the make_prg process recursion tree
"""

from abc import ABC, abstractmethod
from typing import List, Optional, Set, Tuple

from loguru import logger

from make_prg import MSA
from make_prg.from_msa.cluster_sequences import ClusteringResult, kmeans_cluster_seqs
from make_prg.from_msa.interval_partition import IntervalPartitioner, Intervals
from make_prg.update.denovo_variants import UpdateData
from make_prg.update.MLPath import MLPathError
from make_prg.utils.misc import equal_msas
from make_prg.utils.seq_utils import (
    SequenceExpander,
    get_consensus_from_MSA,
    get_number_of_unique_gapped_sequences,
    get_number_of_unique_ungapped_sequences,
    remove_columns_full_of_gaps_from_MSA,
)

SubMSAs = List[MSA]


class RecursiveTreeNode(ABC):
    """
    Abstract base class for the nodes of the recursion tree, abstracting its common
    attributes and methods
    """

    def __init__(
        self,
        nesting_level: int,
        alignment: MSA,
        parent: Optional["RecursiveTreeNode"],
        prg_builder: "PrgBuilder",  # noqa: F821
        children_subalignments: SubMSAs,
    ):
        """
        Builds a new tree node, and its children, recursively
        """
        self.nesting_level: int = nesting_level
        self.alignment: MSA = remove_columns_full_of_gaps_from_MSA(alignment)
        self.parent: Optional["RecursiveTreeNode"] = parent
        self.prg_builder = prg_builder
        self._node_id: int = (
            self.prg_builder.get_next_node_id()
        )  # note node_id is fully protected from writes (see self.__hash__())

        # generate recursion tree
        self._children: List["RecursiveTreeNode"] = self._get_children(
            children_subalignments
        )

        self.log_that_node_was_created()

    @property
    def node_id(self):
        return self._node_id

    def __eq__(self, other: "RecursiveTreeNode") -> bool:
        """
        Compares two RecursiveTreeNodes. Two RecursiveTreeNodes are equal if all their trivial attributes are equal. For
        the non-trivial attributes, we compare them like this:
        1. self._children: we recursively compare these;
        2. self.parent: is enough for us to compare if the parent is set and its id. We should not recursively compare
                        the parent, like the children, because then we would get trapped in an infinite recursive loop.
        3. self.prg_builder: let's compare just the self.prg_builder.locus_name. We should not recursively compare
                             the prg_builders, because then we would get trapped in an infinite recursive loop.
        """
        # first compare trivial attributes
        if (self.nesting_level, self.prg_builder.locus_name, self.node_id) != (
            other.nesting_level,
            other.prg_builder.locus_name,
            other.node_id,
        ):
            return False

        # now compares parent:
        both_parents_are_none = self.parent is None and other.parent is None
        both_parents_are_not_none = self.parent is not None and other.parent is not None
        only_one_parent_is_none = (
            not both_parents_are_none and not both_parents_are_not_none
        )
        if only_one_parent_is_none:
            return False
        different_parents = (
            both_parents_are_not_none and self.parent.node_id != other.parent.node_id
        )
        if different_parents:
            return False

        # now compares the alignment, which requires a special function because Bio.AlignIO.MultipleSeqAlignment
        # does not implement __eq__()
        if not equal_msas(self.alignment, other.alignment):
            return False

        # now recursively compares the children
        different_number_of_children = len(self.children) != len(other.children)
        if different_number_of_children:
            return False
        for first_child, second_child in zip(self.children, other.children):
            if first_child != second_child:
                return False

        return True

    def __hash__(self):
        """
        Properties to satisfy:
            If a == b then hash(a) == hash(b)
            If hash(a) == hash(b), then a might equal b
            If hash(a) != hash(b), then a != b
        Let's then hash on the node_id, which is protected from writes. However, if we hash only on the node_id, we will
        lose performance significantly, because e.g. all first nodes from every PRG will have node_id 0 and will be hashed
        to the same bucket. Thus, we need to mix it with another discriminative property, also protected from writes,
        which could be self.prg_builder.locus_name.
        """
        return hash((self.node_id, self.prg_builder.locus_name))

    @property
    def children(self):
        return self._children

    def _get_children(
        self, children_subalignments: SubMSAs
    ) -> List["RecursiveTreeNode"]:
        return [
            NodeFactory.build(alignment, self.prg_builder, self)
            for alignment in children_subalignments
        ]

    @abstractmethod
    def preorder_traversal_to_build_prg(
        self, prg_as_list: List[str], delim_char: str = " "
    ):
        raise NotImplementedError

    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def is_root(self) -> bool:
        return self.parent is None

    def replace_child(
        self, old_child: "RecursiveTreeNode", new_child: "RecursiveTreeNode"
    ):
        old_child_is_one_of_the_children = old_child in self.children
        assert (
            old_child_is_one_of_the_children
        ), f"Failure to replace a child, {old_child} does not exist"

        old_child_index = self.children.index(old_child)
        self.children[old_child_index] = new_child

    # Note: trivial method, untested
    def log_that_node_was_created(self):
        logger.trace("Created node:\n" + str(self))

    def __repr__(self):
        return (
            f"{self.__class__.__name__}:\n"
            f"Id = {self.node_id}\n"
            f"Nesting level = {self.nesting_level}\n"
            f"Parent = {'None' if self.parent is None else f'Id = {self.parent.node_id}'}\n"
            f"Children = [{', '.join(f'Id = {child.node_id}' for child in self.children)}]\n"
            f"Alignment:\n{format(self.alignment, 'fasta')}"
        )

    def __str__(self):
        return repr(self)


class MultiIntervalNode(RecursiveTreeNode):
    """
    Represents a vertical partition of an MSA
    """

    def __init__(
        self,
        nesting_level: int,
        alignment: MSA,
        parent: Optional["RecursiveTreeNode"],
        prg_builder: "PrgBuilder",  # noqa: F821
        interval_subalignments: SubMSAs,
    ):
        super().__init__(
            nesting_level, alignment, parent, prg_builder, interval_subalignments
        )
        assert not self.is_leaf(), "MultiIntervalNodes should never be leaves"

    def preorder_traversal_to_build_prg(
        self, prg_as_list: List[str], delim_char: str = " "
    ):
        """
        Builds the PRG in prg_as_list. The PRG of a MultiIntervalNode is a concatenation of the PRG of its children
        """
        for child in self.children:
            child.preorder_traversal_to_build_prg(prg_as_list, delim_char)


class MultiClusterNode(RecursiveTreeNode):
    """
    Represents a horizontal partition of an MSA, i.e. sequence clusters
    """

    def __init__(
        self,
        nesting_level: int,
        alignment: MSA,
        parent: Optional["RecursiveTreeNode"],
        prg_builder: "PrgBuilder",  # noqa
        cluster_subalignments: SubMSAs,
    ):
        super().__init__(
            nesting_level, alignment, parent, prg_builder, cluster_subalignments
        )
        assert not self.is_leaf(), "MultiClusterNodes should never be leaves"

    def preorder_traversal_to_build_prg(
        self, prg_as_list: List[str], delim_char: str = " "
    ):
        """
        Builds the PRG in prg_as_list. The PRG of a MultiClusterNode consists of opening a site,
        putting each child as an allele and closing the site
        """
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
    Represents MSAs that are never partitioned.
    These nodes are the only ones that can get updated and indexed.
    """

    def __init__(
        self,
        nesting_level: int,
        alignment: MSA,
        parent: Optional["RecursiveTreeNode"],
        prg_builder: "PrgBuilder",  # noqa: F821
    ):
        super().__init__(nesting_level, alignment, parent, prg_builder, [])
        assert self.is_leaf(), f"Leaf node ({self}) is not a leaf"

        # update-related attributes
        self.new_sequences: Set[str] = set()
        self.indexed_PRG_intervals: Set[Tuple[int, int]] = set()

    def preorder_traversal_to_build_prg(
        self, prg_as_list: List[str], delim_char: str = " ", do_indexing=True
    ):
        """
        Builds the PRG in prg_as_list. The PRG of a leaf node is the sequences themselves it represents
        """
        expanded_sequences = SequenceExpander.get_expanded_sequences_from_MSA(
            self.alignment
        )

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
                    (site_num + 1)
                    if (seq_index < len(expanded_sequences) - 1)
                    else site_num
                )
                start_index = len(prg_as_list)
                prg_as_list.extend(seq)
                end_index = len(prg_as_list)

                prg_as_list.extend(f"{delim_char}{site_num_for_this_seq}{delim_char}")

                if do_indexing:
                    self.prg_builder.update_PRG_index(start_index, end_index, node=self)

    ##################################################################################
    # update methods
    def add_data_to_batch_update(self, update_data: UpdateData):
        """
        Process the given update data and add a new sequence to self.new_sequences
        """
        update_data_PRG_interval = update_data.ml_path_node_key
        update_data_PRG_interval_is_indexed = (
            update_data_PRG_interval in self.indexed_PRG_intervals
        )
        if not update_data_PRG_interval_is_indexed:
            raise UpdateError(
                f"PRG interval {update_data_PRG_interval} not found in indexed "
                f"PRG intervals for node: {self.indexed_PRG_intervals}"
            )

        seq_to_add_as_list = []
        for PRG_interval in sorted(self.indexed_PRG_intervals):
            ML_sequence_should_be_replaced = PRG_interval == update_data_PRG_interval
            if ML_sequence_should_be_replaced:
                seq_to_add_as_list.append(update_data.new_node_sequence)
            else:
                try:
                    ml_path_node = (
                        update_data.ml_path.get_node_given_interval_in_PRG_space(
                            PRG_interval
                        )
                    )
                    seq_to_add_as_list.append(ml_path_node.sequence)
                except MLPathError:
                    pass
        seq_to_add = "".join(seq_to_add_as_list)

        there_has_been_padding = seq_to_add != update_data.new_node_sequence
        if there_has_been_padding:
            logger.trace(
                f"Sequence {update_data.new_node_sequence} padded to {seq_to_add} on "
                f"add_data_to_batch_update() for {self.prg_builder.locus_name}"
            )

        self.new_sequences.add(seq_to_add)

    def add_indexed_PRG_interval(self, interval: Tuple[int, int]):
        self.indexed_PRG_intervals.add(interval)

    def batch_update(self):
        no_update_to_be_done = len(self.new_sequences) == 0
        if no_update_to_be_done:
            return
        self._update_leaf()

    def _update_leaf(self):
        """
        Update this leaf, replacing it by an updated node, which can be of a different subclass.
        Side effects include:
            - Changing the PRG builder root or one of the parent's child;
            - Invalidating the PRG builder index;
        """
        logger.trace(f"Updating subMSA for {self.prg_builder.locus_name}")
        logger.trace(f"Node: {str(self)}")
        logger.trace(f"Sequences added to update: {self.new_sequences}")

        an_aligner_was_given = self.prg_builder.aligner is not None
        assert (
            an_aligner_was_given
        ), "Cannot make updates without a Multiple Sequence Aligner."

        updated_alignment = self.prg_builder.aligner.get_updated_alignment(
            current_alignment=self.alignment, new_sequences=self.new_sequences
        )

        # create a new, updated node, and add it to the tree
        updated_child = NodeFactory.build(
            updated_alignment, self.prg_builder, self.parent
        )

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
    def build(
        alignment: MSA,
        prg_builder: "PrgBuilder",  # noqa: F821
        parent_node: Optional[RecursiveTreeNode] = None,
    ) -> RecursiveTreeNode:
        """
        Builds the correct node given the alignment and other parameters.

        Node building priority depends on the type of node we are building:
        1. If type of node is root
            1.1. Try to build a leaf
            1.2. Force build a multi interval node (even if we have just a single mismatch interval)
        2. If type of node is non-root
            2.1. Try to build a leaf
            2.2. Try to build a multi interval node
            2.3. Try to build a multi cluster node
            2.4. If all fails, force build a leaf

        This can be simplified as:
        1. Try to build a leaf
        2. Try to build a multi interval node (force this build if type of node is root)
        3. Try to build a multi cluster node
        4. Force build a leaf node
        """
        min_match_length = prg_builder.min_match_length
        all_intervals, match_intervals = NodeFactory._get_vertical_partition(
            alignment, min_match_length
        )
        building_a_leaf = NodeFactory._is_single_match_interval(
            all_intervals, match_intervals
        )
        building_multi_interval_node = NodeFactory._is_multi_interval(all_intervals)
        building_the_root = parent_node is None
        nesting_level = 0 if building_the_root else parent_node.nesting_level

        if building_a_leaf:
            return LeafNode(nesting_level, alignment, parent_node, prg_builder)
        elif building_multi_interval_node or building_the_root:
            interval_subalignments = (
                NodeFactory._partition_alignment_into_interval_subalignments(
                    alignment, all_intervals
                )
            )
            return MultiIntervalNode(
                nesting_level,
                alignment,
                parent_node,
                prg_builder,
                interval_subalignments,
            )
        else:  # builds a multi cluster node
            clustering_result = kmeans_cluster_seqs(alignment, min_match_length)
            cluster_further = NodeFactory._infer_if_we_should_cluster_further(
                alignment, clustering_result, nesting_level, prg_builder.max_nesting
            )
            if cluster_further:
                # when building a Multi Cluster node, we open a site, so we go down one nesting level
                nesting_level += 1
                cluster_subalignments = NodeFactory._get_subalignments_by_clustering(
                    alignment, clustering_result
                )
                return MultiClusterNode(
                    nesting_level,
                    alignment,
                    parent_node,
                    prg_builder,
                    cluster_subalignments,
                )
            else:  # can't cluster further, force builds leaf
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

    #####################################################################################################

    #####################################################################################################
    #  interval methods
    @staticmethod
    def _get_vertical_partition(
        alignment: MSA, min_match_length: int
    ) -> Tuple[Intervals, Intervals]:
        consensus = get_consensus_from_MSA(alignment)
        interval_partitioner = IntervalPartitioner(
            consensus, min_match_length, alignment
        )
        (
            match_intervals,
            non_match_intervals,
            all_intervals,
        ) = interval_partitioner.get_intervals()
        return all_intervals, match_intervals

    @staticmethod
    def _is_multi_interval(all_intervals: Intervals) -> bool:
        return len(all_intervals) > 1

    @staticmethod
    def _is_single_match_interval(
        all_intervals: Intervals, match_intervals: Intervals
    ) -> bool:
        return (len(all_intervals) == 1) and (all_intervals[0] in match_intervals)

    @staticmethod
    def _partition_alignment_into_interval_subalignments(
        alignment: MSA, all_intervals: Intervals
    ) -> List[MSA]:
        return [
            alignment[:, interval.start : interval.stop + 1]
            for interval in all_intervals
        ]

    #####################################################################################################

    #####################################################################################################
    #  clustering methods
    @staticmethod
    def _infer_if_we_should_cluster_further(
        alignment: MSA,
        clustering_result: ClusteringResult,
        nesting_level: int,
        max_nesting: int,
    ) -> bool:
        if clustering_result.no_clustering:
            return False

        max_nesting_reached = nesting_level + 1 >= max_nesting
        if max_nesting_reached:
            return False

        alignment_has_issues = NodeFactory._alignment_has_issues(alignment)
        if alignment_has_issues:
            return False

        return True

    @staticmethod
    def _get_subalignments_by_clustering(
        alignment: MSA, clustering_result: ClusteringResult
    ) -> SubMSAs:
        list_sub_alignments = [
            NodeFactory._get_sub_alignment_by_list_id(alignment, clustered_id)
            for clustered_id in clustering_result.clustered_ids
        ]
        return list_sub_alignments

    @staticmethod
    def _get_sub_alignment_by_list_id(alignment: MSA, id_list: List[str]) -> MSA:
        list_records = [record for record in alignment if record.id in id_list]
        sub_alignment = MSA(list_records)
        return sub_alignment

    #####################################################################################################
