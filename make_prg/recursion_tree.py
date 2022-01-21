from typing import List, Set, Optional, Tuple
from loguru import logger
from make_prg.from_msa import MSA
from make_prg.from_msa.cluster_sequences import kmeans_cluster_seqs
from make_prg.utils.seq_utils import (
    SequenceExpander,
    remove_columns_full_of_gaps_from_MSA,
    get_consensus_from_MSA,
    get_number_of_unique_ungapped_sequences,
    get_number_of_unique_gapped_sequences
)
from make_prg.from_msa.interval_partition import IntervalPartitioner, Interval
from make_prg.update.denovo_variants import UpdateData
from make_prg.update.MLPath import MLPathError
from abc import ABC, abstractmethod


class RecursiveTreeNode(ABC):
    def __init__(self, nesting_level: int, alignment: MSA, parent: Optional["RecursiveTreeNode"],
                 prg_builder: "PrgBuilder", force_no_child: bool = False):
        self.nesting_level: int = nesting_level
        self.alignment: MSA = remove_columns_full_of_gaps_from_MSA(alignment)
        self.parent: "RecursiveTreeNode" = parent
        self.prg_builder: "PrgBuilder" = prg_builder
        self.force_no_child = force_no_child
        self.id: int = self.prg_builder.get_next_node_id()

        # generate recursion tree
        self._init_pre_recursion_attributes()
        self._children: List["RecursiveTreeNode"] = self._get_children()

        self.log_that_node_was_created()

    @property
    def children(self):
        return self._children

    @abstractmethod
    def _init_pre_recursion_attributes(self):
        pass

    @abstractmethod
    def _get_children(self) -> List["RecursiveTreeNode"]:
        pass

    @abstractmethod
    def preorder_traversal_to_build_prg(self, prg_as_list: List[str], delim_char: str = " "):
        pass

    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def is_root(self) -> bool:
        return self.parent is None

    def log_that_node_was_created(self):
        logger.trace("Created node:\n" + str(self))

    def __repr__(self):
        return f"Id = {self.id}\n" \
               f"Nesting level = {self.nesting_level}\n" \
               f"Force no child = {self.force_no_child}\n" \
               f"Parent = {'None' if self.parent is None else f'Id = {self.parent.id}'}\n" \
               f"Children = [{', '.join([f'Id = {child.id}' for child in self.children])}]\n" \
               f"Alignment:\n{format(self.alignment, 'fasta')}"

    def __str__(self):
        return repr(self)


class MultiClusterNode(RecursiveTreeNode):
    def __init__(self, nesting_level: int, alignment: MSA, parent: Optional["RecursiveTreeNode"],
                 prg_builder: "PrgBuilder", force_no_child: bool = False):
        super().__init__(nesting_level, alignment, parent, prg_builder, force_no_child)
        assert not self.is_leaf(), "Multicluster nodes should never be leaves"

    def _init_pre_recursion_attributes(self):
        pass  # nothing to init here

    def _get_children(self) -> List["RecursiveTreeNode"]:
        # each child is a PrgBuilderSingleClusterNode for each cluster subalignment
        cluster_subalignments = self._get_subalignments_by_clustering()
        no_clustering_was_done = len(cluster_subalignments) == 1
        children = []
        for alignment in cluster_subalignments:
            child = SingleClusterNode(
                nesting_level=self.nesting_level,
                alignment=alignment,
                parent=self,
                prg_builder=self.prg_builder,
                force_no_child=no_clustering_was_done
            )
            children.append(child)
        return children

    ##################################################################################
    # traversal methods
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
    ##################################################################################

    #####################################################################################################
    #  clustering methods
    def _get_subalignments_by_clustering(self) -> List[MSA]:
        clustering_result = kmeans_cluster_seqs(
            self.alignment,
            self.prg_builder.min_match_length
        )
        list_sub_alignments = [
            self._get_sub_alignment_by_list_id(clustered_id) for clustered_id in clustering_result.clustered_ids
        ]
        return list_sub_alignments

    def _get_sub_alignment_by_list_id(self, id_list: List[str]) -> MSA:
        list_records = [record for record in self.alignment if record.id in id_list]
        sub_alignment = MSA(list_records)
        return sub_alignment
    #####################################################################################################

    def __repr__(self):
        return "MultiClusterNode:\n" + super().__repr__()


class UpdateError(Exception):
    pass


class SingleClusterNode(RecursiveTreeNode):
    def _init_pre_recursion_attributes(self):
        self.consensus: str = get_consensus_from_MSA(self.alignment)
        self.length: int = len(self.consensus)
        interval_partitioner = IntervalPartitioner(self.consensus, self.prg_builder.min_match_length, self.alignment)
        (
            self.match_intervals,
            self.non_match_intervals,
            self.all_intervals,
        ) = interval_partitioner.get_intervals()
        self.new_sequences: Set[str] = set()
        self.indexed_PRG_intervals: Set[Tuple[int, int]] = set()

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

    def _infer_if_this_node_should_have_no_child(self) -> bool:
        if self.force_no_child:
            return True

        single_match_interval = (len(self.all_intervals) == 1) and (
                self.all_intervals[0] in self.match_intervals
        )
        if single_match_interval:
            return True

        max_nesting_level_reached = self.nesting_level == self.prg_builder.max_nesting
        if max_nesting_level_reached:
            return True

        small_variant_site = self.alignment.get_alignment_length() < self.prg_builder.min_match_length
        if small_variant_site:
            return True

        # TODO: this could be computationally optimised (time performance) by receiving this bool variable by
        # TODO: parameter in the constructor, and not recomputing it here. But I think we lose code readability,
        # TODO: maintenability and introduce an annoying dependence... But something to think about if we want to push
        # TODO: performance
        alignment_has_issues = self._alignment_has_issues(self.alignment)
        if alignment_has_issues:
            return True

        return False

    def _infer_if_should_not_cluster(self, interval: Interval, alignment: MSA):
        is_a_match_interval = interval in self.match_intervals
        if is_a_match_interval:
            return True

        alignment_has_issues = self._alignment_has_issues(self.alignment)
        if alignment_has_issues:
            return True

        clustering_result = kmeans_cluster_seqs(alignment, self.prg_builder.min_match_length)
        if clustering_result.no_clustering:
            return True

        return False

    def _get_children(self) -> List["RecursiveTreeNode"]:
        node_has_no_child = self._infer_if_this_node_should_have_no_child()
        if node_has_no_child:
            return list()

        children = []
        for interval in self.all_intervals:
            sub_alignment = self.alignment[:, interval.start : interval.stop + 1]
            sub_alignment_will_not_be_reclustered = self._infer_if_should_not_cluster(interval, sub_alignment)
            if sub_alignment_will_not_be_reclustered:
                subclass = SingleClusterNode
            else:
                subclass = MultiClusterNode
            child = subclass(
                nesting_level=self.nesting_level + 1,
                alignment=sub_alignment,
                parent=self,
                prg_builder=self.prg_builder,
                force_no_child=sub_alignment_will_not_be_reclustered
            )
            children.append(child)

        return children

    def _get_prg(self, prg_as_list: List[str], delim_char: str = " "):
        sequences_can_be_obtained_directly_from_clustering = False
        if not self.is_root():
            clustering_result = kmeans_cluster_seqs(self.alignment, self.prg_builder.min_match_length)
            sequences_can_be_obtained_directly_from_clustering = clustering_result.have_precomputed_sequences

        sequences_of_each_interval = []
        if sequences_can_be_obtained_directly_from_clustering:
            sequences_of_each_interval.append(clustering_result.sequences)
        else:
            for interval in self.all_intervals:
                sub_alignment = self.alignment[:, interval.start:interval.stop + 1]
                seqs = SequenceExpander.get_expanded_sequences_from_MSA(sub_alignment)
                sequences_of_each_interval.append(seqs)

        for sequences_of_this_interval in sequences_of_each_interval:
            single_seq = len(sequences_of_this_interval) == 1
            if single_seq:
                start_index = len(prg_as_list)
                prg_as_list.extend(sequences_of_this_interval[0])
                end_index = len(prg_as_list)
                self.prg_builder.update_PRG_index(start_index, end_index, node=self)
            else:
                # Add the variant seqs to the prg
                site_num = self.prg_builder.get_next_site_num()
                prg_as_list.extend(f"{delim_char}{site_num}{delim_char}")
                for seq_index, seq in enumerate(sequences_of_this_interval):
                    site_num_for_this_seq = (
                        (site_num + 1) if (seq_index < len(sequences_of_this_interval) - 1) else site_num
                    )
                    start_index = len(prg_as_list)
                    prg_as_list.extend(seq)
                    end_index = len(prg_as_list)
                    self.prg_builder.update_PRG_index(
                        start_index, end_index, node=self
                    )
                    prg_as_list.extend(
                        f"{delim_char}{site_num_for_this_seq}{delim_char}"
                    )

    ##################################################################################
    # traversal methods
    def preorder_traversal_to_build_prg(self, prg_as_list: List["str"], delim_char: str = " "):
        if self.is_leaf():
            self._get_prg(prg_as_list, delim_char)
        else:
            for child in self.children:
                child.preorder_traversal_to_build_prg(prg_as_list, delim_char)
    ##################################################################################

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
        logger.trace(f"Updating MSA for {self.prg_builder.locus_name}")
        logger.trace(f"Node: {str(self)}")
        logger.trace(f"Sequences added to update: {self.new_sequences}")

        an_aligner_was_given = self.prg_builder.aligner is not None
        assert an_aligner_was_given, "Cannot make updates without a Multiple Sequence Aligner."

        self.alignment = self.prg_builder.aligner.get_updated_alignment(
            current_alignment=self.alignment,
            new_sequences=self.new_sequences
        )

        # regenerate recursion tree
        self._init_pre_recursion_attributes()
        self._children = self._get_children()
    ##################################################################################

    def __repr__(self):
        return "SingleClusterNode:\n" + RecursiveTreeNode.__repr__(self) + \
               f"Consensus: {self.consensus}\n" \
               f"Match intervals: {self.match_intervals}\n" \
               f"Non-match intervals: {self.non_match_intervals}\n"
