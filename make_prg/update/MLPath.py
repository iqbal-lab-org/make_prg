"""
This module contains classes to represent and operate on maximum likelihood paths
described in denovo_paths.txt files
"""

from typing import Dict, List, Optional, Tuple

from intervaltree.intervaltree import IntervalTree


class MLPathError(Exception):
    pass


class EmptyMLPathSequence(Exception):
    pass


class MLPathNode:
    """
    Represents a maximum likelihood path node described in denovo_paths.txt files, e.g.:
    (2 [117, 118) G)
    """

    def __init__(self, key: Tuple[int, int], sequence: str):
        self.key: Tuple[int, int] = key
        self._set_sequence(sequence)

        # these are set by class MLPath during indexing
        self.start_index_in_linear_path: Optional[int] = None
        self.end_index_in_linear_path: Optional[int] = None

        self._check_is_a_valid_node()

    def _set_sequence(self, sequence: str):
        empty_ML_path_sequence = len(sequence) == 0
        if empty_ML_path_sequence:
            raise EmptyMLPathSequence(
                f"Found a ML path node ({self.key}) with empty sequence"
            )
        self.sequence: str = sequence

    def _check_is_a_valid_node(self):
        interval_size = self.key[1] - self.key[0]
        sequence_size = len(self.sequence)
        valid_node = interval_size == sequence_size
        if not valid_node:
            raise MLPathError(f"{self} is not a valid node")

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.key, self.sequence) == (other.key, other.sequence)
        else:
            return False

    def __hash__(self):
        return hash((self.key, self.sequence))

    def __str__(self):
        return (
            f"PRG key = {self.key}; "
            f"ML seq interval = [{self.start_index_in_linear_path}:"
            f"{self.end_index_in_linear_path}]; "
            f"Seq = {self.sequence}"
        )

    def __repr__(self):
        return f'MLPathNode(key={self.key}, sequence="{self.sequence}")'


class MLPath:
    """
    Represents a maximum likelihood path described in denovo_paths.txt files, e.g.:
    9 nodes
    (0 [0, 110) ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCACCGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC)
    (2 [117, 118) G)
    (3 [121, 171) GAAATGCCCACGACCGGGCTGGATGAGCTGACAGAGGCCGAACGCGAGAC)
    (5 [178, 179) G)
    (6 [182, 301) CTGGCCGAATGGCTGGCCAAGCGCCGGGAAGCCTCGCAGAAGTCGCAGGAGGCCTACACGGCCATGTCTGCGGATCGGTGGCTGGTCACGCTGGCCAAGGCCATCAGGGAAGGGCAGGA)
    (8 [312, 316) ACTG)
    (9 [319, 360) CGCCCCGAACAGGCGGCCGCGATCTGGCACGGCATGGGGGA)
    (11 [369, 370) G)
    (12 [374, 491) GTCGGCAAGGCCTTGCGCAAGGCTGGTCACGCGAAGCCCAAGGCGGTCAGAAAGGGCAAGCCGGTCGATCCGGCTGATCCCAAGGATCAAGGGGAGGGGGCACCAAAGGGGAAATGA)
    """  # noqa

    def __init__(self, ml_path_nodes: List[MLPathNode]):
        if len(ml_path_nodes) == 0:
            raise MLPathError("ML paths cannot be empty")
        self._ml_path_nodes: List[MLPathNode] = ml_path_nodes
        self._ml_path_index_in_linear_path_space: IntervalTree = IntervalTree()
        self._ml_path_index_in_PRG_space: Dict[
            Tuple[int, int], MLPathNode
        ] = {}  # dict because this is exact index
        self._index()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def _index(self):
        start_index_in_linear_path = 0
        for ml_path_node_index, ml_path_node in enumerate(self._ml_path_nodes):
            end_index_in_linear_path = start_index_in_linear_path + len(
                ml_path_node.sequence
            )
            self._ml_path_index_in_linear_path_space.addi(
                start_index_in_linear_path,
                end_index_in_linear_path,
                data=ml_path_node,
            )
            ml_path_node.start_index_in_linear_path = start_index_in_linear_path
            ml_path_node.end_index_in_linear_path = end_index_in_linear_path
            start_index_in_linear_path = end_index_in_linear_path

            self._ml_path_index_in_PRG_space[ml_path_node.key] = ml_path_node

    def get_last_insertion_pos(self) -> int:
        """
        This position is not indexed, but can be where we insert a sequence after the
        last base of the last node
        """
        return self._ml_path_nodes[-1].end_index_in_linear_path

    def get_last_node(self) -> MLPathNode:
        return self._ml_path_nodes[-1]

    def get_node_given_position_in_linear_path_space(self, position: int) -> MLPathNode:
        nodes = self._ml_path_index_in_linear_path_space[position]

        two_or_more_nodes_overlap_this_position = len(nodes) >= 2
        assert not two_or_more_nodes_overlap_this_position, (
            f"2+ nodes overlap at position {position}, "
            f"this ML path ({self}) probably has an indexing bug."
        )
        no_nodes_overlap_this_position = len(nodes) == 0
        if no_nodes_overlap_this_position:
            raise MLPathError(
                f"No nodes overlap this ML path ({self}) at position {position}, "
                f"is the denovo_paths.txt given as input correct?"
            )

        # only one node overlap this position
        node = list(nodes)[0].data
        return node

    def get_node_given_interval_in_PRG_space(
        self, interval: Tuple[int, int]
    ) -> MLPathNode:
        leaf_interval_not_indexed = interval not in self._ml_path_index_in_PRG_space
        if leaf_interval_not_indexed:
            raise MLPathError(
                f"PRG space interval ({interval}) not indexed in this node ({self})"
            )
        return self._ml_path_index_in_PRG_space[interval]

    def __repr__(self):
        return f"MLPath({self._ml_path_nodes})"
