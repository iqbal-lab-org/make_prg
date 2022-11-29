"""
This module contains classes used to parse and represent denovo variants info
contained in a denovo_paths.txt file
"""

import re
import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Deque, Dict, List, Optional, TextIO, Tuple

from loguru import logger

from make_prg.update.MLPath import EmptyMLPathSequence, MLPath, MLPathError, MLPathNode
from make_prg.utils.misc import remove_duplicated_consecutive_elems_from_list
from make_prg.utils.seq_utils import GAP, align


class DenovoError(Exception):
    pass


class TooLongDeletion(Exception):
    pass


class DenovoVariant:
    """
    Represents a denovo variant in a denovo_paths.txt file, e.g.: "44	C	T"
    """

    def __init__(
        self,
        start_index_in_linear_path: int,
        ref: str,
        alt: str,
        ml_path_nodes_it_goes_through: Optional[List[MLPathNode]] = None,
        long_deletion_threshold: int = sys.maxsize,
    ):
        DenovoVariant._param_checking(
            start_index_in_linear_path, ref, alt, long_deletion_threshold
        )
        self.start_index_in_linear_path: int = start_index_in_linear_path
        self.end_index_in_linear_path: int = start_index_in_linear_path + len(ref)
        self.ref: str = ref
        self.alt: str = alt
        self.set_ml_path_nodes_it_goes_through(ml_path_nodes_it_goes_through)
        self.long_deletion_threshold = long_deletion_threshold

    @staticmethod
    def _param_checking(
        start_index_in_linear_path: int,
        ref: str,
        alt: str,
        long_deletion_threshold: int,
    ):
        DenovoVariant._check_sequence_is_composed_of_ACGT_only(ref)
        DenovoVariant._check_sequence_is_composed_of_ACGT_only(alt)
        not_a_variant = ref == alt
        if not_a_variant:
            raise DenovoError(
                f"Found a variant where ref ({ref}) equals alt ({alt}), this is not a "
                f"variant"
            )

        negative_index_for_variant_pos = start_index_in_linear_path < 0
        if negative_index_for_variant_pos:
            raise DenovoError(
                f"Found a negative index for variant pos ({start_index_in_linear_path})"
            )

        deletion_size = len(ref) - len(alt)
        is_a_too_long_deletion = deletion_size >= long_deletion_threshold
        if is_a_too_long_deletion:
            raise TooLongDeletion(
                f"Variant has a too long deletion (delta = {deletion_size}) that "
                f"should be ignored"
            )

    @staticmethod
    def _check_sequence_is_composed_of_ACGT_only(seq: str):
        sequence_is_composed_of_ACGT_only = all([base in "ACGT" for base in seq])
        if not sequence_is_composed_of_ACGT_only:
            raise DenovoError(f"Found a non-ACGT seq ({seq}) in a denovo variant")

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def set_ml_path_nodes_it_goes_through(
        self, ml_path_nodes_it_goes_through: Optional[List[MLPathNode]]
    ):
        """
        self.ml_path_nodes_it_goes_through: a list of MLPathNode, where the i-th
        MLPathNode is the node the i-th base
        of the variant goes through
        """
        if ml_path_nodes_it_goes_through is not None:
            if self.is_strict_insertion_event():
                ml_path_contains_only_the_insertion_point = (
                    len(ml_path_nodes_it_goes_through) == 1
                )
                valid_parameters = ml_path_contains_only_the_insertion_point
            else:
                each_base_is_covered_by_one_node = len(self.ref) == len(
                    ml_path_nodes_it_goes_through
                )
                valid_parameters = each_base_is_covered_by_one_node
            assert valid_parameters, (
                f"Invalid parameters for "
                f"DenovoVariant.set_ml_path_nodes_it_goes_through().\n"
                f"Debug info:\n"
                f"DenovoVariant: {self}\n"
                f"ml_path_nodes_it_goes_through: {ml_path_nodes_it_goes_through}"
            )

        self.ml_path_nodes_it_goes_through: Optional[
            List[MLPathNode]
        ] = ml_path_nodes_it_goes_through

    def get_mutated_sequence(self) -> str:
        """
        We apply this variant to the single MLPathNode in self.ml_path_it_goes_through.
        self.ml_path_it_goes_through needs to have a single node compatible with this
        variant
        """
        ml_path_nodes_it_goes_through_has_a_single_distinct_node = (
            self.ml_path_nodes_it_goes_through is not None
            and len(set(self.ml_path_nodes_it_goes_through)) == 1
        )
        assert ml_path_nodes_it_goes_through_has_a_single_distinct_node, (
            f"Cannot apply variant {self} as it does not go through a single distinct "
            f"node\n"
            f"ML path nodes the variant goes through: "
            f"{self.ml_path_nodes_it_goes_through}"
        )

        node = self.ml_path_nodes_it_goes_through[0]
        node_is_compatible_with_this_variant = (
            node.start_index_in_linear_path <= self.start_index_in_linear_path
            and self.end_index_in_linear_path <= node.end_index_in_linear_path
        )
        assert (
            node_is_compatible_with_this_variant
        ), f"Node {node} is not compatible with variant {self}"

        start_index_inside_node_sequence = (
            self.start_index_in_linear_path - node.start_index_in_linear_path
        )
        end_index_inside_node_sequence = start_index_inside_node_sequence + len(
            self.ref
        )
        ref_wrt_indexes = node.sequence[
            start_index_inside_node_sequence:end_index_inside_node_sequence
        ]
        ref_is_consistent = self.ref == ref_wrt_indexes
        assert ref_is_consistent, (
            f"Ref is not consistent for {self}. Node = {node}. ref_wrt_indexes = "
            f"{ref_wrt_indexes}"
        )

        mutated_sequence = (
            node.sequence[:start_index_inside_node_sequence]
            + self.alt
            + node.sequence[end_index_inside_node_sequence:]
        )
        return mutated_sequence

    def _split_variant_at_boundary_alignment(
        self, ref_alignment: Deque[str], alt_alignment: Deque[str]
    ) -> List["DenovoVariant"]:
        split_variants = []
        current_index_in_linear_path = self.start_index_in_linear_path
        ml_path_node_to_count = Counter(self.ml_path_nodes_it_goes_through)
        deduplicated_ml_path_nodes_it_goes_through = (
            remove_duplicated_consecutive_elems_from_list(
                self.ml_path_nodes_it_goes_through
            )
        )
        for ml_path_node_index, ml_path_node in enumerate(
            deduplicated_ml_path_nodes_it_goes_through
        ):
            sub_ref = []
            sub_alt = []
            nb_of_bases_to_consume = ml_path_node_to_count[ml_path_node]
            current_start_in_linear_path = current_index_in_linear_path

            while nb_of_bases_to_consume > 0:
                ref_base = ref_alignment.popleft()
                if ref_base != GAP:
                    sub_ref.append(ref_base)
                    current_index_in_linear_path += 1
                    nb_of_bases_to_consume -= 1

                alt_base = alt_alignment.popleft()
                if alt_base != GAP:
                    sub_alt.append(alt_base)

            is_last_node = (
                ml_path_node_index
                == len(deduplicated_ml_path_nodes_it_goes_through) - 1
            )
            there_are_remaining_alt_bases = (
                is_last_node and nb_of_bases_to_consume == 0 and len(alt_alignment) > 0
            )
            if there_are_remaining_alt_bases:
                sub_alt.extend(
                    [alt_base for alt_base in alt_alignment if alt_base != GAP]
                )

            sub_ref_seq = "".join(sub_ref)
            sub_alt_seq = "".join(sub_alt)
            sub_ref_and_alt_are_different = sub_ref_seq != sub_alt_seq
            if sub_ref_and_alt_are_different:
                try:
                    split_variant = DenovoVariant(
                        current_start_in_linear_path,
                        sub_ref_seq,
                        sub_alt_seq,
                        long_deletion_threshold=self.long_deletion_threshold,
                    )
                    ml_path_nodes_the_split_variant_goes_through = [ml_path_node] * len(
                        sub_ref_seq
                    )
                    split_variant.set_ml_path_nodes_it_goes_through(
                        ml_path_nodes_the_split_variant_goes_through
                    )
                    split_variants.append(split_variant)
                    logger.debug(f"Split variant to be applied: {split_variant}")
                except TooLongDeletion as error:
                    logger.warning(f"Ignoring split variant: {error}")

        return split_variants

    def split_variant(self) -> List["DenovoVariant"]:
        """
        Split this variant into a list of sub-variants WRT how it is distributed along
        self.ml_path_nodes_it_goes_through
        @return: List of sub-variants that goes through only a single node
        """
        ml_path_nodes_it_goes_through_is_valid = (
            self.ml_path_nodes_it_goes_through is not None
        )
        assert ml_path_nodes_it_goes_through_is_valid, (
            "Error on DenovoVariant.split_variant(): "
            "self.ml_path_nodes_it_goes_through is None"
        )

        nb_of_distinct_ml_path_nodes = len(set(self.ml_path_nodes_it_goes_through))
        variant_goes_through_only_one_leaf = nb_of_distinct_ml_path_nodes == 1
        if variant_goes_through_only_one_leaf:
            split_variants = [self]
            return split_variants

        # here, variant goes through several leaves
        alignment = align(self.ref, self.alt)
        ref_alignment = deque(alignment[0])
        alt_alignment = deque(alignment[1])
        split_variants = self._split_variant_at_boundary_alignment(
            ref_alignment, alt_alignment
        )
        return split_variants

    def is_strict_insertion_event(self) -> bool:
        """
        A strict insertion event is when the ref is empty and the alt has some bases
        """
        return len(self.ref) == 0 and len(self.alt) > 0

    def __str__(self):
        return (
            f"[{self.start_index_in_linear_path}:{self.end_index_in_linear_path}]:"
            f"'{self.ref}'->'{self.alt}'"
        )

    def __repr__(self):
        return (
            "DenovoVariant(start_index_in_linear_path="
            f"{self.start_index_in_linear_path}, "
            f'ref="{self.ref}", '
            f'alt="{self.alt}")'
        )


@dataclass
class UpdateData:
    """
    Represents a trivial minimal class to hold update data when multiprocessing
    """

    ml_path_node_key: Tuple[int, int]
    ml_path: MLPath
    new_node_sequence: str


@dataclass
class DenovoLocusInfo:
    """
    Represents a locus  in a denovo_paths.txt file, e.g.:
    GC00010897
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
    2 denovo variants for this locus
    44	C	T
    422	A	T
    """  # noqa

    sample: str
    locus: str
    ml_path: MLPath
    variants: List[DenovoVariant]

    def _get_ml_path_nodes_spanning_variant(
        self, variant: DenovoVariant
    ) -> List[MLPathNode]:
        """
        Given a variant, return a list of MLPathNodes spanning the variant.
        For each base of the ref of variant, we have one MLPathNode, so if we have a
        ref == "ACGT", the list will
        contain 4 MLPathNodes.
        """
        if variant.is_strict_insertion_event():
            # interval/ref is empty, search for the start index
            try:
                ml_path_node = (
                    self.ml_path.get_node_given_position_in_linear_path_space(
                        variant.start_index_in_linear_path
                    )
                )
            except MLPathError as ml_path_error:
                # this might happen when the insertion is after the last base of the
                # last node
                is_in_the_last_insertion_position = (
                    variant.start_index_in_linear_path
                    == self.ml_path.get_last_insertion_pos()
                )
                if is_in_the_last_insertion_position:
                    ml_path_node = self.ml_path.get_last_node()
                else:
                    # this is indeed an error
                    raise ml_path_error
            ml_path_nodes = [ml_path_node]
        else:
            # we have a real interval / ref is not empty
            ml_path_nodes = []
            for position_in_linear_path in range(
                variant.start_index_in_linear_path, variant.end_index_in_linear_path
            ):
                ml_path_node = (
                    self.ml_path.get_node_given_position_in_linear_path_space(
                        position_in_linear_path
                    )
                )
                ml_path_nodes.append(ml_path_node)

        return ml_path_nodes

    def get_update_data(self) -> List[UpdateData]:
        """
        Get a list of updates to be done to MLPathNodes to update the PRG
        """
        update_data_list = []
        for variant in self.variants:
            ml_path_nodes = self._get_ml_path_nodes_spanning_variant(variant)
            variant.set_ml_path_nodes_it_goes_through(ml_path_nodes)
            split_variants = variant.split_variant()

            for split_variant in split_variants:
                update_data = UpdateData(
                    ml_path_node_key=split_variant.ml_path_nodes_it_goes_through[0].key,
                    ml_path=self.ml_path,
                    new_node_sequence=split_variant.get_mutated_sequence(),
                )
                update_data_list.append(update_data)

        return update_data_list


class DenovoVariantsDB:
    """
    Represents the whole denovo_paths.txt file. Most important member is the attribute
    self.locus_name_to_update_data,
    where users can get the update data given a locus.
    """

    @staticmethod
    def _read_nb_of_samples(filehandler: TextIO) -> int:
        line = filehandler.readline().strip()
        nb_of_samples = int(line.split()[0])
        return nb_of_samples

    @staticmethod
    def _read_sample(filehandler: TextIO) -> str:
        line = filehandler.readline().strip()
        sample = line.split()[1]
        logger.trace(f"Read sample: {sample}")
        return sample

    @staticmethod
    def _read_nb_of_loci_in_sample(filehandler: TextIO) -> int:
        line = filehandler.readline().strip()
        nb_of_loci_in_sample = int(line.split()[0])
        return nb_of_loci_in_sample

    @staticmethod
    def _read_locus(filehandler: TextIO) -> str:
        locus = filehandler.readline().strip()
        logger.trace(f"Read locus: {locus}")
        return locus

    @staticmethod
    def _read_nb_of_nodes_in_ml_path(filehandler: TextIO) -> int:
        line = filehandler.readline().strip()
        nb_of_nodes_in_ml_path = int(line.split()[0])
        return nb_of_nodes_in_ml_path

    @classmethod
    def _read_MLPathNode(cls, filehandler: TextIO) -> MLPathNode:
        line = filehandler.readline().strip()
        try:
            matches = cls.ml_path_regex.search(line)
            start_index = int(matches.group(1))
            end_index = int(matches.group(2))
            sequence = matches.group(3)
        except Exception as exc:
            assert False, (
                f"Failed matching ML path regex to line: {line}\n"
                f"Exception: {str(exc)}"
            )

        ml_path_node = MLPathNode(key=(start_index, end_index), sequence=sequence)
        return ml_path_node

    @classmethod
    def _read_ml_path(cls, filehandler: TextIO) -> MLPath:
        nb_of_nodes_in_ml_path = cls._read_nb_of_nodes_in_ml_path(filehandler)
        ml_path = []
        for _ in range(nb_of_nodes_in_ml_path):
            try:
                ml_path_node = cls._read_MLPathNode(filehandler)
                ml_path.append(ml_path_node)
            except EmptyMLPathSequence:
                # if the ML path node has empty sequence, it has empty interval, so we
                # just ignore it, as it is not amenable to updates
                pass

        logger.trace(f"Read ML path: {ml_path}")
        return MLPath(ml_path)

    @staticmethod
    def _read_nb_of_variants(filehandler: TextIO) -> int:
        line = filehandler.readline().strip()
        nb_of_variants = int(line.split()[0])
        return nb_of_variants

    @staticmethod
    def _read_DenovoVariant(
        filehandler: TextIO, long_deletion_threshold: int = sys.maxsize
    ) -> DenovoVariant:
        line = filehandler.readline().strip("\n")
        line_split = line.split("\t")

        start_index_in_linear_path = int(line_split[0]) - 1
        ref = line_split[1]
        alt = line_split[2]
        denovo_variant = DenovoVariant(
            start_index_in_linear_path=start_index_in_linear_path,
            ref=ref,
            alt=alt,
            long_deletion_threshold=long_deletion_threshold,
        )

        logger.debug(f"Read variant: {denovo_variant}")
        return denovo_variant

    @classmethod
    def _read_variants(
        cls, filehandler: TextIO, long_deletion_threshold: int = sys.maxsize
    ) -> List[DenovoVariant]:
        nb_of_variants = cls._read_nb_of_variants(filehandler)
        variants = []
        for _ in range(nb_of_variants):
            try:
                denovo_variant = cls._read_DenovoVariant(
                    filehandler, long_deletion_threshold
                )
                variants.append(denovo_variant)
            except TooLongDeletion as error:
                logger.warning(f"Ignoring variant: {error}")
        return variants

    def _get_locus_name_to_denovo_loci_core(
        self, filehandler: TextIO
    ) -> Dict[str, List[DenovoLocusInfo]]:
        locus_name_to_denovo_loci = defaultdict(list)
        try:
            nb_of_samples = self._read_nb_of_samples(filehandler)
        except IndexError:
            logger.warning(
                f"File containing denovo paths ({self.filepath}) is empty, is it the "
                f"correct file?"
            )
            return locus_name_to_denovo_loci

        for sample_index in range(nb_of_samples):
            sample = self._read_sample(filehandler)
            nb_of_loci_in_sample = self._read_nb_of_loci_in_sample(filehandler)
            for locus_index in range(nb_of_loci_in_sample):
                locus = self._read_locus(filehandler)
                ml_path = self._read_ml_path(filehandler)
                variants = self._read_variants(
                    filehandler, self.long_deletion_threshold
                )
                denovo_locus = DenovoLocusInfo(sample, locus, ml_path, variants)
                locus_name_to_denovo_loci[locus].append(denovo_locus)

        return locus_name_to_denovo_loci

    def _get_locus_name_to_denovo_loci(self) -> Dict[str, List[DenovoLocusInfo]]:
        with open(self.filepath) as filehandler:
            return self._get_locus_name_to_denovo_loci_core(filehandler)

    @staticmethod
    def _get_locus_name_to_update_data(
        locus_name_to_denovo_loci: Dict[str, List[DenovoLocusInfo]]
    ) -> Dict[str, List[UpdateData]]:
        locus_name_to_update_data = defaultdict(list)
        for locus_name, denovo_loci in locus_name_to_denovo_loci.items():
            for denovo_locus in denovo_loci:
                update_data = denovo_locus.get_update_data()
                locus_name_to_update_data[locus_name].extend(update_data)
        return locus_name_to_update_data

    def __init__(self, filepath: str, long_deletion_threshold: int = sys.maxsize):
        self.filepath: Path = Path(filepath)
        self.long_deletion_threshold = long_deletion_threshold
        locus_name_to_denovo_loci = self._get_locus_name_to_denovo_loci()
        self.locus_name_to_update_data = self._get_locus_name_to_update_data(
            locus_name_to_denovo_loci
        )

    # Example:
    # (0 [0, 110) ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCACCGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC)  # noqa
    ml_path_regex = re.compile(r"\(\d+ \[(\d+), (\d+)\) ([ACGT]*)\)")
