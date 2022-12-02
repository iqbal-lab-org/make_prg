"""
This module contains classes to build PRGs from MSAs and manage sets of PRGs
"""

import pickle
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from zipfile import ZipFile

from make_prg.recursion_tree import LeafNode, NodeFactory, RecursiveTreeNode
from make_prg.utils.io_utils import load_alignment_file, zip_set_of_files
from make_prg.utils.prg_encoder import PRG_Ints, PrgEncoder


class LeafNotFoundException(Exception):
    pass


class PrgBuilder(object):
    """
    Prg builder based from a multiple sequence alignment.
    """

    def __init__(
        self,
        locus_name: str,
        msa_file: Path,
        alignment_format: str,
        max_nesting: int,
        min_match_length: int,
        aligner: Optional["MSAAligner"] = None,  # noqa: F821
    ):
        self._locus_name: str = locus_name  # note locus_name is fully protected from writes (see self.__hash__())
        self.max_nesting: int = max_nesting
        self.min_match_length: int = min_match_length
        self.aligner: Optional["MSAAligner"] = aligner  # noqa: F821
        self.next_node_id: int = 0
        self.site_num: int = 5
        self.prg_index: Dict[Tuple[int, int], LeafNode] = {}

        alignment = load_alignment_file(str(msa_file), alignment_format)
        self.root: RecursiveTreeNode = NodeFactory.build(alignment, self, None)

    @property
    def locus_name(self):
        return self._locus_name

    def __getstate__(self):
        """
        This method is defined so that we don't pickle self.aligner.
        If we pickle self.aligner, then two PRGs built even with the same aligner, but different temp paths won't be
        comparable, and we don't actually care about the aligner/temp path that was used, if the rest is the same.
        """
        state = self.__dict__.copy()
        state["aligner"] = None  # remove aligner from the object to be pickled
        return state

    def __eq__(self, other: "PrgBuilder") -> bool:
        """
        Equality comparison is required because when identical objects are pickled in different environments/machines,
        their byte representation might be different. However, the reconstructed object is still identical.
        Similarly to self.__getstate__(), we don't take self.aligner into account here
        """
        # first compare trivial attributes
        if (
            self.locus_name,
            self.max_nesting,
            self.min_match_length,
            self.next_node_id,
            self.site_num,
        ) != (
            other.locus_name,
            other.max_nesting,
            other.min_match_length,
            other.next_node_id,
            other.site_num,
        ):
            return False

        # now compares attributes that involve RecursiveTreeNode
        # these are costly to compare, so better to do later
        if self.prg_index != other.prg_index or self.root != other.root:
            return False

        return True

    def __hash__(self):
        """
        Properties to satisfy:
            If a == b then hash(a) == hash(b)
            If hash(a) == hash(b), then a might equal b
            If hash(a) != hash(b), then a != b
        Let's then hash on the locus name, which is protected from writes
        """
        return hash(self.locus_name)

    def replace_root(self, new_root: RecursiveTreeNode):
        self.root = new_root

    def build_prg(self) -> str:
        self.site_num = 5
        prg_as_list = []
        self.root.preorder_traversal_to_build_prg(prg_as_list)
        prg = "".join(prg_as_list)
        return prg

    def get_next_site_num(self) -> int:
        site_num = self.site_num
        self.site_num += 2
        return site_num

    def get_next_node_id(self) -> int:
        self.next_node_id += 1
        return self.next_node_id - 1

    def update_PRG_index(self, start_index: int, end_index: int, node: LeafNode):
        interval = (start_index, end_index)
        self.prg_index[interval] = node
        node.add_indexed_PRG_interval(interval)

    def clear_PRG_index(self):
        # first clear the index of all nodes
        for node in self.prg_index.values():
            node.clear_PRG_interval_index()

        self.prg_index.clear()

    def get_node_given_interval(self, interval: Tuple[int, int]) -> LeafNode:
        interval_is_indexed = interval in self.prg_index

        # TODO: move this back to assert once is solved
        # TODO: should it really be an assert?
        # TODO: in the pandora paper data, out of 486k update operations, 12 failed with this error
        # TODO: so, there is an edge-case bug here to be solved in next minor versions
        if not interval_is_indexed:
            raise LeafNotFoundException(
                f"Queried PRG interval {interval} does not exist in PRG index for locus {self.locus_name}.\n"
                f"Indexed PRG intervals: {self.prg_index.keys()}"
            )
        # assert interval in self.prg_index, \
        #     f"Fatal error: Queried interval {interval} does not exist in leaves index for locus {self.locus_name}"

        return self.prg_index[interval]

    def serialize(self, filepath: [Path, str]):
        with open(filepath, "wb") as filehandler:
            pickle.dump(self, filehandler, protocol=4)

    @staticmethod
    def deserialize_from_bytes(array_of_bytes: bytes) -> "PrgBuilder":
        return pickle.loads(array_of_bytes)

    @staticmethod
    def write_prg_as_text(output_prefix: str, prg_string: str):
        sample = Path(output_prefix).name
        prg_filename = Path(output_prefix + ".prg.fa")
        with prg_filename.open("w") as prg:
            print(f">{sample}\n{prg_string}", file=prg)

    @staticmethod
    def write_prg_as_binary(output_prefix: str, prg_string: str):
        prg_ints_fpath = Path(output_prefix + ".bin")
        prg_encoder = PrgEncoder()
        prg_ints: PRG_Ints = prg_encoder.encode(prg_string)
        with prg_ints_fpath.open("wb") as ostream:
            prg_encoder.write(prg_ints, ostream)


class PrgBuilderZipDatabase:
    """
    Represent a collection of PrgBuilders, to be saved to and loaded from a zip file
    """

    def __init__(self, zip_filepath: Path):
        is_a_zip_file = zip_filepath.suffix == ".zip"
        assert (
            is_a_zip_file
        ), "PrgBuilderZipDatabase initialised without a .zip filepath"
        self._zip_filepath: Path = zip_filepath
        self._zip_file: Optional[ZipFile] = None

    def save(self, locus_to_prg_builder_pickle_path: Dict[str, Path]):
        zip_set_of_files(self._zip_filepath, locus_to_prg_builder_pickle_path)

    def load(self):
        self._zip_file = ZipFile(self._zip_filepath)

    def close(self):
        if self._zip_file is not None:
            self._zip_file.close()

    def get_number_of_loci(self) -> int:
        return len(self.get_loci_names())

    def get_loci_names(self) -> List[str]:
        return sorted(self._zip_file.namelist())

    def get_PrgBuilder(self, locus: str) -> PrgBuilder:
        return PrgBuilder.deserialize_from_bytes(self._zip_file.read(locus))

    def __eq__(self, other: "PrgBuilderZipDatabase") -> bool:
        different_loci = self.get_loci_names() != other.get_loci_names()
        if different_loci:
            return False

        for locus in self.get_loci_names():
            prg_builder_1 = self.get_PrgBuilder(locus)
            prg_builder_2 = other.get_PrgBuilder(locus)
            if prg_builder_1 != prg_builder_2:
                return False

        return True
