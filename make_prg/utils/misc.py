import os
from itertools import chain, groupby
from typing import Any, List

from make_prg import MSA


def remove_duplicated_consecutive_elems_from_list(the_list: List[Any]) -> List[Any]:
    return [elem[0] for elem in groupby(the_list)]


def flatten_list(list_of_lists: List[List[Any]]) -> List[Any]:
    return list(chain.from_iterable(list_of_lists))


def equal_msas(msa_1: MSA, msa_2: MSA) -> bool:
    """
    This is required because Bio.AlignIO.MultipleSeqAlignment has no __eq__() defined
    """
    msa_1_as_fasta = format(msa_1, "fasta")
    msa_2_as_fasta = format(msa_2, "fasta")
    return msa_1_as_fasta == msa_2_as_fasta


# Note: not unit tested, just used for in-house debugging
def should_output_debug_graphs() -> bool:
    return "make_prg_output_debug_graphs" in os.environ
