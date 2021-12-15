from typing import List, Any
from itertools import groupby, chain


def remove_duplicated_consecutive_elems_from_list(the_list: List[Any]) -> List[Any]:
    return [elem[0] for elem in groupby(the_list)]


def flatten_list(list_of_lists: List[List[Any]]) -> List[Any]:
    return list(chain.from_iterable(list_of_lists))
