from unittest import TestCase

from make_prg.utils.misc import (
    flatten_list,
    remove_duplicated_consecutive_elems_from_list,
)


class TestMisc(TestCase):
    def test___remove_duplicated_consecutive_elems_from_list___empty_list(self):
        the_list = []

        expected = []
        actual = remove_duplicated_consecutive_elems_from_list(the_list)

        self.assertEqual(expected, actual)

    def test___remove_duplicated_consecutive_elems_from_list___no_duplicates___nothing_removed(
        self,
    ):
        the_list = [1, 2, 3, 4]

        expected = [1, 2, 3, 4]
        actual = remove_duplicated_consecutive_elems_from_list(the_list)

        self.assertEqual(expected, actual)

    def test___remove_duplicated_consecutive_elems_from_list___several_intertwinned_duplicates___nothing_removed(
        self,
    ):
        the_list = [1, 2, 1, 2, 1, 2]

        expected = [1, 2, 1, 2, 1, 2]
        actual = remove_duplicated_consecutive_elems_from_list(the_list)

        self.assertEqual(expected, actual)

    def test___remove_duplicated_consecutive_elems_from_list___several_consecutive_duplicates(
        self,
    ):
        the_list = [1, 2, 2, 2, 3, 3, 1, 1, 1, 1, 1, 2]

        expected = [1, 2, 3, 1, 2]
        actual = remove_duplicated_consecutive_elems_from_list(the_list)

        self.assertEqual(expected, actual)

    def test___flatten_list___empty_list(self):
        the_list = [[]]

        expected = []
        actual = flatten_list(the_list)

        self.assertEqual(expected, actual)

    def test___flatten_list___single_list(self):
        the_list = [[1, 2, 3]]

        expected = [1, 2, 3]
        actual = flatten_list(the_list)

        self.assertEqual(expected, actual)

    def test___flatten_list___three_lists(self):
        the_list = [[1, 2, 3], [4], [5, 6]]

        expected = [1, 2, 3, 4, 5, 6]
        actual = flatten_list(the_list)

        self.assertEqual(expected, actual)
