from unittest import TestCase
from unittest.mock import Mock

from intervaltree import Interval, IntervalTree

from make_prg.update.MLPath import MLPath, MLPathError, MLPathNode


class MLPathTest(TestCase):
    def setUp(self):
        self.ml_path_node_1 = MLPathNode(key=(5, 9), sequence="ACGT")
        self.ml_path_node_2 = MLPathNode(key=(12, 13), sequence="G")
        self.ml_path_node_3 = MLPathNode(key=(20, 22), sequence="AA")
        self.ml_path_nodes = [
            self.ml_path_node_1,
            self.ml_path_node_2,
            self.ml_path_node_3,
        ]
        self.ml_path = MLPath(self.ml_path_nodes)

    def test___constructor___empty_MLPath___raises_MLPathError(self):
        with self.assertRaises(MLPathError):
            MLPath([])

    # TODO: might not be the best practice to check the private attributes in tests
    # TODO: it breaks encapsulation
    def test___constructor___single_MLPath(self):
        ml_path_nodes = [self.ml_path_node_1]
        ml_path = MLPath(ml_path_nodes)

        self.assertEqual(ml_path_nodes, ml_path._ml_path_nodes)
        expected_linear_path_space_index = IntervalTree(
            [Interval(0, 4, self.ml_path_node_1)]
        )
        self.assertEqual(
            expected_linear_path_space_index,
            ml_path._ml_path_index_in_linear_path_space,
        )
        expected_PRG_space_index = {(5, 9): self.ml_path_node_1}
        self.assertEqual(expected_PRG_space_index, ml_path._ml_path_index_in_PRG_space)

    def test___constructor___MLPath_with_two_nodes(self):
        ml_path_nodes = [self.ml_path_node_1, self.ml_path_node_2]
        ml_path = MLPath(ml_path_nodes)

        self.assertEqual(ml_path_nodes, ml_path._ml_path_nodes)
        expected_linear_path_space_index = IntervalTree(
            [Interval(0, 4, self.ml_path_node_1), Interval(4, 5, self.ml_path_node_2)]
        )
        self.assertEqual(
            expected_linear_path_space_index,
            ml_path._ml_path_index_in_linear_path_space,
        )
        expected_PRG_space_index = {
            (5, 9): self.ml_path_node_1,
            (12, 13): self.ml_path_node_2,
        }
        self.assertEqual(expected_PRG_space_index, ml_path._ml_path_index_in_PRG_space)

    def test___constructor___MLPath_with_three_nodes(self):
        self.assertEqual(self.ml_path_nodes, self.ml_path._ml_path_nodes)
        expected_linear_path_space_index = IntervalTree(
            [
                Interval(0, 4, MLPathNode(key=(5, 9), sequence="ACGT")),
                Interval(4, 5, MLPathNode(key=(12, 13), sequence="G")),
                Interval(5, 7, MLPathNode(key=(20, 22), sequence="AA")),
            ]
        )
        self.assertEqual(
            expected_linear_path_space_index,
            self.ml_path._ml_path_index_in_linear_path_space,
        )
        expected_PRG_space_index = {
            (5, 9): self.ml_path_node_1,
            (12, 13): self.ml_path_node_2,
            (20, 22): self.ml_path_node_3,
        }
        self.assertEqual(
            expected_PRG_space_index, self.ml_path._ml_path_index_in_PRG_space
        )

    def test___equality(self):
        ml_path_2 = MLPath(self.ml_path_nodes)
        self.assertEqual(self.ml_path, ml_path_2)

    def test___inequality(self):
        ml_path_2 = MLPath(self.ml_path_nodes)
        ml_path_2._ml_path_nodes = Mock()
        self.assertNotEqual(self.ml_path, ml_path_2)

        ml_path_2 = MLPath(self.ml_path_nodes)
        ml_path_2._ml_path_index_in_linear_path_space = Mock()
        self.assertNotEqual(self.ml_path, ml_path_2)

        ml_path_2 = MLPath(self.ml_path_nodes)
        ml_path_2._ml_path_index_in_PRG_space = Mock()
        self.assertNotEqual(self.ml_path, ml_path_2)

    def test___inequality___other_type(self):
        a_string = "a string"
        self.assertNotEqual(self.ml_path, a_string)

    def test___get_last_insertion_pos(self):
        expected = 7
        actual = self.ml_path.get_last_insertion_pos()
        self.assertEqual(expected, actual)

    def test___get_last_node(self):
        expected = self.ml_path_node_3
        actual = self.ml_path.get_last_node()
        self.assertEqual(expected, actual)

    def test___get_node_given_position_in_linear_path_space___positions_of_ml_path_node_1(
        self,
    ):
        expected = self.ml_path_node_1
        for position in range(0, 4):
            actual = self.ml_path.get_node_given_position_in_linear_path_space(position)
            self.assertEqual(expected, actual)

    def test___get_node_given_position_in_linear_path_space___positions_of_ml_path_node_2(
        self,
    ):
        expected = self.ml_path_node_2
        actual = self.ml_path.get_node_given_position_in_linear_path_space(4)
        self.assertEqual(expected, actual)

    def test___get_node_given_position_in_linear_path_space___positions_of_ml_path_node_3(
        self,
    ):
        expected = self.ml_path_node_3
        for position in range(5, 7):
            actual = self.ml_path.get_node_given_position_in_linear_path_space(position)
            self.assertEqual(expected, actual)

    def test___get_node_given_position_in_linear_path_space___one_position_over_the_last_one___raises_MLPathError(
        self,
    ):
        with self.assertRaises(MLPathError):
            self.ml_path.get_node_given_position_in_linear_path_space(7)

    def test___get_node_given_position_in_linear_path_space___two_positions_over_the_last_one___raises_MLPathError(
        self,
    ):
        with self.assertRaises(MLPathError):
            self.ml_path.get_node_given_position_in_linear_path_space(8)

    def test___get_node_given_interval_in_PRG_space___interval_not_indexed___raises_MLPathError(
        self,
    ):
        with self.assertRaises(MLPathError):
            self.ml_path.get_node_given_interval_in_PRG_space((5, 8))

    def test___get_node_given_interval_in_PRG_space___interval_indexed___first_node(
        self,
    ):
        expected = self.ml_path_node_1
        actual = self.ml_path.get_node_given_interval_in_PRG_space((5, 9))
        self.assertEqual(expected, actual)

    def test___get_node_given_interval_in_PRG_space___interval_indexed___second_node(
        self,
    ):
        expected = self.ml_path_node_2
        actual = self.ml_path.get_node_given_interval_in_PRG_space((12, 13))
        self.assertEqual(expected, actual)

    def test___get_node_given_interval_in_PRG_space___interval_indexed___third_node(
        self,
    ):
        expected = self.ml_path_node_3
        actual = self.ml_path.get_node_given_interval_in_PRG_space((20, 22))
        self.assertEqual(expected, actual)
