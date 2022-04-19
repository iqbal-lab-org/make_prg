from unittest import TestCase
from make_prg.update.MLPath import MLPathNode, EmptyMLPathSequence, MLPathError

class MLPathNodeTest(TestCase):
    def test___constructor___invalid_sequence_raises_EmptyMLPathSequence(self):
        with self.assertRaises(EmptyMLPathSequence):
            MLPathNode((0, 0), sequence="")

    def test___constructor___invalid_node_length_one_base_shorter_raises_MLPathError(self):
        with self.assertRaises(MLPathError):
            MLPathNode((5, 8), sequence="ACGT")

    def test___constructor___invalid_node_length_one_base_longer_raises_MLPathError(self):
        with self.assertRaises(MLPathError):
            MLPathNode((5, 10), sequence="ACGT")

    def test___constructor(self):
        ml_path_node = MLPathNode((5, 9), sequence="ACGT")
        self.assertEqual((5, 9), ml_path_node.key)
        self.assertEqual("ACGT", ml_path_node.sequence)
        self.assertEqual(None, ml_path_node.start_index_in_linear_path)
        self.assertEqual(None, ml_path_node.end_index_in_linear_path)

    def test___equality(self):
        ml_path_node_1 = MLPathNode((5, 9), sequence="ACGT")
        ml_path_node_2 = MLPathNode((5, 9), sequence="ACGT")
        self.assertEqual(ml_path_node_1, ml_path_node_2)

    def test___inequality(self):
        ml_path_node_1 = MLPathNode((5, 9), sequence="ACGT")

        ml_path_node_2 = MLPathNode((6, 10), sequence="ACGT")
        self.assertNotEqual(ml_path_node_1, ml_path_node_2)

        ml_path_node_2 = MLPathNode((5, 9), sequence="ACCT")
        self.assertNotEqual(ml_path_node_1, ml_path_node_2)

    def test___inequality___other_type(self):
        ml_path_node_1 = MLPathNode((5, 9), sequence="ACGT")
        a_string = "a string"
        self.assertNotEqual(ml_path_node_1, a_string)

    def test___str(self):
        ml_path_node = MLPathNode((5, 9), sequence="ACGT")

        expected = "PRG key = (5, 9); ML seq interval = [None:None]; Seq = ACGT"
        actual = str(ml_path_node)

        self.assertEqual(expected, actual)

    def test___repr(self):
        ml_path_node = MLPathNode((5, 9), sequence="ACGT")

        expected = "MLPathNode(key=(5, 9), sequence=\"ACGT\")"
        actual = repr(ml_path_node)

        self.assertEqual(expected, actual)
