from unittest import TestCase
from unittest.mock import Mock
from make_prg.update.denovo_variants import UpdateData

class TestUpdateData(TestCase):
    def test___constructor(self):
        ml_path_node_key = Mock()
        ml_path = Mock()
        new_node_sequence = Mock()

        update_data = UpdateData(ml_path_node_key, ml_path, new_node_sequence)

        self.assertEqual(ml_path_node_key, update_data.ml_path_node_key)
        self.assertEqual(ml_path, update_data.ml_path)
        self.assertEqual(new_node_sequence, update_data.new_node_sequence)

    def test___equality(self):
        ml_path_node_key = Mock()
        ml_path = Mock()
        new_node_sequence = Mock()

        update_data_1 = UpdateData(ml_path_node_key, ml_path, new_node_sequence)
        update_data_2 = UpdateData(ml_path_node_key, ml_path, new_node_sequence)

        self.assertEqual(update_data_1, update_data_2)

    def test___inequality(self):
        ml_path_node_key = Mock()
        ml_path = Mock()
        new_node_sequence = Mock()

        update_data_1 = UpdateData(ml_path_node_key, ml_path, new_node_sequence)

        update_data_2 = UpdateData(Mock(), ml_path, new_node_sequence)
        self.assertNotEqual(update_data_1, update_data_2)

        update_data_2 = UpdateData(ml_path_node_key, Mock(), new_node_sequence)
        self.assertNotEqual(update_data_1, update_data_2)

        update_data_2 = UpdateData(ml_path_node_key, ml_path, Mock())
        self.assertNotEqual(update_data_1, update_data_2)

    def test___inequality___with_other_type(self):
        ml_path_node_key = Mock()
        ml_path = Mock()
        new_node_sequence = Mock()
        update_data_1 = UpdateData(ml_path_node_key, ml_path, new_node_sequence)
        a_string = "a string"

        self.assertNotEqual(update_data_1, a_string)

    def test___repr(self):
        ml_path_node_key = "ml_path_node_key"
        ml_path = "ml_path"
        new_node_sequence = "new_node_sequence"
        update_data_1 = UpdateData(ml_path_node_key, ml_path, new_node_sequence)

        expected = "UpdateData(ml_path_node_key='ml_path_node_key', ml_path='ml_path', new_node_sequence='new_node_sequence')"
        actual = repr(update_data_1)
        self.assertEqual(expected, actual)
