from unittest import TestCase
from unittest.mock import patch, Mock, mock_open
from make_prg.prg_builder import PrgBuilder, LeafNotFoundException, PrgBuilderZipDatabase
from make_prg.recursion_tree import NodeFactory
from pathlib import Path
from tests.test_helpers import sample_prg, are_zip_files_equal
import filecmp
from zipfile import ZipFile

workdir = Path("tests/data/prg_builder/write_prg")


class TestPrgBuilder(TestCase):
    # Note: can't use setUp() as patches are not applied to it
    def setup_prg_builder(self) -> None:
        self.locus = "locus"
        self.aligner = Mock()
        self.max_nesting = 5
        self.min_match_length = 7
        self.prg_builder = PrgBuilder(self.locus, Path("msa_file"), "fasta",
                                      self.max_nesting, self.min_match_length, self.aligner)

    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__, return_value="root_node")
    def test___constructor(self, build_mock, load_alignment_file_mock):
        self.setup_prg_builder()
        self.assertEqual(self.locus, self.prg_builder.locus_name)
        self.assertEqual(self.max_nesting, self.prg_builder.max_nesting)
        self.assertEqual(self.min_match_length, self.prg_builder.min_match_length)
        self.assertEqual(self.aligner, self.prg_builder.aligner)
        self.assertEqual(0, self.prg_builder.next_node_id)
        self.assertEqual(5, self.prg_builder.site_num)
        self.assertEqual({}, self.prg_builder.prg_index)
        self.assertEqual("root_node", self.prg_builder.root)
        load_alignment_file_mock.assert_called_once_with("msa_file", "fasta")
        build_mock.assert_called_once_with("MSA", self.prg_builder, None)

    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__, return_value="root_node")
    def test___replace_root(self, *uninteresting_mocks):
        self.setup_prg_builder()
        self.assertEqual("root_node", self.prg_builder.root)

        new_root = Mock()
        self.prg_builder.replace_root(new_root)

        self.assertEqual(new_root, self.prg_builder.root)

    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__)
    def test___build_prg(self, *uninteresting_mocks):
        self.setup_prg_builder()

        def preorder_traversal_to_build_prg_mock_fun(prg_as_list):
            prg_as_list.extend(["a", "b", "c"])
        preorder_traversal_to_build_prg_mock = Mock(side_effect=preorder_traversal_to_build_prg_mock_fun)
        self.prg_builder.root.preorder_traversal_to_build_prg = preorder_traversal_to_build_prg_mock

        expected = "abc"
        actual = self.prg_builder.build_prg()

        self.assertEqual(expected, actual)
        preorder_traversal_to_build_prg_mock.assert_called_once()

    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__, return_value="root_node")
    def test___get_next_site_num(self, *uninteresting_mocks):
        self.setup_prg_builder()
        self.assertEqual(5, self.prg_builder.get_next_site_num())
        self.assertEqual(7, self.prg_builder.get_next_site_num())
        self.assertEqual(9, self.prg_builder.get_next_site_num())

    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__, return_value="root_node")
    def test___get_next_node_id(self, *uninteresting_mocks):
        self.setup_prg_builder()
        self.assertEqual(0, self.prg_builder.get_next_node_id())
        self.assertEqual(1, self.prg_builder.get_next_node_id())
        self.assertEqual(2, self.prg_builder.get_next_node_id())

    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__, return_value="root_node")
    def test___update_prg_index___single_update(self, *uninteresting_mocks):
        self.setup_prg_builder()
        add_indexed_PRG_interval_mock = Mock()
        node_mock_1 = Mock(add_indexed_PRG_interval=add_indexed_PRG_interval_mock)
        self.prg_builder.update_PRG_index(2, 5, node_mock_1)

        expected = {(2, 5): node_mock_1}
        actual = self.prg_builder.prg_index

        self.assertEqual(expected, actual)
        add_indexed_PRG_interval_mock.assert_called_once_with((2, 5))

    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__, return_value="root_node")
    def test___update_prg_index___multiple_updates(self, *uninteresting_mocks):
        self.setup_prg_builder()
        node_mock_1 = Mock()
        node_mock_2 = Mock()
        node_mock_3 = Mock()
        node_mock_4 = Mock()
        self.prg_builder.update_PRG_index(2, 5, node_mock_1)
        self.prg_builder.update_PRG_index(0, 10, node_mock_2)
        self.prg_builder.update_PRG_index(100, 200, node_mock_3)
        self.prg_builder.update_PRG_index(2, 5, node_mock_4)

        expected = {
            (0, 10): node_mock_2,
            (2, 5): node_mock_4,
            (100, 200): node_mock_3,
        }
        actual = self.prg_builder.prg_index

        self.assertEqual(expected, actual)

    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__, return_value="root_node")
    def test___clear_PRG_index(self, *uninteresting_mocks):
        self.setup_prg_builder()

        clear_PRG_interval_index_mock_1 = Mock()
        node_mock_1 = Mock(clear_PRG_interval_index=clear_PRG_interval_index_mock_1)
        clear_PRG_interval_index_mock_2 = Mock()
        node_mock_2 = Mock(clear_PRG_interval_index=clear_PRG_interval_index_mock_2)
        self.prg_builder.update_PRG_index(2, 5, node_mock_1)
        self.prg_builder.update_PRG_index(0, 10, node_mock_2)

        self.prg_builder.clear_PRG_index()

        clear_PRG_interval_index_mock_1.assert_called_once_with()
        clear_PRG_interval_index_mock_2.assert_called_once_with()
        self.assertEqual({}, self.prg_builder.prg_index)

    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__, return_value="root_node")
    def test___get_node_given_interval(self, *uninteresting_mocks):
        self.setup_prg_builder()
        node_mock_1 = Mock()
        node_mock_2 = Mock()
        node_mock_3 = Mock()
        self.prg_builder.update_PRG_index(2, 5, node_mock_1)
        self.prg_builder.update_PRG_index(3, 8, node_mock_2)
        self.prg_builder.update_PRG_index(5, 9, node_mock_3)

        self.assertEqual(node_mock_1, self.prg_builder.get_node_given_interval((2, 5)))
        self.assertEqual(node_mock_2, self.prg_builder.get_node_given_interval((3, 8)))
        self.assertEqual(node_mock_3, self.prg_builder.get_node_given_interval((5, 9)))
        for i in range(11):
            for j in range(11):
                interval = (i, j)
                if interval != (2, 5) and interval != (3, 8) and interval != (5, 9):
                    with self.assertRaises(LeafNotFoundException):
                        self.prg_builder.get_node_given_interval(interval)

    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__, return_value="root_node")
    @patch("builtins.open", new_callable=mock_open)
    @patch("pickle.dump")
    def test___serialize(self, dump_mock, open_mock, *uninteresting_mocks):
        self.setup_prg_builder()
        self.prg_builder.serialize("filepath")
        open_mock.assert_called_once_with("filepath", 'wb')

        # did not manage to assert equality of arg [1] (which is the open_mock handler)
        # so asserting equality of arg [0]
        dump_mock.assert_called_once()
        self.assertEqual(self.prg_builder, dump_mock.call_args[0][0])

    @patch("pickle.loads")
    def test___array_of_bytes(self, loads_mock, *uninteresting_mocks):
        array_of_bytes = b"test"
        PrgBuilder.deserialize_from_bytes(array_of_bytes)
        loads_mock.assert_called_once_with(array_of_bytes)

    def test___write_prg_as_text(self):
        PrgBuilder.write_prg_as_text(f"{workdir}/sample", sample_prg)
        self.assertTrue(filecmp.cmp(workdir / "sample.prg.fa", workdir / "sample.truth.prg.fa"))

    def test___write_prg_as_binary(self):
        PrgBuilder.write_prg_as_binary(f"{workdir}/sample", sample_prg)
        self.assertTrue(filecmp.cmp(workdir / "sample.bin", workdir / "sample.truth.bin"))


    @patch("make_prg.prg_builder.load_alignment_file", return_value="MSA")
    @patch.object(NodeFactory, NodeFactory.build.__name__, return_value="root_node")
    def test___get_state___aligner_is_not_pickled(self, *uninteresting_mocks):
        self.setup_prg_builder()

        expected = self.prg_builder.__dict__
        expected['aligner'] = None
        actual = self.prg_builder.__getstate__()

        self.assertEqual(expected, actual)


class TestPrgBuilderZipDatabase(TestCase):
    def setUp(self) -> None:
        self.zip_filepath = Path("tests/data/utils/io_utils/zip_set_of_files/files.truth.zip")
        self.prg_builder_zip_db = PrgBuilderZipDatabase(self.zip_filepath)

    def tearDown(self) -> None:
        self.prg_builder_zip_db.close()

    def test___constructor___not_a_zip_file___raises_AssertionError(self):
        with self.assertRaises(AssertionError):
            PrgBuilderZipDatabase(Path("file.txt"))

    def test___constructor(self):
        self.assertTrue(are_zip_files_equal(self.zip_filepath, self.prg_builder_zip_db._zip_filepath))

    @patch("make_prg.prg_builder.zip_set_of_files")
    def test___save(self, zip_set_of_files_mock):
        locus_to_prg_builder_pickle_path_mock = Mock()

        self.prg_builder_zip_db.save(locus_to_prg_builder_pickle_path_mock)

        zip_set_of_files_mock.assert_called_once_with(self.prg_builder_zip_db._zip_filepath,
                                                      locus_to_prg_builder_pickle_path_mock)

    def test___load(self):
        self.assertIsNone(self.prg_builder_zip_db._zip_file)
        self.prg_builder_zip_db.load()
        self.assertIsNotNone(self.prg_builder_zip_db._zip_file)

    @patch.object(ZipFile, ZipFile.close.__name__)
    def test___close___zip_file_not_loaded___close_not_called(self, zipfile_close_mock):
        self.prg_builder_zip_db.close()
        zipfile_close_mock.assert_not_called()

    @patch.object(ZipFile, ZipFile.close.__name__)
    def test___close___zip_file_loaded___close_called(self, zipfile_close_mock):
        self.prg_builder_zip_db.load()
        self.prg_builder_zip_db.close()
        zipfile_close_mock.assert_called_once_with()

    def test___get_number_of_loci(self):
        self.prg_builder_zip_db.load()

        expected = 3
        actual = self.prg_builder_zip_db.get_number_of_loci()

        self.assertEqual(expected, actual)

    def test___get_loci_names(self):
        self.prg_builder_zip_db.load()

        expected = ["f1", "file3", "file_2"]
        actual = self.prg_builder_zip_db.get_loci_names()

        self.assertEqual(expected, actual)

    @patch.object(ZipFile, ZipFile.read.__name__, return_value="read_mock")
    @patch.object(PrgBuilder, PrgBuilder.deserialize_from_bytes.__name__)
    def test___get_PrgBuilder(self, deserialize_from_bytes_mock, zipfile_read_mock):
        self.prg_builder_zip_db.load()
        self.prg_builder_zip_db.get_PrgBuilder("locus")

        zipfile_read_mock.assert_called_once_with("locus")
        deserialize_from_bytes_mock.assert_called_once_with("read_mock")

