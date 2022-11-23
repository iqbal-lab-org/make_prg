from unittest import TestCase
from unittest.mock import patch, Mock
from tests.test_helpers import make_alignment
from make_prg.utils.misc import equal_msas
from make_prg.recursion_tree import RecursiveTreeNode, MultiIntervalNode, MultiClusterNode, LeafNode, NodeFactory, UpdateError
from make_prg.prg_builder import PrgBuilder
from make_prg import MSA
from pathlib import Path
from make_prg.from_msa.interval_partition import IntervalType, Interval
from make_prg.utils.seq_utils import SequenceExpander
from make_prg.update.denovo_variants import UpdateData, MLPathError
from make_prg.from_msa.cluster_sequences import ClusteringResult

# Note we use MultiIntervalNodes here instead as RecursiveTreeNode is an abstract class
@patch.object(PrgBuilder, PrgBuilder.get_next_node_id.__name__, return_value=0)
@patch("make_prg.prg_builder.load_alignment_file")
@patch.object(RecursiveTreeNode, RecursiveTreeNode.log_that_node_was_created.__name__)
@patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA",
       return_value="remove_columns_full_of_gaps_from_MSA_mock")
class TestRecursiveTreeNode(TestCase):
    # Note: can't apply patches to setUp(), so creating this method that is called in every test
    def setup(self) -> None:
        self.alignment_mock = Mock()
        with patch.object(NodeFactory, "build"):
            self.prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)
        self.parent_mock = Mock(node_id=512)

    @patch.object(NodeFactory, NodeFactory.build.__name__, side_effect=["child_build_1"])
    def test____get_children___one_child(self, node_factory_build_mock,
                                              *uninteresting_mocks):
        self.setup()
        node = MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder,
                                 interval_subalignments=["algn1"])
        self.assertEqual(["child_build_1"], node.children)


    @patch.object(NodeFactory, NodeFactory.build.__name__, side_effect=[
                "child_build_1", "child_build_2", "child_build_3"])
    def test____get_children___three_children(self, node_factory_build_mock,
                                              *uninteresting_mocks):
        self.setup()
        node = MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder,
                                 interval_subalignments=["algn1", "algn2", "algn3"])
        self.assertEqual(["child_build_1", "child_build_2", "child_build_3"],
                         node.children)

    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__, return_value=["child"])
    def test___is_root___is_indeed_root(self, *uninteresting_mocks):
        self.setup()
        node = MultiIntervalNode(1, self.alignment_mock, None, self.prg_builder, [])
        self.assertTrue(node.is_root())

    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__, return_value=["child"])
    def test___is_root___is_not_root(self, *uninteresting_mocks):
        self.setup()
        node = MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, [])
        self.assertFalse(node.is_root())

    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__,
                  return_value=["child_1", "child_2", "child_3"])
    def test___replace_child___old_child_not_present___raises_Assertion_Error(self, *uninteresting_mocks):
        self.setup()
        node = MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, [])

        with self.assertRaises(AssertionError):
            node.replace_child("inexistent_child", "new_child")

    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__,
                  return_value=["child_1", "child_2", "child_3"])
    def test___replace_child___replaces_first_child(self, *uninteresting_mocks):
        self.setup()
        node = MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, [])

        node.replace_child("child_1", "new_child")

        self.assertEqual(["new_child", "child_2", "child_3"], node.children)

    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__,
                  return_value=["child_1", "child_2", "child_3"])
    def test___replace_child___replaces_second_child(self, *uninteresting_mocks):
        self.setup()
        node = MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, [])

        node.replace_child("child_2", "new_child")

        self.assertEqual(["child_1", "new_child", "child_3"], node.children)

    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__,
                  return_value=["child_1", "child_2", "child_3"])
    def test___replace_child___replaces_third_child(self, *uninteresting_mocks):
        self.setup()
        node = MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, [])

        node.replace_child("child_3", "new_child")

        self.assertEqual(["child_1", "child_2", "new_child"], node.children)




class Preorder_traversal_to_build_prg_Mock:
    def __init__(self, child_id):
        self.child_id = child_id

    def preorder_traversal_to_build_prg(self, prg_as_list, delim_char):
        prg_as_list.append(f"child_{self.child_id}")




@patch.object(PrgBuilder, PrgBuilder.get_next_node_id.__name__, return_value=0)
@patch.object(NodeFactory, NodeFactory.build.__name__)
@patch("make_prg.prg_builder.load_alignment_file")
@patch.object(RecursiveTreeNode, RecursiveTreeNode.log_that_node_was_created.__name__)
class TestMultiIntervalNode(TestCase):
    # Note: can't apply patches to setUp(), so creating this method that is called in every test
    def setup(self) -> None:
        self.alignment_mock = Mock()
        self.prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)
        self.parent_mock = Mock(node_id=512)

    @patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA", return_value="remove_columns_full_of_gaps_from_MSA_mock")
    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__, return_value=["child_1", "child_2"])
    def test___constructor(self, get_children_mock, remove_columns_full_of_gaps_from_MSA_mock, *uninteresting_mocks):
        self.setup()
        children_mock = Mock()
        node = MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, children_mock)

        self.assertEqual(1, node.nesting_level)
        remove_columns_full_of_gaps_from_MSA_mock.assert_called_once_with(self.alignment_mock)
        self.assertTrue("remove_columns_full_of_gaps_from_MSA_mock", node.alignment)
        self.assertEqual(self.parent_mock, node.parent)
        self.assertEqual(0, node.node_id)
        get_children_mock.assert_called_once_with(children_mock)
        self.assertEqual(["child_1", "child_2"], node.children)
        self.assertEqual(self.prg_builder, node.prg_builder)
        self.assertFalse(node.is_leaf())
        self.assertFalse(node.is_root())

    @patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA",
           return_value="remove_columns_full_of_gaps_from_MSA_mock")
    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__, return_value=[])
    def test___constructor___no_children___raises_AssertionError(self, *uninteresting_mocks):
        self.setup()
        children_mock = Mock()
        with self.assertRaises(AssertionError):
            MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, children_mock)

    @patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA",
           return_value="remove_columns_full_of_gaps_from_MSA_mock")
    @patch.object(PrgBuilder, PrgBuilder.get_next_site_num.__name__, return_value=32)
    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__,
                  return_value=[Preorder_traversal_to_build_prg_Mock(1)])
    def test___preorder_traversal_to_build_prg___single_child(self, *uninteresting_mocks):
        self.setup()
        children_mock = Mock()
        node = MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, children_mock)

        prg_as_list = []
        node.preorder_traversal_to_build_prg(prg_as_list, delim_char="*")

        expected_prg = "child_1"
        actual_prg = "".join(prg_as_list)

        self.assertEqual(expected_prg, actual_prg)

    @patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA",
           return_value="remove_columns_full_of_gaps_from_MSA_mock")
    @patch.object(PrgBuilder, PrgBuilder.get_next_site_num.__name__, return_value=32)
    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__,
                  return_value=[Preorder_traversal_to_build_prg_Mock(1),
                     Preorder_traversal_to_build_prg_Mock(2),
                     Preorder_traversal_to_build_prg_Mock(3)])
    def test___preorder_traversal_to_build_prg___multiple_children(self, *uninteresting_mocks):
        self.setup()
        children_mock = Mock()
        node = MultiIntervalNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, children_mock)

        prg_as_list = []
        node.preorder_traversal_to_build_prg(prg_as_list, delim_char="*")

        expected_prg = "child_1child_2child_3"
        actual_prg = "".join(prg_as_list)

        self.assertEqual(expected_prg, actual_prg)

    @patch.object(MultiIntervalNode, MultiIntervalNode._get_children.__name__,
                  return_value=[Mock(node_id=100), Mock(node_id=200)])
    def test___repr_and_str(self, *uninteresting_mocks):
        self.setup()
        alignment = make_alignment(
            ["AAAT", "C--C", "AATT", "GNGG"], ["s1", "s2", "s3", "s4"]
        )
        node = MultiIntervalNode(1, alignment, self.parent_mock, self.prg_builder, [])

        expected = """MultiIntervalNode:
Id = 0
Nesting level = 1
Parent = Id = 512
Children = [Id = 100, Id = 200]
Alignment:
>s1
AAAT
>s2
C--C
>s3
AATT
>s4
GNGG
"""
        actual = repr(node)
        self.assertEqual(expected, actual)

        actual = str(node)
        self.assertEqual(expected, actual)




@patch.object(PrgBuilder, PrgBuilder.get_next_node_id.__name__, return_value=0)
@patch.object(NodeFactory, NodeFactory.build.__name__)
@patch("make_prg.prg_builder.load_alignment_file")
@patch.object(RecursiveTreeNode, RecursiveTreeNode.log_that_node_was_created.__name__)
class TestMultiClusterNode(TestCase):
    # Note: can't apply patches to setUp(), so creating this method that is called in every test
    def setup(self) -> None:
        self.alignment_mock = Mock()
        self.prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)
        self.parent_mock = Mock(node_id=512)

    @patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA", return_value="remove_columns_full_of_gaps_from_MSA_mock")
    @patch.object(MultiClusterNode, MultiClusterNode._get_children.__name__, return_value=["child_1", "child_2"])
    def test___constructor(self, get_children_mock, remove_columns_full_of_gaps_from_MSA_mock, *uninteresting_mocks):
        self.setup()
        children_mock = Mock()
        node = MultiClusterNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, children_mock)

        self.assertEqual(1, node.nesting_level)
        remove_columns_full_of_gaps_from_MSA_mock.assert_called_once_with(self.alignment_mock)
        self.assertTrue("remove_columns_full_of_gaps_from_MSA_mock", node.alignment)
        self.assertEqual(self.parent_mock, node.parent)
        self.assertEqual(0, node.node_id)
        get_children_mock.assert_called_once_with(children_mock)
        self.assertEqual(["child_1", "child_2"], node.children)
        self.assertEqual(self.prg_builder, node.prg_builder)
        self.assertFalse(node.is_leaf())
        self.assertFalse(node.is_root())

    @patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA",
           return_value="remove_columns_full_of_gaps_from_MSA_mock")
    @patch.object(MultiClusterNode, MultiClusterNode._get_children.__name__, return_value=[])
    def test___constructor___no_children___raises_AssertionError(self, *uninteresting_mocks):
        self.setup()
        children_mock = Mock()
        with self.assertRaises(AssertionError):
            MultiClusterNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, children_mock)

    @patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA",
           return_value="remove_columns_full_of_gaps_from_MSA_mock")
    @patch.object(PrgBuilder, PrgBuilder.get_next_site_num.__name__, return_value=32)
    @patch.object(MultiClusterNode, MultiClusterNode._get_children.__name__,
                  return_value=[Preorder_traversal_to_build_prg_Mock(1)])
    def test___preorder_traversal_to_build_prg___single_child(self, *uninteresting_mocks):
        self.setup()
        children_mock = Mock()
        node = MultiClusterNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, children_mock)

        prg_as_list = []
        node.preorder_traversal_to_build_prg(prg_as_list, delim_char="*")

        expected_prg = "*32*child_1*32*"
        actual_prg = "".join(prg_as_list)

        self.assertEqual(expected_prg, actual_prg)

    @patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA",
           return_value="remove_columns_full_of_gaps_from_MSA_mock")
    @patch.object(PrgBuilder, PrgBuilder.get_next_site_num.__name__, return_value=32)
    @patch.object(MultiClusterNode, MultiClusterNode._get_children.__name__,
                  return_value=[Preorder_traversal_to_build_prg_Mock(1),
                     Preorder_traversal_to_build_prg_Mock(2),
                     Preorder_traversal_to_build_prg_Mock(3)])
    def test___preorder_traversal_to_build_prg___multiple_children(self, *uninteresting_mocks):
        self.setup()
        children_mock = Mock()
        node = MultiClusterNode(1, self.alignment_mock, self.parent_mock, self.prg_builder, children_mock)

        prg_as_list = []
        node.preorder_traversal_to_build_prg(prg_as_list, delim_char="*")

        expected_prg = "*32*child_1*33*child_2*33*child_3*32*"
        actual_prg = "".join(prg_as_list)

        self.assertEqual(expected_prg, actual_prg)

    @patch.object(MultiClusterNode, MultiClusterNode._get_children.__name__,
                  return_value=[Mock(node_id=100), Mock(node_id=200)])
    def test___repr_and_str(self, *uninteresting_mocks):
        self.setup()
        alignment = make_alignment(
            ["AAAT", "C--C", "AATT", "GNGG"], ["s1", "s2", "s3", "s4"]
        )
        node = MultiClusterNode(1, alignment, self.parent_mock, self.prg_builder, [])

        expected = """MultiClusterNode:
Id = 0
Nesting level = 1
Parent = Id = 512
Children = [Id = 100, Id = 200]
Alignment:
>s1
AAAT
>s2
C--C
>s3
AATT
>s4
GNGG
"""
        actual = repr(node)
        self.assertEqual(expected, actual)

        actual = str(node)
        self.assertEqual(expected, actual)



@patch.object(PrgBuilder, PrgBuilder.get_next_node_id.__name__, return_value=0)
@patch.object(NodeFactory, NodeFactory.build.__name__)
@patch("make_prg.prg_builder.load_alignment_file")
@patch.object(RecursiveTreeNode, RecursiveTreeNode.log_that_node_was_created.__name__)
class TestLeafNode(TestCase):
    # Note: can't apply patches to setUp(), so creating this method that is called in every test
    def setup(self) -> None:
        self.alignment = make_alignment(
            ["AAAT", "C--C", "AATT", "GNGG"], ["s1", "s2", "s3", "s4"]
        )
        self.prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)
        self.parent_mock = Mock(node_id=512)
        self.node = LeafNode(1, self.alignment, self.parent_mock, self.prg_builder)

    def test___constructor(self, *uninteresting_mocks):
        self.setup()

        self.assertEqual(1, self.node.nesting_level)
        self.assertTrue(self.alignment, self.node.alignment)
        self.assertEqual(self.parent_mock, self.node.parent)
        self.assertEqual(0, self.node.node_id)
        self.assertEqual([], self.node.children)
        self.assertEqual(self.prg_builder, self.node.prg_builder)
        self.assertEqual(set(), self.node.new_sequences)
        self.assertTrue(self.node.is_leaf())
        self.assertFalse(self.node.is_root())


    @patch.object(SequenceExpander, SequenceExpander.get_expanded_sequences_from_MSA.__name__, return_value = ["ACGT"])
    @patch.object(PrgBuilder, PrgBuilder.update_PRG_index.__name__)
    def test___preorder_traversal_to_build_prg___single_seq(self, update_prg_index_mock, *uninteresting_mocks):
        self.setup()

        expected_prg_as_list = list("ACGT")
        actual_prg_as_list = []
        self.node.preorder_traversal_to_build_prg(actual_prg_as_list, delim_char="*")

        self.assertEqual(expected_prg_as_list, actual_prg_as_list)
        update_prg_index_mock.assert_called_once_with(0, 4, node=self.node)

    @patch.object(SequenceExpander, SequenceExpander.get_expanded_sequences_from_MSA.__name__, return_value=["ACGT"])
    @patch.object(PrgBuilder, PrgBuilder.update_PRG_index.__name__)
    def test___preorder_traversal_to_build_prg___single_seq___no_indexing(self, update_prg_index_mock, *uninteresting_mocks):
        self.setup()

        expected_prg_as_list = list("ACGT")
        actual_prg_as_list = []
        self.node.preorder_traversal_to_build_prg(actual_prg_as_list, delim_char="*", do_indexing=False)

        self.assertEqual(expected_prg_as_list, actual_prg_as_list)
        update_prg_index_mock.assert_not_called()

    @patch.object(SequenceExpander, SequenceExpander.get_expanded_sequences_from_MSA.__name__, return_value=[
        "AA", "C", "GGGG"])
    @patch.object(PrgBuilder, PrgBuilder.get_next_site_num.__name__, return_value=42)
    @patch.object(PrgBuilder, PrgBuilder.update_PRG_index.__name__)
    def test___preorder_traversal_to_build_prg___several_seqs(self, update_prg_index_mock, *uninteresting_mocks):
        self.setup()

        expected_prg_as_list = list("*42*AA*43*C*43*GGGG*42*")
        actual_prg_as_list = []
        self.node.preorder_traversal_to_build_prg(actual_prg_as_list, delim_char="*")

        self.assertEqual(expected_prg_as_list, actual_prg_as_list)
        self.assertEqual(3, update_prg_index_mock.call_count)
        update_prg_index_mock.assert_any_call(4, 6, node=self.node)
        update_prg_index_mock.assert_any_call(10, 11, node=self.node)
        update_prg_index_mock.assert_any_call(15, 19, node=self.node)

    @patch.object(SequenceExpander, SequenceExpander.get_expanded_sequences_from_MSA.__name__, return_value=[
        "AA", "C", "GGGG"])
    @patch.object(PrgBuilder, PrgBuilder.get_next_site_num.__name__, return_value=42)
    @patch.object(PrgBuilder, PrgBuilder.update_PRG_index.__name__)
    def test___preorder_traversal_to_build_prg___several_seqs___no_indexing(self, update_prg_index_mock, *uninteresting_mocks):
        self.setup()

        expected_prg_as_list = list("*42*AA*43*C*43*GGGG*42*")
        actual_prg_as_list = []
        self.node.preorder_traversal_to_build_prg(actual_prg_as_list, delim_char="*", do_indexing=False)

        self.assertEqual(expected_prg_as_list, actual_prg_as_list)
        update_prg_index_mock.assert_not_called()


    def test___add_data_to_batch_update___interval_not_indexed___raises_UpdateError(self, *uninteresting_mocks):
        self.setup()
        ml_path_mock = Mock()
        update_data = UpdateData((0, 4), ml_path_mock, "ACGT")
        self.node.indexed_PRG_intervals = {}
        self.assertEqual(set(), self.node.new_sequences)

        with self.assertRaises(UpdateError):
            self.node.add_data_to_batch_update(update_data)

    def test___add_data_to_batch_update___single_interval(self, *uninteresting_mocks):
        self.setup()
        ml_path_mock = Mock()
        update_data = UpdateData((0, 4), ml_path_mock, "ACGT")
        self.node.indexed_PRG_intervals = {(0, 4)}
        self.assertEqual(set(), self.node.new_sequences)

        self.node.add_data_to_batch_update(update_data)

        expected = {"ACGT"}
        self.assertEqual(expected, self.node.new_sequences)

    def test___add_data_to_batch_update___update_data_interval_flanked(self, *uninteresting_mocks):
        self.setup()

        def get_node_given_interval_in_PRG_space_mock_fun(interval):
            if interval == (0, 4):
                return Mock(sequence="left_flank")
            elif interval == (20, 30):
                return Mock(sequence="right_flank")
            else:
                assert False, "Error on get_node_given_interval_in_PRG_space_mock_fun"

        get_node_given_interval_in_PRG_space_mock = Mock(side_effect=get_node_given_interval_in_PRG_space_mock_fun)
        ml_path_node_mock = Mock(get_node_given_interval_in_PRG_space = get_node_given_interval_in_PRG_space_mock)
        update_data = UpdateData((5, 10), ml_path_node_mock, "ACGT")
        self.node.indexed_PRG_intervals = {(5, 10), (20, 30), (0, 4)}
        self.assertEqual(set(), self.node.new_sequences)

        self.node.add_data_to_batch_update(update_data)

        expected = {"left_flankACGTright_flank"}
        self.assertEqual(expected, self.node.new_sequences)
        get_node_given_interval_in_PRG_space_mock.assert_any_call((0, 4))
        get_node_given_interval_in_PRG_space_mock.assert_any_call((20, 30))

    def test___add_data_to_batch_update___update_data_interval_on_the_left_end(self, *uninteresting_mocks):
        self.setup()

        def get_node_given_interval_in_PRG_space_mock_fun(interval):
            if interval == (5, 10):
                return Mock(sequence="middle_flank")
            elif interval == (20, 30):
                return Mock(sequence="right_flank")
            else:
                assert False, "Error on get_node_given_interval_in_PRG_space_mock_fun"

        get_node_given_interval_in_PRG_space_mock = Mock(side_effect=get_node_given_interval_in_PRG_space_mock_fun)
        ml_path_node_mock = Mock(get_node_given_interval_in_PRG_space = get_node_given_interval_in_PRG_space_mock)
        update_data = UpdateData((0, 4), ml_path_node_mock, "ACGT")
        self.node.indexed_PRG_intervals = {(5, 10), (20, 30), (0, 4)}
        self.assertEqual(set(), self.node.new_sequences)

        self.node.add_data_to_batch_update(update_data)

        expected = {"ACGTmiddle_flankright_flank"}
        self.assertEqual(expected, self.node.new_sequences)
        get_node_given_interval_in_PRG_space_mock.assert_any_call((5, 10))
        get_node_given_interval_in_PRG_space_mock.assert_any_call((20, 30))

    def test___add_data_to_batch_update___update_data_interval_on_the_right_end(self, *uninteresting_mocks):
        self.setup()

        def get_node_given_interval_in_PRG_space_mock_fun(interval):
            if interval == (0, 4):
                return Mock(sequence="left_flank")
            elif interval == (5, 10):
                return Mock(sequence="middle_flank")
            else:
                assert False, "Error on get_node_given_interval_in_PRG_space_mock_fun"

        get_node_given_interval_in_PRG_space_mock = Mock(side_effect=get_node_given_interval_in_PRG_space_mock_fun)
        ml_path_node_mock = Mock(get_node_given_interval_in_PRG_space = get_node_given_interval_in_PRG_space_mock)
        update_data = UpdateData((20, 30), ml_path_node_mock, "ACGT")
        self.node.indexed_PRG_intervals = {(5, 10), (20, 30), (0, 4)}
        self.assertEqual(set(), self.node.new_sequences)

        self.node.add_data_to_batch_update(update_data)

        expected = {"left_flankmiddle_flankACGT"}
        self.assertEqual(expected, self.node.new_sequences)
        get_node_given_interval_in_PRG_space_mock.assert_any_call((0, 4))
        get_node_given_interval_in_PRG_space_mock.assert_any_call((5, 10))

    def test___add_data_to_batch_update___update_data_interval_flanked___left_interval_not_in_ml_path(self, *uninteresting_mocks):
        self.setup()

        def get_node_given_interval_in_PRG_space_mock_fun(interval):
            if interval == (20, 30):
                return Mock(sequence="right_flank")
            else:
                raise MLPathError()

        get_node_given_interval_in_PRG_space_mock = Mock(side_effect=get_node_given_interval_in_PRG_space_mock_fun)
        ml_path_node_mock = Mock(get_node_given_interval_in_PRG_space = get_node_given_interval_in_PRG_space_mock)
        update_data = UpdateData((5, 10), ml_path_node_mock, "ACGT")
        self.node.indexed_PRG_intervals = {(5, 10), (20, 30), (0, 4)}
        self.assertEqual(set(), self.node.new_sequences)

        self.node.add_data_to_batch_update(update_data)

        expected = {"ACGTright_flank"}
        self.assertEqual(expected, self.node.new_sequences)
        get_node_given_interval_in_PRG_space_mock.assert_any_call((0, 4))
        get_node_given_interval_in_PRG_space_mock.assert_any_call((20, 30))

    def test___add_data_to_batch_update___update_data_interval_flanked___right_interval_not_in_ml_path(self, *uninteresting_mocks):
        self.setup()

        def get_node_given_interval_in_PRG_space_mock_fun(interval):
            if interval == (0, 4):
                return Mock(sequence="left_flank")
            else:
                raise MLPathError()

        get_node_given_interval_in_PRG_space_mock = Mock(side_effect=get_node_given_interval_in_PRG_space_mock_fun)
        ml_path_node_mock = Mock(get_node_given_interval_in_PRG_space = get_node_given_interval_in_PRG_space_mock)
        update_data = UpdateData((5, 10), ml_path_node_mock, "ACGT")
        self.node.indexed_PRG_intervals = {(5, 10), (20, 30), (0, 4)}
        self.assertEqual(set(), self.node.new_sequences)

        self.node.add_data_to_batch_update(update_data)

        expected = {"left_flankACGT"}
        self.assertEqual(expected, self.node.new_sequences)
        get_node_given_interval_in_PRG_space_mock.assert_any_call((0, 4))
        get_node_given_interval_in_PRG_space_mock.assert_any_call((20, 30))

    def test___add_indexed_PRG_interval___single_interval(self, *uninteresting_mocks):
        self.setup()

        self.assertEqual(set(), self.node.indexed_PRG_intervals)
        self.node.add_indexed_PRG_interval((0, 4))
        self.assertEqual({(0, 4)}, self.node.indexed_PRG_intervals)

    def test___add_indexed_PRG_interval___multiple_intervals(self, *uninteresting_mocks):
        self.setup()

        self.assertEqual(set(), self.node.indexed_PRG_intervals)
        self.node.add_indexed_PRG_interval((0, 4))
        self.node.add_indexed_PRG_interval((25, 30))
        self.node.add_indexed_PRG_interval((0, 4))
        self.node.add_indexed_PRG_interval((0, 4))
        self.node.add_indexed_PRG_interval((25, 30))
        self.node.add_indexed_PRG_interval((10, 15))
        self.node.add_indexed_PRG_interval((25, 30))
        self.node.add_indexed_PRG_interval((0, 4))
        self.node.add_indexed_PRG_interval((17, 19))
        self.assertEqual({(0, 4), (10, 15), (17, 19), (25, 30)},
                         self.node.indexed_PRG_intervals)

    @patch.object(LeafNode, LeafNode._update_leaf.__name__)
    def test___batch_update___no_new_sequences___no_update_to_be_done(self, update_leaf_mock, *uninteresting_mocks):
        self.setup()
        self.node.batch_update()
        update_leaf_mock.assert_not_called()

    @patch.object(LeafNode, LeafNode._update_leaf.__name__)
    def test___batch_update___one_new_sequence___update_is_called(self, update_leaf_mock, *uninteresting_mocks):
        self.setup()
        self.node.new_sequences.add("ACGT")
        self.node.batch_update()
        update_leaf_mock.assert_called_once_with()

    def test___update_leaf___no_aligner_was_given___raises_AssertionError(self, *uninteresting_mocks):
        self.setup()
        self.node.new_sequences.add("ACGT")
        with self.assertRaises(AssertionError):
            self.node._update_leaf()

    def test___update_leaf___is_not_root(self, *uninteresting_mocks):
        self.setup()

        # prepare aligner mock
        updated_alignment = Mock()
        get_updated_alignment_mock = Mock(return_value=updated_alignment)
        aligner_mock = Mock()
        aligner_mock.get_updated_alignment = get_updated_alignment_mock
        self.prg_builder.aligner = aligner_mock

        # prepare replace_child mock
        replace_child_mock = Mock()
        self.parent_mock.replace_child = replace_child_mock

        # prepare replace_root mock
        replace_root_mock = Mock()
        self.prg_builder.replace_root = replace_root_mock

        # prepare clear_PRG_index mock
        clear_PRG_index_mock = Mock()
        self.prg_builder.clear_PRG_index = clear_PRG_index_mock

        self.node.new_sequences.add("AAAA")
        self.node.new_sequences.add("CC")

        # we repatch some methods
        updated_child = Mock()
        with patch.object(NodeFactory, "build", return_value=updated_child) as node_factory_build_mock:
            self.node._update_leaf()

            # assert_called_once_with() does not work because can't check if MSAs are equal
            get_updated_alignment_mock.assert_called_once()
            kwargs = get_updated_alignment_mock.call_args_list[0][1]
            self.assertTrue(equal_msas(self.alignment, kwargs["current_alignment"]))
            self.assertEqual({"AAAA", "CC"}, kwargs["new_sequences"])

            node_factory_build_mock.assert_called_once_with(
                updated_alignment, self.prg_builder, self.parent_mock
            )

            # this is not a root, child is replaced
            replace_child_mock.assert_called_once_with(self.node, updated_child)
            replace_root_mock.assert_not_called()

            clear_PRG_index_mock.assert_called_once_with()

    def test___update_leaf___is_root(self, *uninteresting_mocks):
        self.setup()
        self.node = LeafNode(1, self.alignment, None, self.prg_builder)

        # prepare aligner mock
        updated_alignment = Mock()
        get_updated_alignment_mock = Mock(return_value=updated_alignment)
        aligner_mock = Mock()
        aligner_mock.get_updated_alignment = get_updated_alignment_mock
        self.prg_builder.aligner = aligner_mock

        # prepare replace_child mock
        replace_child_mock = Mock()
        self.parent_mock.replace_child = replace_child_mock

        # prepare replace_root mock
        replace_root_mock = Mock()
        self.prg_builder.replace_root = replace_root_mock

        # prepare clear_PRG_index mock
        clear_PRG_index_mock = Mock()
        self.prg_builder.clear_PRG_index = clear_PRG_index_mock

        self.node.new_sequences.add("AAAA")
        self.node.new_sequences.add("CC")

        # we repatch some methods
        updated_child = Mock()
        with patch.object(NodeFactory, "build", return_value=updated_child) as node_factory_build_mock:
            self.node._update_leaf()

            # assert_called_once_with() does not work because can't check if MSAs are equal
            get_updated_alignment_mock.assert_called_once()
            kwargs = get_updated_alignment_mock.call_args_list[0][1]
            self.assertTrue(equal_msas(self.alignment, kwargs["current_alignment"]))
            self.assertEqual({"AAAA", "CC"}, kwargs["new_sequences"])

            node_factory_build_mock.assert_called_once_with(
                updated_alignment, self.prg_builder, None
            )

            # this is a root, it is replaced
            replace_child_mock.assert_not_called()
            replace_root_mock.assert_called_once_with(updated_child)

            clear_PRG_index_mock.assert_called_once_with()

    def test___clear_PRG_interval_index(self, *uninteresting_mocks):
        self.setup()
        self.node.indexed_PRG_intervals.add((0, 1))
        self.node.indexed_PRG_intervals.add((10, 20))
        self.node.indexed_PRG_intervals.add((30, 100))
        self.assertEqual(3, len(self.node.indexed_PRG_intervals))

        self.node.clear_PRG_interval_index()
        self.assertEqual(0, len(self.node.indexed_PRG_intervals))


@patch.object(PrgBuilder, PrgBuilder.get_next_node_id.__name__, return_value=0)
@patch("make_prg.prg_builder.load_alignment_file")
@patch.object(RecursiveTreeNode, RecursiveTreeNode.log_that_node_was_created.__name__)
class TestNodeFactory(TestCase):
    def setup(self) -> None:
        self.alignment = make_alignment(["AAAT", "C--C", "AATT", "GNGG"], ["s1", "s2", "s3", "s4"])
        with patch.object(NodeFactory, NodeFactory.build.__name__):
            self.prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)
        self.parent_mock = Mock(node_id=512)

    @patch.object(NodeFactory, NodeFactory._get_vertical_partition.__name__, return_value=(Mock(), Mock()))
    @patch.object(NodeFactory, NodeFactory._is_single_match_interval.__name__, return_value=False)
    @patch.object(NodeFactory, NodeFactory._partition_alignment_into_interval_subalignments.__name__,
                  return_value = "interval_subalignments_mock")
    def test___build___root___multi_interval(self, *uninteresting_mocks):
        self.setup()

        # mock MultiIntervalNode.__init__
        def __init__(node_self, nesting_level, alignment, parent, prg_builder, interval_subalignments):
            self.assertEqual(0, nesting_level)
            self.assertEqual(self.alignment, alignment)
            self.assertEqual(None, parent)
            self.assertEqual(self.prg_builder, prg_builder)
            self.assertEqual("interval_subalignments_mock", interval_subalignments)

        with patch.object(MultiIntervalNode, '__init__', __init__):
            node = NodeFactory.build(self.alignment, self.prg_builder, None)
            self.assertTrue(isinstance(node, MultiIntervalNode))

    @patch.object(NodeFactory, NodeFactory._get_vertical_partition.__name__, return_value=(Mock(), Mock()))
    @patch.object(NodeFactory, NodeFactory._is_single_match_interval.__name__, return_value=True)
    def test___build___root___leaf(self, *uninteresting_mocks):
        self.setup()

        # mock Leaf.__init__
        def __init__(node_self, nesting_level, alignment, parent, prg_builder):
            self.assertEqual(0, nesting_level)
            self.assertEqual(self.alignment, alignment)
            self.assertEqual(None, parent)
            self.assertEqual(self.prg_builder, prg_builder)

        with patch.object(LeafNode, '__init__', __init__):
            node = NodeFactory.build(self.alignment, self.prg_builder, None)
            self.assertTrue(isinstance(node, LeafNode))

    @patch.object(NodeFactory, NodeFactory._get_vertical_partition.__name__, return_value=(Mock(), Mock()))
    @patch.object(NodeFactory, NodeFactory._is_single_match_interval.__name__, return_value=False)
    @patch.object(NodeFactory, NodeFactory._partition_alignment_into_interval_subalignments.__name__,
                  return_value="interval_subalignments_mock")
    def test___build___non_root___parent_is_multi_cluster_creates_multi_interval(self, *uninteresting_mocks):
        self.setup()

        # mock MultiClusterNode.__init__
        def __MultiClusterNode_init__(node_self, nesting_level):
            node_self.nesting_level = nesting_level
        with patch.object(MultiClusterNode, '__init__', __MultiClusterNode_init__):
            outer_parent = MultiClusterNode(4)

        # mock MultiIntervalNode.__init__
        def __MultiIntervalNode_init__(node_self, nesting_level, alignment, parent, prg_builder, interval_subalignments):
            self.assertEqual(4, nesting_level)
            self.assertEqual(id(self.alignment), id(alignment))
            self.assertEqual(id(outer_parent), id(parent))
            self.assertEqual(id(self.prg_builder), id(prg_builder))
            self.assertEqual("interval_subalignments_mock", interval_subalignments)

        with patch.object(MultiIntervalNode, '__init__', __MultiIntervalNode_init__):
            node = NodeFactory.build(self.alignment, self.prg_builder, outer_parent)
            self.assertTrue(isinstance(node, MultiIntervalNode))

    @patch.object(NodeFactory, NodeFactory._get_vertical_partition.__name__, return_value=(Mock(), Mock()))
    @patch.object(NodeFactory, NodeFactory._is_single_match_interval.__name__, return_value=True)
    def test___build___non_root___parent_is_multi_cluster_creates_leaf(self, *uninteresting_mocks):
        self.setup()

        # mock MultiClusterNode.__init__
        def __MultiClusterNode_init__(node_self, nesting_level):
            node_self.nesting_level = nesting_level
        with patch.object(MultiClusterNode, '__init__', __MultiClusterNode_init__):
            outer_parent = MultiClusterNode(4)

        # mock Leaf.__init__
        def __Leaf_init__(node_self, nesting_level, alignment, parent, prg_builder):
            self.assertEqual(4, nesting_level)
            self.assertEqual(id(self.alignment), id(alignment))
            self.assertEqual(id(outer_parent), id(parent))
            self.assertEqual(id(self.prg_builder), id(prg_builder))
        with patch.object(LeafNode, '__init__', __Leaf_init__):
            node = NodeFactory.build(self.alignment, self.prg_builder, outer_parent)
            self.assertTrue(isinstance(node, LeafNode))

    @patch("make_prg.recursion_tree.kmeans_cluster_seqs")
    @patch.object(NodeFactory, NodeFactory._infer_if_we_should_cluster_further.__name__, return_value=True)
    @patch.object(NodeFactory, NodeFactory._get_subalignments_by_clustering.__name__,
                  return_value="get_subalignments_by_clustering_mock")
    def test___build___non_root___parent_is_multi_interval_creates_multi_cluster(self, *uninteresting_mocks):
        self.setup()

        # mock MultiIntervalNode.__init__
        def __MultiIntervalNode_init__(node_self, nesting_level):
            node_self.nesting_level = nesting_level
        with patch.object(MultiIntervalNode, '__init__', __MultiIntervalNode_init__):
            outer_parent = MultiIntervalNode(4)

        # mock MultiClusterNode.__init__
        def __MultiClusterNode_init__(node_self, nesting_level, alignment, parent, prg_builder,
                                       cluster_subalignments):
            self.assertEqual(5, nesting_level)
            self.assertEqual(id(self.alignment), id(alignment))
            self.assertEqual(id(outer_parent), id(parent))
            self.assertEqual(id(self.prg_builder), id(prg_builder))
            self.assertEqual("get_subalignments_by_clustering_mock", cluster_subalignments)

        with patch.object(MultiClusterNode, '__init__', __MultiClusterNode_init__):
            node = NodeFactory.build(self.alignment, self.prg_builder, outer_parent)
            self.assertTrue(isinstance(node, MultiClusterNode))

    @patch("make_prg.recursion_tree.kmeans_cluster_seqs")
    @patch.object(NodeFactory, NodeFactory._infer_if_we_should_cluster_further.__name__, return_value=False)
    def test___build___non_root___parent_is_multi_interval_creates_leaf(self, *uninteresting_mocks):
        self.setup()

        # mock MultiIntervalNode.__init__
        def __MultiIntervalNode_init__(node_self, nesting_level):
            node_self.nesting_level = nesting_level
        with patch.object(MultiIntervalNode, '__init__', __MultiIntervalNode_init__):
            outer_parent = MultiIntervalNode(4)

        # mock LeafNode.__init__
        def __LeafNode_init__(node_self, nesting_level, alignment, parent, prg_builder):
            self.assertEqual(4, nesting_level)
            self.assertEqual(id(self.alignment), id(alignment))
            self.assertEqual(id(outer_parent), id(parent))
            self.assertEqual(id(self.prg_builder), id(prg_builder))

        with patch.object(LeafNode, '__init__', __LeafNode_init__):
            node = NodeFactory.build(self.alignment, self.prg_builder, outer_parent)
            self.assertTrue(isinstance(node, LeafNode))

    def test___build___non_root___parent_is_leaf_ValueError_is_raised(self, *uninteresting_mocks):
        self.setup()

        # mock LeafNode.__init__
        def __LeafNode_init__(node_self, nesting_level):
            node_self.nesting_level = nesting_level
        with patch.object(LeafNode, '__init__', __LeafNode_init__):
            outer_parent = LeafNode(4)

        with self.assertRaises(ValueError):
            NodeFactory.build(self.alignment, self.prg_builder, outer_parent)

    @patch("make_prg.recursion_tree.get_number_of_unique_ungapped_sequences", return_value=1)
    def test___alignment_has_issues___too_few_unique_sequences_1(self,
                 get_number_of_unique_ungapped_sequences_mock, *uninteresting_mocks):
        self.setup()
        self.assertTrue(NodeFactory._alignment_has_issues(self.alignment))
        get_number_of_unique_ungapped_sequences_mock.assert_called_once_with(self.alignment)

    @patch("make_prg.recursion_tree.get_number_of_unique_ungapped_sequences", return_value=2)
    def test___alignment_has_issues___too_few_unique_sequences_2(self,
                 get_number_of_unique_ungapped_sequences_mock, *uninteresting_mocks):
        self.setup()
        self.assertTrue(NodeFactory._alignment_has_issues(self.alignment))
        get_number_of_unique_ungapped_sequences_mock.assert_called_once_with(self.alignment)

    @patch("make_prg.recursion_tree.get_number_of_unique_ungapped_sequences", return_value=3)
    @patch("make_prg.recursion_tree.get_number_of_unique_gapped_sequences", return_value=4)
    def test___alignment_has_issues___alignment_has_ambiguity___ungapped_smaller_gapped(self,
         get_number_of_unique_gapped_sequences_mock, get_number_of_unique_ungapped_sequences_mock, *uninteresting_mocks):
        self.setup()
        self.assertTrue(NodeFactory._alignment_has_issues(self.alignment))
        get_number_of_unique_gapped_sequences_mock.assert_called_once_with(self.alignment)
        get_number_of_unique_ungapped_sequences_mock.assert_called_once_with(self.alignment)

    @patch("make_prg.recursion_tree.get_number_of_unique_ungapped_sequences", return_value=4)
    @patch("make_prg.recursion_tree.get_number_of_unique_gapped_sequences", return_value=4)
    def test___alignment_has_issues___alignment_has_ambiguity___ungapped_equals_gapped___no_issues(self,
        get_number_of_unique_gapped_sequences_mock, get_number_of_unique_ungapped_sequences_mock, *uninteresting_mocks):
        self.setup()
        self.assertFalse(NodeFactory._alignment_has_issues(self.alignment))
        get_number_of_unique_gapped_sequences_mock.assert_called_once_with(self.alignment)
        get_number_of_unique_ungapped_sequences_mock.assert_called_once_with(self.alignment)

    @patch("make_prg.recursion_tree.get_number_of_unique_ungapped_sequences", return_value=5)
    @patch("make_prg.recursion_tree.get_number_of_unique_gapped_sequences", return_value=4)
    def test___alignment_has_issues___alignment_has_ambiguity___ungapped_larger_gapped___raises_AssertionError(self,
        *uninteresting_mocks):
        self.setup()
        with self.assertRaises(AssertionError):
            NodeFactory._alignment_has_issues(self.alignment)

    def test___get_vertical_partition(self, *uninteresting_mocks):
        msa = make_alignment(["AAAAATTTTTGGGGG", "AAAAACCCCCGGGGG"])

        all_intervals, match_intervals = NodeFactory._get_vertical_partition(msa, 3)

        first_interval = Interval(IntervalType.Match, 0, 4)
        second_interval = Interval(IntervalType.NonMatch, 5, 9)
        third_interval = Interval(IntervalType.Match, 10, 14)
        expected_all_intervals = [first_interval, second_interval, third_interval]
        expected_match_intervals = [first_interval, third_interval]
        self.assertEqual(expected_all_intervals, all_intervals)
        self.assertEqual(expected_match_intervals, match_intervals)

    def test___is_single_match_interval___single_match_interval(self, *uninteresting_mocks):
        interval = Interval(IntervalType.Match, 3, 10)
        all_intervals = [interval]
        match_intervals = [interval]
        self.assertTrue(NodeFactory._is_single_match_interval(all_intervals, match_intervals))

    def test___is_single_match_interval___no_intervals(self, *uninteresting_mocks):
        self.assertFalse(NodeFactory._is_single_match_interval([], []))

    def test___is_single_match_interval___two_match_intervals(self, *uninteresting_mocks):
        interval_1 = Interval(IntervalType.Match, 3, 10)
        interval_2 = Interval(IntervalType.Match, 30, 100)
        all_intervals = [interval_1, interval_2]
        match_intervals = [interval_1, interval_2]
        self.assertFalse(NodeFactory._is_single_match_interval(all_intervals, match_intervals))

    def test___is_single_match_interval___single_mismatch_interval(self, *uninteresting_mocks):
        interval_1 = Interval(IntervalType.NonMatch, 3, 10)
        all_intervals = [interval_1]
        match_intervals = []
        self.assertFalse(NodeFactory._is_single_match_interval(all_intervals, match_intervals))

    def test___is_single_match_interval___two_mismatch_intervals(self, *uninteresting_mocks):
        interval_1 = Interval(IntervalType.NonMatch, 3, 10)
        interval_2 = Interval(IntervalType.NonMatch, 30, 100)
        all_intervals = [interval_1, interval_2]
        match_intervals = []
        self.assertFalse(NodeFactory._is_single_match_interval(all_intervals, match_intervals))

    def test___is_single_match_interval___single_match_interval_but_two_intervals(self, *uninteresting_mocks):
        interval_1 = Interval(IntervalType.Match, 3, 10)
        interval_2 = Interval(IntervalType.NonMatch, 30, 100)
        all_intervals = [interval_1, interval_2]
        match_intervals = [interval_1]
        self.assertFalse(NodeFactory._is_single_match_interval(all_intervals, match_intervals))


    def test___partition_alignment_into_interval_subalignments(self, *uninteresting_mocks):
        msa = make_alignment(["AAAAATTTTTGGGGG", "AAAAACCCCCGGGGG"])
        first_interval = Interval(IntervalType.Match, 0, 4)
        second_interval = Interval(IntervalType.NonMatch, 5, 9)
        third_interval = Interval(IntervalType.Match, 10, 14)
        all_intervals = [first_interval, second_interval, third_interval]

        interval_subalignments = NodeFactory._partition_alignment_into_interval_subalignments(msa, all_intervals)

        expected_interval_subalignments = [
            make_alignment(["AAAAA", "AAAAA"]),
            make_alignment(["TTTTT", "CCCCC"]),
            make_alignment(["GGGGG", "GGGGG"])
        ]
        self.assertEqual(len(expected_interval_subalignments), len(interval_subalignments))
        self.assertTrue(all(map(equal_msas, expected_interval_subalignments, interval_subalignments)))

    def test___infer_if_we_should_cluster_further___no_clustering(self, *uninteresting_mocks):
        clustering_result_mock = Mock(no_clustering=True)
        self.assertFalse(NodeFactory._infer_if_we_should_cluster_further(Mock(), clustering_result_mock, 3, 5))

    def test___infer_if_we_should_cluster_further___max_nesting_reached___smaller(self, *uninteresting_mocks):
        clustering_result_mock = Mock(no_clustering=False)
        self.assertFalse(NodeFactory._infer_if_we_should_cluster_further(Mock(), clustering_result_mock, 4, 5))

    def test___infer_if_we_should_cluster_further___max_nesting_reached___equal(self, *uninteresting_mocks):
        clustering_result_mock = Mock(no_clustering=False)
        self.assertFalse(NodeFactory._infer_if_we_should_cluster_further(Mock(), clustering_result_mock, 5, 5))

    def test___infer_if_we_should_cluster_further___max_nesting_reached___larger(self, *uninteresting_mocks):
        clustering_result_mock = Mock(no_clustering=False)
        self.assertFalse(NodeFactory._infer_if_we_should_cluster_further(Mock(), clustering_result_mock, 6, 5))

    @patch.object(NodeFactory, NodeFactory._alignment_has_issues.__name__, return_value=True)
    def test___infer_if_we_should_cluster_further___alignment_has_issues(self,
         alignment_has_issues_mock, *uninteresting_mocks):
        clustering_result_mock = Mock(no_clustering=False)
        alignment_mock = Mock()
        self.assertFalse(NodeFactory._infer_if_we_should_cluster_further(alignment_mock, clustering_result_mock, 3, 5))
        alignment_has_issues_mock.assert_called_once_with(alignment_mock)

    @patch.object(NodeFactory, NodeFactory._alignment_has_issues.__name__, return_value=False)
    def test___infer_if_we_should_cluster_further___no_issues___ok_to_cluster(self,
                                                                              alignment_has_issues_mock,
                                                                              *uninteresting_mocks):
        clustering_result_mock = Mock(no_clustering=False)
        alignment_mock = Mock()
        self.assertTrue(NodeFactory._infer_if_we_should_cluster_further(alignment_mock, clustering_result_mock, 3, 5))
        alignment_has_issues_mock.assert_called_once_with(alignment_mock)

    def test___get_subalignments_by_clustering(self, *uninteresting_mocks):
        alignment = make_alignment(
            ["AAAT", "C--C", "AATT", "GNGG", "CCCC", "TTTT", "AAAA"], ["s0", "s1", "s2", "s3", "s4", "s5", "s6"]
        )
        clustering_result = ClusteringResult([["s5", "s0", "s1", "s6"], ["s3"], ["s2", "s4"]])

        expected = [
            make_alignment(
                ["AAAT", "C--C", "TTTT", "AAAA"], ["s0", "s1", "s5", "s6"]),
            make_alignment(
                ["GNGG"], ["s3"]),
            make_alignment(
                ["AATT", "CCCC"],
                ["s2", "s4"]),
        ]
        actual = NodeFactory._get_subalignments_by_clustering(alignment, clustering_result)

        self.assertEqual(3, len(actual))
        for i in range(3):
            self.assertTrue(equal_msas(expected[i], actual[i]))

    def test___get_sub_alignment_by_list_id___GivenOrderedIds_SubalignmentInSequenceOrder(self, *uninteresting_mocks):
        self.setup()
        expected = MSA([self.alignment[0], self.alignment[2]])
        actual = NodeFactory._get_sub_alignment_by_list_id(self.alignment, ["s1", "s3"])
        self.assertTrue(equal_msas(expected, actual))

    def test___get_sub_alignment_by_list_id___GivenUnorderedIds_SubalignmentStillInSequenceOrder(self, *uninteresting_mocks):
        """
        Sequences given rearranged are still output in input order
        """
        self.setup()
        expected = MSA([self.alignment[0], self.alignment[2]])
        actual = NodeFactory._get_sub_alignment_by_list_id(self.alignment, ["s3", "s1"])
        self.assertTrue(equal_msas(expected, actual))
