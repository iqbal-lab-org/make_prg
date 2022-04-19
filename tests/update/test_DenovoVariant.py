from unittest import TestCase
from unittest.mock import Mock
from make_prg.update.denovo_variants import DenovoVariant, DenovoError
from make_prg.update.MLPath import MLPathNode

class DenovoVariantTest(TestCase):
    def setUp(self) -> None:
        self.sample_denovo_variant = DenovoVariant(5, "AACC", "GT")

    def test___ref_not_composed_of_ACGT_only___DenovoError_raised(self):
        with self.assertRaises(DenovoError):
            DenovoVariant(0, ref="ACGTB", alt="")

    def test___alt_not_composed_of_ACGT_only___DenovoError_raised(self):
        with self.assertRaises(DenovoError):
            DenovoVariant(0, ref="", alt="ACGTB")

    def test___ref_and_alt_are_identical___DenovoError_raised(self):
        seq="ACGT"
        with self.assertRaises(DenovoError):
            DenovoVariant(0, ref=seq, alt=seq)

    def test___negative_index_for_variant_pos___DenovoError_raised(self):
        with self.assertRaises(DenovoError):
            DenovoVariant(-1, "A", "C")

    def test___constructor___variant_correctly_built(self):
        ml_path_nodes_it_goes_through = [Mock(), Mock(), Mock(), Mock()]
        denovo_variant = DenovoVariant(5, "AACC", "GT", ml_path_nodes_it_goes_through)
        self.assertEqual(5, denovo_variant.start_index_in_linear_path)
        self.assertEqual(9, denovo_variant.end_index_in_linear_path)
        self.assertEqual("AACC", denovo_variant.ref)
        self.assertEqual("GT", denovo_variant.alt)
        self.assertEqual(ml_path_nodes_it_goes_through, denovo_variant.ml_path_nodes_it_goes_through)

    def test___equality___same_variant(self):
        ml_path_nodes_it_goes_through = [Mock()]
        denovo_variant_1 = DenovoVariant(1, "A", "C", ml_path_nodes_it_goes_through)
        denovo_variant_2 = DenovoVariant(1, "A", "C", ml_path_nodes_it_goes_through)
        self.assertEqual(denovo_variant_1, denovo_variant_2)

    def test___equality___different_variants(self):
        ml_path_nodes_it_goes_through_1 = [Mock()]
        denovo_variant_1 = DenovoVariant(1, "A", "C", ml_path_nodes_it_goes_through_1)

        denovo_variant_2 = DenovoVariant(2, "A", "C", ml_path_nodes_it_goes_through_1)
        self.assertNotEqual(denovo_variant_1, denovo_variant_2)

        denovo_variant_2 = DenovoVariant(1, "T", "C", ml_path_nodes_it_goes_through_1)
        self.assertNotEqual(denovo_variant_1, denovo_variant_2)

        denovo_variant_2 = DenovoVariant(1, "A", "G", ml_path_nodes_it_goes_through_1)
        self.assertNotEqual(denovo_variant_1, denovo_variant_2)

        ml_path_nodes_it_goes_through_2 = [Mock()]
        denovo_variant_2 = DenovoVariant(1, "A", "C", ml_path_nodes_it_goes_through_2)
        self.assertNotEqual(denovo_variant_1, denovo_variant_2)

    def test___equality___comparing_variant_with_other_type(self):
        denovo_variant_1 = DenovoVariant(1, "A", "C")
        other = "a string"
        self.assertNotEqual(denovo_variant_1, other)

    def test___set_ml_path_nodes_it_goes_through___ml_path_nodes_it_goes_through_is_None(self):
        self.sample_denovo_variant.set_ml_path_nodes_it_goes_through(None)
        self.assertIsNone(self.sample_denovo_variant.ml_path_nodes_it_goes_through)

    def test___set_ml_path_nodes_it_goes_through___empty_path_error_out(self):
        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.set_ml_path_nodes_it_goes_through([])

    def test___set_ml_path_nodes_it_goes_through___one_node_less_than_expected_in_path_error_out(self):
        dummy_node = MLPathNode(key=(0, 1), sequence="A")
        ml_path_nodes_it_goes_through = [dummy_node, dummy_node, dummy_node]
        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.set_ml_path_nodes_it_goes_through(ml_path_nodes_it_goes_through)

    def test___set_ml_path_nodes_it_goes_through___one_node_more_than_expected_in_path_error_out(self):
        dummy_node = MLPathNode(key=(0, 1), sequence="A")
        ml_path_nodes_it_goes_through = [dummy_node, dummy_node, dummy_node, dummy_node, dummy_node]
        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.set_ml_path_nodes_it_goes_through(ml_path_nodes_it_goes_through)

    def test___set_ml_path_nodes_it_goes_through___each_base_is_covered_by_one_node(self):
        dummy_node = MLPathNode(key=(0, 1), sequence="A")
        ml_path_nodes_it_goes_through = [dummy_node, dummy_node, dummy_node, dummy_node]

        self.sample_denovo_variant.set_ml_path_nodes_it_goes_through(ml_path_nodes_it_goes_through)

        self.assertEqual(ml_path_nodes_it_goes_through, self.sample_denovo_variant.ml_path_nodes_it_goes_through)

    def test___set_ml_path_nodes_it_goes_through___strict_insertion___no_ml_path_nodes_going_through_it_error_out(self):
        denovo_variant = DenovoVariant(5, "", "ACGT")
        with self.assertRaises(AssertionError):
            denovo_variant.set_ml_path_nodes_it_goes_through([])

    def test___set_ml_path_nodes_it_goes_through___strict_insertion___two_ml_path_nodes_going_through_it_error_out(self):
        denovo_variant = DenovoVariant(5, "", "ACGT")
        dummy_node = MLPathNode(key=(0, 1), sequence="A")
        ml_path_nodes_it_goes_through = [dummy_node, dummy_node]
        with self.assertRaises(AssertionError):
            denovo_variant.set_ml_path_nodes_it_goes_through(ml_path_nodes_it_goes_through)

    def test___set_ml_path_nodes_it_goes_through___strict_insertion_is_ok(self):
        denovo_variant = DenovoVariant(5, "", "ACGT")
        dummy_node = MLPathNode(key=(0, 1), sequence="A")
        ml_path_nodes_it_goes_through = [dummy_node]

        denovo_variant.set_ml_path_nodes_it_goes_through(ml_path_nodes_it_goes_through)

        self.assertEqual(ml_path_nodes_it_goes_through, denovo_variant.ml_path_nodes_it_goes_through)

    def test___get_mutated_sequence___ml_path_is_none___AssertionError_raised(self):
        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.get_mutated_sequence()

    def test___get_mutated_sequence___ml_path_is_empty___AssertionError_raised(self):
        self.sample_denovo_variant.ml_path_nodes_it_goes_through = []
        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.get_mutated_sequence()

    def test___get_mutated_sequence___ml_path_has_two_distinct_nodes___AssertionError_raised(self):
        self.sample_denovo_variant.ml_path_nodes_it_goes_through = [Mock(), Mock()]
        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.get_mutated_sequence()

    def test___get_mutated_sequence___bad_start_index_before_node_start_AssertionError_raised(self):
        ml_path_node = MLPathNode(key=(6, 9), sequence="ACG")
        ml_path_node.start_index_in_linear_path = 6
        ml_path_node.end_index_in_linear_path = 9
        self.sample_denovo_variant.set_ml_path_nodes_it_goes_through([ml_path_node]*4)

        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.get_mutated_sequence()

    def test___get_mutated_sequence___bad_end_index_after_node_end_AssertionError_raised(self):
        ml_path_node = MLPathNode(key=(5, 8), sequence="ACG")
        ml_path_node.start_index_in_linear_path = 5
        ml_path_node.end_index_in_linear_path = 8
        self.sample_denovo_variant.set_ml_path_nodes_it_goes_through([ml_path_node]*4)

        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.get_mutated_sequence()

    def test___get_mutated_sequence___node_seq_not_consistent_with_ref_AssertionError_raised(self):
        ml_path_node = MLPathNode(key=(3, 13), sequence="GTAAGCGTCA")
        ml_path_node.start_index_in_linear_path = 3
        ml_path_node.end_index_in_linear_path = 13
        self.sample_denovo_variant.set_ml_path_nodes_it_goes_through([ml_path_node]*4)

        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.get_mutated_sequence()

    def test___get_mutated_sequence(self):
        ml_path_node = MLPathNode(key=(3, 13), sequence="GTAACCGTCA")
        ml_path_node.start_index_in_linear_path = 3
        ml_path_node.end_index_in_linear_path = 13
        self.sample_denovo_variant.set_ml_path_nodes_it_goes_through([ml_path_node]*4)

        expected = "GTGTGTCA"
        actual = self.sample_denovo_variant.get_mutated_sequence()

        self.assertEqual(expected, actual)

    def test___split_variant___ml_path_nodes_it_goes_through_is_None___raises_Assertion_Error(self):
        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.split_variant()

    def test___split_variant___variant_goes_through_single_node___no_split(self):
        dummy_node = MLPathNode(key=(0, 1), sequence="A")
        # all 4 bases goes through this single node
        self.sample_denovo_variant.ml_path_nodes_it_goes_through = [dummy_node, dummy_node, dummy_node, dummy_node]

        expected = [self.sample_denovo_variant]
        actual = self.sample_denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_var_goes_through_two_nodes_split_in_middle(self):
        denovo_variant = DenovoVariant(5, "ACGT", "CG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AC", "C", [dummy_node_1, dummy_node_1]),
                    DenovoVariant(7, "GT", "G", [dummy_node_2, dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_variant_goes_through_two_nodes_split_in_middle_with_mismatch(self):
        denovo_variant = DenovoVariant(5, "ACGTA", "CAT")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "ACG", "CA", [dummy_node_1, dummy_node_1, dummy_node_1]),
                    DenovoVariant(8, "TA", "T", [dummy_node_2, dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_var_goes_through_two_nodes_split_in_start(self):
        denovo_variant = DenovoVariant(5, "ACGT", "GCG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AC", "GC", [dummy_node_1, dummy_node_1]),
                    DenovoVariant(7, "GT", "G", [dummy_node_2, dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_var_goes_through_two_nodes_split_in_start___one_variant_removed(self):
        denovo_variant = DenovoVariant(5, "ACGT", "AC")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(7, "GT", "", [dummy_node_2, dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_var_goes_through_two_nodes_split_in_end(self):
        denovo_variant = DenovoVariant(5, "ACGT", "CCT")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AC", "C", [dummy_node_1, dummy_node_1]),
                    DenovoVariant(7, "GT", "CT", [dummy_node_2, dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_var_goes_through_two_nodes_split_in_end___one_variant_removed(self):
        denovo_variant = DenovoVariant(5, "ACGT", "GT")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AC", "", [dummy_node_1, dummy_node_1])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___phased_snps_var_goes_through_two_nodes_split_in_middle(self):
        denovo_variant = DenovoVariant(5, "AGCT", "TGCA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AG", "TG", [dummy_node_1, dummy_node_1]),
                    DenovoVariant(7, "CT", "CA", [dummy_node_2, dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___phased_snps_var_goes_through_two_nodes_split_in_start(self):
        denovo_variant = DenovoVariant(5, "ACGT", "TCTT")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "A", "T", [dummy_node_1]),
                    DenovoVariant(6, "CGT", "CTT", [dummy_node_2, dummy_node_2, dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___phased_snps_var_goes_through_two_nodes_split_in_end(self):
        denovo_variant = DenovoVariant(5, "ACGT", "AGGA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "ACG", "AGG", [dummy_node_1, dummy_node_1, dummy_node_1]),
                    DenovoVariant(8, "T", "A", [dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___insertion_var_goes_through_two_nodes_split_in_middle(self):
        denovo_variant = DenovoVariant(5, "GC", "TGCA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "G", "TG", [dummy_node_1]),
                    DenovoVariant(6, "C", "CA", [dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___insertion_var_goes_through_two_nodes_split_in_start___one_variant_removed(self):
        denovo_variant = DenovoVariant(5, "TG", "CATG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "T", "CAT", [dummy_node_1])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___insertion_var_goes_through_two_nodes_split_in_end___one_variant_removed(self):
        denovo_variant = DenovoVariant(5, "TG", "TGCA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(6, "G", "GCA", [dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)


    def test___split_variant___long_insertion_var_before_and_after_goes_through_two_nodes(self):
        denovo_variant = DenovoVariant(5, "TG", "AAAATGCCCC")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "T", "AAAAT", [dummy_node_1]),
                    DenovoVariant(6, "G", "GCCCC", [dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___long_insertion_var_between_goes_through_two_nodes(self):
        denovo_variant = DenovoVariant(5, "TG", "TAAAAG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(6, "G", "AAAAG", [dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___long_insertion_var_before_after_and_between_goes_through_two_nodes(self):
        denovo_variant = DenovoVariant(5, "TG", "AAAATCCCCGAAAA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "T", "AAAAT", [dummy_node_1]),
                    DenovoVariant(6, "G", "CCCCGAAAA", [dummy_node_2])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___del_snp_ins_goes_through_three_nodes(self):
        denovo_variant = DenovoVariant(5, "AAAAAACGGGGGCTTTT", "AATGGGGGATTCCTTCC")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        dummy_node_3 = MLPathNode(key=(2, 3), sequence="T")
        # all 4 bases goes through this single node
        denovo_variant.ml_path_nodes_it_goes_through = [
            dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_1,
            dummy_node_2, dummy_node_2, dummy_node_2, dummy_node_2, dummy_node_2, dummy_node_2, dummy_node_2,
            dummy_node_3, dummy_node_3, dummy_node_3, dummy_node_3]

        expected = [DenovoVariant(5, "AAAAAA", "AA", [dummy_node_1]*6),
                    DenovoVariant(11, "CGGGGGC", "TGGGGGA", [dummy_node_2]*7),
                    DenovoVariant(18, "TTTT", "TTCCTTCC", [dummy_node_3]*4)]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___single_SNPs_through_several_nodes(self):
        denovo_variant = DenovoVariant(5, "ACGTTCGT", "TCCAACGG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        dummy_node_3 = MLPathNode(key=(2, 3), sequence="T")
        dummy_node_4 = MLPathNode(key=(3, 4), sequence="A")
        dummy_node_5 = MLPathNode(key=(4, 5), sequence="C")
        dummy_node_6 = MLPathNode(key=(5, 6), sequence="T")
        dummy_node_7 = MLPathNode(key=(6, 7), sequence="A")
        dummy_node_8 = MLPathNode(key=(7, 8), sequence="C")
        denovo_variant.ml_path_nodes_it_goes_through = [
            dummy_node_1, dummy_node_2, dummy_node_3, dummy_node_4, dummy_node_5, dummy_node_6, dummy_node_7,
            dummy_node_8]

        expected = [DenovoVariant(5, "A", "T", [dummy_node_1]),
                    DenovoVariant(7, "G", "C", [dummy_node_3]),
                    DenovoVariant(8, "T", "A", [dummy_node_4]),
                    DenovoVariant(9, "T", "A", [dummy_node_5]),
                    DenovoVariant(12, "T", "G", [dummy_node_8])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___dels_through_several_nodes(self):
        denovo_variant = DenovoVariant(5, "ACGTTCGT", "AGTG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        dummy_node_3 = MLPathNode(key=(2, 3), sequence="T")
        dummy_node_4 = MLPathNode(key=(3, 4), sequence="A")
        dummy_node_5 = MLPathNode(key=(4, 5), sequence="C")
        dummy_node_6 = MLPathNode(key=(5, 6), sequence="T")
        dummy_node_7 = MLPathNode(key=(6, 7), sequence="A")
        dummy_node_8 = MLPathNode(key=(7, 8), sequence="C")
        denovo_variant.ml_path_nodes_it_goes_through = [
            dummy_node_1, dummy_node_2, dummy_node_3, dummy_node_4, dummy_node_5, dummy_node_6, dummy_node_7,
            dummy_node_8]

        expected = [DenovoVariant(6, "C", "", [dummy_node_2]),
                    DenovoVariant(9, "T", "", [dummy_node_5]),
                    DenovoVariant(10, "C", "", [dummy_node_6]),
                    DenovoVariant(12, "T", "", [dummy_node_8])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___strict_insertion(self):
        denovo_variant = DenovoVariant(5, "", "ACGT")
        dummy_node = MLPathNode(key=(0, 1), sequence="A")
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node]

        expected = [denovo_variant]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___split_variant___middle_ml_path_node_is_removed(self):
        denovo_variant = DenovoVariant(156, "TAG", "CTA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        dummy_node_3 = MLPathNode(key=(2, 3), sequence="T")
        denovo_variant.ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2, dummy_node_3]

        expected = [DenovoVariant(156, "T", "CT", [dummy_node_1]),
                    DenovoVariant(158, "G", "", [dummy_node_3])]
        actual = denovo_variant.split_variant()

        self.assertEqual(expected, actual)

    def test___is_strict_insertion_event___true(self):
        denovo_variant = DenovoVariant(5, "", "AGTG")
        self.assertTrue(denovo_variant.is_strict_insertion_event())

    def test___is_strict_insertion_event___false(self):
        denovo_variant = DenovoVariant(5, "A", "AGTG")
        self.assertFalse(denovo_variant.is_strict_insertion_event())

    def test___str(self):
        expected = "[5:9]:'AACC'->'GT'"
        actual = str(self.sample_denovo_variant)
        self.assertEqual(expected, actual)

    def test___repr(self):
        expected = "DenovoVariant(start_index_in_linear_path=5, ref=\"AACC\", alt=\"GT\")"
        actual = repr(self.sample_denovo_variant)
        self.assertEqual(expected, actual)
