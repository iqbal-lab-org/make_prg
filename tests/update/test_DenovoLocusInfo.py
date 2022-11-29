from unittest import TestCase

from make_prg.update.denovo_variants import DenovoLocusInfo, DenovoVariant, UpdateData
from make_prg.update.MLPath import MLPath, MLPathError, MLPathNode


class DenovoLocusInfoTest(TestCase):
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

    def test___get_ml_path_nodes_spanning_variant___strict_insertion_event(self):
        denovo_variant = DenovoVariant(2, "", "ACGT")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])

        expected = [self.ml_path_node_1]
        actual = denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

        self.assertEqual(expected, actual)

    def test___get_ml_path_nodes_spanning_variant___variant_spanning_one_base_of_first_node(
        self,
    ):
        denovo_variant = DenovoVariant(0, "A", "ACGT")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])

        expected = [self.ml_path_node_1]
        actual = denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

        self.assertEqual(expected, actual)

    def test___get_ml_path_nodes_spanning_variant___variant_spanning_all_bases_of_first_node(
        self,
    ):
        denovo_variant = DenovoVariant(0, "ACGT", "GGGG")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])

        expected = [
            self.ml_path_node_1,
            self.ml_path_node_1,
            self.ml_path_node_1,
            self.ml_path_node_1,
        ]
        actual = denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

        self.assertEqual(expected, actual)

    def test___get_ml_path_nodes_spanning_variant___variant_spanning_all_bases_of_second_node(
        self,
    ):
        denovo_variant = DenovoVariant(4, "G", "AAA")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])

        expected = [self.ml_path_node_2]
        actual = denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

        self.assertEqual(expected, actual)

    def test___get_ml_path_nodes_spanning_variant___variant_spanning_last_base_of_third_node(
        self,
    ):
        denovo_variant = DenovoVariant(6, "A", "T")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])

        expected = [self.ml_path_node_3]
        actual = denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

        self.assertEqual(expected, actual)

    def test___get_ml_path_nodes_spanning_variant___variant_spanning_all_bases_of_third_node(
        self,
    ):
        denovo_variant = DenovoVariant(5, "AA", "T")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])

        expected = [self.ml_path_node_3, self.ml_path_node_3]
        actual = denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

        self.assertEqual(expected, actual)

    def test___get_ml_path_nodes_spanning_variant___variant_spanning_all_nodes(self):
        denovo_variant = DenovoVariant(0, "ACGTGAA", "T")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])

        expected = [
            self.ml_path_node_1,
            self.ml_path_node_1,
            self.ml_path_node_1,
            self.ml_path_node_1,
            self.ml_path_node_2,
            self.ml_path_node_3,
            self.ml_path_node_3,
        ]
        actual = denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

        self.assertEqual(expected, actual)

    def test___get_ml_path_nodes_spanning_variant___variant_spanning_all_nodes_on_border(
        self,
    ):
        denovo_variant = DenovoVariant(3, "TGA", "")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])

        expected = [self.ml_path_node_1, self.ml_path_node_2, self.ml_path_node_3]
        actual = denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

        self.assertEqual(expected, actual)

    def test___get_ml_path_nodes_spanning_variant___variant_spanning_last_base_and_an_additional_base_of_third_node___raises_MLPathError(
        self,
    ):
        denovo_variant = DenovoVariant(6, "AA", "T")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])
        with self.assertRaises(MLPathError):
            denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

    def test___get_ml_path_nodes_spanning_variant___variant_is_a_strict_insertion_at_the_last_insertion_position(
        self,
    ):
        denovo_variant = DenovoVariant(7, "", "AAA")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])

        expected = [self.ml_path_node_3]
        actual = denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

        self.assertEqual(expected, actual)

    def test___get_ml_path_nodes_spanning_variant___variant_is_a_strict_insertion_one_base_after_the_last_insertion_position___raises_MLPathError(
        self,
    ):
        denovo_variant = DenovoVariant(8, "", "AAA")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])
        with self.assertRaises(MLPathError):
            denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

    def test___get_ml_path_nodes_spanning_variant___variant_is_a_non_strict_insertion_at_the_last_insertion_position___raises_MLPathError(
        self,
    ):
        denovo_variant = DenovoVariant(7, "A", "AAA")
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])
        with self.assertRaises(MLPathError):
            denovo_locus_info._get_ml_path_nodes_spanning_variant(denovo_variant)

    def test___get_update_data___empty_variants___returns_empty_update_data(self):
        denovo_locus_info = DenovoLocusInfo("sample", "locus", self.ml_path, [])
        expected = []
        actual = denovo_locus_info.get_update_data()
        self.assertEqual(expected, actual)

    def test___get_update_data___several_variants_in_first_node(self):
        denovo_variant = DenovoVariant(0, "ACGT", "TT")
        denovo_locus_info = DenovoLocusInfo(
            "sample", "locus", self.ml_path, [denovo_variant]
        )

        expected = [UpdateData(self.ml_path_node_1.key, self.ml_path, "TT")]
        actual = denovo_locus_info.get_update_data()

        self.assertEqual(expected, actual)

    def test___get_update_data___indels_in_last_node(self):
        denovo_variant = DenovoVariant(5, "AA", "ATTTTT")
        denovo_locus_info = DenovoLocusInfo(
            "sample", "locus", self.ml_path, [denovo_variant]
        )

        expected = [UpdateData(self.ml_path_node_3.key, self.ml_path, "ATTTTT")]
        actual = denovo_locus_info.get_update_data()

        self.assertEqual(expected, actual)

    def test___get_update_data___one_SNP_in_each_node(self):
        denovo_variant = DenovoVariant(0, "ACGTGAA", "AGGTCAT")
        denovo_locus_info = DenovoLocusInfo(
            "sample", "locus", self.ml_path, [denovo_variant]
        )

        expected = [
            UpdateData(self.ml_path_node_1.key, self.ml_path, "AGGT"),
            UpdateData(self.ml_path_node_2.key, self.ml_path, "C"),
            UpdateData(self.ml_path_node_3.key, self.ml_path, "AT"),
        ]
        actual = denovo_locus_info.get_update_data()

        self.assertEqual(expected, actual)

    def test___get_update_data___one_SNP_in_each_node_border_only(self):
        denovo_variant = DenovoVariant(3, "TGA", "CCC")
        denovo_locus_info = DenovoLocusInfo(
            "sample", "locus", self.ml_path, [denovo_variant]
        )

        expected = [
            UpdateData(self.ml_path_node_1.key, self.ml_path, "ACGC"),
            UpdateData(self.ml_path_node_2.key, self.ml_path, "C"),
            UpdateData(self.ml_path_node_3.key, self.ml_path, "CA"),
        ]
        actual = denovo_locus_info.get_update_data()

        self.assertEqual(expected, actual)
