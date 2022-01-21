from unittest import TestCase
from make_prg.update.denovo_variants import DenovoVariantsDB, DenovoVariant, UpdateData
from make_prg.update.MLPath import MLPathNode, MLPathError, MLPath
from io import StringIO
from pathlib import Path


class DenovoVariantsDBTest(TestCase):
    def test___read_nb_of_samples(self):
        handler = StringIO("2 samples\n")
        expected = 2
        actual = DenovoVariantsDB._read_nb_of_samples(handler)
        self.assertEqual(expected, actual)

    def test___read_sample(self):
        handler = StringIO("Sample toy_sample_1\n")
        expected = "toy_sample_1"
        actual = DenovoVariantsDB._read_sample(handler)
        self.assertEqual(expected, actual)

    def test___read_nb_of_loci_in_sample(self):
        handler = StringIO("23 loci with denovo variants\n")
        expected = 23
        actual = DenovoVariantsDB._read_nb_of_loci_in_sample(handler)
        self.assertEqual(expected, actual)

    def test___read_locus(self):
        handler = StringIO("GC00010897\n")
        expected = "GC00010897"
        actual = DenovoVariantsDB._read_locus(handler)
        self.assertEqual(expected, actual)

    def test___read_nb_of_nodes_in_ml_path(self):
        handler = StringIO("9 nodes \n")
        expected = 9
        actual = DenovoVariantsDB._read_nb_of_nodes_in_ml_path(handler)
        self.assertEqual(expected, actual)

    def test___read_MLPathNode(self):
        handler = StringIO("(3 [121, 171) GAAATGCCCACGACCGGGCTGGATGAGCTGACAGAGGCCGAACGCGAGAC)\n")
        expected = MLPathNode(key=(121, 171), sequence="GAAATGCCCACGACCGGGCTGGATGAGCTGACAGAGGCCGAACGCGAGAC")
        actual = DenovoVariantsDB._read_MLPathNode(handler)
        self.assertEqual(expected, actual)

    def test___read_MLPathNode___badly_formed_line(self):
        handler = StringIO("(3 [121, 171))\n")
        with self.assertRaises(AssertionError):
            DenovoVariantsDB._read_MLPathNode(handler)

    def test___read_ml_path___empty_path___raises_MLPathError(self):
        handler = StringIO("0 nodes \n")
        with self.assertRaises(MLPathError):
            DenovoVariantsDB._read_ml_path(handler)

    def test___read_ml_path___empty_path_with_1_empty_node___raises_MLPathError(self):
        handler = StringIO("1 nodes \n"
                           "(0 [1, 1) )\n")
        with self.assertRaises(MLPathError):
            DenovoVariantsDB._read_ml_path(handler)

    def test___read_ml_path___path_with_1_empty_node_1_non_empty_node(self):
        handler = StringIO("2 nodes \n"
                           "(0 [1, 1) )\n"
                           "(1 [3, 10) ACGTACG)\n")

        expected = MLPath([MLPathNode(key=(3, 10), sequence="ACGTACG")])
        actual = DenovoVariantsDB._read_ml_path(handler)

        self.assertEqual(expected, actual)

    def test___read_ml_path___empty_path_with_2_empty_nodes_2_non_empty_node(self):
        handler = StringIO("4 nodes \n"
                           "(0 [1, 1) )\n"
                           "(1 [3, 10) ACGTACG)\n"
                           "(2 [15, 20) GGCCA)\n"
                           "(3 [30, 30) )\n")

        expected = MLPath([
            MLPathNode(key=(3, 10), sequence="ACGTACG"),
            MLPathNode(key=(15, 20), sequence="GGCCA"),
        ])
        actual = DenovoVariantsDB._read_ml_path(handler)

        self.assertEqual(expected, actual)

    def test___read_nb_of_variants(self):
        handler = StringIO("42 denovo variants for this locus\n")

        expected = 42
        actual = DenovoVariantsDB._read_nb_of_variants(handler)

        self.assertEqual(expected, actual)

    def test___read_DenovoVariant___simple_SNP(self):
        handler = StringIO("49\tA\tG\n")

        expected = DenovoVariant(48, "A", "G")
        actual = DenovoVariantsDB._read_DenovoVariant(handler)

        self.assertEqual(expected, actual)

    def test___read_DenovoVariant___simple_insertion(self):
        handler = StringIO("49\t\tG\n")

        expected = DenovoVariant(48, "", "G")
        actual = DenovoVariantsDB._read_DenovoVariant(handler)

        self.assertEqual(expected, actual)

    def test___read_DenovoVariant___simple_deletion(self):
        handler = StringIO("49\tA\t\n")

        expected = DenovoVariant(48, "A", "")
        actual = DenovoVariantsDB._read_DenovoVariant(handler)

        self.assertEqual(expected, actual)

    def test___read_DenovoVariant___larger_variation(self):
        handler = StringIO("49\tACTGACTG\tGGAGCT\n")

        expected = DenovoVariant(48, "ACTGACTG", "GGAGCT")
        actual = DenovoVariantsDB._read_DenovoVariant(handler)

        self.assertEqual(expected, actual)

    def test___read_variants___no_variants(self):
        handler = StringIO("0 denovo variants for this locus\n")

        expected = []
        actual = DenovoVariantsDB._read_variants(handler)

        self.assertEqual(expected, actual)

    def test___read_variants___one_variant(self):
        handler = StringIO("1 denovo variants for this locus\n"
                           "49\tA\tG\n")

        expected = [DenovoVariant(48, "A", "G")]
        actual = DenovoVariantsDB._read_variants(handler)

        self.assertEqual(expected, actual)

    def test___read_variants___three_variants(self):
        handler = StringIO("3 denovo variants for this locus\n"
                           "49\tA\tG\n"
                           "314\tGT\tAC\n"
                           "500\t\tA\n")

        expected = [DenovoVariant(48, "A", "G"), DenovoVariant(313, "GT", "AC"), DenovoVariant(499, "", "A")]
        actual = DenovoVariantsDB._read_variants(handler)

        self.assertEqual(expected, actual)

    def test___big_bang(self):
        denovo_paths_filepath = "tests/data/update/denovo_paths.txt"

        denovo_variants_DB = DenovoVariantsDB(denovo_paths_filepath)

        self.assertEqual(Path(denovo_paths_filepath), denovo_variants_DB.filepath)
        expected_locus_name_to_update_data = {
            'GC00010897': [
                UpdateData(ml_path_node_key=(0, 110),
                           ml_path=denovo_variants_DB.locus_name_to_update_data['GC00010897'][0].ml_path,
                           new_node_sequence="ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCATCGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC"),
                UpdateData(ml_path_node_key=(374, 491),
                           ml_path=denovo_variants_DB.locus_name_to_update_data['GC00010897'][1].ml_path,
                           new_node_sequence="GTCGGCAAGGCCTTGCGCAAGGCTGGTCACGCGAAGCCCAAGGCGGTCAGAAAGGGCAAGCCGGTCGATCCGGCTGATCCCAAGGATCAAGGGGTGGGGGCACCAAAGGGGAAATGA")
            ],
            'GC00006032': [
                UpdateData(ml_path_node_key=(0, 145),
                           ml_path=denovo_variants_DB.locus_name_to_update_data['GC00006032'][0].ml_path,
                           new_node_sequence="TTGAGTAAAACAATCCCCCGCGCTTATATAAGCGCGTTGATATTTTTAGTTATTAACAAGCAACATCATGCTAATACAGACATACAAGGAGATCATCTCTCTTTGCCTGTTTTTTATTATTTCAGGAGTGTAAACACATTTTCCG")
            ],
            'Cluster_1011': [
                UpdateData(ml_path_node_key=(930, 946),
                           ml_path=denovo_variants_DB.locus_name_to_update_data['Cluster_1011'][0].ml_path,
                           new_node_sequence="TTTTTGACCATTTCCA"),

                UpdateData(ml_path_node_key=(955, 956),
                           ml_path=denovo_variants_DB.locus_name_to_update_data['Cluster_1011'][1].ml_path,
                           new_node_sequence="C")
            ]}
        self.assertEqual(expected_locus_name_to_update_data, denovo_variants_DB.locus_name_to_update_data)

    def test___big_bang___empty(self):
        denovo_paths_filepath = "tests/data/update/empty_denovo_paths.txt"

        denovo_variants_DB = DenovoVariantsDB(denovo_paths_filepath)

        self.assertEqual(Path(denovo_paths_filepath), denovo_variants_DB.filepath)
        self.assertEqual({}, denovo_variants_DB.locus_name_to_update_data)
