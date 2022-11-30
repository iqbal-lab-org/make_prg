from argparse import Namespace
from pathlib import Path
from unittest import TestCase

import pytest

from make_prg.subcommands import output_type, update
from tests.test_helpers import are_dir_trees_equal, remove_dir_if_exists

data_dir = Path("tests/integration_tests/data")


# NOTE: for these tests we need mafft in $PATH
@pytest.mark.forked
class Test_Update_Integration_Full_Builds(TestCase):
    def prepare_options(self, test_name: str, update_DS: Path):
        output_folder = data_dir / "output_update" / test_name
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / test_name)

        options = Namespace(
            denovo_paths=str(data_dir / test_name / "denovo_paths.txt"),
            update_DS=update_DS,
            output_prefix=output_prefix,
            long_deletion_threshold=1000000,
            mafft="mafft",
            log=None,
            output_type=output_type.OutputType("a"),
            force=False,
            output_graphs=False,
            threads=1,
            verbose=False,
        )

        return options

    def test___update___match_update_simple(self):
        options = self.prepare_options(
            test_name="match_update_simple",
            update_DS=data_dir / "truth_output/match/match.update_DS.zip",
        )
        update.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/match_update_simple",
                data_dir / "output_update/match_update_simple",
            )
        )

    def test___update___match_update_simple_with_long_deletion(self):
        options = self.prepare_options(
            test_name="match_update_simple_with_long_deletion",
            update_DS=data_dir / "truth_output/match/match.update_DS.zip",
        )
        options.long_deletion_threshold = 10
        update.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/match_update_simple_with_long_deletion",
                data_dir / "output_update/match_update_simple_with_long_deletion",
            )
        )

    def test___update___match_update_complex___several_samples_and_variants(self):
        options = self.prepare_options(
            test_name="match_update_complex",
            update_DS=data_dir / "truth_output/match/match.update_DS.zip",
        )
        update.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/match_update_complex",
                data_dir / "output_update/match_update_complex",
            )
        )

    def test___update___match_nonmatch_match___locus_with_different_ML_path_for_each_sample(
        self,
    ):
        options = self.prepare_options(
            test_name="match.nonmatch.match_update",
            update_DS=data_dir
            / "truth_output/match.nonmatch.match/match.nonmatch.match.update_DS.zip",
        )
        update.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/match.nonmatch.match_update",
                data_dir / "output_update/match.nonmatch.match_update",
            )
        )

    def test___update___nested_snp_backgrounds(self):
        options = self.prepare_options(
            test_name="nested_snps_seq_backgrounds_update",
            update_DS=data_dir
            / "truth_output/nested_snps_seq_backgrounds/nested_snps_seq_backgrounds.update_DS.zip",
        )
        update.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/nested_snps_seq_backgrounds_update",
                data_dir / "output_update/nested_snps_seq_backgrounds_update",
            )
        )

    def test___update___strict_insertions_and_deletions(self):
        """
        Strict insertions happen "before" the specified position, e.g.:
        This is the variant:
        10		A
        This is the sequence
        ACGTGTTTT G TAACTGTG...
                  ^ this is position 10 (1-based), the insertion point
        Where we will end up inserting will be:
        ACGTGTTTT AG TAACTGTG...
                  ^ insertion
        """
        options = self.prepare_options(
            test_name="strict_insertions_and_deletions_update",
            update_DS=data_dir / "truth_output/match/match.update_DS.zip",
        )
        update.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/strict_insertions_and_deletions_update",
                data_dir / "output_update/strict_insertions_and_deletions_update",
            )
        )

    def test___update___output_files_already_exist(self):
        output_folder = data_dir / "output_update" / "match_update_simple"
        output_prefix = str(output_folder / "match_update_simple")

        options = Namespace(
            denovo_paths=str(data_dir / "match_update_simple" / "denovo_paths.txt"),
            update_DS=data_dir / "truth_output/match/match.update_DS.zip",
            output_prefix=output_prefix,
            mafft="mafft",
            log=None,
            output_type=output_type.OutputType("a"),
            force=False,
            output_graphs=False,
            threads=1,
            verbose=False,
        )

        with self.assertRaises(RuntimeError):
            update.run(options)

    def test___update___output_files_already_exist___force_overwrite(self):
        output_folder = data_dir / "output_update" / "match_update_simple"
        assert (
            output_folder.exists()
        ), f"Error, {output_folder} should exist already when running this test"
        output_prefix = str(output_folder / "match_update_simple")

        options = Namespace(
            denovo_paths=str(data_dir / "match_update_simple" / "denovo_paths.txt"),
            update_DS=data_dir / "truth_output/match/match.update_DS.zip",
            output_prefix=output_prefix,
            long_deletion_threshold=1000000,
            mafft="mafft",
            log=None,
            output_type=output_type.OutputType("a"),
            force=True,
            output_graphs=False,
            threads=1,
            verbose=False,
        )

        update.run(options)

        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/match_update_simple",
                data_dir / "output_update/match_update_simple",
            )
        )

    def test___update___match_update_simple___produce_only_PRGs(self):
        options = self.prepare_options(
            test_name="match_update_simple_prg_only",
            update_DS=data_dir / "truth_output/match/match.update_DS.zip",
        )
        options.output_type = output_type.OutputType("p")
        update.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/match_update_simple_prg_only",
                data_dir / "output_update/match_update_simple_prg_only",
            )
        )

    def test___update___match_update_simple___produce_only_GFAs(self):
        options = self.prepare_options(
            test_name="match_update_simple_gfa_only",
            update_DS=data_dir / "truth_output/match/match.update_DS.zip",
        )
        options.output_type = output_type.OutputType("g")
        update.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/match_update_simple_gfa_only",
                data_dir / "output_update/match_update_simple_gfa_only",
            )
        )

    def test___update___match_update_simple___produce_only_binary_PRGs(self):
        options = self.prepare_options(
            test_name="match_update_simple_bin_only",
            update_DS=data_dir / "truth_output/match/match.update_DS.zip",
        )
        options.output_type = output_type.OutputType("b")
        update.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/match_update_simple_bin_only",
                data_dir / "output_update/match_update_simple_bin_only",
            )
        )

    def test___update_sample_example(self):
        options = self.prepare_options(
            test_name="sample_example_update",
            update_DS=data_dir
            / "truth_output/sample_example/sample_example.update_DS.zip",
        )
        update.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output_update/sample_example_update",
                data_dir / "output_update/sample_example_update",
            )
        )
