import os
from argparse import Namespace
from pathlib import Path
from unittest import TestCase

import pytest

from make_prg.subcommands import from_msa, output_type
from make_prg.subcommands.from_msa import EmptyMSAError
from tests.test_helpers import are_dir_trees_equal, remove_dir_if_exists

data_dir = Path("tests/integration_tests/data")


class Test_From_MSA_Integration_Full_Builds(TestCase):
    def prepare_options(self, test_name):
        input_data = str(data_dir / f"{test_name}.fa")
        output_folder = data_dir / "output" / test_name
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / test_name)

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_type=output_type.OutputType("a"),
            force=False,
            threads=1,
            verbose=False,
        )

        return options

    def test___match(self):
        options = self.prepare_options("match")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/match", data_dir / "output/match"
            )
        )

    def test___match_nonmatch(self):
        options = self.prepare_options("match.nonmatch")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/match.nonmatch",
                data_dir / "output/match.nonmatch",
            )
        )

    def test___match_nonmatch_match(self):
        options = self.prepare_options("match.nonmatch.match")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/match.nonmatch.match",
                data_dir / "output/match.nonmatch.match",
            )
        )

    def test___match_nonmatch_shortmatch(self):
        options = self.prepare_options("match.nonmatch.shortmatch")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/match.nonmatch.shortmatch",
                data_dir / "output/match.nonmatch.shortmatch",
            )
        )

    def test___match_staggereddash(self):
        options = self.prepare_options("match.staggereddash")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/match.staggereddash",
                data_dir / "output/match.staggereddash",
            )
        )

    def test___nonmatch(self):
        options = self.prepare_options("nonmatch")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/nonmatch", data_dir / "output/nonmatch"
            )
        )

    def test___nonmatch_match(self):
        options = self.prepare_options("nonmatch.match")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/nonmatch.match",
                data_dir / "output/nonmatch.match",
            )
        )

    def test___nonmatch_shortmatch(self):
        options = self.prepare_options("nonmatch.shortmatch")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/nonmatch.shortmatch",
                data_dir / "output/nonmatch.shortmatch",
            )
        )

    def test___shortmatch_nonmatch(self):
        options = self.prepare_options("shortmatch.nonmatch")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/shortmatch.nonmatch",
                data_dir / "output/shortmatch.nonmatch",
            )
        )

    def test___shortmatch_nonmatch_match(self):
        options = self.prepare_options("shortmatch.nonmatch.match")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/shortmatch.nonmatch.match",
                data_dir / "output/shortmatch.nonmatch.match",
            )
        )

    def test___contains_n(self):
        options = self.prepare_options("contains_n")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/contains_n", data_dir / "output/contains_n"
            )
        )

    def test___contains_n_and_RYKMSW(self):
        options = self.prepare_options("contains_n_and_RYKMSW")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/contains_n_and_RYKMSW",
                data_dir / "output/contains_n_and_RYKMSW",
            )
        )

    def test___contains_n_no_variants(self):
        options = self.prepare_options("contains_n_no_variants")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/contains_n_no_variants",
                data_dir / "output/contains_n_no_variants",
            )
        )

    def test___contains_RYKMSW(self):
        options = self.prepare_options("contains_RYKMSW")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/contains_RYKMSW",
                data_dir / "output/contains_RYKMSW",
            )
        )

    def test___a_column_full_of_Ns(self):
        options = self.prepare_options("a_column_full_of_Ns")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/a_column_full_of_Ns",
                data_dir / "output/a_column_full_of_Ns",
            )
        )

    def test___fails___unexpected_char_in_MSA(self):
        options = self.prepare_options("fails_2")
        from_msa.run(options)
        output_dir_is_empty = len(os.listdir(data_dir / "output/fails_2")) == 0
        self.assertTrue(output_dir_is_empty)

    def test___nested_snp_backgrounds(self):
        options = self.prepare_options("nested_snps_seq_backgrounds")
        options.min_match_length = 3
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/nested_snps_seq_backgrounds",
                data_dir / "output/nested_snps_seq_backgrounds",
            )
        )

    def test___nested_snp_backgrounds_more_seqs(self):
        options = self.prepare_options("nested_snps_seq_backgrounds_more_seqs")
        options.min_match_length = 3
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/nested_snps_seq_backgrounds_more_seqs",
                data_dir / "output/nested_snps_seq_backgrounds_more_seqs",
            )
        )

    def test___nested_snps_under_del(self):
        options = self.prepare_options("nested_snps_deletion")
        options.min_match_length = 1
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/nested_snps_deletion",
                data_dir / "output/nested_snps_deletion",
            )
        )

    def test___input_file_does_not_exist___raises_FileNotFoundError(self):
        options = self.prepare_options("unexistent_file")
        with self.assertRaises(FileNotFoundError):
            from_msa.run(options)

    # This tests the multi loci input and output
    def test___several_alignments(self):
        input_data = str(data_dir / "several")
        output_folder = data_dir / "output" / "several"
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / "several")

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            force=False,
            output_type=output_type.OutputType("a"),
            threads=1,
            verbose=False,
        )

        from_msa.run(options)

        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/several", data_dir / "output/several"
            )
        )

    def test___several_alignments___one_empty_MSA___raises_EmptyMSAError(self):
        input_data = str(data_dir / "several_empty")
        output_folder = data_dir / "output" / "several_empty"
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / "several_empty")

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            force=False,
            output_type=output_type.OutputType("a"),
            threads=1,
            verbose=False,
        )

        with self.assertRaises(EmptyMSAError):
            from_msa.run(options)

    def test___input_dir_has_no_files___raises_FileNotFoundError(self):
        unexistent_folder = data_dir / "unexistent_folder"
        os.makedirs(unexistent_folder, exist_ok=True)
        input_data = str(unexistent_folder)
        output_folder = data_dir / "output" / "unexistent_folder"
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / "unexistent_folder")

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            threads=1,
            verbose=False,
        )

        with self.assertRaises(FileNotFoundError):
            from_msa.run(options)

    def test___output_files_already_exist(self):
        input_data = str(data_dir / "match.fa")
        output_prefix = str(data_dir / "truth_output/match/match")

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_type=output_type.OutputType("a"),
            force=False,
            threads=1,
            verbose=False,
        )

        with self.assertRaises(RuntimeError):
            from_msa.run(options)

    def prepare_files_for_testing_force_overwrite(self, test_case):
        output_dir = data_dir / f"output/{test_case}"
        assert (
            output_dir.exists()
        ), f"Error, {output_dir} should exist already when running this test"

        input_data = str(data_dir / f"{test_case}.fa")
        output_prefix = str(output_dir / test_case)

        return input_data, output_prefix

    @pytest.mark.run(after="test___match")
    def test___output_files_already_exist___force_overwrite(self):
        input_data, output_prefix = self.prepare_files_for_testing_force_overwrite(
            "match"
        )

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_type=output_type.OutputType("a"),
            force=True,
            threads=1,
            verbose=False,
        )

        from_msa.run(options)

        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/match", data_dir / "output/match"
            )
        )

    def test___match___produce_only_PRGs(self):
        options = self.prepare_options("match_prg_only")
        options.output_type = output_type.OutputType("p")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/match_prg_only",
                data_dir / "output/match_prg_only",
            )
        )

    @pytest.mark.run(after="test___match___produce_only_PRGs")
    def test___match___produce_only_PRGs___output_files_already_exist(self):
        input_data, output_prefix = self.prepare_files_for_testing_force_overwrite(
            "match_prg_only"
        )

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_type=output_type.OutputType("p"),
            force=False,
            threads=1,
            verbose=False,
        )

        with self.assertRaises(RuntimeError):
            from_msa.run(options)

    def test___match___produce_only_GFAs(self):
        options = self.prepare_options("match_gfa_only")
        options.output_type = output_type.OutputType("g")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/match_gfa_only",
                data_dir / "output/match_gfa_only",
            )
        )

    @pytest.mark.run(after="test___match___produce_only_GFAs")
    def test___match___produce_only_GFAs___output_files_already_exist(self):
        input_data, output_prefix = self.prepare_files_for_testing_force_overwrite(
            "match_gfa_only"
        )

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_type=output_type.OutputType("g"),
            force=False,
            threads=1,
            verbose=False,
        )

        with self.assertRaises(RuntimeError):
            from_msa.run(options)

    def test___match___produce_only_binary_PRGs(self):
        options = self.prepare_options("match_bin_only")
        options.output_type = output_type.OutputType("b")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/match_bin_only",
                data_dir / "output/match_bin_only",
            )
        )

    @pytest.mark.run(after="test___match___produce_only_binary_PRGs")
    def test___match___produce_only_binary___output_files_already_exist(self):
        input_data, output_prefix = self.prepare_files_for_testing_force_overwrite(
            "match_bin_only"
        )

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_type=output_type.OutputType("b"),
            force=False,
            threads=1,
            verbose=False,
        )

        with self.assertRaises(RuntimeError):
            from_msa.run(options)

    def test___match___compressed_input_file(self):
        input_data = str(data_dir / "match.fa.gz")
        output_folder = data_dir / "output/match_compressed"
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / "match_compressed")

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_type=output_type.OutputType("a"),
            force=False,
            threads=1,
            verbose=False,
        )

        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/match_compressed",
                data_dir / "output/match_compressed",
            )
        )

    def test___several___compressed_and_uncompressed_input_files(self):
        input_data = str(data_dir / "several_compressed")
        output_folder = data_dir / "output" / "several_compressed"
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / "several_compressed")

        options = Namespace(
            input=input_data,
            suffix="",
            output_prefix=output_prefix,
            alignment_format="fasta",
            log=None,
            max_nesting=5,
            min_match_length=7,
            force=False,
            output_type=output_type.OutputType("a"),
            threads=1,
            verbose=False,
        )

        from_msa.run(options)

        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/several_compressed",
                data_dir / "output/several_compressed",
            )
        )

    def test___sample_example(self):
        options = self.prepare_options("sample_example")
        options.input = str(data_dir / "sample_example")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/sample_example",
                data_dir / "output/sample_example",
            )
        )

    def test___amira(self):
        options = self.prepare_options("amira_MSAs")
        options.input = str(data_dir / "amira_MSAs")
        from_msa.run(options)
        self.assertTrue(
            are_dir_trees_equal(
                data_dir / "truth_output/amira_MSAs",
                data_dir / "output/amira_MSAs",
            )
        )
