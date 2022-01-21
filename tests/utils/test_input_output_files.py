from unittest import TestCase
from pathlib import Path
from make_prg.utils.input_output_files import InputOutputFiles, InputOutputFilesFromMSA
import filecmp

workdir = Path("tests/data/utils/input_output_files/")


class Test_concatenate_text_files(TestCase):
    def test___concatenate_text_files(self):
        input_filepaths = [workdir / f"concatenate_text_files/file{file_number}.txt" for file_number in range(1, 4)]
        output_filepath = workdir / f"concatenate_text_files/concatenated.txt"

        InputOutputFiles.concatenate_text_files(input_filepaths, output_filepath)

        self.assertTrue(filecmp.cmp(output_filepath, workdir / f"concatenate_text_files/concatenated.truth.txt"))

class Test_remove_known_input_extensions(TestCase):
    def test___not_a_fasta_nor_gz_file___no_removes(self):
        filepath = "fake_dir/test.txt"

        expected = "fake_dir/test.txt"
        actual = InputOutputFilesFromMSA.remove_known_input_extensions(filepath)

        self.assertEqual(expected, actual)

    def test___fa_file___one_remove(self):
        filepath = "fake_dir/test.fa"

        expected = "fake_dir/test"
        actual = InputOutputFilesFromMSA.remove_known_input_extensions(filepath)

        self.assertEqual(expected, actual)

    def test___fasta_file___one_remove(self):
        filepath = "fake_dir/test.fasta"

        expected = "fake_dir/test"
        actual = InputOutputFilesFromMSA.remove_known_input_extensions(filepath)

        self.assertEqual(expected, actual)

    def test___gz_file_only___no_removes(self):
        filepath = "fake_dir/test.gz"

        expected = "fake_dir/test.gz"
        actual = InputOutputFilesFromMSA.remove_known_input_extensions(filepath)

        self.assertEqual(expected, actual)

    def test___fa_gz_file___two_removes(self):
        filepath = "fake_dir/test.fa.gz"

        expected = "fake_dir/test"
        actual = InputOutputFilesFromMSA.remove_known_input_extensions(filepath)

        self.assertEqual(expected, actual)

    def test___fasta_gz_file___two_removes(self):
        filepath = "fake_dir/test.fasta.gz"

        expected = "fake_dir/test"
        actual = InputOutputFilesFromMSA.remove_known_input_extensions(filepath)

        self.assertEqual(expected, actual)

    def test___fa_gz_file_with_dot___two_removes(self):
        filepath = "fake_dir/match.nonmatch.fa.gz"

        expected = "fake_dir/match.nonmatch"
        actual = InputOutputFilesFromMSA.remove_known_input_extensions(filepath)

        self.assertEqual(expected, actual)
