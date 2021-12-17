from unittest import TestCase
from unittest.mock import patch, Mock
from make_prg.utils.io_utils import *
from pathlib import Path
from tests.test_helpers import make_alignment, equal_msas, are_zip_files_equal
import filecmp

workdir = Path("tests/data/utils/io_utils/")


class Test_load_alignment_file(TestCase):
    def test___load_alignment_file___load_fasta(self):
        expected = make_alignment(
            ["ACGT", "ACGT", "A--C", "A--G", "C---"],
            ["seq1", "seq2", "seq3", "seq4", "seq5"])
        actual = load_alignment_file(workdir / "load_alignment_file/sample_msa.fa", "fasta")

        self.assertTrue(equal_msas(expected, actual))

    def test___load_alignment_file___load_fasta_gz(self):
        expected = make_alignment(
            ["ACGT", "ACGT", "A--C", "A--G", "C---"],
            ["seq1", "seq2", "seq3", "seq4", "seq5"])
        actual = load_alignment_file(workdir / "load_alignment_file/sample_msa.fa.gz", "fasta")

        self.assertTrue(equal_msas(expected, actual))


class Test_concatenate_text_files(TestCase):
    def test___concatenate_text_files(self):
        input_filepaths = [workdir / f"concatenate_text_files/file{file_number}.txt" for file_number in range(1, 4)]
        output_filepath = workdir / f"concatenate_text_files/concatenated.txt"

        concatenate_text_files(input_filepaths, output_filepath)

        self.assertTrue(filecmp.cmp(output_filepath, workdir / f"concatenate_text_files/concatenated.truth.txt"))


class Test_get_temp_dir_for_multiprocess(TestCase):
    @patch("os.makedirs")
    def test___get_temp_dir_for_multiprocess(self, makedirs_mock):
        root_temp_dir = Path("root_temp_dir")

        expected = root_temp_dir / "mp_temp"
        actual = get_temp_dir_for_multiprocess(root_temp_dir)

        self.assertEqual(expected, actual)
        makedirs_mock.assert_called_once_with(expected, exist_ok=True)


class Test_zip_set_of_files(TestCase):
    def test___zip_set_of_files___not_a_zip_filepath_given___raises_AssertionError(self):
        with self.assertRaises(AssertionError):
            zip_set_of_files(Path("file.txt"), {})

    def test___zip_set_of_files(self):
        zip_set_of_files(workdir / "zip_set_of_files/files.zip", {
            "f1": workdir / f"concatenate_text_files/file1.txt",
            "file_2": workdir / f"concatenate_text_files/file2.txt",
            "file3": workdir / f"concatenate_text_files/file3.txt"
        })

        self.assertTrue(are_zip_files_equal(workdir / "zip_set_of_files/files.zip",
                                            workdir / "zip_set_of_files/files.truth.zip"))