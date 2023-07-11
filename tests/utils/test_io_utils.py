from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

from make_prg.utils.io_utils import (
    get_temp_dir_for_multiprocess,
    load_alignment_file,
    zip_set_of_files,
)
from make_prg.utils.misc import equal_msas
from make_prg.utils.seq_utils import get_majority_consensus_from_MSA
from tests.test_helpers import are_zip_files_equal, make_alignment
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from io import StringIO


workdir = Path("tests/data/utils/io_utils/")


class Test_load_alignment_file(TestCase):
    def setUp(self):
        # Define a couple of MSAs for testing purposes
        self.msa1 = make_alignment(
            ["ATGCN", "ATGCG", "ATGCT"],
            ["seq1", "seq2", "seq3"],
        )
        self.sample_msa = make_alignment(
            ["ACGT", "ACGT", "A--C", "A--G", "C---"],
            ["seq1", "seq2", "seq3", "seq4", "seq5"],
        )

    def test___load_alignment_file___load_fasta(self):
        actual = load_alignment_file(
            workdir / "load_alignment_file/sample_msa.fa", "fasta"
        )

        self.assertTrue(equal_msas(self.sample_msa, actual))

    def test___load_alignment_file___load_fasta_gz(self):
        actual = load_alignment_file(
            workdir / "load_alignment_file/sample_msa.fa.gz", "fasta"
        )

        self.assertTrue(equal_msas(self.sample_msa, actual))

    def test___load_alignment_file___alignment_with_N(self):
        # Convert the MSA to a string in FASTA format
        msa_string = StringIO()
        AlignIO.write(self.msa1, msa_string, 'fasta')
        msa_string.seek(0)

        # Load the alignment
        alignment = load_alignment_file(msa_string, 'fasta')

        # Check the type of the returned object
        self.assertIsInstance(alignment, MultipleSeqAlignment)

        # Check that the sequences are upper case
        for record in alignment:
            self.assertEqual(record.seq, record.seq.upper())

        # Check that 'N' has been replaced by the consensus
        consensus = get_majority_consensus_from_MSA(alignment)
        for record in alignment:
            for i, nucleotide in enumerate(str(record.seq)):
                if nucleotide == 'N':
                    self.assertEqual(nucleotide, consensus[i])

    def test___load_alignment_file___invalid_format(self):
        # Convert the MSA to a string in a format not supported by AlignIO
        msa_string = StringIO()
        msa_string.write(str(self.msa1))
        msa_string.seek(0)

        # Attempt to load the alignment
        with self.assertRaises(ValueError):
            load_alignment_file(msa_string, 'invalid_format')

    def test___load_alignment_file___file_not_exists(self):
        with self.assertRaises(FileNotFoundError):
            load_alignment_file('nonexistent_file.fa', 'fasta')

    def test___load_alignment_file___empty_alignment(self):
        with self.assertRaises(ValueError):
            msa_string = StringIO()
            AlignIO.write(MultipleSeqAlignment([]), msa_string, 'fasta')
            msa_string.seek(0)
            load_alignment_file(msa_string, 'fasta')

    def test___load_alignment_file___all_Ns(self):
        msa_string = StringIO()
        AlignIO.write(make_alignment(["NNNN", "NNNN"], ["seq1", "seq2"]), msa_string, 'fasta')
        msa_string.seek(0)
        alignment = load_alignment_file(msa_string, 'fasta')
        for record in alignment:
            self.assertTrue(len(record.seq)==4 and "N" not in record.seq)
            for char in record.seq:
                self.assertIn(char, "ACGT")

    def test___load_alignment_file___large_alignment(self):
        msa_string = StringIO()
        large_alignment = make_alignment(["ACGT" * 250] * 100, [f'seq{i}' for i in range(100)])
        AlignIO.write(large_alignment, msa_string, 'fasta')
        msa_string.seek(0)
        alignment = load_alignment_file(msa_string, 'fasta')
        self.assertEqual(len(alignment), 100)


class Test_get_temp_dir_for_multiprocess(TestCase):
    @patch("os.makedirs")
    def test___get_temp_dir_for_multiprocess(self, makedirs_mock):
        root_temp_dir = Path("root_temp_dir")

        expected = root_temp_dir / "mp_temp"
        actual = get_temp_dir_for_multiprocess(root_temp_dir)

        self.assertEqual(expected, actual)
        makedirs_mock.assert_called_once_with(expected, exist_ok=True)


class Test_zip_set_of_files(TestCase):
    def test___zip_set_of_files___not_a_zip_filepath_given___raises_AssertionError(
        self,
    ):
        with self.assertRaises(AssertionError):
            zip_set_of_files(Path("file.txt"), {})

    def test___zip_set_of_files(self):
        concatenated_files_dir = Path("tests/data/utils/input_output_files/")
        zip_set_of_files(
            workdir / "zip_set_of_files/files.zip",
            {
                "f1": concatenated_files_dir / "concatenate_text_files/file1.txt",
                "file_2": concatenated_files_dir / "concatenate_text_files/file2.txt",
                "file3": concatenated_files_dir / "concatenate_text_files/file3.txt",
            },
        )

        self.assertTrue(
            are_zip_files_equal(
                workdir / "zip_set_of_files/files.zip",
                workdir / "zip_set_of_files/files.truth.zip",
            )
        )
