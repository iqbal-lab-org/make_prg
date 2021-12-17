from unittest import TestCase
from unittest.mock import patch, Mock, mock_open
import os
from make_prg.utils.msa_aligner import MAFFT, NotAValidExecutableError, ExecutionError
from pathlib import Path
import subprocess
import filecmp
from Bio import SeqIO
from tests.test_helpers import make_alignment


class TestMAFFT(TestCase):
    def setUp(self) -> None:
        self.executable = "fake_exec"
        self.tmpdir = Path("fake_path")

    @patch("shutil.which", return_value=None)
    def test___constructor___invalid_executable(self, shutil_which_mock):
        with self.assertRaises(NotAValidExecutableError):
            MAFFT(self.executable)

        shutil_which_mock.assert_called_once_with(self.executable, mode=os.X_OK)

    @patch("shutil.which", return_value=1)
    @patch.object(Path, Path.mkdir.__name__)
    def test___constructor___valid(self, mkdir_mock, shutil_which_mock):
        mafft = MAFFT(self.executable, self.tmpdir)

        self.assertEqual(self.executable, mafft._executable)
        self.assertEqual(self.tmpdir, mafft._tmpdir)
        mkdir_mock.assert_called_once_with(parents=True, exist_ok=True)
        shutil_which_mock.assert_called_once_with(self.executable, mode=os.X_OK)

    @patch("shutil.which", return_value=1)
    @patch.object(Path, Path.mkdir.__name__)
    @patch("subprocess.Popen")
    def test___run_aligner___ok_run(self, popen_mock, *uninteresting_mocks):
        wait_mock = Mock()
        wait_mock.wait.return_value = 0
        popen_mock.return_value = wait_mock
        args = "fake_arg_1 fake_arg_2 fake_arg_3"
        env = {"fake_key_1": "fake_value_1", "fake_key_2": "fake_value_2"}

        mafft = MAFFT(self.executable, self.tmpdir)
        mafft._run_aligner(args, env)

        popen_mock.assert_called_once_with(args, stderr=subprocess.PIPE, encoding="utf-8", shell=True, env=env)
        wait_mock.wait.assert_called_once_with()

    @patch("shutil.which", return_value=1)
    @patch.object(Path, Path.mkdir.__name__)
    @patch("subprocess.Popen")
    def test___run_aligner___fail_run(self, popen_mock, *uninteresting_mocks):
        wait_mock = Mock()
        wait_mock.wait.return_value = 1  # this mocks an exit code of 1 (i.e. fail) when invoking MAFFT
        popen_mock.return_value = wait_mock
        args = "fake_arg_1 fake_arg_2 fake_arg_3"
        env = {"fake_key_1": "fake_value_1", "fake_key_2": "fake_value_2"}

        mafft = MAFFT(self.executable, self.tmpdir)
        with self.assertRaises(ExecutionError):
            mafft._run_aligner(args, env)

        popen_mock.assert_called_once_with(args, stderr=subprocess.PIPE, encoding="utf-8", shell=True, env=env)
        wait_mock.wait.assert_called_once_with()

    @patch("shutil.which", return_value=1)
    @patch.object(Path, Path.mkdir.__name__)
    def test___create_new_sequences_file(self, *uninteresting_mocks):
        directory = Path("tests/data/utils/create_new_sequences_file")
        new_sequences = {"ACGT", "AAAA", "TTTT"}

        mafft = MAFFT(self.executable, self.tmpdir)
        mafft._create_new_sequences_file(directory, new_sequences)

        self.assertTrue(filecmp.cmp(directory/"new_sequences.fa", directory/"new_sequences.truth.fa"))

    @patch("shutil.which", return_value=1)
    @patch.object(Path, Path.mkdir.__name__)
    def test___get_aligner_name(self, *uninteresting_mocks):
        mafft = MAFFT(self.executable, self.tmpdir)

        expected = "MAFFT"
        actual = mafft.get_aligner_name()

        self.assertEqual(expected, actual)

    @patch("shutil.which", return_value=1)
    @patch.object(Path, Path.mkdir.__name__)
    @patch("shutil.rmtree")
    def test___cleanup_run(self, shutil_rmtree_mock, *uninteresting_mocks):
        run_tmpdir = Path("run_tmpdir")

        mafft = MAFFT(self.executable, self.tmpdir)
        mafft._cleanup_run(run_tmpdir)

        shutil_rmtree_mock.assert_called_once_with(run_tmpdir)

    @patch("shutil.which", return_value=1)
    @patch.object(MAFFT, MAFFT._create_new_sequences_file.__name__, return_value=Path("_create_new_sequences_file"))
    @patch("os.environ", {"fake_key_1": "fake_value_1"})
    @patch.object(Path, Path.mkdir.__name__)
    @patch("make_prg.utils.msa_aligner.create_temp_dir", return_value=Path("fake_tmp_dir"))
    @patch("builtins.open", new_callable=mock_open)
    @patch.object(SeqIO, SeqIO.write.__name__)
    @patch.object(MAFFT, MAFFT._run_aligner.__name__)
    @patch("make_prg.utils.msa_aligner.load_alignment_file", return_value="updated_alignment")
    @patch.object(MAFFT, MAFFT._cleanup_run.__name__)
    def test___get_updated_alignment(self, cleanup_run_mock, load_alignment_file_mock, run_aligner_mock, write_mock,
                                     open_mock, create_temp_dir_mock, *uninteresting_mocks):
        mafft = MAFFT(self.executable, self.tmpdir)
        alignment = make_alignment([
            "AAAA",
            "CCCC",
            "GGGG",
            "TTTT",
        ])
        new_seqs = {
            "AATA",
            "CCGC",
            "GGCG",
            "TTAT"

        }

        expected_updated_alignment = "updated_alignment"
        actual_updated_alignment = mafft.get_updated_alignment(alignment, new_seqs)

        self.assertEqual(expected_updated_alignment, actual_updated_alignment)

        # TODO: ensure calls are in the correct order
        create_temp_dir_mock.assert_called_once_with(mafft._tmpdir)
        open_mock.assert_called_once_with(Path('fake_tmp_dir/previous_msa.fa'), 'w')

        write_mock.assert_called_once()
        # did not manage to assert equality of arg [1] (which is the open_mock handler)
        # so asserting equality of args [0] and [2]
        self.assertEqual(alignment, write_mock.call_args[0][0])
        self.assertEqual("fasta", write_mock.call_args[0][2])

        run_aligner_mock.assert_called_once_with(
            'fake_exec --auto --quiet --thread 1 --add _create_new_sequences_file fake_tmp_dir/previous_msa.fa > fake_tmp_dir/updated_msa.fa',
            {"fake_key_1": "fake_value_1", "TMPDIR": "fake_tmp_dir"})
        load_alignment_file_mock.assert_called_once_with("fake_tmp_dir/updated_msa.fa", "fasta")
        cleanup_run_mock.assert_called_once_with(Path("fake_tmp_dir"))