import os
import shutil
import subprocess
import time
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Set

from Bio import SeqIO
from Bio.AlignIO import MultipleSeqAlignment
from loguru import logger

from make_prg.utils.io_utils import create_temp_dir, load_alignment_file


class NotAValidExecutableError(Exception):
    pass


class ExecutionError(Exception):
    pass


class MSAAligner(ABC):
    def _set_executable(self, executable: str):
        is_valid_executable = shutil.which(executable, mode=os.X_OK) is not None
        if not is_valid_executable:
            raise NotAValidExecutableError(
                f"Given MSA executable {executable} does not work or is invalid"
            )
        self._executable: str = executable

    def _set_tmpdir(self, tmpdir: Path):
        tmpdir.mkdir(parents=True, exist_ok=True)
        self._tmpdir: Path = tmpdir

    # tmpdir == Path("..") means the root temp dir is ".." and we will build per-run temp dir alongside "."
    def __init__(self, executable: str, tmpdir: Path = Path("..")):
        self._set_executable(executable)
        self._set_tmpdir(tmpdir)

    @abstractmethod
    def get_updated_alignment(
        self, current_alignment: MultipleSeqAlignment, new_sequences: Set[str]
    ) -> MultipleSeqAlignment:
        raise NotImplementedError

    @classmethod
    def get_aligner_name(cls) -> str:
        raise NotImplementedError

    def _run_aligner(self, args: str, env=None):
        start = time.time()
        process = subprocess.Popen(
            args,
            stderr=subprocess.PIPE,
            encoding="utf-8",
            shell=True,
            env=env,
        )
        exit_code = process.wait()

        if exit_code != 0:
            raise ExecutionError(
                f"Failed to execute {self.__class__.get_aligner_name()} for arguments "
                f"{args} due to the following error:\n{process.stderr.read()}"
            )
        stop = time.time()
        runtime = stop - start
        logger.debug(
            f"{self.__class__.get_aligner_name()} runtime for arguments {args} in "
            f"seconds: {runtime:.3f}"
        )

    @staticmethod
    def _create_new_sequences_file(directory: Path, new_sequences: Set[str]) -> Path:
        # sorted() is done so that we have a deterministic order of new_sequences and tests run correctly
        new_sequences = sorted(list(new_sequences))
        new_sequences_filepath = directory / "new_sequences.fa"
        with open(new_sequences_filepath, "w") as new_sequences_handler:
            for index_new_seq, new_seq in enumerate(new_sequences):
                print(
                    f">Denovo_path_{index_new_seq}",
                    file=new_sequences_handler,
                )
                print(new_seq, file=new_sequences_handler)
        return new_sequences_filepath


class MAFFT(MSAAligner):
    def __init__(self, tmpdir: Path = Path("..")):
        source_file_dir = os.path.dirname(os.path.abspath(__file__))
        super().__init__(f"{source_file_dir}/mafft-linux64/mafft.bat", tmpdir)

    @classmethod
    def get_aligner_name(cls) -> str:
        return "MAFFT"

    def _cleanup_run(self, run_tmpdir: Path):
        shutil.rmtree(run_tmpdir)

    def get_updated_alignment(
        self, current_alignment: MultipleSeqAlignment, new_sequences: Set[str]
    ) -> MultipleSeqAlignment:
        # setup
        run_tmpdir = create_temp_dir(self._tmpdir)

        current_msa_filepath = run_tmpdir / "previous_msa.fa"
        with open(current_msa_filepath, "w") as current_msa_handler:
            SeqIO.write(current_alignment, current_msa_handler, "fasta")

        new_sequences_filepath = self._create_new_sequences_file(
            run_tmpdir, new_sequences
        )
        new_msa = run_tmpdir / "updated_msa.fa"

        # run
        args = " ".join(
            [
                self._executable,
                "--auto",
                "--quiet",
                "--thread",
                "1",
                "--add",
                str(new_sequences_filepath),
                str(current_msa_filepath),
                ">",
                str(new_msa),
            ]
        )
        env = os.environ
        env["TMPDIR"] = str(run_tmpdir)
        self._run_aligner(args, env)

        # load the updated alignment
        updated_alignment = load_alignment_file(str(new_msa), "fasta")

        # clean up
        self._cleanup_run(run_tmpdir)

        return updated_alignment
