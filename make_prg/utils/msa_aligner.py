from abc import ABC, abstractmethod
import shutil
from pathlib import Path
from typing import Set
import time
import subprocess
from loguru import logger
import os
from Bio.AlignIO import MultipleSeqAlignment
from Bio import SeqIO
from make_prg.utils.io_utils import load_alignment_file, create_temp_dir


class NotAValidExecutableError(Exception):
    pass


class ExecutionError(Exception):
    pass


class MSAAligner(ABC):
    def _set_executable(self, executable: str):
        is_valid_executable = shutil.which(executable, mode=os.X_OK) is not None
        if not is_valid_executable:
            raise NotAValidExecutableError(f"Given MSA executable {executable} does not work or is invalid")
        self._executable: str = executable

    def _set_tmpdir(self, tmpdir: Path):
        tmpdir.mkdir(parents=True, exist_ok=True)
        self._tmpdir: Path = tmpdir

    def __init__(self, executable: str, tmpdir: Path = Path("..")):
        self._set_executable(executable)
        self._set_tmpdir(tmpdir)

    @abstractmethod
    def get_updated_alignment(self, current_alignment: MultipleSeqAlignment, new_sequences: Set[str]) -> MultipleSeqAlignment:
        pass

    @classmethod
    def get_aligner_name(cls) -> str:
        pass

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
                f"Failed to execute {self.__class__.get_aligner_name()} for arguments {args} due to the following "
                f"error:\n{process.stderr.read()}"
            )
        stop = time.time()
        runtime = stop - start
        logger.debug(f"{self.__class__.get_aligner_name()} runtime for arguments {args} in seconds: {runtime:.3f}")

    @staticmethod
    def _create_new_sequences_file(directory: Path, new_sequences: Set[str]) -> Path:
        # this is just done so that we have a deterministic order of new_sequences and tests run correctly
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
    @classmethod
    def get_aligner_name(cls) -> str:
        return "MAFFT"

    def _cleanup_run(self, run_tmpdir: Path):
        shutil.rmtree(run_tmpdir)

    def get_updated_alignment(self, current_alignment: MultipleSeqAlignment, new_sequences: Set[str]) -> MultipleSeqAlignment:
        # setup
        run_tmpdir = create_temp_dir(self._tmpdir)

        current_msa_filepath = run_tmpdir / "previous_msa.fa"
        with open(current_msa_filepath, "w") as current_msa_handler:
            SeqIO.write(current_alignment, current_msa_handler, "fasta")

        new_sequences_filepath = self._create_new_sequences_file(run_tmpdir, new_sequences)
        new_msa = run_tmpdir / f"updated_msa.fa"

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
