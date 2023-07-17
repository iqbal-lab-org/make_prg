import gzip
import os
import tempfile
from io import StringIO
from pathlib import Path
from typing import Dict, Union
from zipfile import ZipFile

from Bio import AlignIO
from Bio.Seq import Seq

from make_prg import MSA
from make_prg.subcommands.output_type import OutputType
from make_prg.utils.seq_utils import get_majority_consensus_from_MSA


def load_alignment_file(
    msa_file: Union[str, Path, StringIO], alignment_format: str
) -> MSA:
    if isinstance(msa_file, StringIO):
        alignment = AlignIO.read(msa_file, alignment_format)
    else:
        msa_file = str(msa_file)
        if msa_file.endswith(".gz"):
            with gzip.open(msa_file, "rt") as handle:
                alignment = AlignIO.read(handle, alignment_format)
        else:
            with open(msa_file, "r") as handle:
                alignment = AlignIO.read(handle, alignment_format)

    # upper case seqs
    for record in alignment:
        record.seq = record.seq.upper()

    # Compute the consensus sequence
    consensus = get_majority_consensus_from_MSA(alignment)

    # Replace 'N' with the consensus sequence in each record
    for record in alignment:
        record.seq = Seq(
            "".join(
                [
                    consensus[i] if nucleotide == "N" else nucleotide
                    for i, nucleotide in enumerate(str(record.seq))
                ]
            )
        )

    return alignment


# Note: not unit tested
# From https://gist.github.com/jacobtomlinson/9031697
def remove_empty_folders(path: str, remove_root: bool = True):
    if not os.path.isdir(path):
        return

    # remove empty subfolders
    files = os.listdir(path)
    for f in files:
        fullpath = os.path.join(path, f)
        if os.path.isdir(fullpath):
            remove_empty_folders(fullpath)

    # if folder empty, delete it
    files = os.listdir(path)
    if len(files) == 0 and remove_root:
        os.rmdir(path)


# Note: not unit tested
def output_files_already_exist(output_type: OutputType, output_prefix: str):
    files_to_check = []
    if output_type.prg:
        files_to_check.extend(
            [Path(output_prefix + ".prg.fa"), Path(output_prefix + ".update_DS.zip")]
        )
    if output_type.gfa:
        files_to_check.extend(
            [Path(output_prefix + ".prg.gfa"), Path(output_prefix + ".prg.gfa.zip")]
        )
    if output_type.binary:
        files_to_check.extend(
            [Path(output_prefix + ".prg.bin"), Path(output_prefix + ".prg.bin.zip")]
        )

    for file_to_check in files_to_check:
        if file_to_check.exists():
            return True

    return False


def create_temp_dir(output_dir: Path) -> Path:
    temp_dir = tempfile.mkdtemp(dir=output_dir)
    return Path(temp_dir)


def get_temp_dir_for_multiprocess(root_temp_dir: Path):
    temp_dir = root_temp_dir / "mp_temp"
    os.makedirs(temp_dir, exist_ok=True)
    return temp_dir


def zip_set_of_files(zip_filepath: Path, filename_to_filepath: Dict[str, Path]):
    is_a_zip_file = zip_filepath.suffix == ".zip"
    assert is_a_zip_file, "zip_set_of_files() was not given a .zip filepath"
    with ZipFile(zip_filepath, "w") as zip_file:
        for filename, filepath in filename_to_filepath.items():
            zip_file.write(filepath, filename)
