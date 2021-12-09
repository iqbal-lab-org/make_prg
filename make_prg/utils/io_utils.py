import gzip
import fileinput
import shutil
import uuid
from Bio import AlignIO
import os
from loguru import logger
from make_prg.from_msa import MSA
from typing import Dict, Optional, List, Tuple
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile


def load_alignment_file(msa_file: [str, Path], alignment_format: str) -> MSA:
    msa_file = str(msa_file)
    if msa_file.endswith(".gz"):
        handle = gzip.open(msa_file, "rt")
        alignment = AlignIO.read(handle, alignment_format)
        handle.close()
    else:
        alignment = AlignIO.read(msa_file, alignment_format)
    for record in alignment:
        record.seq = record.seq.upper()
    return alignment


def concatenate_text_files(input_filepaths, output_filepath):
    output_filepath_parent_dir = Path(output_filepath).parent
    os.makedirs(output_filepath_parent_dir, exist_ok=True)

    with open(output_filepath, "w") as fout, fileinput.input(input_filepaths) as fin:
        for line in fin:
            fout.write(line)


# Note: not unit tested
# From https://gist.github.com/jacobtomlinson/9031697
def remove_empty_folders(path: str, remove_root: bool = True):
    if not os.path.isdir(path):
        return

    # remove empty subfolders
    files = os.listdir(path)
    if len(files):
        for f in files:
            fullpath = os.path.join(path, f)
            if os.path.isdir(fullpath):
                remove_empty_folders(fullpath)

    # if folder empty, delete it
    files = os.listdir(path)
    if len(files) == 0 and remove_root:
        os.rmdir(path)


# Note: not unit tested
def output_files_already_exist(output_prefix: str):
    return (
        Path(output_prefix + ".prg.fa").exists()
        or Path(output_prefix + ".update_DS.zip").exists()
        or Path(output_prefix + ".prg.bin").exists()
        or Path(output_prefix + ".prg.bin.zip").exists()
        or Path(output_prefix + ".prg.gfa").exists()
        or Path(output_prefix + ".prg.gfa.zip").exists()
    )


def create_temp_dir(output_prefix: str) -> Path:
    temp_dir = Path(output_prefix + "_" + str(uuid.uuid4()))
    temp_dir.mkdir(parents=True, exist_ok=False)
    return temp_dir


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


# Note: this whole class is not unit tested
class SetOutputFiles:
    """
    Class that represents all files that were generated in processing a PRG (PRG fasta, GFA, pickles, etc)
    """
    def __init__(self, PRG: Optional[Path] = None, binary_PRG: Optional[Path] = None,
                 gfa: Optional[Path] = None, pickle: Optional[Path] = None,
                 stats: Optional[Path] = None):
        self.PRG: Optional[Path] = PRG
        self.binary_PRG: Optional[Path] = binary_PRG
        self.gfa: Optional[Path] = gfa
        self.pickle: Optional[Path] = pickle
        self.stats: Optional[Path] = stats

    @staticmethod
    def clear_temp_extensions(filename: str) -> str:
        while True:
            file_should_be_cleared = filename.endswith(".fa") or filename.endswith(".prg") or \
                                     filename.endswith(".pickle") or filename.endswith(".stats") or \
                                     filename.endswith(".bin") or filename.endswith(".gfa")
            if file_should_be_cleared:
                filename = Path(filename).with_suffix("").name
            else:
                return filename

    @staticmethod
    def _delete_file(filepath: Optional[Path]):
        if filepath is not None:
            filepath.unlink(missing_ok=True)

    def delete_files(self):
        self._delete_file(self.PRG)
        self._delete_file(self.binary_PRG)
        self._delete_file(self.gfa)
        self._delete_file(self.pickle)
        self._delete_file(self.stats)

    @staticmethod
    def get_locus_to_set_of_output_files(temp_root: Path) -> Dict[str, "SetOutputFiles"]:
        locus_to_set_of_output_files = defaultdict(SetOutputFiles)
        for file in temp_root.iterdir():
            if file.is_file():
                locus_name = SetOutputFiles.clear_temp_extensions(file.name)
                if file.name.endswith(".prg.fa"):
                    locus_to_set_of_output_files[locus_name].PRG = file
                elif file.name.endswith(".bin"):
                    locus_to_set_of_output_files[locus_name].binary_PRG = file
                elif file.name.endswith(".gfa"):
                    locus_to_set_of_output_files[locus_name].gfa = file
                elif file.name.endswith(".pickle"):
                    locus_to_set_of_output_files[locus_name].pickle = file
                elif file.name.endswith(".stats"):
                    locus_to_set_of_output_files[locus_name].stats = file
        return locus_to_set_of_output_files


# Note: not unit tested
def get_stats_on_variants(stats_files: List[Path]) -> Tuple[int, int]:
    nb_of_variants_successfully_applied = 0
    nb_of_variants_that_failed_to_be_applied = 0
    for stat_file in stats_files:
        with open(stat_file) as stat_file_fh:
            line_split = stat_file_fh.readline().strip().split()
            nb_of_variants_successfully_applied_for_this_locus = int(line_split[1])
            nb_of_variants_successfully_applied += (
                nb_of_variants_successfully_applied_for_this_locus
            )
            nb_of_variants_that_failed_to_be_applied_for_this_locus = int(line_split[2])
            nb_of_variants_that_failed_to_be_applied += (
                nb_of_variants_that_failed_to_be_applied_for_this_locus
            )
    return nb_of_variants_successfully_applied, nb_of_variants_that_failed_to_be_applied


# this avoids circular import issues, probably has a better way through this though
from make_prg import prg_builder


# Note: not unit tested
def create_final_files(temp_dir: Path, output_prefix: str, is_a_single_MSA: bool, output_stats: bool = False):
    logger.info("Concatenating files from several threads into single final files...")

    logger.info("Creating FASTA file of PRGs...")
    locus_to_set_of_output_files = SetOutputFiles.get_locus_to_set_of_output_files(temp_dir / "mp_temp")

    prg_files = [output_files.PRG for output_files in locus_to_set_of_output_files.values()]
    prg_files = sorted(prg_files)  # guarantees reproducibility
    concatenate_text_files(prg_files, output_prefix + ".prg.fa")

    # zip all PRG Builders
    logger.info("Creating update data structures...")
    prg_builder_zip_db = prg_builder.PrgBuilderZipDatabase(Path(f"{output_prefix}.update_DS.zip"))
    locus_to_prg_builder_pickle_path = {locus: output_files.pickle
                                        for locus, output_files in locus_to_set_of_output_files.items()}
    prg_builder_zip_db.save(locus_to_prg_builder_pickle_path)

    # zip all encoded PRGs
    logger.info("Creating encoded PRGs...")
    filename_to_encoded_PRG_paths = {output_files.binary_PRG.name: output_files.binary_PRG
                                     for output_files in locus_to_set_of_output_files.values()}
    if is_a_single_MSA:
        shutil.copy(list(filename_to_encoded_PRG_paths.values())[0], f"{output_prefix}.prg.bin")
    else:
        zip_set_of_files(Path(f"{output_prefix}.prg.bin.zip"), filename_to_encoded_PRG_paths)

    # zip all GFAs
    logger.info("Creating GFAs...")
    filename_to_gfa_paths = {output_files.gfa.name: output_files.gfa
                             for output_files in locus_to_set_of_output_files.values()}
    if is_a_single_MSA:
        shutil.copy(list(filename_to_gfa_paths.values())[0], f"{output_prefix}.prg.gfa")
    else:
        zip_set_of_files(Path(f"{output_prefix}.prg.gfa.zip"), filename_to_gfa_paths)

    # sum up stats files and output stats
    if output_stats:
        logger.info("Computing stats on updates...")
        stats_files = [output_files.stats for output_files in locus_to_set_of_output_files.values()]
        (
            nb_of_variants_successfully_applied,
            nb_of_variants_that_failed_to_be_applied,
        ) = get_stats_on_variants(stats_files)
        logger.success(
            f"Number of variants successfully applied: {nb_of_variants_successfully_applied}"
        )
        logger.warning(
            f"Number of variants that failed to be applied: {nb_of_variants_that_failed_to_be_applied}"
        )

    # cleanup
    for output_files in locus_to_set_of_output_files.values():
        output_files.delete_files()
    remove_empty_folders(str(temp_dir))
