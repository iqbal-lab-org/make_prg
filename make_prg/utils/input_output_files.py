import fileinput
import os
import re
import shutil
from abc import ABC
from pathlib import Path
from typing import List, Optional, Tuple

from loguru import logger

from make_prg.prg_builder import PrgBuilderZipDatabase
from make_prg.subcommands.output_type import OutputType
from make_prg.utils.io_utils import zip_set_of_files


# Note: several methods in this class are just input/output gathering, so they are not
# unit tested but they are integration tested
class InputOutputFiles(ABC):
    """
    Class that represents the set of input files consumed and output files created for
    each make_prg subcommand.
    This is an abstract class grouping common attributes and methods.
    Each make_prg subcommand should have a child concrete class.
    """

    def __init__(
        self,
        locus_name: str,
        output_type: OutputType,
        root_temp_dir: Path,
        output_stats: bool,
    ):
        self.locus_name: str = locus_name
        self.output_type: OutputType = output_type
        self.temp_prefix: str = str(root_temp_dir / self.locus_name)

        self._set_output_file(output_type.prg, "prg", ".fa")
        self._set_output_file(output_type.prg, "pickle")
        self._set_output_file(output_type.gfa, "gfa")
        self._set_output_file(output_type.binary, "bin")
        self._set_output_file(output_stats, "stats")

    def _set_output_file(self, condition, extension, post_extension=""):
        self.__dict__[extension] = (
            Path(f"{self.temp_prefix}.{extension}{post_extension}")
            if condition
            else None
        )

    def run_was_successful(self) -> bool:
        if self.output_type.prg:
            return self.prg.exists()
        if self.output_type.gfa:
            return self.gfa.exists()
        if self.output_type.binary:
            return self.bin.exists()
        return False

    @staticmethod
    def get_successfull_runs(
        list_of_InputOutputFiles: List["InputOutputFiles"],
    ) -> List["InputOutputFiles"]:
        list_of_InputOutputFiles_that_succeeded = []

        for input_output_files in list_of_InputOutputFiles:
            if input_output_files.run_was_successful():
                list_of_InputOutputFiles_that_succeeded.append(input_output_files)
            else:
                input_output_files.delete_files()

        return list_of_InputOutputFiles_that_succeeded

    @staticmethod
    def create_final_files(
        list_of_InputOutputFiles: List["InputOutputFiles"], output_prefix: str
    ):
        logger.info(
            "Concatenating files from several threads into single final files..."
        )

        is_a_single_MSA = len(list_of_InputOutputFiles) == 1
        output_type = list_of_InputOutputFiles[0].output_type

        if output_type.prg:
            logger.info("Creating FASTA file of PRGs...")
            prg_files = [
                input_output_file.prg for input_output_file in list_of_InputOutputFiles
            ]
            prg_files = sorted(prg_files)  # guarantees reproducibility
            InputOutputFiles.concatenate_text_files(
                prg_files, output_prefix + ".prg.fa"
            )

            # zip all PRG Builders
            logger.info("Creating update data structures...")
            prg_builder_zip_db = PrgBuilderZipDatabase(
                Path(f"{output_prefix}.update_DS.zip")
            )
            locus_to_prg_builder_pickle_path = {
                input_output_file.locus_name: input_output_file.pickle
                for input_output_file in list_of_InputOutputFiles
            }
            prg_builder_zip_db.save(locus_to_prg_builder_pickle_path)

        # zip all encoded PRGs
        if output_type.binary:
            logger.info("Creating binary PRGs...")
            filename_to_encoded_PRG_paths = {
                input_output_file.bin.name: input_output_file.bin
                for input_output_file in list_of_InputOutputFiles
            }

            if is_a_single_MSA:
                shutil.copy(
                    list(filename_to_encoded_PRG_paths.values())[0],
                    f"{output_prefix}.prg.bin",
                )
            else:
                zip_set_of_files(
                    Path(f"{output_prefix}.prg.bin.zip"), filename_to_encoded_PRG_paths
                )

        # zip all GFAs
        if output_type.gfa:
            logger.info("Creating GFAs...")
            filename_to_gfa_paths = {
                input_output_file.gfa.name: input_output_file.gfa
                for input_output_file in list_of_InputOutputFiles
            }
            if is_a_single_MSA:
                shutil.copy(
                    list(filename_to_gfa_paths.values())[0], f"{output_prefix}.prg.gfa"
                )
            else:
                zip_set_of_files(
                    Path(f"{output_prefix}.prg.gfa.zip"), filename_to_gfa_paths
                )

        # sum up stats files and output stats
        output_stats = list_of_InputOutputFiles[0].stats is not None
        if output_stats:
            logger.info("Computing stats on updates...")
            stats_files = [
                input_output_file.stats
                for input_output_file in list_of_InputOutputFiles
            ]
            (
                nb_of_variants_successfully_applied,
                nb_of_variants_that_failed_to_be_applied,
            ) = InputOutputFiles.get_stats_on_variants(stats_files)
            logger.success(
                "Number of variants successfully applied: "
                f"{nb_of_variants_successfully_applied}"
            )
            logger.warning(
                "Number of variants that failed to be applied: "
                f"{nb_of_variants_that_failed_to_be_applied}"
            )

        # cleanup
        for input_output_file in list_of_InputOutputFiles:
            input_output_file.delete_files()

    #################################################################################
    # Input/output files helper methods
    @staticmethod
    def concatenate_text_files(input_filepaths, output_filepath):
        output_filepath_parent_dir = Path(output_filepath).parent
        os.makedirs(output_filepath_parent_dir, exist_ok=True)

        with open(output_filepath, "w") as fout, fileinput.input(
            input_filepaths
        ) as fin:
            for line in fin:
                fout.write(line)

    @staticmethod
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
                nb_of_variants_that_failed_to_be_applied_for_this_locus = int(
                    line_split[2]
                )
                nb_of_variants_that_failed_to_be_applied += (
                    nb_of_variants_that_failed_to_be_applied_for_this_locus
                )
        return (
            nb_of_variants_successfully_applied,
            nb_of_variants_that_failed_to_be_applied,
        )

    @staticmethod
    def _delete_file(filepath: Optional[Path]):
        if filepath is not None and filepath.exists():
            filepath.unlink()

    def delete_files(self):
        self._delete_file(self.prg)
        self._delete_file(self.bin)
        self._delete_file(self.gfa)
        self._delete_file(self.pickle)
        self._delete_file(self.stats)

    #################################################################################


class InputOutputFilesFromMSA(InputOutputFiles):
    def __init__(
        self, input_filepath: Path, output_type: OutputType, root_temp_dir: Path
    ):
        self.input_filepath: Path = input_filepath
        locus_name: str = InputOutputFilesFromMSA.remove_known_input_extensions(
            str(input_filepath.name)
        )
        super().__init__(locus_name, output_type, root_temp_dir, output_stats=False)

    @staticmethod
    def get_list_of_InputOutputFilesFromMSA(
        input_files: List[Path], output_type: OutputType, root_temp_dir: Path
    ) -> List["InputOutputFilesFromMSA"]:
        return [
            InputOutputFilesFromMSA(input_file, output_type, root_temp_dir)
            for input_file in input_files
        ]

    @staticmethod
    def remove_known_input_extensions(input_filepath: str) -> str:
        return re.sub(r"\.(fa|fasta)(\.gz)?$", "", input_filepath)


class InputOutputFilesUpdate(InputOutputFiles):
    def __init__(self, locus_name: str, output_type: OutputType, root_temp_dir: Path):
        super().__init__(locus_name, output_type, root_temp_dir, output_stats=True)

    @staticmethod
    def get_list_of_InputOutputFilesUpdate(
        loci: List[str], output_type: OutputType, root_temp_dir: Path
    ) -> List["InputOutputFilesUpdate"]:
        return [
            InputOutputFilesUpdate(locus, output_type, root_temp_dir) for locus in loci
        ]
