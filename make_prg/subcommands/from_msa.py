import multiprocessing
from pathlib import Path
from typing import List

from loguru import logger

from make_prg import prg_builder
from make_prg.from_msa import MIN_MATCH_LEN, NESTING_LVL
from make_prg.utils import gfa, io_utils, seq_utils
from make_prg.utils.input_output_files import InputOutputFilesFromMSA
from make_prg.utils.misc import should_output_debug_graphs


class EmptyMSAError(Exception):
    pass


def register_parser(subparsers):
    subparser_msa = subparsers.add_parser(
        "from_msa",
        usage="make_prg from_msa",
        help="Make PRG from multiple sequence alignment",
    )
    subparser_msa.add_argument(
        "-i",
        "--input",
        action="store",
        type=str,
        required=True,
        help="Multiple sequence alignment file or a directory containing such files",
    )
    subparser_msa.add_argument(
        "-s",
        "--suffix",
        action="store",
        type=str,
        default="",
        help=(
            "If the input parameter (-i, --input) is a directory, then filter for "
            "files with this suffix. If this parameter is not given, all files in the "
            "input directory is considered."
        ),
    )
    subparser_msa.add_argument(
        "-o",
        "--output-prefix",
        dest="output_prefix",
        action="store",
        type=str,
        required=True,
        help="Prefix for the output files",
    )
    subparser_msa.add_argument(
        "-f",
        "--alignment-format",
        dest="alignment_format",
        action="store",
        default="fasta",
        help=(
            "Alignment format of MSA, must be a biopython AlignIO input "
            "alignment_format. See http://biopython.org/wiki/AlignIO. "
            "Default: %(default)s"
        ),
    )
    subparser_msa.add_argument(
        "-N",
        "--max-nesting",
        dest="max_nesting",
        action="store",
        type=int,
        default=NESTING_LVL,
        help="Maximum number of levels to use for nesting. Default: %(default)d",
    )
    subparser_msa.add_argument(
        "-L",
        "--min-match-length",
        dest="min_match_length",
        action="store",
        type=int,
        default=MIN_MATCH_LEN,
        help=(
            "Minimum number of consecutive characters which must be identical for a "
            "match. Default: %(default)d"
        ),
    )

    subparser_msa.set_defaults(func=run)

    return subparser_msa


def get_all_input_files(input_path: str, suffix: str) -> List[Path]:
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"{input_path} does not exist")

    if input_path.is_file():
        all_files = [input_path]
    else:
        all_files = [
            path.resolve()
            for path in input_path.iterdir()
            if path.is_file() and path.name.endswith(suffix)
        ]
    return all_files


def process_MSA(options, input_and_output_files: InputOutputFilesFromMSA):
    locus_name = input_and_output_files.locus_name
    prefix = input_and_output_files.temp_prefix
    logger.info(f"Generating PRG for {locus_name}...")

    try:
        builder = prg_builder.PrgBuilder(
            locus_name=locus_name,
            msa_file=input_and_output_files.input_filepath,
            alignment_format=options.alignment_format,
            max_nesting=options.max_nesting,
            min_match_length=options.min_match_length,
        )

        logger.info(f"Writing output files of locus {locus_name}")
        prg = builder.build_prg()

        if options.output_type.prg:
            builder.write_prg_as_text(prefix, prg)
            builder.serialize(f"{prefix}.pickle")

        if options.output_type.binary:
            builder.write_prg_as_binary(prefix, prg)

        if options.output_type.gfa:
            gfa.GFA_Output.write_gfa(prefix, prg)

        if should_output_debug_graphs():
            from make_prg.utils.recursive_tree_drawer import RecursiveTreeDrawer

            RecursiveTreeDrawer.output_debug_graphs(
                builder, Path(options.output_prefix + "_debug_graphs")
            )

    except ValueError as value_error:
        if "No records found in handle" in value_error.args[0]:
            raise EmptyMSAError(f"No records found in MSA of locus {locus_name}")
        else:
            raise value_error
    except seq_utils.SequenceCurationError as sequence_curation_error:
        logger.warning(
            f"Skipping building PRG for {locus_name}. Error: "
            f"{str(sequence_curation_error)}"
        )


def run(cl_options):
    options = cl_options

    logger.info("Getting input files...")
    input_files = get_all_input_files(options.input, options.suffix)

    there_is_no_input_files = len(input_files) == 0
    if there_is_no_input_files:
        raise FileNotFoundError(f"No input files found in {options.input}")

    if not options.force and io_utils.output_files_already_exist(
        options.output_type, options.output_prefix
    ):
        raise RuntimeError("One or more output files already exists, aborting run...")

    output_dir = Path(options.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    root_temp_dir = io_utils.create_temp_dir(output_dir)
    mp_temp_dir = io_utils.get_temp_dir_for_multiprocess(root_temp_dir)
    input_and_output_files = (
        InputOutputFilesFromMSA.get_list_of_InputOutputFilesFromMSA(
            input_files, options.output_type, mp_temp_dir
        )
    )
    args = [(options, iof) for iof in input_and_output_files]

    logger.info(f"Using {options.threads} threads to generate PRGs...")
    with multiprocessing.Pool(options.threads, maxtasksperchild=1) as pool:
        pool.starmap(process_MSA, args, chunksize=1)
    logger.success("All PRGs generated!")

    successful_input_and_output_files = InputOutputFilesFromMSA.get_successfull_runs(
        input_and_output_files
    )
    all_runs_failed = len(successful_input_and_output_files) == 0
    if all_runs_failed:
        logger.error("No PRGs were built, please check errors")
    else:
        InputOutputFilesFromMSA.create_final_files(
            successful_input_and_output_files, options.output_prefix
        )

    io_utils.remove_empty_folders(str(root_temp_dir))
    logger.success("All done!")
