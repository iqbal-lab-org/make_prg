from typing import List
import multiprocessing
from pathlib import Path
from loguru import logger
from make_prg import prg_builder
from make_prg.from_msa import NESTING_LVL, MIN_MATCH_LEN
from make_prg.utils import io_utils, gfa


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
        help=(
            "Multiple sequence alignment file or a directory containing such files"
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
            "alignment_format. See http://biopython.org/wiki/AlignIO. Default: %(default)s"
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


def get_all_input_files(input_path: str) -> List[Path]:
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"{input_path} does not exist")

    if input_path.is_file():
        all_files = [input_path]
    else:
        all_files = [
            path.resolve() for path in input_path.iterdir() if path.is_file()
        ]
    return all_files


def process_MSA(msa_filepath: Path):
    global options
    logger.info(f"Generating PRG for {msa_filepath}...")
    msa_name = msa_filepath.name
    locus_name = msa_filepath.with_suffix("").name

    temp_dir = io_utils.get_temp_dir_for_multiprocess(options.temp_dir)
    prefix = str(temp_dir / msa_name)

    try:
        builder = prg_builder.PrgBuilder(
            locus_name=locus_name,
            msa_file=msa_filepath,
            alignment_format=options.alignment_format,
            max_nesting=options.max_nesting,
            min_match_length=options.min_match_length
        )

        logger.info(f"Writing output files of locus {msa_name}")
        prg = builder.build_prg()

        if options.output_type.prg or options.output_type.binary:
            builder.write_prg(prefix, prg)

        if options.output_type.gfa:
            gfa.GFA_Output.write_gfa(prefix, prg)

        if options.output_type.prg:
            builder.serialize(f"{prefix}.pickle")

        if options.output_graphs:
            builder.output_debug_graphs(Path(options.output_prefix + "_debug_graphs"))

    except ValueError as value_error:
        if "No records found in handle" in value_error.args[0]:
            logger.warning(f"No records found in MSA {msa_filepath}, skipping...")
        else:
            raise value_error


def run(cl_options):
    global options
    options = cl_options

    logger.info("Getting input files...")
    input_files = get_all_input_files(options.input)

    there_is_no_input_files = len(input_files) == 0
    if there_is_no_input_files:
        raise FileNotFoundError(f"No input files found in {options.input}")

    if io_utils.output_files_already_exist(options.output_prefix):
        raise RuntimeError("One or more output files already exists, aborting run...")

    output_dir = Path(options.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = io_utils.create_temp_dir(output_dir)
    options.temp_dir = temp_dir

    logger.info(f"Using {options.threads} threads to generate PRGs...")
    with multiprocessing.Pool(options.threads, maxtasksperchild=1) as pool:
        pool.map(process_MSA, input_files, chunksize=1)
    logger.success(f"All PRGs generated!")

    is_a_single_file = len(input_files) == 1
    io_utils.create_final_files(temp_dir, options.output_prefix, options.output_type,
                                is_a_single_MSA=is_a_single_file)
    logger.success("All done!")
