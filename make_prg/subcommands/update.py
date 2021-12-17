import multiprocessing
import os
from pathlib import Path
from loguru import logger
from make_prg.utils import io_utils, gfa
from make_prg.update.denovo_variants import DenovoVariantsDB
from make_prg.prg_builder import PrgBuilderZipDatabase, LeafNotFoundException
from make_prg.utils.msa_aligner import MAFFT, MSAAligner
from dataclasses import dataclass


def register_parser(subparsers):
    subparser_update_prg = subparsers.add_parser(
        "update",
        usage="make_prg update",
        help="Update PRGs given new sequences.",
    )

    def check_if_is_update_file(argument: str):
        argument = Path(argument)
        is_zip_file = argument.suffix == ".zip"
        is_update_DS = argument.suffix == ".update_DS"
        if not is_zip_file and not is_update_DS:
            subparser_update_prg.error(f"{argument} is not a update_DS nor a zip file.")
        return argument
    subparser_update_prg.add_argument(
        "-u",
        "--update_DS",
        action="store",
        type=lambda argument: check_if_is_update_file(argument),
        required=True,
        help=(
            "Filepath to the update data structures (a *.update_DS.zip file created "
            "from make_prg from_msa)"
        ),

    )
    subparser_update_prg.add_argument(
        "-d",
        "--denovo_paths",
        action="store",
        type=str,
        required=True,
        help=(
            "Filepath containing denovo sequences. Should point to a denovo_paths.txt file"
        ),
    )
    subparser_update_prg.add_argument(
        "-o",
        "--output_prefix",
        action="store",
        type=str,
        required=True,
        help="Prefix for the output files",
    )
    subparser_update_prg.add_argument(
        "-t",
        "--threads",
        action="store",
        type=int,
        default=1,
        help="Number of threads Default: %(default)d",
    )
    subparser_update_prg.add_argument(
        "--mafft",
        help="Path to MAFFT executable. By default, it is assumed to be on $PATH",
        default="mafft",
    )
    subparser_update_prg.add_argument(
        "--output_graphs",
        dest="output_graphs",
        action="store_true",
        default=False,
        help="Outputs the recursive tree graphical representation (for development use only)",
    )
    subparser_update_prg.set_defaults(func=run)

    return subparser_update_prg


@dataclass
class UpdateSharedData:
    denovo_variants_db: DenovoVariantsDB
    aligner: MSAAligner


def update(locus_name: str):
    global options
    prg_builder_zip_db = PrgBuilderZipDatabase(options.update_DS)
    prg_builder_zip_db.load()
    prg_builder_for_locus = prg_builder_zip_db.get_PrgBuilder(locus_name)

    global update_shared_data
    prg_builder_for_locus.aligner = update_shared_data.aligner
    update_data_list = update_shared_data.denovo_variants_db.locus_name_to_update_data.get(locus_name, [])
    nb_of_variants_sucessfully_updated = 0
    nb_of_variants_with_failed_update = 0

    we_have_variants = len(update_data_list) > 0
    if we_have_variants:
        logger.debug(f"Updating {locus_name} ...")

        leaves_to_update = set()
        for update_data in update_data_list:
            try:
                prg_builder_tree_node = prg_builder_for_locus.get_node_given_interval(
                    update_data.ml_path_node_key
                )
                prg_builder_tree_node.add_data_to_batch_update(
                    update_data
                )
                leaves_to_update.add(prg_builder_tree_node)
                nb_of_variants_sucessfully_updated += 1
            except LeafNotFoundException as exc:
                logger.warning(f"Failed finding leaf: {exc}")
                nb_of_variants_with_failed_update += 1

        # update the modified leaves
        # Note: the sorted() is needed for determinism so that we can compare indexes and test outputs
        for leaf in sorted(leaves_to_update, key=lambda node: node.id):
            leaf.batch_update()
        logger.debug(
            f"Updated {locus_name}: {nb_of_variants_sucessfully_updated} denovo sequences added!"
        )
    else:
        logger.debug(f"{locus_name} has no new variants, no update needed")

    # regenerate PRG
    temp_dir = io_utils.get_temp_dir_for_multiprocess(options.temp_dir)
    locus_prefix = temp_dir / locus_name

    logger.info(f"Writing output files of locus {locus_name}")
    prg = prg_builder_for_locus.build_prg()
    prg_builder_for_locus.write_prg(str(locus_prefix), prg)
    gfa.GFA_Output.write_gfa(str(locus_prefix), prg)
    prg_builder_for_locus.serialize(f"{locus_prefix}.pickle")
    with open(f"{locus_prefix}.stats", "w") as stats_filehandler:
        print(
            f"{locus_name} {nb_of_variants_sucessfully_updated} {nb_of_variants_with_failed_update}",
            file=stats_filehandler,
        )
    if options.output_graphs:
        prg_builder_for_locus.output_debug_graphs(Path(options.output_prefix + "_debug_graphs"))


def run(cl_options):
    global options
    options = cl_options
    global update_shared_data

    if io_utils.output_files_already_exist(options.output_prefix):
        raise RuntimeError("One or more output files already exists, aborting run...")

    # read input data
    logger.info("Checking Multiple Sequence Aligner...")
    output_dir = Path(options.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = io_utils.create_temp_dir(output_dir)
    options.temp_dir = temp_dir
    msa_temp_path = temp_dir / "msa_temp"
    mafft_aligner = MAFFT(executable=options.mafft, tmpdir=msa_temp_path)

    prg_builder_zip_db = None
    try:
        logger.info("Reading update data structures...")
        prg_builder_zip_db = PrgBuilderZipDatabase(options.update_DS)
        prg_builder_zip_db.load()
        logger.info(f"Reading {options.denovo_paths}...")
        denovo_variants_db = DenovoVariantsDB(options.denovo_paths)
        update_shared_data = UpdateSharedData(denovo_variants_db, mafft_aligner)

        output_dir = Path(options.output_prefix).parent
        output_dir.mkdir(exist_ok=True)

        # update all PRGs with denovo sequences
        logger.info(f"Using {options.threads} threads to update PRGs...")
        multithreaded_input = prg_builder_zip_db.get_loci_names()
        with multiprocessing.Pool(options.threads, maxtasksperchild=1) as pool:
            pool.map(update, multithreaded_input, chunksize=1)
        logger.success(f"All PRGs updated!")

        is_a_single_MSA = prg_builder_zip_db.get_number_of_loci() == 1
        io_utils.create_final_files(temp_dir, options.output_prefix, is_a_single_MSA=is_a_single_MSA, output_stats=True)
        logger.success("All done!")
    finally:
        if prg_builder_zip_db is not None:
            prg_builder_zip_db.close()
