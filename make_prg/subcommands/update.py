import multiprocessing
from pathlib import Path

from loguru import logger

from make_prg.prg_builder import LeafNotFoundException, PrgBuilderZipDatabase
from make_prg.update.denovo_variants import DenovoVariantsDB
from make_prg.update.update_shared_data import SingletonUpdateSharedData
from make_prg.utils import gfa, io_utils
from make_prg.utils.input_output_files import InputOutputFilesUpdate
from make_prg.utils.misc import should_output_debug_graphs
from make_prg.utils.msa_aligner import MAFFT


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
        "--update-DS",
        dest="update_DS",
        action="store",
        type=lambda argument: check_if_is_update_file(argument),
        required=True,
        help=(
            "Filepath to the update data structures (a *.update_DS.zip file created "
            "from make_prg from_msa or update)"
        ),
    )
    subparser_update_prg.add_argument(
        "-o",
        "--output-prefix",
        dest="output_prefix",
        action="store",
        type=str,
        required=True,
        help="Prefix for the output files",
    )
    subparser_update_prg.add_argument(
        "-d",
        "--denovo-paths",
        dest="denovo_paths",
        action="store",
        type=str,
        required=True,
        help=(
            "Filepath containing denovo sequences. Should point to a denovo_paths.txt "
            "file"
        ),
    )
    subparser_update_prg.add_argument(
        "-D",
        "--deletion-threshold",
        dest="long_deletion_threshold",
        action="store",
        type=int,
        default=10,
        help=(
            "Ignores long deletions of the given size or longer. If long deletions "
            "should not be ignored, "
            "put a large value. Default: %(default)d"
        ),
    )
    subparser_update_prg.set_defaults(func=run)

    return subparser_update_prg


def update(options, input_and_output_files: InputOutputFilesUpdate):
    # TODO: handle failed runs here?

    locus_name = input_and_output_files.locus_name
    prg_builder_zip_db = PrgBuilderZipDatabase(options.update_DS)
    prg_builder_zip_db.load()
    prg_builder_for_locus = prg_builder_zip_db.get_PrgBuilder(locus_name)

    # retrieves singleton data
    update_shared_data = SingletonUpdateSharedData()
    prg_builder_for_locus.aligner = update_shared_data.data.aligner
    update_data_list = (
        update_shared_data.data.denovo_variants_db.locus_name_to_update_data.get(
            locus_name, []
        )
    )
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
                prg_builder_tree_node.add_data_to_batch_update(update_data)
                leaves_to_update.add(prg_builder_tree_node)
                nb_of_variants_sucessfully_updated += 1
            except LeafNotFoundException as exc:
                logger.warning(f"Failed finding leaf: {exc}")
                nb_of_variants_with_failed_update += 1

        # update the modified leaves
        # Note: the sorted() is needed for determinism so that we can compare indexes
        # and test outputs
        for leaf in sorted(leaves_to_update, key=lambda node: node.node_id):
            leaf.batch_update()
        logger.debug(
            f"Updated {locus_name}: {nb_of_variants_sucessfully_updated} denovo "
            "sequences added!"
        )
    else:
        logger.debug(f"{locus_name} has no new variants, no update needed")

    # regenerate PRG
    logger.info(f"Writing output files of locus {locus_name}")
    temp_prefix = input_and_output_files.temp_prefix
    prg = prg_builder_for_locus.build_prg()

    if options.output_type.prg:
        prg_builder_for_locus.write_prg_as_text(str(temp_prefix), prg)
        prg_builder_for_locus.serialize(f"{temp_prefix}.pickle")

    if options.output_type.binary:
        prg_builder_for_locus.write_prg_as_binary(str(temp_prefix), prg)

    if options.output_type.gfa:
        gfa.GFA_Output.write_gfa(str(temp_prefix), prg)

    if should_output_debug_graphs():
        from make_prg.utils.recursive_tree_drawer import RecursiveTreeDrawer

        RecursiveTreeDrawer.output_debug_graphs(
            prg_builder_for_locus, Path(options.output_prefix + "_debug_graphs")
        )

    with open(f"{temp_prefix}.stats", "w") as stats_filehandler:
        print(
            f"{locus_name} {nb_of_variants_sucessfully_updated} "
            f"{nb_of_variants_with_failed_update}",
            file=stats_filehandler,
        )


def run(cl_options):
    options = cl_options

    if not options.force and io_utils.output_files_already_exist(
        options.output_type, options.output_prefix
    ):
        raise RuntimeError("One or more output files already exists, aborting run...")

    # read input data
    logger.debug("Checking Multiple Sequence Aligner...")
    output_dir = Path(options.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    root_temp_dir = io_utils.create_temp_dir(output_dir)
    msa_temp_path = root_temp_dir / "msa_temp"
    mafft_aligner = MAFFT(tmpdir=msa_temp_path)

    prg_builder_zip_db = None
    try:
        logger.info("Reading update data structures...")
        prg_builder_zip_db = PrgBuilderZipDatabase(options.update_DS)
        prg_builder_zip_db.load()
        logger.info(f"Reading {options.denovo_paths}...")
        denovo_variants_db = DenovoVariantsDB(
            options.denovo_paths, options.long_deletion_threshold
        )

        # initialises a singleton storing the denovo_variants_db and the aligner for the update operation
        SingletonUpdateSharedData(denovo_variants_db, mafft_aligner)

        output_dir = Path(options.output_prefix).parent
        output_dir.mkdir(exist_ok=True)

        mp_temp_dir = io_utils.get_temp_dir_for_multiprocess(root_temp_dir)
        loci = prg_builder_zip_db.get_loci_names()
        input_and_output_files = (
            InputOutputFilesUpdate.get_list_of_InputOutputFilesUpdate(
                loci, options.output_type, mp_temp_dir
            )
        )
        args = [(options, iof) for iof in input_and_output_files]

        # update all PRGs with denovo sequences
        logger.info(f"Using {options.threads} threads to update PRGs...")

        with multiprocessing.Pool(options.threads, maxtasksperchild=1) as pool:
            pool.starmap(update, args, chunksize=1)
        logger.success("All PRGs updated!")

        InputOutputFilesUpdate.create_final_files(
            input_and_output_files, options.output_prefix
        )
        io_utils.remove_empty_folders(str(root_temp_dir))
        logger.success("All done!")
    finally:
        if prg_builder_zip_db is not None:
            prg_builder_zip_db.close()
