import logging
import os

from make_prg import make_prg_from_msa, utils


def run(options):
    if options.output_prefix is None:
        prefix = options.MSA
    else:
        if os.path.isdir(options.output_prefix):
            prefix = os.path.join(options.output_prefix, os.path.basename(options.MSA))
        else:
            prefix = options.output_prefix
    prefix += ".max_nest%d.min_match%d" % (
        options.max_nesting,
        options.min_match_length,
    )

    if options.verbose:
        log_level = logging.DEBUG
        msg = "Using debug logging"
    else:
        log_level = logging.INFO
        msg = "Using info logging"

    log_file = f"{prefix}.log"
    if os.path.exists(log_file):
        os.unlink(log_file)
    logging.basicConfig(
        filename=log_file,
        level=log_level,
        format="%(asctime)s %(message)s",
        datefmt="%d/%m/%Y %I:%M:%S",
    )
    logging.info(msg)
    logging.info(
        "Input parameters max_nesting: %d, min_match_length: %d",
        options.max_nesting,
        options.min_match_length,
    )

    if os.path.isfile("%s.prg" % prefix) and options.no_overwrite:
        prg_file = "%s.prg" % prefix
        logging.info(f"Re-using existing prg file {prg_file}")
        aseq = make_prg_from_msa.AlignedSeq(
            options.MSA,
            alignment_format=options.alignment_format,
            max_nesting=options.max_nesting,
            min_match_length=options.min_match_length,
            prg_file=prg_file,
        )
    else:
        aseq = make_prg_from_msa.AlignedSeq(
            options.MSA,
            alignment_format=options.alignment_format,
            max_nesting=options.max_nesting,
            min_match_length=options.min_match_length,
        )
        logging.info("Write PRG file to %s.prg", prefix)
        utils.write_prg(prefix, aseq.prg)
        m = aseq.max_nesting_level_reached
        logging.info("Max_nesting_reached\t%d", m)

    logging.info("Write GFA file to %s.gfa", prefix)
    utils.write_gfa("%s.gfa" % prefix, aseq.prg)

    with open("summary.tsv", "a") as s:
        s.write(
            "%s\t%d\t%d\t%f\n"
            % (
                options.MSA,
                aseq.site - 2,
                aseq.max_nesting_level_reached,
                aseq.prop_in_match_intervals,
            )
        )
