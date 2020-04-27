import logging
import os
from pathlib import Path

from make_prg import make_prg_from_msa, io_utils


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
        logging.info(f"Write PRG file to {prefix}.prg")
        io_utils.write_prg(prefix, aseq.prg)
        m = aseq.max_nesting_level_reached
        logging.info(f"Max_nesting_reached\t{m}")

    logging.info(f"Write GFA file to {prefix}.gfa")
    io_utils.write_gfa(f"{prefix}.gfa", aseq.prg)

    summary_file = Path(prefix).parent / "summary.tsv"
    with summary_file.open("a") as s:
        s.write(
            f"{options.MSA}\t{aseq.site - 2}\t"
            f"{aseq.max_nesting_level_reached}\t{aseq.prop_in_match_intervals}\n"
        )
