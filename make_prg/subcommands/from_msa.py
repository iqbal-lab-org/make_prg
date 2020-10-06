import logging
import os
from pathlib import Path

from make_prg.from_msa import aligned_seq, NESTING_LVL, MIN_MATCH_LEN
from make_prg import io_utils


def register_parser(subparsers):
    subparser_msa = subparsers.add_parser(
        "msa",
        usage="make_prg msa [options] <MSA input file>",
        help="Make PRG from multiple sequence alignment",
    )

    subparser_msa.add_argument(
        "MSA",
        action="store",
        type=str,
        help=(
            "Input file: a multiple sequence alignment in supported alignment_format. "
            "If not in aligned fasta alignment_format, use -f to input the "
            "alignment_format type"
        ),
    )
    subparser_msa.add_argument(
        "-f",
        "--alignment_format",
        dest="alignment_format",
        action="store",
        default="fasta",
        help=(
            "Alignment format of MSA, must be a biopython AlignIO input "
            "alignment_format. See http://biopython.org/wiki/AlignIO. Default: fasta"
        ),
    )
    subparser_msa.add_argument(
        "--max_nesting",
        dest="max_nesting",
        action="store",
        type=int,
        default=NESTING_LVL,
        help="Maximum number of levels to use for nesting. Default: {}".format(
            NESTING_LVL
        ),
    )
    subparser_msa.add_argument(
        "--min_match_length",
        dest="min_match_length",
        action="store",
        type=int,
        default=MIN_MATCH_LEN,
        help=(
            "Minimum number of consecutive characters which must be identical for a "
            "match. Default: {}".format(MIN_MATCH_LEN)
        ),
    )
    subparser_msa.add_argument(
        "-p", "--prefix", dest="output_prefix", action="store", help="Output prefix"
    )
    subparser_msa.add_argument(
        "--no_overwrite",
        dest="no_overwrite",
        action="store_true",
        help="Do not overwrite pre-existing prg file with same name",
    )
    subparser_msa.set_defaults(func=run)

    return subparser_msa


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

    # Set up file logging
    log_file = f"{prefix}.log"
    if os.path.exists(log_file):
        os.unlink(log_file)
    formatter = logging.Formatter(
        fmt="%(levelname)s %(asctime)s %(message)s", datefmt="%d/%m/%Y %I:%M:%S"
    )
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)

    logging.info(
        "Input parameters max_nesting: %d, min_match_length: %d",
        options.max_nesting,
        options.min_match_length,
    )

    if os.path.isfile("%s.prg" % prefix) and options.no_overwrite:
        prg_file = "%s.prg" % prefix
        logging.info(f"Re-using existing prg file {prg_file}")
        aseq = aligned_seq.AlignedSeq(
            options.MSA,
            alignment_format=options.alignment_format,
            max_nesting=options.max_nesting,
            min_match_length=options.min_match_length,
            prg_file=prg_file,
        )
    else:
        aseq = aligned_seq.AlignedSeq(
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
