import logging
from pathlib import Path

from make_prg.from_msa import prg_builder, NESTING_LVL, MIN_MATCH_LEN
from make_prg import io_utils


def register_parser(subparsers):
    subparser_msa = subparsers.add_parser(
        "from_msa",
        usage="make_prg from_msa [options] <MSA input file>",
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
        "-N",
        "--max_nesting",
        dest="max_nesting",
        action="store",
        type=int,
        default=NESTING_LVL,
        help="Maximum number of levels to use for nesting. Default: %(default)s",
    )
    subparser_msa.add_argument(
        "-L",
        "--min_match_length",
        dest="min_match_length",
        action="store",
        type=int,
        default=MIN_MATCH_LEN,
        help=(
            "Minimum number of consecutive characters which must be identical for a "
            "match. Default: %(default)s"
        ),
    )
    subparser_msa.add_argument(
        "-o",
        "--outdir",
        dest="output_dir",
        action="store",
        default=".",
        help="Output directory. Default: Current working directory",
    )
    subparser_msa.add_argument(
        "-n",
        "--prg_name",
        dest="prg_name",
        action="store",
        help="Prg file name. Default: MSA file name",
    )
    subparser_msa.add_argument(
        "-S",
        "--seqid",
        help="Sequence identifier to use for the output sequence/PRG. Default is the file name",
    )
    subparser_msa.add_argument(
        "--no_overwrite",
        dest="no_overwrite",
        action="store_true",
        help="By default, a prg file which already exists is overwritten and the prg is recomputed. Use this option to keep an existing prg file instead.",
    )
    subparser_msa.add_argument("--log", help="Path to write log to. Default is stderr")
    subparser_msa.set_defaults(func=run)

    return subparser_msa


def run(options):
    MSA_file = Path(options.MSA).resolve()
    if not MSA_file.exists():
        raise ValueError(f"File not found: {options.MSA}")
    output_dir = Path(options.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if options.prg_name is None:
        options.prg_name = MSA_file.stem
    ofile_prefix = output_dir / options.prg_name

    # Set up file logging
    formatter = logging.Formatter(
        fmt="%(levelname)s %(asctime)s %(message)s", datefmt="%d/%m/%Y %I:%M:%S"
    )
    handler = (
        logging.FileHandler(options.log) if options.log else logging.StreamHandler()
    )
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)

    logging.info(
        "Input parameters max_nesting: %d, min_match_length: %d",
        options.max_nesting,
        options.min_match_length,
    )

    prg_fname = ofile_prefix.with_suffix(".prg")
    if prg_fname.exists() and options.no_overwrite:
        logging.info(f"Re-using existing prg file {prg_fname}")
        aseq = prg_builder.PrgBuilder(
            options.MSA,
            alignment_format=options.alignment_format,
            max_nesting=options.max_nesting,
            min_match_length=options.min_match_length,
            prg_file=prg_fname,
        )
    else:
        aseq = prg_builder.PrgBuilder(
            options.MSA,
            alignment_format=options.alignment_format,
            max_nesting=options.max_nesting,
            min_match_length=options.min_match_length,
        )
        logging.info(f"Write PRG file to {prg_fname}")
        io_utils.write_prg(prg_fname, aseq.prg, options)
        m = aseq.max_nesting_level_reached
        logging.info(f"Max_nesting_reached\t{m}")

    gfa_fname = ofile_prefix.with_suffix(".gfa")
    logging.info(f"Write GFA file to {gfa_fname}")
    io_utils.write_gfa(gfa_fname, aseq.prg)

    summary_fname = output_dir / "summary.tsv"
    with summary_fname.open("a") as s:
        s.write(
            f"{options.MSA}\t{aseq.site - 2}\t"
            f"{aseq.max_nesting_level_reached}\t{aseq.prop_in_match_intervals}\n"
        )
