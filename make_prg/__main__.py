import argparse

import make_prg

NESTING_LVL = 5
MIN_MATCH_LEN = 7


def main():
    parser = argparse.ArgumentParser(
        prog="make_prg",
        usage="make_prg <subcommand> <options>",
        description="script to run make_prg subcommands",
    )

    parser.add_argument("--version", action="version", version=make_prg.__version__)
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )

    # _____________________________ prg_from_msa ______________________________#
    subparser_prg_from_msa = subparsers.add_parser(
        "prg_from_msa",
        usage="make_prg prg_from_msa [options] <MSA input file>",
        help="Make PRG from multiple sequence alignment",
    )

    subparser_prg_from_msa.add_argument(
        "MSA",
        action="store",
        type=str,
        help=(
            "Input file: a multiple sequence alignment in supported alignment_format. "
            "If not in aligned fasta alignment_format, use -f to input the "
            "alignment_format type"
        ),
    )
    subparser_prg_from_msa.add_argument(
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
    subparser_prg_from_msa.add_argument(
        "--max_nesting",
        dest="max_nesting",
        action="store",
        type=int,
        default=NESTING_LVL,
        help="Maximum number of levels to use for nesting. Default: {}".format(
            NESTING_LVL
        ),
    )
    subparser_prg_from_msa.add_argument(
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
    subparser_prg_from_msa.add_argument(
        "-p", "--prefix", dest="output_prefix", action="store", help="Output prefix"
    )
    subparser_prg_from_msa.add_argument(
        "--staggered_start",
        dest="staggered_start",
        action="store_true",
        help="Assume dashes at start of alignment are due to missing data, rather than indels",
    )
    subparser_prg_from_msa.add_argument(
        "--no_overwrite",
        dest="no_overwrite",
        action="store_true",
        help="Do not overwrite pre-existing prg file with same name",
    )
    subparser_prg_from_msa.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )
    subparser_prg_from_msa.set_defaults(func=make_prg.subcommands.prg_from_msa.run)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
