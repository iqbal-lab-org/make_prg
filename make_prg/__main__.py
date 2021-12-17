import argparse
import sys

from loguru import logger

from make_prg import __version__
from make_prg.subcommands import from_msa, update
from typing import List


def setup_common_last_options(parsers: List[argparse.ArgumentParser]):
    for par in parsers:
        par.add_argument(
            "-t",
            "--threads",
            action="store",
            type=int,
            default=1,
            help="Number of threads. Default: %(default)d",
        )

        par.add_argument(
            "-g",
            "--output-graphs",
            dest="output_graphs",
            action="store_true",
            default=False,
            help="Outputs the recursive tree graphical representation (for development use only)"
        )

        par.add_argument(
            "-v",
            "--verbose",
            action='count',
            default=0,
            help="Increase output verbosity (-v for debug, -vv for trace - trace is for developers only)",
        )

        par.add_argument("--log", help="Path to write log to. Default is stderr")


def setup_logger(args: argparse.Namespace):
    if "verbose" in args:  # args.verbose does not exist if make_prg is called with no args
        log_levels = ["INFO", "DEBUG", "TRACE"]
        log_level = log_levels[min(args.verbose, len(log_levels) - 1)]
        log_file = args.log or sys.stderr
        handlers = [
            dict(sink=log_file, enqueue=True, level=log_level),
        ]
        logger.configure(handlers=handlers)


def main():
    parser = argparse.ArgumentParser(
        prog="make_prg",
        usage="make_prg <subcommand> <options>",
        description="Subcommand entrypoint",
    )

    parser.add_argument("-V", "--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )

    msa_parser = from_msa.register_parser(subparsers)
    update_parser = update.register_parser(subparsers)
    setup_common_last_options([msa_parser, update_parser])

    args = parser.parse_args()

    setup_logger(args)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
