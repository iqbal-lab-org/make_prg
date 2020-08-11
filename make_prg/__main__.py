import argparse

from make_prg import __version__
from make_prg.subcommands import prg_from_msa


def main():
    parser = argparse.ArgumentParser(
        prog="make_prg",
        usage="make_prg <subcommand> <options>",
        description="script to run make_prg subcommands",
    )

    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )

    prg_from_msa.register_parser(subparsers)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
