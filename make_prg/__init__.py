from importlib import metadata

__version__ = metadata.version("make_prg")

__all__ = ["from_msa", "subcommands"]

from Bio.AlignIO import MultipleSeqAlignment

MSA = MultipleSeqAlignment
