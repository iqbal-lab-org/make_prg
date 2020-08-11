from Bio.AlignIO import MultipleSeqAlignment

MSA = MultipleSeqAlignment

from pkg_resources import get_distribution

try:
    __version__ = get_distribution("make_prg").version
except:
    __version__ = "local"

__all__ = ["make_prg_from_msa", "subcommands", "io_utils", "seq_utils"]

from make_prg import *
