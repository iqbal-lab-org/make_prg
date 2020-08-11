# ___Constants/Aliases___ #
from Bio.AlignIO import MultipleSeqAlignment

MSA = MultipleSeqAlignment
NESTING_LVL = 5
MIN_MATCH_LEN = 7

# ___Version___ #
from pkg_resources import get_distribution

try:
    __version__ = get_distribution("make_prg").version
except:
    __version__ = "local"

__all__ = ["make_prg_from_msa", "subcommands", "io_utils", "utils"]

from make_prg import *
