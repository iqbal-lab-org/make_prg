# ___Version___ #
from pkg_resources import get_distribution

try:
    __version__ = get_distribution("make_prg").version
except:
    __version__ = "1.0.0"

__all__ = ["from_msa", "subcommands"]
