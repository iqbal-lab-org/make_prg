from dataclasses import dataclass
from typing import Optional

from loguru import logger

from make_prg.update.denovo_variants import DenovoVariantsDB
from make_prg.utils.msa_aligner import MSAAligner


@dataclass
class UpdateSharedData:
    """
    Data that is shared between child processes
    """

    denovo_variants_db: DenovoVariantsDB
    aligner: MSAAligner


class SingletonUpdateSharedData(object):
    """
    This class follows the singleton pattern and represents read-only data that is shared between child processes
    when running make_prg update. As python does not implement multithreading, we have to use multiprocessing. However,
    we can still efficiently share in-memory data between a parent process and a list of child processes with the copy-on-write
    optimisation present in Unix systems (https://en.wikipedia.org/wiki/Copy-on-write). As long as the data does not
    change, the child processes will refer to the parent process data in memory. If we give this large data as input to
    the multiprocessing task, it will be pickled/unpickled, severely reducing efficiency. The programmer is responsible
    to not modify the data in the multiprocessing parts.
    Adapted from https://www.geeksforgeeks.org/singleton-pattern-in-python-a-complete-guide/
    """

    # create the singleton
    def __new__(
        cls,
        denovo_variants_db: Optional[DenovoVariantsDB] = None,
        aligner: Optional[MSAAligner] = None,
    ):
        if not hasattr(cls, "instance"):
            logger.trace("Creating SingletonUpdateSharedData")
            cls.instance = super(SingletonUpdateSharedData, cls).__new__(cls)
            cls.instance.data = None
        return cls.instance

    # the first time constructing an object of this class, we need valid denovo_variants_db and aligner to correctly
    # initialise the singleton. For the subsequent constructions, these parameters are ignored
    def __init__(
        self,
        denovo_variants_db: Optional[DenovoVariantsDB] = None,
        aligner: Optional[MSAAligner] = None,
    ):
        singleton_has_not_been_initialised = self.__class__.instance.data is None
        if singleton_has_not_been_initialised:
            logger.trace("Initialising SingletonUpdateSharedData")
            self.__class__.instance.data = UpdateSharedData(denovo_variants_db, aligner)
