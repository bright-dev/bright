"""A class to setup, run, and parse Serpent."""

import isoname
import metasci.nuke as msn
from metasci import SafeRemove

from char import reactor
from char import ENDF_FLAG
from char import GroupStructure
from char import FMdic, dicFM

from char import FineTimeIndex, FineTime
from char import CoarseTimeIndex, CoarseTime
from char import CoreLoad_zzaaam, CoreLoad_LLAAAM, CoreLoad_MCNP
from char import CoreTran_zzaaam, CoreTran_LLAAAM, CoreTran_MCNP
from char import metastabletrak, metastableMCNP
from char import mat_number, number_mat
from char import InitialFuelStream
from char import kParticles, kCycles, kCyclesSkip

from n_code import NCode


class NCodeSerpent(NCode):
    """A Serpent neutronics code wrapper class."""

    def __init__(self):
        self.name    = "Serpent"
        self.run_str = "sss-dev"

        # Remote file lists
        self.place_remote_files = [reactor]
        self.fetch_remote_files = [reactor]

