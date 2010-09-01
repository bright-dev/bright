"""A class to setup, run, and parse ORIGEN."""

import isoname
import metasci.nuke as msn
from metasci import SafeRemove

from char import reactor

from char import FineTimeIndex, FineTime
from char import CoarseTimeIndex, CoarseTime
from char import CoreLoad_zzaaam, CoreLoad_LLAAAM, CoreLoad_MCNP
from char import CoreTran_zzaaam, CoreTran_LLAAAM, CoreTran_MCNP

from n_code import NCode


class NCodeORIGEN(NCode):
    """An ORIGEN neutronics code wrapper class."""

    def __init__(self):
        self.name    = "ORIGEN"
        self.run_str = "o2"

        # Remote file lists
        self.place_remote_files = []
        self.fetch_remote_files = []

