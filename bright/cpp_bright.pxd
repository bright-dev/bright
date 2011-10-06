"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

from pyne cimport std

cdef extern from "bright.h" namespace "bright":
    std.string BRIGHT_DATA

    void bright_start() except +

    set[int] track_isos
    vector[int] track_isos_order

    void load_track_isos_hdf5(std.string, std.string, bint) except +
    void load_track_isos_text(std.string, bint) except +

    void sort_track_isos()

    int verbosity
    bint write_hdf5
    bint write_text

    std.string output_filename


