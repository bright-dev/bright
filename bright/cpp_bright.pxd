"""Cython header for bright library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

from pyne cimport std

cdef extern from "bright.h" namespace "bright":
    std.string BRIGHT_DATA

    void bright_start() except +

    set[int] track_nucs
    vector[int] track_nucs_order

    void load_track_nucs_hdf5(std.string, std.string, bint) except +
    void load_track_nucs_text(std.string, bint) except +

    void sort_track_nucs()

    int verbosity
    bint write_hdf5
    bint write_text

    std.string output_filename


