"""Cython header for utils library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector
from libcpp.string cimport string as std_string

cdef extern from "utils.h" namespace "bright":
    std_string BRIGHT_DATA

    void bright_start() except +

    set[int] track_nucs
    vector[int] track_nucs_order

    void load_track_nucs_hdf5(std_string, std_string, bint) except +
    void load_track_nucs_text(std_string, bint) except +

    void sort_track_nucs()

    int verbosity
    bint write_hdf5
    bint write_text

    std_string output_filename
