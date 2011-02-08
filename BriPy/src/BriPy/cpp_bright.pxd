"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std


cdef extern from "../FCComp.h" namespace "FCComps":
    set[int] isos2track
    void load_isos2track_hdf5(std.string, std.string, bint)
    void load_isos2track_text(std.string, bint)

    int verbosity
    bint write_hdf5
    bint write_text

    std.string output_filename


cdef extern from "../FCComp.h":
    cdef cppclass FCComp:
        FCComp()
        FCComp(std.string)
        FCComp(set[std.string], std.string)

        std.string name 
