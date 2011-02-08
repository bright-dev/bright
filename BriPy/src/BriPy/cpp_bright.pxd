"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std


cdef extern from "../FCComp.h" namespace "FCComps":
    set[int] isos2track

    int verbosity
