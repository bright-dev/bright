"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std


cdef extern from "../MassStream.h":
    cdef cppclass MassStream:
        # Constuctors
        MassStream()
        MassStream(map[int, double], float, std.string)
        MassStream(char *, float, std.string)

        # Attributes
        map[int, double] comp
        double mass
        std.string name
