"""Cython header for FR1G library."""
from pyne cimport std

from bright.cpp_reactor1g cimport Reactor1G
from bright.cpp_reactor_parameters cimport ReactorParameters

cdef extern from "fast_reactor1g.h" namespace "bright":

    cdef cppclass FastReactor1G(Reactor1G):
        # Constructors        
        FastReactor1G() except +
        FastReactor1G(std.string, std.string) except +
        FastReactor1G(ReactorParameters, std.string) except +
        FastReactor1G(std.string, ReactorParameters, std.string) except +

        # Methods
        void calc_params() except +


