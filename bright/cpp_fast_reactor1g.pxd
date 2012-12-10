"""Cython header for FR1G library."""
from libcpp.string cimport string as std_string

from bright.cpp_reactor1g cimport Reactor1G
from bright.cpp_reactor_parameters cimport ReactorParameters

cdef extern from "fast_reactor1g.h" namespace "bright":

    cdef cppclass FastReactor1G(Reactor1G):
        # Constructors        
        FastReactor1G() except +
        FastReactor1G(std_string, std_string) except +
        FastReactor1G(ReactorParameters, std_string) except +
        FastReactor1G(std_string, ReactorParameters, std_string) except +

        # Methods
        void calc_params() except +


