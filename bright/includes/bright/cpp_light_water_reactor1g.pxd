"""Cython header for LWR1G library."""
from pyne cimport std

from bright.cpp_reactor1g cimport Reactor1G
from bright.cpp_reactor_parameters cimport ReactorParameters


cdef extern from "light_water_reactor1g.h" namespace "bright":

    cdef cppclass LightWaterReactor1G(Reactor1G):
        # Constructors        
        LightWaterReactor1G() except +
        LightWaterReactor1G(std.string, std.string) except +
        LightWaterReactor1G(ReactorParameters, std.string) except +
        LightWaterReactor1G(std.string, ReactorParameters, std.string) except +

        # Methods
        void calc_params() except +




