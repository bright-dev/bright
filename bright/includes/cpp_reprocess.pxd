"""Cython header for reprocess library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

from pyne cimport std
from pyne cimport cpp_material
from pyne cimport material

from bright.cpp_fccomp cimport FCComp

cdef extern from "reprocess.h" namespace "bright":

    cdef cppclass Reprocess(FCComp):
        # Constructors
        Reprocess() except +
        Reprocess(map[int, double], std.string) except +

        # Attributes
        map[int, double] sepeff

        # Methods
        void initialize(map[int, double]) except +
        void calc_params() except +
        cpp_material.Material calc() except +
        cpp_material.Material calc(map[int, double]) except +
        cpp_material.Material calc(cpp_material.Material) except +


