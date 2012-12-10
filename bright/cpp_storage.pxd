"""Cython header for  library."""
from libcpp.map cimport map
from libcpp.string cimport string as std_string

from pyne cimport cpp_material
from pyne cimport material

from bright.cpp_fccomp cimport FCComp


cdef extern from "storage.h" namespace "bright":

    cdef cppclass Storage(FCComp):
        # Constructors
        Storage() except +
        Storage(std_string) except +

        # Attributes
        double decay_time

        # Methods
        void calc_params() except +
        cpp_material.Material calc() except +
        cpp_material.Material calc(map[int, double]) except +
        cpp_material.Material calc(cpp_material.Material) except +
        cpp_material.Material calc(double) except +
        cpp_material.Material calc(map[int, double], double) except +
        cpp_material.Material calc(cpp_material.Material, double) except +


