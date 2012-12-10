"""Cython header for fccomp library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector
from libcpp.string cimport string as std_string

from pyne cimport cpp_material
from pyne cimport material


cdef extern from "fccomp.h" namespace "bright":
    cdef cppclass FCComp:
        # Constructors
        FCComp() except +
        FCComp(std_string) except +
        FCComp(set[std_string], std_string) except +

        # Attributes
        std_string name 
        std_string natural_name

        cpp_material.Material mat_feed
        cpp_material.Material mat_prod

        map[std_string, double] params_prior_calc
        map[std_string, double] params_after_calc

        int pass_num
        set[std_string] track_params

        # Methods
        void calc_params() except +
        void write_mat_pass() except +
        void write_params_pass() except +
        void write_text() except +
        void write_hdf5() except +
        void write() except +
        cpp_material.Material calc() except +


