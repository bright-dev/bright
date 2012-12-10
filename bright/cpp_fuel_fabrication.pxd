"""Cython header for Fuel Fabrication library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string as std_string

from pyne cimport cpp_material
from pyne cimport material

from bright.cpp_fccomp cimport FCComp
from bright.cpp_reactor1g cimport Reactor1G


cdef extern from "fuel_fabrication.h" namespace "bright":

    cdef cppclass FuelFabrication(FCComp):
        # Constructors        
        FuelFabrication() except +
        FuelFabrication(std_string) except +
        FuelFabrication(set[std_string], std_string) except +
        FuelFabrication(map[std_string, material.matp], map[std_string, double], Reactor1G, std_string) except +
        FuelFabrication(map[std_string, material.matp], map[std_string, double], Reactor1G, set[std_string], std_string) except +

        # Attributes
        map[std_string, material.matp] materials
        map[std_string, double] mass_weights_in
        map[std_string, double] mass_weights_out
        map[std_string, double] deltaRs

        Reactor1G reactor

        # Methods
        void initialize(map[std_string, material.matp], map[std_string, double], Reactor1G) except +
        void calc_params() except +

        void calc_deltaRs() except +
        cpp_material.Material calc_core_input() except +
        void calc_mass_ratios() except +

        cpp_material.Material calc() except +
        cpp_material.Material calc(map[std_string, material.matp], map[std_string, double], Reactor1G) except +


