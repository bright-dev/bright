"""Cython header for  library."""
from libcpp.map cimport map
from libcpp.vector cimport vector

from pyne cimport std


cdef extern from "reactor_parameters.h" namespace "bright":

    cdef cppclass ReactorParameters:
        # Constructors        
        ReactorParameters() except +

        # Attributes
        int batches
        double flux

        map[std.string, double] fuel_form
        map[std.string, double] cladding_form
        map[std.string, double] coolant_form

        double fuel_density
        double cladding_density
        double coolant_density

        double pnl
        double BUt
        double specific_power
        int burn_regions
        vector[double] burn_times

        bint use_disadvantage_factor
        std.string lattice_type
        bint rescale_hydrogen
        std.string burnup_via_constant
        double branch_ratio_cutoff

        double fuel_radius
        double void_radius
        double clad_radius
        double unit_cell_pitch

        double open_slots
        double total_slots

    ReactorParameters fill_lwr_defaults() except +

    ReactorParameters fill_fr_defaults() except +

