"""Python wrapper for  ."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.vector cimport vector as cpp_vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free

from pyne cimport std
from pyne cimport stlconverters as conv

cimport numpy as np
import numpy as np

cimport cpp_reactor_parameters



###############################
### ReactorParameters Class ###
###############################

cdef class ReactorParameters:
    """This data structure is a set of physical reactor parameters. It may be used to instantiate new reactor objects **OR**
    to define default settings for a reactor type.  The data stored in this class is copied over to 
    a reactor instance in the :meth:`Reactor1G.initialize` method.  However, the attributes of this objects 
    take on more natural names than their :class:`Reactor1G` analogies.  This is because it is this 
    object that Bright users will more often be interacting with. 

    Attributes:
        * batches (int): This is the total number of batches in the fuel management scheme. 
          This is typically indexed by b.
        * flux (float): The nominal flux value that the library for this reactor 
          type was generated with.  Used to correctly weight batch-specific fluxes.
        * fuel_form (dict): This is the chemical form of fuel as dictionary.  Keys are 
          strings that represent isotopes (mixed form) while values represent the 
          corresponding mass weights.  The heavy metal concentration by the key "IHM".  
          This will automatically fill in the nuclides in mat_feed for the "IHM" weight.  
          For example, LWRs typically use a UOX fuel form::

            ReactorParameters.fuel_form = {"IHM": 1.0, "O16": 2.0}

        * coolant_form (dict): This is the chemical form of coolant as dictionary.  
          This uses the same notation as fuel_form except that "IHM" is no longer 
          a valid key.  The term 'coolant' is used in preference over the term 
          'moderator' because not all reactors moderate neutrons.  For example, 
          LWRs often cool the reactor core with borated water::

            ReactorParamters.coolant_form = {}

            ReactorParamters.coolant_form["H1"]  = 2.0
            ReactorParamters.coolant_form["O16"] = 1.0
            ReactorParamters.coolant_form["B10"] = 0.199 * 550 * 10.0**-6
            ReactorParamters.coolant_form["B11"] = 0.801 * 550 * 10.0**-6

        * fuel_density (float): The fuel region density.  A float in units of [g/cm^3].
        * coolant_density (float): The coolant region density.  A float in units of [g/cm^3].
        * pnl (float): The reactor's non-leakage probability.  This is often 
          used as a calibration parameter.
        * BUt (float): The reactor's target discharge burnup.  This is given 
          in units of [MWd/kgIHM].  Often the actual discharge burnup BUd does not 
          quite hit this value, but comes acceptably close.
        * use_disadvantage_factor (bool): Determines whether the thermal disadvantage 
          factor is employed or not.  LWRs typically set this as True while FRs 
          have a False value.
        * lattice_type (str): A flag that represents what lattice type the fuel 
          assemblies are arranged in.  Currently accepted values are "Planar", 
          "Spherical", and "Cylindrical".
        * rescale_hydrogen (bool): This determines whether the reactor should 
          rescale the Hydrogen-1 destruction rate in the coolant as a
          function of fluence.  The scaling factor is calculated via the 
          following equation

            .. math:: f(F) = 1.36927 - 0.01119 \cdot BU(F)

          This is typically not done for fast reactors but is a useful correction 
          for LWRs.
        * burnup_via_constant (str): Flag for constant "flux" or "power" calculations.
        * branch_ratio_cutoff (float): the cutoff value below which the bateman equations
          are not solved.
        * radius (float): The radius of the fuel region.  In units of [cm].
        * pitch (float): The pitch or length of the unit fuel pin cell.  In units of [cm].
        * open_slots (float): The number of slots in a fuel assembly that are open.  
          Thus this is the number of slots that do not contain a fuel pin and are instead 
          filled in by coolant. 
        * total_slots (float): The total number of fuel pin slots in a fuel assembly.  
          For a 17x17 bundle, S_T is 289.0. 
    """

    def __cinit__(self):
        self.rp_pointer = new cpp_reactor_parameters.ReactorParameters()

    def __dealloc__(self):
        del self.rp_pointer


    #
    # Class Attributes
    #

    property batches:
        def __get__(self):
            return self.rp_pointer.batches

        def __set__(self, int value):
            self.rp_pointer.batches = value


    property flux:
        def __get__(self):
            return self.rp_pointer.flux

        def __set__(self, double value):
            self.rp_pointer.flux = value





    property fuel_form:
        def __get__(self):
            return conv.map_to_dict_str_dbl(self.rp_pointer.fuel_form)

        def __set__(self, dict value):
            self.rp_pointer.fuel_form = conv.dict_to_map_str_dbl(value)


    property cladding_form:
        def __get__(self):
            return conv.map_to_dict_str_dbl(self.rp_pointer.cladding_form)

        def __set__(self, dict value):
            self.rp_pointer.cladding_form = conv.dict_to_map_str_dbl(value)


    property coolant_form:
        def __get__(self):
            return conv.map_to_dict_str_dbl(self.rp_pointer.coolant_form)

        def __set__(self, dict value):
            self.rp_pointer.coolant_form = conv.dict_to_map_str_dbl(value)





    property fuel_density:
        def __get__(self):
            return self.rp_pointer.fuel_density

        def __set__(self, double value):
            self.rp_pointer.fuel_density = value


    property cladding_density:
        def __get__(self):
            return self.rp_pointer.cladding_density

        def __set__(self, double value):
            self.rp_pointer.cladding_density = value


    property coolant_density:
        def __get__(self):
            return self.rp_pointer.coolant_density

        def __set__(self, double value):
            self.rp_pointer.coolant_density = value





    property pnl:
        def __get__(self):
            return self.rp_pointer.pnl

        def __set__(self, double value):
            self.rp_pointer.pnl = value


    property BUt:
        def __get__(self):
            return self.rp_pointer.BUt

        def __set__(self, double value):
            self.rp_pointer.BUt = value


    property specific_power:
        def __get__(self):
            return self.rp_pointer.specific_power

        def __set__(self, double value):
            self.rp_pointer.specific_power = value


    property burn_regions:
        def __get__(self):
            return self.rp_pointer.burn_regions

        def __set__(self, int value):
            self.rp_pointer.burn_regions = value


    property burn_times:
        def __get__(self):
            return conv.vector_to_array_1d_dbl(self.rp_pointer.burn_times)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            self.rp_pointer.burn_times = conv.array_to_vector_1d_dbl(value)





    property use_disadvantage_factor:
        def __get__(self):
            return self.rp_pointer.use_disadvantage_factor

        def __set__(self, bint value):
            self.rp_pointer.use_disadvantage_factor = value


    property lattice_type:
        def __get__(self):
            cdef std.string value = self.rp_pointer.lattice_type
            return value.c_str()

        def __set__(self, char * value):
            self.rp_pointer.lattice_type = std.string(value)


    property rescale_hydrogen:
        def __get__(self):
            return self.rp_pointer.rescale_hydrogen

        def __set__(self, bint value):
            self.rp_pointer.rescale_hydrogen = value


    property burnup_via_constant:
        def __get__(self):
            cdef std.string value = self.rp_pointer.burnup_via_constant
            return value.c_str()

        def __set__(self, char * value):
            self.rp_pointer.burnup_via_constant = std.string(value)


    property branch_ratio_cutoff:
        def __get__(self):
            return self.rp_pointer.branch_ratio_cutoff

        def __set__(self, double value):
            self.rp_pointer.branch_ratio_cutoff = value







    property fuel_radius:
        def __get__(self):
            return self.rp_pointer.fuel_radius

        def __set__(self, double value):
            self.rp_pointer.fuel_radius = value


    property void_radius:
        def __get__(self):
            return self.rp_pointer.void_radius

        def __set__(self, double value):
            self.rp_pointer.void_radius = value


    property clad_radius:
        def __get__(self):
            return self.rp_pointer.clad_radius

        def __set__(self, double value):
            self.rp_pointer.clad_radius = value


    property unit_cell_pitch:
        def __get__(self):
            return self.rp_pointer.unit_cell_pitch

        def __set__(self, double value):
            self.rp_pointer.unit_cell_pitch = value





    property open_slots:
        def __get__(self):
            return self.rp_pointer.open_slots

        def __set__(self, double value):
            self.rp_pointer.open_slots = value


    property total_slots:
        def __get__(self):
            return self.rp_pointer.total_slots

        def __set__(self, double value):
            self.rp_pointer.total_slots = value




##################################
### Reactor Parameter Defaults ###
##################################


def lwr_defaults():
    cdef cpp_reactor_parameters.ReactorParameters cpp_lwrd = cpp_reactor_parameters.fill_lwr_defaults()
    cdef ReactorParameters lwrd = ReactorParameters()
    lwrd.rp_pointer[0] = cpp_lwrd
    return lwrd



def fr_defaults():
    cdef cpp_reactor_parameters.ReactorParameters cpp_frd = cpp_reactor_parameters.fill_fr_defaults()
    cdef ReactorParameters frd = ReactorParameters()
    frd.rp_pointer[0] = cpp_frd
    return frd


