"""Python wrapper for reactor parameters."""
# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.map cimport map as cpp_map
from libcpp.vector cimport vector as cpp_vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free

from pyne cimport std

from pyne cimport stlconverters as conv
from pyne import stlconverters as conv

cimport numpy as np
import numpy as np

cimport cpp_reactor_parameters



###############################
### ReactorParameters Class ###
###############################

cdef class ReactorParameters:
    """This data structure is a set of physical reactor parameters. It may be used to instantiate new reactor 
    objects **OR** to define default settings for a reactor type.  The data stored in this class is copied over 
    to a reactor instance in the initialize() method.  However, the attributes of this objects 
    take on more natural names than their reactor attribute analogies.  This is because it is this 
    object that Bright users will more often be interacting with. 
    """

    def __cinit__(self):
        self.rp_pointer = new cpp_reactor_parameters.ReactorParameters()
        self._fuel_form = None
        self._cladding_form = None
        self._coolant_form = None

    def __dealloc__(self):
        del self.rp_pointer


    #
    # Class Attributes
    #

    property batches:
        """This is the total number of batches (int) in the fuel management scheme. 
        This is typically indexed by b."""
        def __get__(self):
            return self.rp_pointer.batches

        def __set__(self, int value):
            self.rp_pointer.batches = value


    property flux:
        """The nominal flux value (float) that the library for this reactor type was generated with.  
        Often used to correctly weight batch-specific fluxes."""
        def __get__(self):
            return self.rp_pointer.flux

        def __set__(self, double value):
            self.rp_pointer.flux = value





    property fuel_form:
        """This is the chemical form of fuel as a dictionary or other mapping.  Keys are 
        often strings that represent isotopes while values represent the corresponding 
        mass weights.  The heavy metal concentration by the key "IHM".  
        This will automatically fill in the nuclides in mat_feed for the "IHM" weight.  
        For example, LWRs typically use a UOX fuel form::

            ReactorParameters.fuel_form = {"IHM": 1.0, "O16": 2.0}

        """
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._fuel_form is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).fuel_form
                self._fuel_form = proxy

            return self._fuel_form

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).fuel_form = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).fuel_form = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).fuel_form = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._fuel_form = None


    property cladding_form:
        """This is the chemical form of cladding as a dictionary or other mapping.  
        This uses the same notation as fuel_form except that "IHM" is no longer 
        a valid key.  Cladding is often made from some zircalloy.
        """
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._cladding_form is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).cladding_form
                self._cladding_form = proxy

            return self._cladding_form

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).cladding_form = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).cladding_form = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).cladding_form = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._cladding_form = None


    property coolant_form:
        """This is the chemical form of coolant as a dictionary or other mapping.  
        This uses the same notation as fuel_form except that "IHM" is no longer 
        a valid key.  The term 'coolant' is used in preference over the term 
        'moderator' because not all reactors moderate neutrons.  For example, 
        LWRs often cool the reactor core with borated water::

            ReactorParamters.coolant_form = {}

            ReactorParamters.coolant_form["H1"]  = 2.0
            ReactorParamters.coolant_form["O16"] = 1.0
            ReactorParamters.coolant_form["B10"] = 0.199 * 550 * 10.0**-6
            ReactorParamters.coolant_form["B11"] = 0.801 * 550 * 10.0**-6

        """
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._coolant_form is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).coolant_form
                self._coolant_form = proxy

            return self._coolant_form

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).coolant_form = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).coolant_form = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_reactor_parameters.ReactorParameters *> self.rp_pointer).coolant_form = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._coolant_form = None





    property fuel_density:
        """The fuel region density.  A float in units of [g/cm^3]."""
        def __get__(self):
            return self.rp_pointer.fuel_density

        def __set__(self, double value):
            self.rp_pointer.fuel_density = value


    property cladding_density:
        """The cladding region density.  A float in units of [g/cm^3]."""
        def __get__(self):
            return self.rp_pointer.cladding_density

        def __set__(self, double value):
            self.rp_pointer.cladding_density = value


    property coolant_density:
        """The coolant region density.  A float in units of [g/cm^3]."""
        def __get__(self):
            return self.rp_pointer.coolant_density

        def __set__(self, double value):
            self.rp_pointer.coolant_density = value





    property pnl:
        """The reactor's non-leakage probability (float).  This is often used as a calibration parameter."""
        def __get__(self):
            return self.rp_pointer.pnl

        def __set__(self, double value):
            self.rp_pointer.pnl = value


    property BUt:
        """The reactor's target discharge burnup (float).  This is given 
        in units of [MWd/kgIHM].  Often the actual discharge burnup BUd does not 
        quite hit this value, but comes acceptably close."""
        def __get__(self):
            return self.rp_pointer.BUt

        def __set__(self, double value):
            self.rp_pointer.BUt = value


    property specific_power:
        """The specific power of the fuel (float) in units of [MW/kgIHM]"""
        def __get__(self):
            return self.rp_pointer.specific_power

        def __set__(self, double value):
            self.rp_pointer.specific_power = value


    property burn_regions:
        """Number of annular burn regions (int)."""
        def __get__(self):
            return self.rp_pointer.burn_regions

        def __set__(self, int value):
            self.rp_pointer.burn_regions = value


    property burn_times:
        """A non-negative, monotonically increasing numpy float array (C++ vector<double>) of burnup times [days]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl(self.rp_pointer.burn_times)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            self.rp_pointer.burn_times = conv.array_to_vector_1d_dbl(value)





    property use_disadvantage_factor:
        """Boolaean to determine whether the thermal disadvantage factor is employed or not.  
        LWRs typically set this as True while FRs have a False value."""
        def __get__(self):
            return self.rp_pointer.use_disadvantage_factor

        def __set__(self, bint value):
            self.rp_pointer.use_disadvantage_factor = value


    property lattice_type:
        """Flag (str) that represents what lattice type the fuel assemblies are arranged in.  
        Currently accepted values are "Planar", "Spherical", and "Cylindrical"."""
        def __get__(self):
            cdef std.string value = self.rp_pointer.lattice_type
            return value.c_str()

        def __set__(self, char * value):
            self.rp_pointer.lattice_type = std.string(value)


    property rescale_hydrogen:
        """Boolean to determine whether the reactor should rescale the Hydrogen-1 destruction 
        rate in the coolant as a function of fluence.  The scaling factor is calculated via the 
        following equation

            .. math:: f(F) = 1.36927 - 0.01119 \cdot BU(F)

        This is typically not done for fast reactors but is a useful correction for LWRs."""
        def __get__(self):
            return self.rp_pointer.rescale_hydrogen

        def __set__(self, bint value):
            self.rp_pointer.rescale_hydrogen = value


    property burnup_via_constant:
        """Flag (str) for constant "flux" or "power" calculations."""
        def __get__(self):
            cdef std.string value = self.rp_pointer.burnup_via_constant
            return value.c_str()

        def __set__(self, char * value):
            self.rp_pointer.burnup_via_constant = std.string(value)


    property branch_ratio_cutoff:
        """The cutoff value (float) below which the bateman equations are not solved."""
        def __get__(self):
            return self.rp_pointer.branch_ratio_cutoff

        def __set__(self, double value):
            self.rp_pointer.branch_ratio_cutoff = value







    property fuel_radius:
        """The radius (float) of the fuel region [cm]."""
        def __get__(self):
            return self.rp_pointer.fuel_radius

        def __set__(self, double value):
            self.rp_pointer.fuel_radius = value


    property void_radius:
        """The radius (float) of the void region [cm]."""
        def __get__(self):
            return self.rp_pointer.void_radius

        def __set__(self, double value):
            self.rp_pointer.void_radius = value


    property clad_radius:
        """The radius (float) of the cladding region [cm]."""
        def __get__(self):
            return self.rp_pointer.clad_radius

        def __set__(self, double value):
            self.rp_pointer.clad_radius = value


    property unit_cell_pitch:
        """The pitch or length (float) of the unit fuel pin cell [cm]."""
        def __get__(self):
            return self.rp_pointer.unit_cell_pitch

        def __set__(self, double value):
            self.rp_pointer.unit_cell_pitch = value





    property open_slots:
        """The number of slots (float) in a fuel assembly that are open.  Thus this is the 
        number of slots that do not contain a fuel pin and are instead filled in by coolant."""
        def __get__(self):
            return self.rp_pointer.open_slots

        def __set__(self, double value):
            self.rp_pointer.open_slots = value


    property total_slots:
        """The total number of fuel pin slots (float) in a fuel assembly.  For a 17x17 bundle 
        this is 289.0."""
        def __get__(self):
            return self.rp_pointer.total_slots

        def __set__(self, double value):
            self.rp_pointer.total_slots = value




##################################
### Reactor Parameter Defaults ###
##################################


def lwr_defaults():
    """This function returns a copy of the LWR default presets. These are applicable to most cases.
    However, if you want to use your own LWR parameters, it is recommended you use this function
    and then only change the necessary attributes.  

    Returns
    -------
    lwrd : ReactorParameters 
        Light water reactor default parameters.

    Warnings
    --------
    Note that the target burnup default value is zero.  Generally, at least this value should be overridden.
    """
    cdef cpp_reactor_parameters.ReactorParameters cpp_lwrd = cpp_reactor_parameters.fill_lwr_defaults()
    cdef ReactorParameters lwrd = ReactorParameters()
    lwrd.rp_pointer[0] = cpp_lwrd
    return lwrd



def fr_defaults():
    """This function returns a copy of the FR default presets. These are applicable to most cases.
    However, if you want to use your own FR parameters, it is recommended you use this function
    and then only change the necessary attributes.  

    Returns
    -------
    frd : ReactorParameters 
        Fast reactor default parameters.

    Warnings
    --------
    Note that the target burnup default value is zero.  Generally, at least this value should be overridden.
    """
    cdef cpp_reactor_parameters.ReactorParameters cpp_frd = cpp_reactor_parameters.fill_fr_defaults()
    cdef ReactorParameters frd = ReactorParameters()
    frd.rp_pointer[0] = cpp_frd
    return frd


