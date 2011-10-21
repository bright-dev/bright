"""Python wrapper for RMG."""
# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free

cimport numpy as np
import numpy as np

from pyne cimport std
from pyne cimport nucname

from pyne cimport stlconverters as conv
from pyne import stlconverters as conv

cimport pyne.cpp_material
cimport pyne.material
import pyne.material

cimport cpp_fccomp
cimport cpp_reactor_parameters
cimport cpp_fluence_point
cimport cpp_reactormg

from bright.reactor_parameters cimport ReactorParameters
from bright.reactor_parameters import ReactorParameters

from bright.fluence_point cimport FluencePoint
from bright.fluence_point import FluencePoint

cimport fccomp
import fccomp


#######################
### ReactorMG Class ###
#######################

cdef class ReactorMG(fccomp.FCComp):
    """Multi-Group Reactor Fuel Cycle Component Class.  Daughter of FCComp class.

    Parameters
    ----------
    reactor_parameters : ReactorParameters or None, optional 
        A special data structure that contains information on how to setup and run the reactor.
    track_params : set of str or None, optional 
        A set of strings that represents what parameter data the reactor should store and set.  
        Different reactor types may have different characteristic parameters that are of interest.
    name : str, optional 
        The name of the reactor fuel cycle component instance.

    """

    def __cinit__(self, *args, **kwargs):
        # Set property defaults
        self._chemical_form_fuel = None
        self._chemical_form_clad = None
        self._chemical_form_cool = None

        self._I = None
        self._J = None
        self._K = None

        self._K_ind = None

    def __init__(self, reactor_parameters=None, track_params=None, char * name="", *args, **kwargs):
        cdef ReactorParameters rp
        cdef std.string cpp_name = std.string(name)

        if (reactor_parameters is None) and (track_params is None):
            self._inst = new cpp_reactormg.ReactorMG(cpp_name)

        elif (reactor_parameters is None) and isinstance(track_params, set):
            self._inst = new cpp_reactormg.ReactorMG(conv.py_to_cpp_set_str(track_params), cpp_name)

        elif isinstance(reactor_parameters, ReactorParameters) and (track_params is None):
            rp = reactor_parameters
            self._inst = new cpp_reactormg.ReactorMG(<cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0], cpp_name)

        elif isinstance(reactor_parameters, ReactorParameters) and isinstance(track_params, set):
            rp = reactor_parameters
            self._inst = new cpp_reactormg.ReactorMG(<cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0], conv.py_to_cpp_set_str(track_params), cpp_name)

        else:
            if reactor_parameters is not None:
                raise TypeError("The reactor_parameters keyword must be an instance of the ReactorParameters class or None.  Got " + str(type(reactor_parameters)))

            if track_params is not None:
                raise TypeError("The track_params keyword must be a set of strings or None.  Got " + str(type(track_params)))


    #
    # Class Attributes
    #

    # ReactorMG attributes

    property B:
        """This integer is the total number of batches in the fuel management scheme.  
        B is typically indexed by b."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).B

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).B = value


    property flux:
        """The nominal flux value (float) that the library for this reactor type was generated with.  
        Used to correctly weight batch-specific fluxes."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).flux

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).flux = value





    property chemical_form_fuel:
        """This is the chemical form of fuel as a dictionary or other mapping.  Keys are 
        often strings that represent isotopes while values represent the corresponding 
        mass weights.  The heavy metal concentration by the key "IHM".  
        This will automatically fill in the nuclides in mat_feed for the "IHM" weight.  
        For example, LWRs typically use a UOX fuel form::

            ReactorMG.chemical_form_fuel = {"IHM": 1.0, "O16": 2.0}

        """
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._chemical_form_fuel is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_reactormg.ReactorMG *> self._inst).chemical_form_fuel
                self._chemical_form_fuel = proxy

            return self._chemical_form_fuel

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_fuel = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_fuel = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_fuel = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._chemical_form_fuel = None


    property chemical_form_clad:
        """This is the chemical form of cladding as a dictionary or other mapping.  
        This uses the same notation as fuel_form except that "IHM" is no longer 
        a valid key.  Cladding is often made from some zircalloy.
        """
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._chemical_form_clad is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_reactormg.ReactorMG *> self._inst).chemical_form_clad
                self._chemical_form_clad = proxy

            return self._chemical_form_clad

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_clad = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_clad = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_clad = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._chemical_form_clad = None


    property chemical_form_cool:
        """This is the chemical form of coolant as a dictionary or other mapping.  
        This uses the same notation as fuel_form except that "IHM" is no longer 
        a valid key.  The term 'coolant' is used in preference over the term 
        'moderator' because not all reactors moderate neutrons."""
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._chemical_form_cool is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_reactormg.ReactorMG *> self._inst).chemical_form_cool
                self._chemical_form_cool = proxy

            return self._chemical_form_cool

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_cool = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_cool = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_cool = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._chemical_form_cool = None




    property rho_fuel:
        """The fuel region density.  A float in units of [g/cm^3]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).rho_fuel

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).rho_fuel = value


    property rho_clad:
        """The cladding region density.  A float in units of [g/cm^3]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).rho_clad

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).rho_clad = value


    property rho_cool:
        """The coolant region density.  A float in units of [g/cm^3]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).rho_cool

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).rho_cool = value





    property P_NL:
        """The reactor's non-leakage probability (float).  This is often used as a calibration parameter."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).P_NL

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).P_NL = value


    property target_BU:
        """The reactor's target discharge burnup (float).  This is given 
        in units of [MWd/kgIHM].  Often the actual discharge burnup BUd does not 
        quite hit this value, but comes acceptably close."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).target_BU

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).target_BU = value


    property specific_power:
        """The specific power of the fuel (float) in units of [MW/kgIHM]"""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).specific_power

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).specific_power = value


    property burn_regions:
        """Number of annular burn regions (int)."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).burn_regions

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).burn_regions = value


    property S:
        """Number of burnup time steps."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).S

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).S = value


    property burn_time:
        """Curent burnup time [d]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).burn_time

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).burn_time = value


    property bt_s:
        """Curent burnup time index.  burn_time == burn_times[bt_s]"""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).bt_s

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).bt_s = value


    property burn_times:
        """A non-negative, monotonically increasing numpy float array (C++ vector<double>) of burnup times [days]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).burn_times)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).burn_times = conv.array_to_vector_1d_dbl(value)





    property use_zeta:
        """Boolaean to determine whether the thermal disadvantage factor is employed or not.  
        LWRs typically set this as True while FRs have a False value."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).use_zeta

        def __set__(self, bint value):
            (<cpp_reactormg.ReactorMG *> self._inst).use_zeta = value


    property lattice_flag:
        """Flag (str) that represents what lattice type the fuel assemblies are arranged in.  
        Currently accepted values are "Planar", "Spherical", and "Cylindrical"."""
        def __get__(self):
            cdef std.string value = (<cpp_reactormg.ReactorMG *> self._inst).lattice_flag
            return value.c_str()

        def __set__(self, char * value):
            (<cpp_reactormg.ReactorMG *> self._inst).lattice_flag = std.string(value)


    property rescale_hydrogen_xs:
        """Boolean to determine whether the reactor should rescale the Hydrogen-1 destruction 
        rate in the coolant as a function of fluence.  The scaling factor is calculated via the 
        following equation

            .. math:: f(F) = 1.36927 - 0.01119 \cdot BU(F)

        This is typically not done for fast reactors but is a useful correction for LWRs."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).rescale_hydrogen_xs

        def __set__(self, bint value):
            (<cpp_reactormg.ReactorMG *> self._inst).rescale_hydrogen_xs = value


    property burnup_via_constant:
        """Flag (str) for constant "flux" or "power" calculations."""
        def __get__(self):
            cdef std.string value = (<cpp_reactormg.ReactorMG *> self._inst).burnup_via_constant
            return value.c_str()

        def __set__(self, char * value):
            (<cpp_reactormg.ReactorMG *> self._inst).burnup_via_constant = std.string(value)


    property branch_ratio_cutoff:
        """The cutoff value (float) below which the bateman equations are not solved."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).branch_ratio_cutoff

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).branch_ratio_cutoff = value





    property r_fuel:
        """The radius (float) of the fuel region [cm]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).r_fuel

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).r_fuel = value


    property r_void:
        """The radius (float) of the void region [cm]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).r_void

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).r_void = value


    property r_clad:
        """The radius (float) of the cladding region [cm]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).r_clad

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).r_clad = value


    property pitch:
        """The pitch or length (float) of the unit fuel pin cell [cm]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).pitch

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).pitch = value





    property S_O:
        """The number of slots (float) in a fuel assembly that are open.  Thus this is the 
        number of slots that do not contain a fuel pin and are instead filled in by coolant."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).S_O

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).S_O = value


    property S_T:
        """The total number of fuel pin slots (float) in a fuel assembly.  For a 17x17 bundle 
        this is 289.0."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).S_T

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).S_T = value


    property V_fuel:
        """The relative fuel region volume."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).V_fuel

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).V_fuel = value


    property V_clad:
        """The relative cladding region volume."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).V_clad

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).V_clad = value


    property V_cool:
        """The relative coolant region volume."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).V_cool

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).V_cool = value







    property libfile:
        """The path (str) to the reactor data library; usually something like "lwr_mg.h5"."""
        def __get__(self):
            cdef std.string value = (<cpp_reactormg.ReactorMG *> self._inst).libfile
            return value.c_str()

        def __set__(self, char * value):
            (<cpp_reactormg.ReactorMG *> self._inst).libfile = std.string(value)






    property I:
        """Set of nuclides that may be in mat_feed.  Indexed by i."""
        def __get__(self):
            cdef conv._SetInt proxy

            if self._I is None:
                proxy = conv.SetInt(False, False)
                proxy.set_ptr = &(<cpp_reactormg.ReactorMG *> self._inst).I
                self._I = proxy

            return self._I

        def __set__(self, value):
            cdef cpp_set[int] s

            if isinstance(value, conv._SetInt):
                (<cpp_reactormg.ReactorMG *> self._inst).I = deref((<conv._SetInt> value).set_ptr)
            elif hasattr(value, '__len__'):
                s = cpp_set[int]()
                for nuc in value:
                    s.insert(pyne.nucname.zzaaam(nuc))
                (<cpp_reactormg.ReactorMG *> self._inst).I = s
            else:
                raise TypeError('{0} cannot be converted to a C++ set.'.format(type(value)))

            self._I = None


    property J:
        """Set of nuclides that may be in mat_prod.  Indexed by j."""
        def __get__(self):
            cdef conv._SetInt proxy

            if self._J is None:
                proxy = conv.SetInt(False, False)
                proxy.set_ptr = &(<cpp_reactormg.ReactorMG *> self._inst).J
                self._J = proxy

            return self._J

        def __set__(self, value):
            cdef cpp_set[int] s

            if isinstance(value, conv._SetInt):
                (<cpp_reactormg.ReactorMG *> self._inst).J = deref((<conv._SetInt> value).set_ptr)
            elif hasattr(value, '__len__'):
                s = cpp_set[int]()
                for nuc in value:
                    s.insert(pyne.nucname.zzaaam(nuc))
                (<cpp_reactormg.ReactorMG *> self._inst).J = s
            else:
                raise TypeError('{0} cannot be converted to a C++ set.'.format(type(value)))

            self._J = None


    property K:
        """Set of nuclides that is the union of all nucs in mat_feed and all nucs in nuc_data.
        Indexed by k."""
        def __get__(self):
            cdef conv._SetInt proxy

            if self._K is None:
                proxy = conv.SetInt(False, False)
                proxy.set_ptr = &(<cpp_reactormg.ReactorMG *> self._inst).K
                self._K = proxy

            return self._K

        def __set__(self, value):
            cdef cpp_set[int] s

            if isinstance(value, conv._SetInt):
                (<cpp_reactormg.ReactorMG *> self._inst).K = deref((<conv._SetInt> value).set_ptr)
            elif hasattr(value, '__len__'):
                s = cpp_set[int]()
                for nuc in value:
                    s.insert(pyne.nucname.zzaaam(nuc))
                (<cpp_reactormg.ReactorMG *> self._inst).K = s
            else:
                raise TypeError('{0} cannot be converted to a C++ set.'.format(type(value)))

            self._K = None







    property K_num:
        """Size of K set."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).K_num

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).K_num = value


    property K_ord:
        """Lowest-to-highest order of K."""
        def __get__(self):
            return conv.vector_to_array_1d_int((<cpp_reactormg.ReactorMG *> self._inst).K_ord)

        def __set__(self, set value):
            (<cpp_reactormg.ReactorMG *> self._inst).K_ord = conv.array_to_vector_1d_int(value)


    property K_ind:
        """Lowest-to-highest map of J into matrix position."""
        def __get__(self):
            cdef conv._MapIntInt proxy

            if self._K_ind is None:
                proxy = conv.MapIntInt(False, False)
                proxy.map_ptr = &(<cpp_reactormg.ReactorMG *> self._inst).K_ind
                self._K_ind = proxy

            return self._K_ind

        def __set__(self, value):
            cdef cpp_pair[int, int] item
            cdef cpp_map[int, int]  m

            if isinstance(value, conv._MapIntInt):
                (<cpp_reactormg.ReactorMG *> self._inst).K_ind = deref((<conv._MapIntInt> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[int, int]()
                for k, v in value.items():
                    item = cpp_pair[int, int](k, v)
                    m.insert(item)
                (<cpp_reactormg.ReactorMG *> self._inst).K_ind = m
            elif hasattr(value, '__len__'):
                m = cpp_map[int, int]()
                for i in value:
                    item = cpp_pair[int, int](i[0], i[1])
                    m.insert(item)
                (<cpp_reactormg.ReactorMG *> self._inst).K_ind = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._K_ind = None




    """\
    property decay_matrix:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).decay_matrix)

        def __set__(self, value):
            (<cpp_reactormg.ReactorMG *> self._inst).decay_matrix = conv.array_to_vector_2d_dbl(value)


    property thermal_yield_matrix:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).thermal_yield_matrix)

        def __set__(self, value):
            (<cpp_reactormg.ReactorMG *> self._inst).thermal_yield_matrix = conv.array_to_vector_2d_dbl(value)


    property fast_yield_matrix:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).fast_yield_matrix)

        def __set__(self, value):
            (<cpp_reactormg.ReactorMG *> self._inst).fast_yield_matrix = conv.array_to_vector_2d_dbl(value)


    property fission_product_yield_matrix:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).fission_product_yield_matrix)

        def __set__(self, value):
            (<cpp_reactormg.ReactorMG *> self._inst).fission_product_yield_matrix = conv.array_to_vector_3d_dbl(value)

    """






    # Perturbation table goes here


    property nperturbations:
        """Number of rows in the pertubtaion table.  Indexed by p."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).nperturbations

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).nperturbations = value


    property perturbed_fields:
        """Mapping of the form {field_name: [min, max, delta]}."""
        def __get__(self):
            return conv.map_to_dict_str_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).perturbed_fields)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).perturbed_fields = conv.dict_to_map_str_array_to_vector_1d_dbl(value)








    property G:
        """Number of energy bins.  Indexed by g for incident neutrons and h for exiting neutrons."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).G

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).G = value


    property E_g:
        """Energy bin boundaries [MeV].  Vector of doubles."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).E_g)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).E_g = conv.array_to_vector_1d_dbl(value)


    property phi_g:
        """Group fluxes [n/s/cm^2] (ie flux per energy group)."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).phi_g)

        def __set__(self, np.ndarray[np.float64_t, ndim=2] value):
            (<cpp_reactormg.ReactorMG *> self._inst).phi_g = conv.array_to_vector_2d_dbl(value)


    property phi:
        """Total flux [n/s/cm^2]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).phi)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).phi = conv.array_to_vector_1d_dbl(value)


    property Phi:
        """Fluence [n/kb]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Phi)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).Phi = conv.array_to_vector_1d_dbl(value)


    property time0:
        """Time steps used in data library [d]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).time0)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).time0 = conv.array_to_vector_1d_dbl(value)


    property BU0:
        """Burnup vector used in data library [MWd/kgIHM]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).BU0)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).BU0 = conv.array_to_vector_1d_dbl(value)






    property Ti0:
        """Data library's transmutation vector [kg_i]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Ti0)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Ti0 = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property sigma_t_pg:
        """Total cross section from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_t_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_t_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)
    
    
    property sigma_a_pg:
        """Absorption cross section from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_a_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_a_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)
    
    
    property nubar_sigma_f_pg:
        """Neutrons per fission times fission cross section from data library [n barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_sigma_f_pg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_sigma_f_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property chi_pg:
        """Fission energy spectrum from data library [MeV]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_s_pgh:
        """Group to group scattering cross section from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_s_pgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_s_pgh = conv.dict_to_map_int_array_to_vector_3d_dbl(value)


    property sigma_f_pg:
        """Fission cross section from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_f_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_f_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_gamma_pg:
        """Capture cross section from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)
    

    property sigma_2n_pg:
        """(n, 2n) cross section from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_3n_pg:
        """(n, 3n) cross section from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_3n_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_3n_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_alpha_pg:
        """(n, alpha) cross section from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_alpha_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_alpha_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_proton_pg:
        """(n, proton) cross section from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_proton_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_proton_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_gamma_x_pg:
        """Capture cross section (excited) from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_x_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_x_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)
    

    property sigma_2n_x_pg:
        """(n, 2n *) cross section from data library [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_x_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_x_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)






    property A_HM_t:
        """Atomic weight of heavy metal."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).A_HM_t)

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).A_HM_t = conv.array_to_vector_1d_dbl(value)


    property MW_fuel_t:
        """Fuel Molecular Weight [amu]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).MW_fuel_t)

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).MW_fuel_t = conv.array_to_vector_1d_dbl(value)


    property MW_clad_t:
        """Cladding Molecular Weight [amu]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).MW_clad_t)

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).MW_clad_t = conv.array_to_vector_1d_dbl(value)


    property MW_cool_t:
        """Coolant Molecular Weight [amu]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).MW_cool_t)

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).MW_cool_t = conv.array_to_vector_1d_dbl(value)


    property n_fuel_it:
        """Fuel Atom Number Weight."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).n_fuel_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).n_fuel_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property n_clad_it:
        """Cladding Atom Number Weight."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).n_clad_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).n_clad_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)

    
    property n_cool_it:
        """Coolant Atom Number Weight."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).n_cool_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).n_cool_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property m_fuel_it:
        """Fuel Mass Weight."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).m_fuel_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).m_fuel_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property m_clad_it:
        """Cladding Mass Weight."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).m_clad_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).m_clad_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)

    
    property m_cool_it:
        """Coolant Mass Weight."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).m_cool_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).m_cool_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property N_fuel_it:
        """Fuel Number Density [atoms/cm^3]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).N_fuel_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).N_fuel_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property N_clad_it:
        """Cladding Number Density [atoms/cm^3]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).N_clad_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).N_clad_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)

    
    property N_cool_it:
        """Coolant Number Density [atoms/cm^3]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).N_cool_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).N_cool_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)









    property phi_tg:
        """Group fluxes as a function of time [n/s/cm^2]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).phi_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).phi_tg = conv.array_to_vector_2d_dbl(value)


    property phi_t:
        """Total flux as a function of time [n/s/cm^2]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).phi_t)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).phi_t = conv.array_to_vector_1d_dbl(value)


    property Phi_t:
        """Fluence as a function of time [n/kb]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Phi_t)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).Phi_t = conv.array_to_vector_1d_dbl(value)


    property BU_t:
        """Burnup as a function of time [MWd/kgIHM]."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).BU_t)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).BU_t = conv.array_to_vector_1d_dbl(value)


    property zeta_tg:
        """Disadvantage factors per group as a function of time."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).zeta_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).zeta_tg = conv.array_to_vector_2d_dbl(value)







    property T_it:
        """Transformation Matrix [kg_i/kgIHM]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).T_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).T_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property sigma_t_itg:
        """Total cross section as a function of nuclide and burn_time [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_t_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_t_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_a_itg:
        """Absorption cross section as a function of nuclide and burn_time [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_a_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_a_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property nubar_sigma_f_itg:
        """Neutrons per fission times Fission cross section as a function of nuclide and burn_time [n barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_sigma_f_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_sigma_f_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property chi_itg:
        """Fission neutron energy spectrum as a function of nuclide and burn_time [MeV]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_s_itgh:
        """Group to group scattering cross section as a function of nuclide and burn_time [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_s_itgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_s_itgh = conv.dict_to_map_int_array_to_vector_3d_dbl(value)


    property sigma_f_itg:
        """Fission cross section as a function of nuclide and burn_time [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_f_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_f_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_gamma_itg:
        """Capture cross section as a function of nuclide and burn_time [barns]"""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_2n_itg:
        """(n, 2n) cross section as a function of nuclide and burn_time [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_3n_itg:
        """(n, 3n) cross section as a function of nuclide and burn_time [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_3n_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_3n_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_alpha_itg:
        """(n, alpha) cross section as a function of nuclide and burn_time [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_alpha_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_alpha_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_proton_itg:
        """(n, proton) cross section as a function of nuclide and burn_time [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_proton_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_proton_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_gamma_x_itg:
        """Capture cross section (excited) as a function of nuclide and burn_time [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_x_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_x_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_2n_x_itg:
        """(n, 2n *) cross section as a function of nuclide and burn_time [barns]."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_x_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_x_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)









    property Sigma_t_fuel_tg:
        """Fuel-averaged macroscopic total cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_a_fuel_tg:
        """Fuel-averaged macroscopic absorption cross section as a function of time and energy group [1/cm]"""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property nubar_Sigma_f_fuel_tg:
        """Fuel-averaged nubar times the macroscopic fission cross-section as a function of time and energy group [n/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property chi_fuel_tg:
        """Fuel-averaged fission neutron energy spectrum as a function of time and energy group [MeV]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_s_fuel_tgh:
        """Fuel-averaged macroscopic scattering kernel cross-section as a function of time [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_fuel_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_fuel_tgh = conv.array_to_vector_3d_dbl(value)


    property Sigma_f_fuel_tg:
        """Fuel-averaged macroscopic fission cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_fuel_tg:
        """Fuel-averaged macroscopic capture cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_fuel_tg:
        """Fuel-averaged macroscopic (n, 2n) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_3n_fuel_tg:
        """Fuel-averaged macroscopic (n, 3n) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_alpha_fuel_tg:
        """Fuel-averaged macroscopic (n, alpha) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_proton_fuel_tg:
        """Fuel-averaged macroscopic proton cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_x_fuel_tg:
        """Fuel-averaged macroscopic capture (excited) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_x_fuel_tg:
        """Fuel-averaged macroscopic (n, 2n *) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property kappa_fuel_tg:
        """Inverse of the fuel diffusion coefficent [1/cm^2]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).kappa_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).kappa_fuel_tg = conv.array_to_vector_2d_dbl(value)









    property Sigma_t_clad_tg:
        """Cladding-averaged macroscopic total cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_a_clad_tg:
        """Cladding-averaged macroscopic absorption cross section as a function of time and energy group [1/cm]"""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_clad_tg = conv.array_to_vector_2d_dbl(value)


    property nubar_Sigma_f_clad_tg:
        """Cladding-averaged nubar times the macroscopic fission cross-section as a function of time and energy group [n/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_clad_tg = conv.array_to_vector_2d_dbl(value)


    property chi_clad_tg:
        """Cladding-averaged fission neutron energy spectrum as a function of time and energy group [MeV]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_s_clad_tgh:
        """Cladding-averaged macroscopic scattering kernel cross-section as a function of time [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_clad_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_clad_tgh = conv.array_to_vector_3d_dbl(value)


    property Sigma_f_clad_tg:
        """Cladding-averaged macroscopic fission cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_clad_tg:
        """Cladding-averaged macroscopic capture cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_clad_tg:
        """Cladding-averaged macroscopic (n, 2n) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_3n_clad_tg:
        """Cldding-averaged macroscopic (n, 3n) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_alpha_clad_tg:
        """Cladding-averaged macroscopic (n, alpha) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_proton_clad_tg:
        """Cladding-averaged macroscopic proton cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_x_clad_tg:
        """Cladding-averaged macroscopic capture (excited) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_x_clad_tg:
        """Cladding-averaged macroscopic (n, 2n *) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_clad_tg = conv.array_to_vector_2d_dbl(value)


    property kappa_clad_tg:
        """Inverse of the cladding diffusion coefficent [1/cm^2]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).kappa_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).kappa_clad_tg = conv.array_to_vector_2d_dbl(value)










    property Sigma_t_cool_tg:
        """Coolant-averaged macroscopic total cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_a_cool_tg:
        """Coolant-averaged macroscopic absorption cross section as a function of time and energy group [1/cm]"""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_cool_tg = conv.array_to_vector_2d_dbl(value)


    property nubar_Sigma_f_cool_tg:
        """Coolant-averaged nubar times the macroscopic fission cross-section as a function of time and energy group [n/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_cool_tg = conv.array_to_vector_2d_dbl(value)


    property chi_cool_tg:
        """Coolant-averaged fission neutron energy spectrum as a function of time and energy group [MeV]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_s_cool_tgh:
        """Coolant-averaged macroscopic scattering kernel cross-section as a function of time [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_cool_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_cool_tgh = conv.array_to_vector_3d_dbl(value)


    property Sigma_f_cool_tg:
        """Coolant-averaged macroscopic fission cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_cool_tg:
        """Coolant-averaged macroscopic capture cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_cool_tg:
        """Coolant-averaged macroscopic (n, 2n) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_3n_cool_tg:
        """Coolant-averaged macroscopic (n, 3n) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_alpha_cool_tg:
        """Coolant-averaged macroscopic (n, alpha) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_proton_cool_tg:
        """Coolant-averaged macroscopic proton cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_x_cool_tg:
        """Coolant-averaged macroscopic capture (excited) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_x_cool_tg:
        """Coolant-averaged macroscopic (n, 2n *) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_cool_tg = conv.array_to_vector_2d_dbl(value)


    property kappa_cool_tg:
        """Inverse of the coolant diffusion coefficent [1/cm^2]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).kappa_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).kappa_cool_tg = conv.array_to_vector_2d_dbl(value)












    property Sigma_t_tg:
        """Core-averaged macroscopic total cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_a_tg:
        """Core-averaged macroscopic absorption cross section as a function of time and energy group [1/cm]"""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_tg = conv.array_to_vector_2d_dbl(value)


    property nubar_Sigma_f_tg:
        """Core-averaged nubar times the macroscopic fission cross-section as a function of time and energy group [n/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_tg = conv.array_to_vector_2d_dbl(value)


    property chi_tg:
        """Core-averaged fission neutron energy spectrum as a function of time and energy group [MeV]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_s_tgh:
        """Core-averaged macroscopic scattering kernel cross-section as a function of time [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_tgh = conv.array_to_vector_3d_dbl(value)


    property Sigma_f_tg:
        """Core-averaged macroscopic fission cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_tg:
        """Core-averaged macroscopic capture cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_tg:
        """Core-averaged macroscopic (n, 2n) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_3n_tg:
        """Core-averaged macroscopic (n, 3n) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_alpha_tg:
        """Core-averaged macroscopic (n, alpha) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_proton_tg:
        """Core-averaged macroscopic proton cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_x_tg:
        """Core-averaged macroscopic capture (excited) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_x_tg:
        """Core-averaged macroscopic (n, 2n *) cross-section as a function of time and energy group [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_tg = conv.array_to_vector_2d_dbl(value)









    property A_tgh:
        """Absorprion matrix, as a function of time [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).A_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).A_tgh = conv.array_to_vector_3d_dbl(value)


    property F_tgh:
        """Fission Matrix, as a function of time [1/cm]."""
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).F_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).F_tgh = conv.array_to_vector_3d_dbl(value)


    property A_inv_tgh:
        """Inverse of absorprion matrix, as a function of time [cm]."""
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).A_inv_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).A_inv_tgh = conv.array_to_vector_3d_dbl(value)


    property A_inv_F_tgh:
        """Inverse of absorprion matrix mult by the fission matrix, as a function of time."""
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).A_inv_F_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).A_inv_F_tgh = conv.array_to_vector_3d_dbl(value)

    """\
    property T_int_tij:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).T_int_tij)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).T_int_tij = conv.array_to_vector_3d_dbl(value)


    property M_tij:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).M_tij)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).M_tij = conv.array_to_vector_3d_dbl(value)

    """










    property nearest_neighbors:
        """Attribute that denotes the indices of the perturbation table which are closest to 
        the current state of the reactor (ie densities, burn_time, etc)."""
        def __get__(self):
            return conv.vector_to_array_1d_int((<cpp_reactormg.ReactorMG *> self._inst).nearest_neighbors)

        def __set__(self, np.ndarray[np.int32_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).nearest_neighbors = conv.array_to_vector_1d_int(value)






    property k_t:
        """Multiplication factor of the core as a function of time."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).k_t)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).k_t = conv.array_to_vector_1d_dbl(value)





    property td_n:
        """Lower index of discharge time."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).td_n

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).td_n = value


    property td:
        """Discharge time [days]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).td

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).td = value


    property BUd:
        """Discharge Burnup [MWd/kgIHM]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).BUd

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).BUd = value


    property Phid:
        """Discharge Fluence [n/kb]."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).Phid

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).Phid = value


    property k:
        """Core multiplication factor."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).k

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).k = value






    property mat_feed_u:
        """The input uranium material, mat_feed.sub_u()."""
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_u
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_u = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_feed_tru:
        """The input transuranic material, mat_feed.sub_tru()."""
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_tru
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_tru = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_feed_lan:
        """The input lanthanide material, mat_feed.sub_lan()."""
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_lan
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_lan = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_feed_act:
        """The input actininide material, mat_feed.sub_act()."""
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_act
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_act = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_u:
        """The output urnaium material, mat_prod.sub_u()."""
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_u
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_u = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_tru:
        """The output transuranic material, mat_prod.sub_tru()."""
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_tru
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_tru = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_lan:
        """The output lanthanide material, mat_prod.sub_lan()."""
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_lan
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_lan = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_act:
        """The output actininide material, mat_prod.sub_act()."""
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_act
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_act = <pyne.cpp_material.Material> mat.mat_pointer[0]





    property deltaR:
        """The :math:`\delta R` value of the core with the current mat_feed.  This is equal 
        to the production rate minus the destruction rate at the target burnup::

            deltaR = batch_average(target_BU, "P") - batch_average(target_BU, "D")

        This is computed via the calc_deltaR() method."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).deltaR

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).deltaR = value


    property tru_cr:
        """The transuranic conversion ratio of the reactor (float).  
        This is set via the calc_tru_cr() method."""
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).tru_cr

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).tru_cr = value





    #
    # Class Methods
    # 

    def initialize(self, ReactorParameters reactor_parameters):
        """The initialize() method for reactors copies all of the reactor specific parameters to this instance.
        Additionally, it calculates and sets the volumes VF and VC.

        Parameters
        ----------
        reactor_parameters : ReactorParameters 
            A special data structure that contains information on how to setup and run the reactor.

        """
        cdef ReactorParameters rp = reactor_parameters
        (<cpp_reactormg.ReactorMG *> self._inst).initialize(<cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0])


    def loadlib(self, char * libfile="reactor.h5"):
        """This method finds the HDF5 library for this reactor and extracts the necessary information from it.
        This method is typically called by the constructor of the child reactor type object.  It must be 
        called before attempting to do any real computation.

        Parameters
        ----------
        libfile : str, optional 
            Path to the reactor library.

        """
        (<cpp_reactormg.ReactorMG *> self._inst).loadlib(std.string(libfile))


    def interpolate_cross_sections(self):
        """This method iterpolates the isotopic, time-dependent cross-sections based on the current 
        state of the burn_time, bt_s, and nearest_neighbors attributes.  It is prudent to call 
        the calc_nearest_neighbors() method before this one.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).interpolate_cross_sections()


    def calc_mass_weights(self):
        """Calculates the mass weights for this time step.  Needed for fold_mass_weights() method.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_mass_weights()


    def fold_mass_weights(self):
        """This method performs the all-important task of doing the isotopically-weighted linear combination of raw data. 
        In a very real sense this is what makes this reactor *this specific reactor*.  The weights are taken 
        as the values of mat_feed.  The raw data must have previously been read in from loadlib().  

        Warnings
        --------
        Anytime any reactor parameter whatsoever (mat_feed, P_NL, *etc*) is altered in any way, 
        the fold_mass_weights() function must be called to reset all of the resultant data.
        If you are unsure, please call this function anyway to be safe.  There is little 
        harm in calling it twice by accident.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).fold_mass_weights()


    def assemble_multigroup_matrices(self):
        """Folds mass weight in with cross-sections for current time step.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).assemble_multigroup_matrices()


    def assemble_transmutation_matrices(self):
        """Calculates the transmutation matrices for the current time step
        """
        (<cpp_reactormg.ReactorMG *> self._inst).assemble_transmutation_matrices()


    def calc_criticality(self):
        """Performs the criticality calculation to find k for this time step.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_criticality()


    def calc_transmutation(self):
        """Burns up the fuel using a matrix exponential method for this time step.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_transmutation()








    def init_core(self):
        """This method generates a time-dependent parameters from an reactor's initial conditions.
        This includes all burnup and criticality calculations.  These time-dependent data
        are then used to determine discharge compositions and other parameters.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).init_core()


    def burnup_core(self):
        """This method generates a time-dependent parameters from an reactor's initial conditions.
        This includes all burnup and criticality calculations.  These time-dependent data
        are then used to determine discharge compositions and other parameters.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).burnup_core()






    def calc_nearest_neighbors(self):
        """Calculates a sorted array that indexes the nearest neighbors of the 
        perturbations based off of the current state of the reactor.  The results may
        be found in the neareest_neighbors attribute.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_nearest_neighbors()







    def calc_T_itd(self):
        """This function evaluates transmutation matrix at the discharge time td.
        The resultant isotopic dictionary is then converted into the mat_prod mass stream
        for this pass through the reactor.  Thus if ever you need to calculate mat_prod
        without going through calc(), use this function.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_T_itd()






    def calc_mat_prod(self):
        """This is a convenience function that wraps the transmutation matrix methods.  
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_mat_prod()

    def calc_sub_mats(self):
        """This sets possibly relevant reactor input and output substreams.  Specifically, it calculates the 
        attributes:

            * mat_feed_u
            * mat_feed_tru
            * mat_feed_lan
            * mat_feed_act
            * mat_prod_u
            * mat_prod_tru
            * mat_prod_lan
            * mat_prod_act

        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_sub_mats()


    def calc_tru_cr(self):
        """This calculates and sets the transuranic conversion ratio tru_cr through the equation:

        .. math:: \\mbox{tru\_cr} = \\frac{\\mbox{mat\_feed\_tru.mass} - \\mbox{mat\_prod\_tru.mass}}{\\frac{\\mbox{BUd}}{935.0}}

        Returns
        -------
        tru_cr : float 
            The value of the transuranic conversion ratio just calculated.

        """
        return (<cpp_reactormg.ReactorMG *> self._inst).calc_tru_cr()







    def fluence_at_BU(self, double burnup):
        """This function takes a burnup value  and returns a special fluence point object.  
        The fluence point is an amalgamation of data where the at which the burnup occurs.
        This object instance FP contains three pieces of information::
    
            FP.f    #Index immediately lower than where BU achieved (int)
            FP.F    #Fluence value itself (float)
            FP.m    #Slope dBU/dF between points f and f+1 (double)

        Parameters
        ----------
        burnup : float 
            Burnup [MWd/kgIHM] at which to calculate the corresponding fluence.

        Returns
        -------
        fp : FluencePoint 
            A class containing fluence information.

        """
        cdef FluencePoint fp = FluencePoint()
        fp.fp_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).fluence_at_BU(burnup)
        return fp


    def batch_average_k(self, double BUd):
        """Finds the batch average k(F) when at discharge burnup BUd.
        This function is typically iterated over until a BUd is found such that k(F) = 1.0 + err.

        Parameters
        ----------
        BUd : float 
            The discharge burnup [MWd/kgIHM] to obtain a batch-averaged value for.

        Returns
        -------
        k : float 
            The batch averaged multiplication factor.

        """
        cdef double k = (<cpp_reactormg.ReactorMG *> self._inst).batch_average_k(BUd)
        return k


    def BUd_bisection_method(self):
        """Calculates the maximum discharge burnup via the Bisection Method for a given mat_feed
        in this reactor.  This iterates over values of BUd to find a batch averaged multiplication factor 
        that is closest to 1.0.

        Other root finding methods for determining maximum discharge burnup are certainly possible.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).BUd_bisection_method()


    def run_P_NL(self, double pnl):
        """Performs a reactor run for a specific non-leakage probability value.
        This requires that mat_feed be (meaningfully) set and is for use with calibrate_P_NL_to_BUd().

        This function amounts to the following code::

            self.P_NL = pnl
            self.fold_mass_weights()
            self.BUd_bisection_method()

        Parameters
        ----------
        pnl : float 
            The new non-leakage probability for the reactor.

        """
        (<cpp_reactormg.ReactorMG *> self._inst).run_P_NL(pnl)
    

    def calibrate_P_NL_to_BUd(self):
        """Often times the non-leakage probability of a reactor is not known, though the input isotopics 
        and the target discharge burnup are.  This function handles that situation by
        calibrating the non-leakage probability of this reactor P_NL to hit its target burnup target_BU.
        Such a calibration proceeds by bisection method as well.  This function is extremely useful for 
        benchmarking calculations.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calibrate_P_NL_to_BUd()




    def calc(self, input=None):
        """Since many other methods provide the computational heavy-lifting of reactor calculations, 
        the calc() method is relatively simple::

            self.mat_feed = input
            self.fold_mass_weights()
            self.BUd_bisection_method()
            self.calc_mat_prod()
            return self.mat_prod

        As you can see, all this function does is set burn an input stream to its maximum 
        discharge burnup and then reports on the output isotopics.

        Parameters
        ----------
        input : dict or Material or None, optional 
            If input is present, it set as the component's mat_feed.  If input is a nuclide 
            dictionary (zzaaam keys, float values), this dictionary is first converted into 
            a Material before being set as mat_feed.

        Returns
        -------
        output : Material 
            mat_prod

        """
        cdef pyne.material._Material in_mat 
        cdef pyne.material._Material output = pyne.material.Material()

        if input is None:
            output.mat_pointer[0] = (<cpp_reactormg.FCComp *> self._inst).calc()
        elif isinstance(input, dict):
            output.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).calc(conv.dict_to_map_int_dbl(input))
        elif isinstance(input, pyne.material._Material):
            in_mat = input
            output.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).calc(<pyne.cpp_material.Material> in_mat.mat_pointer[0])

        return output
