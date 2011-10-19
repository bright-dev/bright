"""Python wrapper for reactor1g."""
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
cimport cpp_reactor1g

from bright.reactor_parameters cimport ReactorParameters
from bright.reactor_parameters import ReactorParameters

from bright.fluence_point cimport FluencePoint
from bright.fluence_point import FluencePoint

cimport fccomp
import fccomp


#######################
### Reactor1G Class ###
#######################


cdef class Reactor1G(fccomp.FCComp):
    """One-Group Reactor Fuel Cycle Component Class.  Daughter of FCComp class.

    Parameters
    ----------
    reactor_parameters : ReactorParameters, optional 
        A special data structure that contains information on how to setup and run the reactor.
    track_params : string set, optional 
        A set of strings that represents what parameter data the reactor should store and set.  
        Different reactor types may have different characteristic parameters of interest.
    name : str, optional 
        The name of the reactor fuel cycle component instance.

    Notes
    -----
    Some data members and functions have names that end in '_F_'.  This indicates that these are a 
    function of fluence, the time integral of the flux.  The '_Fd_' suffix implies that the data is 
    evaluated at the discharge fluence.

    """

    def __cinit__(self, reactor_parameters=None, track_params=None, char * name="", *args, **kwargs):
        cdef ReactorParameters rp
        cdef std.string cpp_name = std.string(name)

        if (reactor_parameters is None) and (track_params is None):
            self._inst = new cpp_reactor1g.Reactor1G(cpp_name)

        elif (reactor_parameters is None) and isinstance(track_params, set):
            self._inst = new cpp_reactor1g.Reactor1G(conv.py_to_cpp_set_str(track_params), cpp_name)

        elif isinstance(reactor_parameters, ReactorParameters) and (track_params is None):
            rp = reactor_parameters
            self._inst = new cpp_reactor1g.Reactor1G(<cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0], cpp_name)

        elif isinstance(reactor_parameters, ReactorParameters) and isinstance(track_params, set):
            rp = reactor_parameters
            self._inst = new cpp_reactor1g.Reactor1G(<cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0], conv.py_to_cpp_set_str(track_params), cpp_name)

        else:
            if reactor_parameters is not None:
                raise TypeError("The reactor_parameters keyword must be an instance of the ReactorParameters class or None.  Got " + str(type(reactor_parameters)))

            if track_params is not None:
                raise TypeError("The track_params keyword must be a set of strings or None.  Got " + str(type(track_params)))

        # Set property defaults
        self._fuel_chemical_form = None
        self._coolant_chemical_form = None

        self._niF = None
        self._niC = None
        self._miF = None
        self._miC = None
        self._NiF = None
        self._NiC = None

    #
    # Class Attributes
    #

    # Reactor1G attributes

    property B:
        """This integer is the total number of batches in the fuel management scheme.  
        B is typically indexed by b."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).B

        def __set__(self, int value):
            (<cpp_reactor1g.Reactor1G *> self._inst).B = value


    property phi:
        """The nominal flux value (float) that the library for this reactor type was generated with.  
        Used to correctly weight batch-specific fluxes."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).phi

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).phi = value


    property fuel_chemical_form:
        """This is the chemical form of fuel as a dictionary or other mapping.  Keys are 
        often strings that represent isotopes while values represent the corresponding 
        mass weights.  The heavy metal concentration by the key "IHM".  
        This will automatically fill in the nuclides in mat_feed for the "IHM" weight.  
        For example, LWRs typically use a UOX fuel form::

            ReactorParameters.fuel_form = {"IHM": 1.0, "O16": 2.0}

        """
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._fuel_chemical_form is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor1g.Reactor1G *> self._inst).fuel_chemical_form
                self._fuel_chemical_form = proxy

            return self._fuel_chemical_form

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_reactor1g.Reactor1G *> self._inst).fuel_chemical_form = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).fuel_chemical_form = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).fuel_chemical_form = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._fuel_chemical_form = None


    property coolant_chemical_form:
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

            if self._coolant_chemical_form is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor1g.Reactor1G *> self._inst).coolant_chemical_form
                self._coolant_chemical_form = proxy

            return self._coolant_chemical_form

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_reactor1g.Reactor1G *> self._inst).coolant_chemical_form = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).coolant_chemical_form = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).coolant_chemical_form = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._coolant_chemical_form = None


    property rhoF:
        """The fuel region density.  A float in units of [g/cm^3]."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).rhoF

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).rhoF = value


    property rhoC:
        """The coolant region density.  A float in units of [g/cm^3]."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).rhoC

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).rhoC = value


    property P_NL:
        """The reactor's non-leakage probability (float).  This is often used as a calibration parameter."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).P_NL

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).P_NL = value


    property target_BU:
        """The reactor's target discharge burnup (float).  This is given 
        in units of [MWd/kgIHM].  Often the actual discharge burnup BUd does not 
        quite hit this value, but comes acceptably close."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).target_BU

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).target_BU = value


    property use_zeta:
        """Boolaean to determine whether the thermal disadvantage factor is employed or not.  
        LWRs typically set this as True while FRs have a False value."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).use_zeta

        def __set__(self, bint value):
            (<cpp_reactor1g.Reactor1G *> self._inst).use_zeta = value


    property lattice_flag:
        """Flag (str) that represents what lattice type the fuel assemblies are arranged in.  
        Currently accepted values are "Planar", "Spherical", and "Cylindrical"."""
        def __get__(self):
            cdef std.string value = (<cpp_reactor1g.Reactor1G *> self._inst).lattice_flag
            return value.c_str()

        def __set__(self, char * value):
            (<cpp_reactor1g.Reactor1G *> self._inst).lattice_flag = std.string(value)


    property rescale_hydrogen_xs:
        """Boolean to determine whether the reactor should rescale the Hydrogen-1 destruction 
        rate in the coolant as a function of fluence.  The scaling factor is calculated via the 
        following equation

            .. math:: f(F) = 1.36927 - 0.01119 \cdot BU(F)

        This is typically not done for fast reactors but is a useful correction for LWRs."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).rescale_hydrogen_xs

        def __set__(self, bint value):
            (<cpp_reactor1g.Reactor1G *> self._inst).rescale_hydrogen_xs = value





    property r:
        """The radius (float) of the fuel region [cm]."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).r

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).r = value


    property l:
        """The pitch or length (float) of the unit fuel pin cell [cm]."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).l

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).l = value


    property S_O:
        """The number of slots (float) in a fuel assembly that are open.  Thus this is the 
        number of slots that do not contain a fuel pin and are instead filled in by coolant."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).S_O

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).S_O = value


    property S_T:
        """The total number of fuel pin slots (float) in a fuel assembly.  For a 17x17 bundle 
        this is 289.0."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).S_T

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).S_T = value


    property VF:
        """The relative fuel region volume."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).VF

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).VF = value


    property VC:
        """The relative coolant region volume."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).VC

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).VC = value





    property libfile:
        """The path (str) to the reactor data library; usually something like "LWR.h5" or "FR.h5"."""
        def __get__(self):
            cdef std.string value = (<cpp_reactor1g.Reactor1G *> self._inst).libfile
            return value.c_str()

        def __set__(self, char * value):
            (<cpp_reactor1g.Reactor1G *> self._inst).libfile = std.string(value)


    property F:
        """The fluence points that the reactor library is based on.  This is a vector of floats that 
        have units [n/kb].  This is read in from libfile."""
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).F)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).F = conv.array_to_vector_1d_dbl(value)


    property BUi_F_:
        """The burnup of each initial isotope in the core as a function of fluence.  This is a dictionary whose
        keys are initial nuclides and whose values are vectors of floats.  This data has units of [MWd/kgIHM] 
        and is read in from libfile."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).BUi_F_)

        def __set__(self, dict value):
            (<cpp_reactor1g.Reactor1G *> self._inst).BUi_F_ = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property pi_F_:
        """The neutron production rate of each initial isotope in the core as a function of fluence.  
        This is a dictionary whose keys are initial nuclides and whose values are vectors of floats.  
        This data has units of [neutrons/seconds] (abbr [n/s]) is read in from libfile."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).pi_F_)

        def __set__(self, dict value):
            (<cpp_reactor1g.Reactor1G *> self._inst).pi_F_ = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property di_F_:
        """The neutron destruction rate of each initial isotope in the core as a function of fluence. 
         This is a dictionary whose keys are initial nuclides and whose values are vectors of floats.  
        This data has units of [n/s] and is read in from libfile."""
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).di_F_)

        def __set__(self, dict value):
            (<cpp_reactor1g.Reactor1G *> self._inst).di_F_ = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property Tij_F_:
        """The transmutation matrix of each initial isotope in the core into daughter nuclides as a 
        function of fluence. This is a dictionary whose keys are initial nuclides and whose values 
        also dictionaries.  These nested dictionaries have keys that are final nuclides and whose 
        values are vectors of floats.  This data has units of [kg_i/kgIHM] or kilogram of each initial 
        nuclide per kg of total initial heavy metal. This matrix is read in from libfile."""
        def __get__(self):
            return conv.map_to_dict_int_int_vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).Tij_F_)

        def __set__(self, dict value):
            (<cpp_reactor1g.Reactor1G *> self._inst).Tij_F_ = conv.dict_to_map_int_int_array_to_vector_1d_dbl(value)





    property A_IHM:
        """The atomic weight of the initial heavy metal (float)."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).A_IHM

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).A_IHM = value


    property MWF:
        """The molecular weight of the fuel (float)."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).MWF

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).MWF = value


    property MWC:
        """The molecular weight of the coolant (float)."""
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).MWC

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).MWC = value

    
    property niF:
        """Atomic number weight of the fuel as a function of initial nuclide. 
        Map with zzaaam-integer keys and float values."""
        def __get__(self):
            cdef conv._MapIntDouble proxy

            if self._niF is None:
                proxy = conv.MapIntDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor1g.Reactor1G *> self._inst).niF
                self._niF = proxy

            return self._niF

        def __set__(self, value):
            cdef cpp_pair[int, double] item
            cdef cpp_map[int, double]  m

            if isinstance(value, conv._MapIntDouble):
                (<cpp_reactor1g.Reactor1G *> self._inst).niF = deref((<conv._MapIntDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[int, double]()
                for k, v in value.items():
                    item = cpp_pair[int, double](k, v)
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).niF = m
            elif hasattr(value, '__len__'):
                m = cpp_map[int, double]()
                for i in value:
                    item = cpp_pair[int, double](i[0], i[1])
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).niF = m
            else:
                raise TypeError('{0} cannot be converted to a C++ map.'.format(type(value)))

            self._niF = None

    
    property niC:
        """Atomic number weight of the coolant as a function of initial nuclide.  
        Map with zzaaam-integer keys and float values."""
        def __get__(self):
            cdef conv._MapIntDouble proxy

            if self._niC is None:
                proxy = conv.MapIntDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor1g.Reactor1G *> self._inst).niC
                self._niC = proxy

            return self._niC

        def __set__(self, value):
            cdef cpp_pair[int, double] item
            cdef cpp_map[int, double]  m

            if isinstance(value, conv._MapIntDouble):
                (<cpp_reactor1g.Reactor1G *> self._inst).niC = deref((<conv._MapIntDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[int, double]()
                for k, v in value.items():
                    item = cpp_pair[int, double](k, v)
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).niC = m
            elif hasattr(value, '__len__'):
                m = cpp_map[int, double]()
                for i in value:
                    item = cpp_pair[int, double](i[0], i[1])
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).niC = m
            else:
                raise TypeError('{0} cannot be converted to a C++ map.'.format(type(value)))

            self._niC = None

    
    property miF:
        """Mass weight of the fuel as a function of initial nuclide.  
        Map with zzaaam-integer keys and float values."""
        def __get__(self):
            cdef conv._MapIntDouble proxy

            if self._miF is None:
                proxy = conv.MapIntDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor1g.Reactor1G *> self._inst).miF
                self._miF = proxy

            return self._miF

        def __set__(self, value):
            cdef cpp_pair[int, double] item
            cdef cpp_map[int, double]  m

            if isinstance(value, conv._MapIntDouble):
                (<cpp_reactor1g.Reactor1G *> self._inst).miF = deref((<conv._MapIntDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[int, double]()
                for k, v in value.items():
                    item = cpp_pair[int, double](k, v)
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).miF = m
            elif hasattr(value, '__len__'):
                m = cpp_map[int, double]()
                for i in value:
                    item = cpp_pair[int, double](i[0], i[1])
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).miF = m
            else:
                raise TypeError('{0} cannot be converted to a C++ map.'.format(type(value)))

            self._miF = None

    
    property miC:
        """Mass weight of the coolant as a function of initial nuclide.  
        Map with zzaaam-integer keys and float values."""
        def __get__(self):
            cdef conv._MapIntDouble proxy

            if self._miC is None:
                proxy = conv.MapIntDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor1g.Reactor1G *> self._inst).miC
                self._miC = proxy

            return self._miC

        def __set__(self, value):
            cdef cpp_pair[int, double] item
            cdef cpp_map[int, double]  m

            if isinstance(value, conv._MapIntDouble):
                (<cpp_reactor1g.Reactor1G *> self._inst).miC = deref((<conv._MapIntDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[int, double]()
                for k, v in value.items():
                    item = cpp_pair[int, double](k, v)
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).miC = m
            elif hasattr(value, '__len__'):
                m = cpp_map[int, double]()
                for i in value:
                    item = cpp_pair[int, double](i[0], i[1])
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).miC = m
            else:
                raise TypeError('{0} cannot be converted to a C++ map.'.format(type(value)))

            self._miC = None

    
    property NiF:
        """Number density of the fuel as a function of initial nuclide.  
        Map with zzaaam-integer keys and float values."""
        def __get__(self):
            cdef conv._MapIntDouble proxy

            if self._NiF is None:
                proxy = conv.MapIntDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor1g.Reactor1G *> self._inst).NiF
                self._NiF = proxy

            return self._NiF

        def __set__(self, value):
            cdef cpp_pair[int, double] item
            cdef cpp_map[int, double]  m

            if isinstance(value, conv._MapIntDouble):
                (<cpp_reactor1g.Reactor1G *> self._inst).NiF = deref((<conv._MapIntDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[int, double]()
                for k, v in value.items():
                    item = cpp_pair[int, double](k, v)
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).NiF = m
            elif hasattr(value, '__len__'):
                m = cpp_map[int, double]()
                for i in value:
                    item = cpp_pair[int, double](i[0], i[1])
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).NiF = m
            else:
                raise TypeError('{0} cannot be converted to a C++ map.'.format(type(value)))

            self._NiF = None

    
    property NiC:
        """Number density of the coolant as a function of initial nuclide.  
        Map with zzaaam-integer keys and float values."""
        def __get__(self):
            cdef conv._MapIntDouble proxy

            if self._NiC is None:
                proxy = conv.MapIntDouble(False, False)
                proxy.map_ptr = &(<cpp_reactor1g.Reactor1G *> self._inst).NiC
                self._NiC = proxy

            return self._NiC

        def __set__(self, value):
            cdef cpp_pair[int, double] item
            cdef cpp_map[int, double]  m

            if isinstance(value, conv._MapIntDouble):
                (<cpp_reactor1g.Reactor1G *> self._inst).NiC = deref((<conv._MapIntDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[int, double]()
                for k, v in value.items():
                    item = cpp_pair[int, double](k, v)
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).NiC = m
            elif hasattr(value, '__len__'):
                m = cpp_map[int, double]()
                for i in value:
                    item = cpp_pair[int, double](i[0], i[1])
                    m.insert(item)
                (<cpp_reactor1g.Reactor1G *> self._inst).NiC = m
            else:
                raise TypeError('{0} cannot be converted to a C++ map.'.format(type(value)))

            self._NiC = None





    property dF_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).dF_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).dF_F_ = conv.array_to_vector_1d_dbl(value)


    property dC_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).dC_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).dC_F_ = conv.array_to_vector_1d_dbl(value)


    property BU_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).BU_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).BU_F_ = conv.array_to_vector_1d_dbl(value)


    property P_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).P_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).P_F_ = conv.array_to_vector_1d_dbl(value)


    property D_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).D_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).D_F_ = conv.array_to_vector_1d_dbl(value)


    property k_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).k_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).k_F_ = conv.array_to_vector_1d_dbl(value)


    property Mj_F_:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).Mj_F_)

        def __set__(self, dict value):
            (<cpp_reactor1g.Reactor1G *> self._inst).Mj_F_ = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property zeta_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).zeta_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).zeta_F_ = conv.array_to_vector_1d_dbl(value)





    property fd:
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).fd

        def __set__(self, int value):
            (<cpp_reactor1g.Reactor1G *> self._inst).fd = value


    property Fd:
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).Fd

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).Fd = value


    property BUd:
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).BUd

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).BUd = value


    property k:
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).k

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).k = value





    property mat_feed_u:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).mat_feed_u
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactor1g.Reactor1G *> self._inst).mat_feed_u = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_feed_tru:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).mat_feed_tru
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactor1g.Reactor1G *> self._inst).mat_feed_tru = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_feed_lan:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).mat_feed_lan
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactor1g.Reactor1G *> self._inst).mat_feed_lan = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_feed_act:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).mat_feed_act
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactor1g.Reactor1G *> self._inst).mat_feed_act = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_u:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).mat_prod_u
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactor1g.Reactor1G *> self._inst).mat_prod_u = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_tru:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).mat_prod_tru
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactor1g.Reactor1G *> self._inst).mat_prod_tru = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_lan:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).mat_prod_lan
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactor1g.Reactor1G *> self._inst).mat_prod_lan = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_act:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).mat_prod_act
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactor1g.Reactor1G *> self._inst).mat_prod_act = <pyne.cpp_material.Material> mat.mat_pointer[0]




    property deltaR:
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).deltaR

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).deltaR = value


    property tru_cr:
        def __get__(self):
            return (<cpp_reactor1g.Reactor1G *> self._inst).tru_cr

        def __set__(self, double value):
            (<cpp_reactor1g.Reactor1G *> self._inst).tru_cr = value





    property SigmaFa_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).SigmaFa_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).SigmaFa_F_ = conv.array_to_vector_1d_dbl(value)


    property SigmaFtr_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).SigmaFtr_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).SigmaFtr_F_ = conv.array_to_vector_1d_dbl(value)


    property kappaF_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).kappaF_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).kappaF_F_ = conv.array_to_vector_1d_dbl(value)





    property SigmaCa_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).SigmaCa_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).SigmaCa_F_ = conv.array_to_vector_1d_dbl(value)


    property SigmaCtr_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).SigmaCtr_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).SigmaCtr_F_ = conv.array_to_vector_1d_dbl(value)


    property kappaC_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).kappaC_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).kappaC_F_ = conv.array_to_vector_1d_dbl(value)





    property lattice_E_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).lattice_E_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).lattice_E_F_ = conv.array_to_vector_1d_dbl(value)


    property lattice_F_F_:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactor1g.Reactor1G *> self._inst).lattice_F_F_)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactor1g.Reactor1G *> self._inst).lattice_F_F_ = conv.array_to_vector_1d_dbl(value)


    #
    # Class Methods
    # 

    def initialize(self, ReactorParameters reactor_parameters):
        """The initialize() method for reactors copies all of the reactor specific parameters to this instance.
        Additionally, it calculates and sets the volumes VF and VC.

        Args:
            * reactor_parameters (ReactorParameters): A special data structure that contains information
              on how to setup and run the reactor.
        """
        cdef ReactorParameters rp = reactor_parameters
        (<cpp_reactor1g.Reactor1G *> self._inst).initialize(<cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0])


    def loadlib(self, char * libfile="reactor.h5"):
        """This method finds the HDF5 library for this reactor and extracts the necessary information from it.
        This method is typically called by the constructor of the child reactor type object.  It must be 
        called before attempting to do any real computation.

        Args: 
            * libfile (str): Path to the reactor library.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).loadlib(std.string(libfile))


    def fold_mass_weights(self):
        """This method performs the all-important task of doing the isotopically-weighted linear combination of raw data. 
        In a very real sense this is what makes this reactor *this specific reactor*.  The weights are taken 
        as the values of mat_feed.  The raw data must have previously been read in from loadlib().  

        .. warning::

            Anytime any reactor parameter whatsoever (mat_feed, P_NL, *etc*) is altered in any way, 
            the fold_mass_weights() function must be called to reset all of the resultant data.
            If you are unsure, please call this function anyway to be safe.  There is little 
            harm in calling it twice by accident.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).fold_mass_weights()





    def calc_Mj_F_(self):
        """This function calculates and sets the Mj_F_ attribute from mat_feed and the 
        raw reactor data Tij_F_.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).calc_Mj_F_()


    def calc_Mj_Fd_(self):
        """This function evaluates Mj_F_ calculated from calc_Mj_F_() at the discharge fluence Fd.
        The resultant isotopic dictionary is then converted into the mat_prod mass stream
        for this pass through the reactor.  Thus if ever you need to calculate mat_prod
        without going through calc(), use this function.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).calc_Mj_Fd_()





    def calc_mat_prod(self):
        """This is a convenience function that wraps the transmutation matrix methods.  
        It is equivalent to::

            #Wrapper to calculate discharge isotopics.
            calc_Mj_F_()
            calc_Mj_Fd_()

        """
        (<cpp_reactor1g.Reactor1G *> self._inst).calc_mat_prod()

    def calcSubStreams(self):
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
        (<cpp_reactor1g.Reactor1G *> self._inst).calcSubStreams()


    def calc_tru_cr(self):
        """This calculates and sets the transuranic conversion ratio tru_cr through the equation:

        .. math:: \mbox{tru_cr} = \frac{\mbox{mat_feed_tru.mass} - \mbox{mat_prod_tru.mass}}{\frac{\mbox{BUd}}{935.0}}

        Returns:
            * tru_cr (float): The value of the transuranic conversion ratio just calculated.
        """
        return (<cpp_reactor1g.Reactor1G *> self._inst).calc_tru_cr()





    def calc_deltaR(self, input=None):
        """Calculates and sets the deltaR value of the reactor.  
        This is equal to the production rate minus the destruction rate at the target burnup::

            deltaR = batch_average(target_BU, "P") - batch_average(target_BU, "D")

        Args:
            * input (dict or Material): If input is present, it set as the component's 
              mat_feed.  If input is a isotopic dictionary (zzaaam keys, float values), this
              dictionary is first converted into a Material before being set as mat_feed.

        Returns:
            * deltaR (float): deltaR.
        """
        cdef pyne.material._Material in_mat 
        cdef double deltaR 

        if input is None:
            deltaR = (<cpp_reactor1g.Reactor1G *> self._inst).calc_deltaR()
        elif isinstance(input, dict):
            deltaR = (<cpp_reactor1g.Reactor1G *> self._inst).calc_deltaR(conv.dict_to_map_int_dbl(input))
        elif isinstance(input, pyne.material._Material):
            in_mat = input
            deltaR = (<cpp_reactor1g.Reactor1G *> self._inst).calc_deltaR(<pyne.cpp_material.Material> in_mat.mat_pointer[0])

        return deltaR





    def fluence_at_BU(self, double burnup):
        """This function takes a burnup value  and returns a special fluence point object.  
        The fluence point is an amalgamation of data where the at which the burnup occurs.
        This object instance FP contains three pieces of information::
    
            FP.f    #Index immediately lower than where BU achieved (int)
            FP.F    #Fluence value itself (float)
            FP.m    #Slope dBU/dF between points f and f+1 (double)

        Args:
            * burnup (float): Burnup [MWd/kgIHM] at which to calculate the corresponding fluence.

        Returns:
            * fp (FluencePoint): A class containing fluence information.
        """
        cdef FluencePoint fp = FluencePoint()
        fp.fp_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).fluence_at_BU(burnup)
        return fp


    def batch_average(self, double BUd, char * PDk_flag="K"):
        """Finds the batch-averaged P(F), D(F), or k(F) when at discharge burnup BUd.
        This function is typically iterated over until a BUd is found such that k(F) = 1.0 + err.

        Args:
            * BUd (float): The discharge burnup [MWd/kgIHM] to obtain a batch-averaged value for.

        Keyword Args:
            * PDk_flag (string): Flag that determines whether the neutron production rate "P" [n/s], 
              the neutron destruction rate "D" [n/s], or the multiplication factor "K" is reported in the output.

        Returns:
            * PDk (float): the batch averaged neutron production rate,
        """
        cdef double PDk = (<cpp_reactor1g.Reactor1G *> self._inst).batch_average(BUd, std.string(PDk_flag))
        return PDk


    def batch_average_k(self, double BUd):
        """Convenience function that calls batch_average(BUd, "K").

        Args:
            * BUd (float): The discharge burnup [MWd/kgIHM] to obtain a batch-averaged value for.

        Returns:
            * k (float): the batch averaged multiplication factor.
        """
        cdef double PDk = (<cpp_reactor1g.Reactor1G *> self._inst).batch_average_k(BUd)
        return PDk


    def BUd_bisection_method(self):
        """Calculates the maximum discharge burnup via the Bisection Method for a given mat_feed
        in this reactor.  This iterates over values of BUd to find a batch averaged multiplication factor 
        that is closest to 1.0.

        Other root finding methods for determining maximum discharge burnup are certainly possible.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).BUd_bisection_method()


    def run_P_NL(self, double pnl):
        """Performs a reactor run for a specific non-leakage probability value.
        This requires that mat_feed be (meaningfully) set and is for use with calibrate_P_NL_to_BUd().

        This function amounts to the following code::

            self.P_NL = pnl
            self.fold_mass_weights()
            self.BUd_bisection_method()

        Args:
            * pnl (float): The new non-leakage probability for the reactor.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).run_P_NL(pnl)
    

    def calibrate_P_NL_to_BUd(self):
        """Often times the non-leakage probability of a reactor is not known, though the input isotopics 
        and the target discharge burnup are.  This function handles that situation by
        calibrating the non-leakage probability of this reactor P_NL to hit its target burnup target_BU.
        Such a calibration proceeds by bisection method as well.  This function is extremely useful for 
        benchmarking calculations.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).calibrate_P_NL_to_BUd()




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

        Args:
            * input (dict or Material): If input is present, it set as the component's 
              mat_feed.  If input is a isotopic dictionary (zzaaam keys, float values), this
              dictionary is first converted into a Material before being set as mat_feed.

        Returns:
            * output (Material): mat_prod.
        """
        cdef pyne.material._Material in_mat 
        cdef pyne.material._Material output = pyne.material.Material()

        if input is None:
            output.mat_pointer[0] = (<cpp_fccomp.FCComp *> self._inst).calc()
        elif isinstance(input, dict):
            output.mat_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).calc(conv.dict_to_map_int_dbl(input))
        elif isinstance(input, pyne.material._Material):
            in_mat = input
            output.mat_pointer[0] = (<cpp_reactor1g.Reactor1G *> self._inst).calc(<pyne.cpp_material.Material> in_mat.mat_pointer[0])

        return output






    def lattice_E_planar(self, double a, double b):
        """Calculates the lattice function E(F) for planar geometry.  Sets value as lattice_E_F_

        Args:
            * a (float): Fuel region radius equivalent [cm].
            * b (float): Unit fuel cell pitch length equivalent [cm].
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).lattice_E_planar(a, b)


    def lattice_F_planar(self, double a, double b):
        """Calculates the lattice function F(F) for planar geometry.  Sets value as lattice_F_F_

        Args:
            * a (float): Fuel region radius equivalent [cm].
            * b (float): Unit fuel cell pitch length equivalent [cm].
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).lattice_F_planar(a, b)


    def lattice_E_spherical(self, double a, double b):
        """Calculates the lattice function E(F) for spherical geometry.  Sets value as lattice_E_F_

        Args:
            * a (float): Fuel region radius equivalent [cm].
            * b (float): Unit fuel cell pitch length equivalent [cm].
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).lattice_E_spherical(a, b)


    def lattice_F_spherical(self, double a, double b):
        """Calculates the lattice function F(F) for spherical geometry.  Sets value as lattice_F_F_

        Args:
            * a (float): Fuel region radius equivalent [cm].
            * b (float): Unit fuel cell pitch length equivalent [cm].
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).lattice_F_spherical(a, b)


    def lattice_E_cylindrical(self, double a, double b):
        """Calculates the lattice function E(F) for cylindrical geometry.  Sets value as lattice_E_F_

        Args:
            * a (float): Fuel region radius equivalent [cm].
            * b (float): Unit fuel cell pitch length equivalent [cm].
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).lattice_E_cylindrical(a, b)


    def lattice_F_cylindrical(self, double a, double b):
        """Calculates the lattice function F(F) for cylindrical geometry.  Sets value as lattice_F_F_

        Args:
            * a (float): Fuel region radius equivalent [cm].
            * b (float): Unit fuel cell pitch length equivalent [cm].
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).lattice_F_cylindrical(a, b)





    def calc_zeta(self):
        """This calculates the thermal disadvantage factor for the geometry specified by Lattice.  The results
        are set to zeta_F_.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).calc_zeta()


    def calc_zeta_planar(self):
        """This calculates the thermal disadvantage factor for a planar geometry.  The results
        are set to zeta_F_.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).calc_zeta_planar()


    def calc_zeta_spherical(self):
        """This calculates the thermal disadvantage factor for a spherical geometry.  The results
        are set to zeta_F_.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).calc_zeta_spherical()


    def calc_zeta_cylindrical(self):
        """This calculates the thermal disadvantage factor for a clyindrical geometry.  The results
        are set to zeta_F_.
        """
        (<cpp_reactor1g.Reactor1G *> self._inst).calc_zeta_cylindrical()



