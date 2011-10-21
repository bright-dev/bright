"""Python wrapper for fuel fabrication."""
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
cimport cpp_fuel_fabrication
cimport cpp_reactor1g

from reactor1g cimport Reactor1G
from reactor1g import Reactor1G

cimport fccomp
import fccomp


#############################
### FuelFabrication Class ###
#############################


cdef class FuelFabrication(fccomp.FCComp):
    """Fuel Fabrication Fuel Cycle Component Class.  Daughter of FCComp class.

    Parameters
    ----------
    materials : dict or map or None, optional 
        A dictionary whose keys are string labels (eg, "U-235", "TRU", "My Fuel") 
        and whose values are mass streams.  For example::

            materials = {
                "U235": Material({922350: 1.0}, 1.0, "U-235"),
                "U236": Material({922360: 1.0}, 1.0, "U-236"),
                "U238": Material({922380: 1.0}, 1.0, "U-238"),
                }

        would be valid for a light water reactor.
    mass_weights_in : dict or map or None, optional 
        A dictionary whose keys are the same as for materials and whose values are the 
        associated weight (float) for that stream.  If a material should be allowed to 
        vary (ie optimized over), specify its weight as a negative number. For instance::

            mass_weights_in = {
                "U235": -1.0,
                "U236": 0.005,
                "U238": -1.0,        
                }

        would be valid for a light water reactor with half a percent of U-236 always present.
    reactor : Reactor1G or None, optional 
        An instance of a Reactor1G class to fabricate fuel for.
    track_params : list of str or None, optional 
        Additional parameters to track, if any.        
    name : str, optional 
        The name of the fuel fabrication fuel cycle component instance.

    """

    def __cinit__(self, *args, **kwargs):
        # Set property defaults
        self._materials = None
        self._mass_weights_in = None
        self._mass_weights_out = None
        self._deltaRs = None


    def __init__(self, materials=None, mass_weights_in=None, reactor=None, track_params=None, char * name="", *args, **kwargs):
        cdef std.string cpp_name = std.string(name)
        cdef Reactor1G r1g 

        if (materials is None) and (mass_weights_in is None) and (reactor is None) and (track_params is None):
            self._inst = new cpp_fuel_fabrication.FuelFabrication(std.string(name))

        elif (materials is None) and (mass_weights_in is None) and (reactor is None) and isinstance(track_params, set):
            self._inst = new cpp_fuel_fabrication.FuelFabrication(conv.py_to_cpp_set_str(track_params), cpp_name)

        elif isinstance(materials, dict) and isinstance(mass_weights_in, dict) and isinstance(reactor, Reactor1G) and (track_params is None):
            r1g = reactor
            self._inst = new cpp_fuel_fabrication.FuelFabrication(
                                pyne.material.dict_to_map_str_matp(materials), 
                                conv.dict_to_map_str_dbl(mass_weights_in), 
                                deref(<cpp_reactor1g.Reactor1G *> r1g._inst), 
                                std.string(name))

        elif isinstance(materials, dict) and isinstance(mass_weights_in, dict) and isinstance(reactor, Reactor1G) and isinstance(track_params, set):
            r1g = reactor
            self._inst = new cpp_fuel_fabrication.FuelFabrication(
                                pyne.material.dict_to_map_str_matp(materials), 
                                conv.dict_to_map_str_dbl(mass_weights_in), 
                                deref(<cpp_reactor1g.Reactor1G *> r1g._inst), 
                                conv.py_to_cpp_set_str(track_params),
                                std.string(name))

        else:
            if materials is not None:
                raise TypeError("The materials keyword must be a dictionary of (string, Material) items or None.  Got " + str(type(materials)))

            if mass_weights_in is not None:
                raise TypeError("The mass_weights_in keyword must be a dictionary of (string, float) items or None.  Got " + str(type(mass_weights_in)))

            if reactor is not None:
                raise TypeError("The reactor keyword must be a Reactor1G instance or None.  Got " + str(type(reactor)))

            if track_params is not None:
                raise TypeError("The track_params keyword must be a set of strings or None.  Got " + str(type(track_params)))


    #
    # Class Attributes
    #

    # FuelFabrication attributes

    property materials:
        """A mapping of materials which are mixed to create a valid fuel form for reactor."""
        def __get__(self):
            cdef pyne.material._MapStrMaterial proxy

            if self._materials is None:
                proxy = pyne.material.MapStrMaterial(False, False)
                proxy.map_ptr = &(<cpp_fuel_fabrication.FuelFabrication *> self._inst).materials
                self._materials = proxy

            return self._materials

        def __set__(self, value):
            cdef std.string s
            cdef pyne.material._MapStrMaterial proxy

            if isinstance(value, pyne.material._MapStrMaterial):
                (<cpp_fuel_fabrication.FuelFabrication *> self._inst).materials = deref((<pyne.material._MapStrMaterial> value).map_ptr)
                self._materials = None
            elif hasattr(value, 'items'):
                proxy = pyne.material.MapStrMaterial()
                for k, v in value.items():
                    proxy[k] = v
                self._materials = proxy
            elif hasattr(value, '__len__'):
                proxy = pyne.material.MapStrMaterial()
                for i in value:
                    proxy[i[0]] = i[1]
                self._materials = proxy
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

                self._materials = proxy


    property mass_weights_in:
        """A mapping representing the initial specification of the mass weights.  
        Exactly two of these should have negative values."""
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._mass_weights_in is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_fuel_fabrication.FuelFabrication *> self._inst).mass_weights_in
                self._mass_weights_in = proxy

            return self._mass_weights_in

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_fuel_fabrication.FuelFabrication *> self._inst).mass_weights_in = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_fuel_fabrication.FuelFabrication *> self._inst).mass_weights_in = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_fuel_fabrication.FuelFabrication *> self._inst).mass_weights_in = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._mass_weights_in = None


    property mass_weights_out:
        """A mapping representing the mass weights that are calculated to generate a valid fuel 
        from the materials  provided."""
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._mass_weights_out is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_fuel_fabrication.FuelFabrication *> self._inst).mass_weights_out
                self._mass_weights_out = proxy

            return self._mass_weights_out

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_fuel_fabrication.FuelFabrication *> self._inst).mass_weights_out = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_fuel_fabrication.FuelFabrication *> self._inst).mass_weights_out = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_fuel_fabrication.FuelFabrication *> self._inst).mass_weights_out = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._mass_weights_out = None


    property deltaRs:
        """A mapping representing the deltaR values of each of the materials."""
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._deltaRs is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_fuel_fabrication.FuelFabrication *> self._inst).deltaRs
                self._deltaRs = proxy

            return self._deltaRs

        def __set__(self, dict value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_fuel_fabrication.FuelFabrication *> self._inst).deltaRs = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_fuel_fabrication.FuelFabrication *> self._inst).deltaRs = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_fuel_fabrication.FuelFabrication *> self._inst).deltaRs = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._deltaRs = None


    property reactor:
        """An instance of a Reactor1G class that the mat_prod is valid as a fuel for."""
        def __get__(self):
            cdef Reactor1G value = Reactor1G()
            cdef cpp_reactor1g.Reactor1G cpp_value = (<cpp_fuel_fabrication.FuelFabrication *> self._inst).reactor
            (<Reactor1G> value)._free_inst = False
            (<Reactor1G> value)._inst = &cpp_value
            return value

        def __set__(self, Reactor1G value):
            (<cpp_fuel_fabrication.FuelFabrication *> self._inst).reactor = deref(<cpp_reactor1g.Reactor1G *> value._inst)


    #
    # Class Methods
    # 

    def initialize(self, dict materials, dict mass_weights_in, Reactor1G reactor):
        """The initialize() method takes the appropriate materials, input mass weights,
        and a reactor and resets the current FuelFabrication instance.

        Parameters
        ----------
        materials : dict
            A dictionary whose keys are string labels and whose values are mass streams.  
        mass_weights_in : dict 
            A dictionary whose keys are the same as for materials and whose values are the 
            associated weight (float) for that stream.
        reactor : Reactor1G 
            An instance of a Reactor1G class to fabricate fuel for.

        """
        cdef Reactor1G r1g = reactor
        (<cpp_fuel_fabrication.FuelFabrication *> self._inst).initialize(pyne.material.dict_to_map_str_matp(materials), 
                                   conv.dict_to_map_str_dbl(mass_weights_in), 
                                   deref(<cpp_reactor1g.Reactor1G *> r1g._inst))


    def calc_params(self):
        """Here the parameters for FuelFabrication are set.  For example::

            self.params_prior_calc["Weight_U235"]  = self.mass_weights_in["U235"]
            self.Paramsout["Weight_U235"] = self.mass_weights_out["U235"]

            self.params_prior_calc["deltaR_U235"]  = self.deltaRs["U235"]
            self.Paramsout["deltaR_U235"] = self.deltaRs["U235"]

            self.params_prior_calc["Weight_U238"]  = self.mass_weights_in["U238"]
            self.Paramsout["Weight_U238"] = self.mass_weights_out["U238"]

            self.params_prior_calc["deltaR_U238"]  = self.deltaRs["U238"]
            self.Paramsout["deltaR_U238"] = self.deltaRs["U238"]
        """
        (<cpp_fccomp.FCComp *> self._inst).calc_params()



    def calc_deltaRs(self):
        """Computes deltaRs for each mass stream."""
        (<cpp_fuel_fabrication.FuelFabrication *> self._inst).calc_deltaRs()


    def calc_core_input(self):
        """Computes the core input material that becomes mat_prod based on materials and 
        mass_weights_out.

        Returns
        -------
        core_input : Material 
            mat_prod.
        """
        cdef pyne.material._Material pymat = pyne.pyne.material.Material()
        pymat.mat_pointer[0] = (<cpp_fuel_fabrication.FuelFabrication *> self._inst).calc_core_input()
        return pymat


    def calc_mass_ratios(self):
        """Calculates mass_weights_out by varying the values of the two parameter which had 
        negative values in mass_weights_in.  Therefore, this is the portion of the code that 
        performs the optimization calculation.
        """
        (<cpp_fuel_fabrication.FuelFabrication *> self._inst).calc_mass_ratios()





    def calc(self, materials=None, mass_weights_in=None, reactor=None):
        """This method performs an optimization calculation on all input materials to determine
        the mass ratios that generate the correct fuel form for the reactor.  It then compiles 
        the fuel and returns the resultant material. 

        Parameters
        ----------
        materials : dict or None, optional
            A dictionary whose keys are string labels and whose values are materials.  
        mass_weights_in : dict or None, optional
            A dictionary whose keys are the same as for materials and whose values are 
            the associated weight (float) for that material.
        reactor : Reactor1G or None, optional 
            An instance of a Reactor1G class to fabricate fuel for.

        Returns
        -------
        core_input : Material 
            mat_prod

        """
        cdef Reactor1G r1g 
        cdef pyne.material._Material core_input = pyne.material.Material()

        if (materials is None) and (mass_weights_in is None) and (reactor is None):
            core_input.mat_pointer[0] = (<cpp_fccomp.FCComp *> self._inst).calc()

        elif isinstance(materials, dict) and isinstance(mass_weights_in, dict) and isinstance(reactor, Reactor1G):
            r1g = reactor
            core_input.mat_pointer[0] = (<cpp_fuel_fabrication.FuelFabrication *> self._inst).calc(
                                                            pyne.material.dict_to_map_str_matp(materials), 
                                                            conv.dict_to_map_str_dbl(mass_weights_in), 
                                                            deref(<cpp_reactor1g.Reactor1G *> r1g._inst))

        else:
            if materials is not None:
                raise TypeError("The materials keyword must be a dictionary of (string, Material) items or None.  Got " + str(type(materials)))

            if mass_weights_in is not None:
                raise TypeError("The mass_weights_in keyword must be a dictionary of (string, float) items or None.  Got " + str(type(mass_weights_in)))

            if reactor is not None:
                raise TypeError("The reactor keyword must be a Reactor1G instance or None.  Got " + str(type(reactor)))

        return core_input



