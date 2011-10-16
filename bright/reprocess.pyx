"""Python wrapper for Reprocess."""
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

from pyne import nucname

cimport pyne.cpp_material
cimport pyne.material
import pyne.material

cimport cpp_fccomp
cimport cpp_reprocess

cimport fccomp
import fccomp

from reprocess cimport convert_sepeff_key, Reprocess

#######################
### Reprocess Class ###
#######################

cdef int convert_sepeff_key(object key):
    cdef int k
    if isinstance(key, int):
        k = key
    elif isinstance(key, basestring):
        if key in nucname.name_zz:
            k = nucname.name_zz[key]
        else:
            k = nucname.zzaaam(key)
    else:
        raise TypeError("Separation keys must be strings or integers.")
    return k
    

cdef class Reprocess(fccomp.FCComp):
    """Reprocess Fuel Cycle Component Class.  Daughter of bright.FCComp class.

    Args:
        * sepeff (dict): A dictionary containing the separation efficiencies (float) to initialize
          the instance with.  The keys of this dictionary must be strings.  However, the strings may 
          represent either elements or isotopes or both::

                #ssed = string dictionary of separation efficiencies.  
                #Of form {zz: 0.99}, eg 
                ssed = {"92": 0.999, "94": 0.99} 
                #of form {LL: 0.99}, eg 
                ssed = {"U": 0.999, "PU": 0.99} 
                #or of form {mixed: 0.99}, eg 
                ssed = {"U235": 0.9, "922350": 0.999, "94239": 0.99}

        * name (str): The name of the reprocessing fuel cycle component instance.

    Note that this automatically calls the public initialize C function.

    .. note::
       The C++ version of the code also allows you to initialize from an int-keyed dictionary (map).
       However, due to a from_python C++ signature ambiguity, you cannot do use this directly in Python.
       Separation efficiencies must therefore be automatically initialized through string dictionaries.
       If you need to initialize via an int dictionary in python, you can always init with an empty
       string dictionary and then manually initialize with an int one.  For example::

            R = Reprocess({}, name)
            R.initialize( {92: 0.99, 942390: 0.9} )

    """

    def __cinit__(self, sepeff=None, char * name=""):
        if sepeff is None:
            sepeff = {}

        cdef dict sepdict = {}
        for key, val in sepeff.items():
            sepdict[convert_sepeff_key(key)] = val
        self._inst = new cpp_reprocess.Reprocess(conv.dict_to_map_int_dbl(sepdict), std.string(name))

        self._sepeff = None

    #
    # Class Attributes
    #

    # Reprocess attributes

    property sepeff:
        def __get__(self):
            cdef conv._MapIntDouble proxy

            if self._params_after_calc is None:
                proxy = conv.MapIntDouble(False, False)
                proxy.map_ptr = &(<cpp_reprocess.Reprocess *> self._inst).sepeff
                self._sepeff = proxy

            return self._sepeff
            #return conv.map_to_dict_int_dbl((<cpp_reprocess.Reprocess *> self._inst).sepeff)

        def __set__(self, value):
            cdef cpp_pair[int, double] item
            cdef cpp_map[int, double]  m

            if isinstance(value, conv._MapIntDouble):
                (<cpp_reprocess.Reprocess *> self._inst).sepeff = deref((<conv._MapIntDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[int, double]()
                for k, v in value.items():
                    item = cpp_pair[int, double](convert_sepeff_key(k), v)
                    m.insert(item)
                (<cpp_reprocess.Reprocess *> self._inst).sepeff = m
            elif hasattr(value, '__len__'):
                m = cpp_map[int, double]()
                for i in value:
                    item = cpp_pair[int, double](convert_sepeff_key(i[0]), i[1])
                    m.insert(item)
                (<cpp_reprocess.Reprocess *> self._inst).sepeff = m
            else:
                raise TypeError('{0} cannot be converted to a C++ map.'.format(type(value)))

            self._sepeff = None

            #value = self._cpp_sepeff(value)
            #(<cpp_reprocess.Reprocess *> self._inst).sepeff = conv.dict_to_map_int_dbl(value)



    #
    # Class Methods
    # 

    def initialize(self, sepeff):
        """The initialize() function calculates the sepeff from an integer-keyed dictionary
        of separation efficiencies.  The difference is that sepdict may contain either elemental or
        isotopic keys and need not contain every isotope tracked.  On the other hand, sepeff
        must have only zzaaam keys that match exactly the isotopes in bright.track_nucs.

        Args:
            * sepeff (dict): Integer valued dictionary of SE to be converted to sepeff.
        """
        cdef dict sepdict = {}
        for key, val in sepeff.items():
            sepdict[convert_sepeff_key(key)] = val
        (<cpp_reprocess.Reprocess *> self._inst).initialize(conv.dict_to_map_int_dbl(sepdict))


    def calc_params(self):
        """Here the parameters for Reprocess are set.  For reprocessing, this amounts to just
        a "Mass" parameter::

            self.params_prior_calc["Mass"]  = self.mat_feed.mass
            self.params_after_calc["Mass"] = self.mat_prod.mass

        """
        (<cpp_fccomp.FCComp *> self._inst).calc_params()


    def calc(self, input=None):
        """This method performs the relatively simply task of multiplying the current input stream by 
        the SE to form a new output stream::

            incomp  = self.mat_feed.mult_by_mass()
            outcomp = {}
            for iso in incomp.keys():
                outcomp[iso] = incomp[iso] * sepeff[iso]
            self.mat_prod = Material(outcomp)
            return self.mat_prod

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
            output.mat_pointer[0] = (<cpp_reprocess.Reprocess *> self._inst).calc(conv.dict_to_map_int_dbl(input))
        elif isinstance(input, pyne.material._Material):
            in_mat = input
            output.mat_pointer[0] = (<cpp_reprocess.Reprocess *> self._inst).calc(<pyne.cpp_material.Material> in_mat.mat_pointer[0])

        return output


