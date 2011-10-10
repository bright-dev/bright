"""Python wrapper for  ."""
# Cython imports
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

cimport pyne.cpp_material
cimport pyne.material
import pyne.material

cimport cpp_fccomp
cimport cpp_storage

cimport fccomp
import fccomp



#####################
### Storage Class ###
#####################


cdef class Storage(fccomp.FCComp):
    """Storage Fuel Cycle Component Class.  Daughter of bright.FCComp class.

    Args:
        * name (str): The name of the storage fuel cycle component instance.
    """

    #cdef cpp_fccomp.Storage * s_pointer

    def __cinit__(self, char * name=""):
        self._inst = new cpp_storage.Storage(std.string(name))


    #
    # Class Attributes
    #

    # Stroage attributes

    property decay_time:
        def __get__(self):
            return (<cpp_storage.Storage *> self._inst).decay_time

        def __set__(self, value):
            (<cpp_storage.Storage *> self._inst).decay_time = <double> value



    #
    # Class Methods
    # 

    def calc_params(self):
        """Here the parameters for Storage are set.  For storage, this amounts to just
        a "Mass" parameter::

            self.params_prior_calc["Mass"]  = self.mat_feed.mass
            self.params_after_calc["Mass"] = self.mat_prod.mass
        """
        (<cpp_fccomp.FCComp *> self._inst).calc_params()


    def calc(self, input=None, decay_time=None):
        """As usual, calc sets up the Storage component's input stream and calculates the corresponding 
        output Material.  Here, this amounts to calling bateman() for every nuclide in 
        mat_feed, for each chain that ends with a nuclide in track_nucs.

        This method is public and accessible from Python.

        Args:
            * input (dict or Material): If input is present, it set as the component's 
              mat_feed.  If input is a isotopic dictionary (zzaaam keys, float values), this
              dictionary is first converted into a Material before being set as mat_feed.
            * decay_time (float): decay_time is set to the time value here prior to any other calculations.  This
              time has units of seconds.

        Returns:
            * output (Material): mat_prod.
        """
        cdef pyne.material._Material in_mat 
        cdef pyne.material._Material output = pyne.material.Material()

        if decay_time is None:
            if input is None:
                output.mat_pointer[0] = (<cpp_fccomp.FCComp *> self._inst).calc()
            elif isinstance(input, dict):
                output.mat_pointer[0] = (<cpp_storage.Storage *> self._inst).calc(conv.dict_to_map_int_dbl(input))
            elif isinstance(input, pyne.material._Material):
                in_mat = input
                output.mat_pointer[0] = (<cpp_storage.Storage *> self._inst).calc(<pyne.cpp_material.Material> in_mat.mat_pointer[0])
            else:
                raise TypeError("'input' must be a Material, dict, or None.")
        else:
            if input is None:
                output.mat_pointer[0] = (<cpp_storage.Storage *> self._inst).calc(<double> decay_time)
            elif isinstance(input, dict):
                output.mat_pointer[0] = (<cpp_storage.Storage *> self._inst).calc(conv.dict_to_map_int_dbl(input), <double> decay_time)
            elif isinstance(input, pyne.material._Material):
                in_mat = input
                output.mat_pointer[0] = (<cpp_storage.Storage *> self._inst).calc(<pyne.cpp_material.Material> in_mat.mat_pointer[0], <double> decay_time)
            else:
                raise TypeError("'input' must be a Material, dict, or None.")

        return output


