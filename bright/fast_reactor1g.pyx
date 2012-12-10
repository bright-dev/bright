"""Python wrapper for  ."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free
from libcpp.string cimport string as std_string

cimport numpy as np
import numpy as np

from pyne cimport nucname
from pyne cimport stlconverters as conv

cimport pyne.cpp_material
cimport pyne.material
import pyne.material

cimport cpp_fccomp
cimport cpp_reactor_parameters
cimport cpp_reactor1g
cimport cpp_fast_reactor1g

from bright.reactor_parameters cimport ReactorParameters
from bright.reactor_parameters import ReactorParameters

cimport reactor1g
import reactor1g



#######################
### Fast Reactor 1G ###
#######################


cdef class FastReactor1G(reactor1g.Reactor1G):
    """A One-Group Fast Reactor Fuel Cycle Component.  This is a daughter class of Reactor1G and 
    a granddaughter of FCComp.

    Parameters
    ----------
    libfile : str, optional
        The path the the FR HDF5 data library.  This value is set to Reactor1G.libfile and 
        used by Reactor1G.loadlib().
    reactor_parameters : ReactorParameters, optional
        The physical reactor parameter data to initialize this FR instance with.  If this argument 
        is not provided, default values are taken.
    name : str, optional
        The name of this FR instance.

    """

    def __cinit__(self, *args, **kwargs):
        pass


    def __init__(self, libfile=None, reactor_parameters=None, char * name="", *args, **kwargs):
        cdef ReactorParameters rp
        cdef std_string cpp_name = std_string(name)

        if (libfile is None) and (reactor_parameters is None):
            self._inst = new cpp_fast_reactor1g.FastReactor1G()
            (<cpp_fast_reactor1g.FastReactor1G *> self._inst).name = cpp_name

        elif isinstance(libfile, basestring) and (reactor_parameters is None):
            self._inst = new cpp_fast_reactor1g.FastReactor1G(std_string(<char *> libfile), cpp_name)

        elif (libfile is None) and isinstance(reactor_parameters, ReactorParameters):
            rp = reactor_parameters
            self._inst = new cpp_fast_reactor1g.FastReactor1G(<cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0], cpp_name)

        elif isinstance(libfile, basestring) and isinstance(reactor_parameters, ReactorParameters):
            rp = reactor_parameters
            self._inst = new cpp_fast_reactor1g.FastReactor1G(std_string(<char *> libfile), <cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0], cpp_name)

        else:
            if libfile is not None:
                raise TypeError("The libfile keyword must be a string or None.  Got " + str(type(libfile)))

            if reactor_parameters is not None:
                raise TypeError("The reactor_parameters keyword must be an instance of the ReactorParameters class or None.  Got " + str(type(reactor_parameters)))

    #
    # Class Methods
    # 

    # FastReactor1G Methods

    def calc_params(self):
        """Along with its own parameter set to track, the FR model implements its own function to set these
        parameters.  This function is equivalent to the following::

            self.params_prior_calc["BUd"]  = 0.0
            self.params_after_calc["BUd"] = self.BUd

            self.params_prior_calc["TRUCR"]  = 0.0
            self.params_after_calc["TRUCR"] = self.calc_tru_cr()

            self.params_prior_calc["P_NL"]  = 0.0
            self.params_after_calc["P_NL"] = self.P_NL

            self.params_prior_calc["U"]  = self.mat_feed_u.mass
            self.params_after_calc["U"] = self.mat_prod_u.mass

            self.params_prior_calc["TRU"]  = self.mat_feed_tru.mass
            self.params_after_calc["TRU"] = self.mat_prod_tru.mass

            self.params_prior_calc["ACT"]  = self.mat_feed_act.mass
            self.params_after_calc["ACT"] = self.mat_prod_act.mass

            self.params_prior_calc["LAN"]  = self.mat_feed_lan.mass
            self.params_after_calc["LAN"] = self.mat_prod_lan.mass

            self.params_prior_calc["FP"]  = 1.0 - self.mat_feed_act.mass  - self.mat_feed_lan.mass
            self.params_after_calc["FP"] = 1.0 - self.mat_prod_act.mass - self.mat_prod_lan.mass

        """
        (<cpp_fccomp.FCComp *> self._inst).calc_params()


