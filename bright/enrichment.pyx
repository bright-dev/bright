"""Python wrapper for enrichment."""
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
cimport cpp_enrichment

cimport fccomp
import fccomp



########################
### Enrichment Class ###
########################


cdef class EnrichmentParameters:
    """This class is a collection of values that mirror the attributes in 
    Enrichment that are required for the cascade model to run.
    In C-code this a simple `struct.  Like ReactorParameters, this class 
    takes no arguments on initialization.  An empty ErichmentParameters
    instance has all values (weakly) set to zero.
    """

    def __cinit__(self):
        self.ep_pointer = new cpp_enrichment.EnrichmentParameters()

    def __dealloc__(self):
        del self.ep_pointer


    #
    # Class Attributes
    #

    property alpha_0:
        def __get__(self):
            return self.ep_pointer.alpha_0

        def __set__(self, value):
            self.ep_pointer.alpha_0 = <double> value


    property Mstar_0:
        def __get__(self):
            return self.ep_pointer.Mstar_0

        def __set__(self, value):
            self.ep_pointer.Mstar_0 = <double> value


    property j:
        def __get__(self):
            return self.ep_pointer.j

        def __set__(self, value):
            self.ep_pointer.j = <int> value


    property k:
        def __get__(self):
            return self.ep_pointer.k

        def __set__(self, value):
            self.ep_pointer.k = <int> value


    property N0:
        def __get__(self):
            return self.ep_pointer.N0

        def __set__(self, value):
            self.ep_pointer.N0 = <double> value


    property M0:
        def __get__(self):
            return self.ep_pointer.M0

        def __set__(self, value):
            self.ep_pointer.M0 = <double> value


    property xP_j:
        def __get__(self):
            return self.ep_pointer.xP_j

        def __set__(self, value):
            self.ep_pointer.xP_j = <double> value


    property xW_j:
        def __get__(self):
            return self.ep_pointer.xW_j

        def __set__(self, value):
            self.ep_pointer.xW_j = <double> value



def uranium_enrichment_defaults():
    cdef cpp_enrichment.EnrichmentParameters cpp_ued = cpp_enrichment.fillUraniumEnrichmentDefaults()
    cdef EnrichmentParameters ued = EnrichmentParameters()
    ued.ep_pointer[0] = cpp_ued
    return ued



cdef class Enrichment(fccomp.FCComp):
    """Enrichment Fuel Cycle Component Class.  Daughter of bright.FCComp class.

    Args:
        * enrich_params (EnrichmentParameters): This specifies how the enrichment 
          cascade should be set up.  It is a EnrichmentParameters
          helper object.  If enrich_params is not specified, then the cascade 
          is initialized with UraniumEnrichmentDefaults.
        * name (str): The name of the enrichment fuel cycle component instance.

    Note that this automatically calls the public initialize() C function.
    """

    def __cinit__(self, enrich_params=None, char * name=""):
        cdef EnrichmentParameters enr_par

        if enrich_params is None:
            self.e_pointer = new cpp_enrichment.Enrichment(std.string(name))
        elif isinstance(enrich_params, EnrichmentParameters):
            enr_par = enrich_params
            self.e_pointer = new cpp_enrichment.Enrichment(<cpp_enrichment.EnrichmentParameters> enr_par.ep_pointer[0], std.string(name))

    def __dealloc__(self):
        del self.e_pointer


    #
    # Class Attributes
    #

    # Enrichment Attributes

    property alpha_0:
        def __get__(self):
            return self.e_pointer.alpha_0

        def __set__(self, value):
            self.e_pointer.alpha_0 = <double> value


    property Mstar_0:
        def __get__(self):
            return self.e_pointer.Mstar_0

        def __set__(self, value):
            self.e_pointer.Mstar_0 = <double> value


    property Mstar:
        def __get__(self):
            return self.e_pointer.Mstar

        def __set__(self, value):
            self.e_pointer.Mstar = <double> value


    property mat_tail:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = self.e_pointer.mat_tail
            return pymat

        def __set__(self, pyne.material._Material mat):
            self.e_pointer.mat_tail = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property j:
        def __get__(self):
            return self.e_pointer.j

        def __set__(self, value):
            self.e_pointer.j = <int> value


    property k:
        def __get__(self):
            return self.e_pointer.k

        def __set__(self, value):
            self.e_pointer.k = <int> value


    property xP_j:
        def __get__(self):
            return self.e_pointer.xP_j

        def __set__(self, value):
            self.e_pointer.xP_j = <double> value


    property xW_j:
        def __get__(self):
            return self.e_pointer.xW_j

        def __set__(self, value):
            self.e_pointer.xW_j = <double> value


    property N:
        def __get__(self):
            return self.e_pointer.N

        def __set__(self, value):
            self.e_pointer.N = <double> value


    property M:
        def __get__(self):
            return self.e_pointer.M

        def __set__(self, value):
            self.e_pointer.M = <double> value


    property N0:
        def __get__(self):
            return self.e_pointer.N0

        def __set__(self, value):
            self.e_pointer.N0 = <double> value


    property M0:
        def __get__(self):
            return self.e_pointer.M0

        def __set__(self, value):
            self.e_pointer.M0 = <double> value


    property TotalPerFeed:
        def __get__(self):
            return self.e_pointer.TotalPerFeed

        def __set__(self, value):
            self.e_pointer.TotalPerFeed = <double> value


    property SWUperFeed:
        def __get__(self):
            return self.e_pointer.SWUperFeed

        def __set__(self, value):
            self.e_pointer.SWUperFeed = <double> value


    property SWUperProduct:
        def __get__(self):
            return self.e_pointer.SWUperProduct

        def __set__(self, value):
            self.e_pointer.SWUperProduct = <double> value


    # FCComps inherited attributes

    property name:
        def __get__(self):
            cdef std.string n = self.e_pointer.name
            return n.c_str()

        def __set__(self, char * n):
            self.e_pointer.name = std.string(n)


    property natural_name:
        def __get__(self):
            cdef std.string n = self.e_pointer.natural_name
            return n.c_str()

        def __set__(self, char * n):
            self.e_pointer.natural_name = std.string(n)


    property mat_feed:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = self.e_pointer.mat_feed
            return pymat

        def __set__(self, pyne.material._Material mat):
            self.e_pointer.mat_feed = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = self.e_pointer.mat_prod
            return pymat

        def __set__(self, pyne.material._Material mat):
            self.e_pointer.mat_prod = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property params_prior_calc:
        def __get__(self):
            return conv.map_to_dict_str_dbl(self.e_pointer.params_prior_calc)

        def __set__(self, dict pi):
            self.e_pointer.params_prior_calc = conv.dict_to_map_str_dbl(pi)


    property params_after_calc:
        def __get__(self):
            return conv.map_to_dict_str_dbl(self.e_pointer.params_after_calc)

        def __set__(self, dict po):
            self.e_pointer.params_after_calc = conv.dict_to_map_str_dbl(po)


    property pass_num:
        def __get__(self):
            return self.e_pointer.pass_num

        def __set__(self, int pn):
            self.e_pointer.pass_num = pn


    property track_params:
        def __get__(self):
            return conv.cpp_to_py_set_str(self.e_pointer.track_params)

        def __set__(self, set p2t):
            self.e_pointer.track_params = conv.py_to_cpp_set_str(p2t)



    #
    # Class Methods
    # 

    def initialize(self, EnrichmentParameters enrich_params):
        """The initialize function takes an enrichment parameter object and sets
        the corresponding Enrichment attributes to the same value.

        Args:
            * enrich_params (EnrichmentParameters): A class containing the values to
              (re-)initialize an Enrichment cascade with.
        """
        cdef EnrichmentParameters enr_par = enrich_params
        self.e_pointer.initialize(<cpp_enrichment.EnrichmentParameters> enr_par.ep_pointer[0])


    def calc_params(self):
        """Here the parameters for Enrichment are set::

            self.params_prior_calc["MassFeed"]  = self.mat_feed.mass
            self.params_after_calc["MassFeed"] = 0.0

            self.params_prior_calc["MassProduct"]  = 0.0
            self.params_after_calc["MassProduct"] = self.mat_prod.mass

            self.params_prior_calc["MassTails"]  = 0.0
            self.params_after_calc["MassTails"] = self.mat_tail.mass

            self.params_prior_calc["N"]  = self.N
            self.params_after_calc["N"] = self.N

            self.params_prior_calc["M"]  = self.M
            self.params_after_calc["M"] = self.M

            self.params_prior_calc["Mstar"]  = self.Mstar
            self.params_after_calc["Mstar"] = self.Mstar

            self.params_prior_calc["TotalPerFeed"]  = self.TotalPerFeed
            self.params_after_calc["TotalPerFeed"] = self.TotalPerFeed

            self.params_prior_calc["SWUperFeed"]  = self.SWUperFeed
            self.params_after_calc["SWUperFeed"] = 0.0

            self.params_prior_calc["SWUperProduct"]  = 0.0
            self.params_after_calc["SWUperProduct"] = self.SWUperProduct

        """
        (<cpp_enrichment.FCComp *> self.e_pointer).calc_params()


    def calc(self, input=None):
        """This method performs an optimization calculation on M* and solves for 
        appropriate values for all Enrichment attributes.  This includes the 
        product and waste streams flowing out of the the cascade as well.

        Args:
            * input (dict or MassStream or None): If input is present, it is set as the component's 
            mat_feed.  If input is a isotopic dictionary (zzaaam keys, float values), this dictionary 
            is first converted into a MassStream before being set as mat_feed.

        Returns:
            * output (MassStream): mat_prod.

        """
        cdef pyne.material._Material in_mat 
        cdef pyne.material._Material output = pyne.material.Material()

        if input is None:
            output.mat_pointer[0] = (<cpp_enrichment.FCComp *> self.e_pointer).calc()
        elif isinstance(input, dict):
            output.mat_pointer[0] = self.e_pointer.calc(conv.dict_to_map_int_dbl(input))
        elif isinstance(input, pyne.material._Material):
            in_mat = input
            output.mat_pointer[0] = self.e_pointer.calc(<pyne.cpp_material.Material> in_mat.mat_pointer[0])

        return output


    def PoverF(self, double x_F, double x_P, double x_W):
        """Solves for the product over feed enrichment ratio.

        .. math::

            \\frac{p}{f} = \\frac{(x_F - x_W)}{(x_P - x_W)}

        Args:
            * x_F (float): Feed enrichment.
            * x_P (float): Product enrichment.
            * x_W (float): Waste enrichment.

        Returns:
            * pfratio (float): As calculated above.
        """
        return self.e_pointer.PoverF(x_F, x_P, x_W)


    def WoverF(self, double x_F, double x_P, double x_W):
        """Solves for the waste over feed enrichment ratio.

        .. math::

            \\frac{p}{f} = \\frac{(x_F - x_P)}{(x_W - x_P)}

        Args:
            * x_F (float): Feed enrichment.
            * x_P (float): Product enrichment.
            * x_W (float): Waste enrichment.

        Returns:
            * wfratio (float): As calculated above.
        """
        return self.e_pointer.WoverF(x_F, x_P, x_W)

