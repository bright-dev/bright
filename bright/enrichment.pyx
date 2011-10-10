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
            self._inst = new cpp_enrichment.Enrichment(std.string(name))
        elif isinstance(enrich_params, EnrichmentParameters):
            enr_par = enrich_params
            self._inst = new cpp_enrichment.Enrichment(<cpp_enrichment.EnrichmentParameters> enr_par.ep_pointer[0], std.string(name))


    #
    # Class Attributes
    #

    # Enrichment Attributes

    property alpha_0:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).alpha_0

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).alpha_0 = <double> value


    property Mstar_0:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).Mstar_0

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).Mstar_0 = <double> value


    property Mstar:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).Mstar

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).Mstar = <double> value


    property mat_tail:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_enrichment.Enrichment *> self._inst).mat_tail
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_enrichment.Enrichment *> self._inst).mat_tail = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property j:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).j

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).j = <int> value


    property k:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).k

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).k = <int> value


    property xP_j:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).xP_j

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).xP_j = <double> value


    property xW_j:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).xW_j

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).xW_j = <double> value


    property N:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).N

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).N = <double> value


    property M:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).M

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).M = <double> value


    property N0:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).N0

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).N0 = <double> value


    property M0:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).M0

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).M0 = <double> value


    property TotalPerFeed:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).TotalPerFeed

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).TotalPerFeed = <double> value


    property SWUperFeed:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).SWUperFeed

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).SWUperFeed = <double> value


    property SWUperProduct:
        def __get__(self):
            return (<cpp_enrichment.Enrichment *> self._inst).SWUperProduct

        def __set__(self, value):
            (<cpp_enrichment.Enrichment *> self._inst).SWUperProduct = <double> value



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
        (<cpp_enrichment.Enrichment *> self._inst).initialize(<cpp_enrichment.EnrichmentParameters> enr_par.ep_pointer[0])


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
        (<cpp_enrichment.FCComp *> self._inst).calc_params()


    def calc(self, input=None):
        """This method performs an optimization calculation on M* and solves for 
        appropriate values for all Enrichment attributes.  This includes the 
        product and waste streams flowing out of the the cascade as well.

        Args:
            * input (dict or Material or None): If input is present, it is set as the component's 
            mat_feed.  If input is a isotopic dictionary (zzaaam keys, float values), this dictionary 
            is first converted into a Material before being set as mat_feed.

        Returns:
            * output (Material): mat_prod.

        """
        cdef pyne.material._Material in_mat 
        cdef pyne.material._Material output = pyne.material.Material()

        if input is None:
            output.mat_pointer[0] = (<cpp_enrichment.FCComp *> self._inst).calc()
        elif isinstance(input, dict):
            output.mat_pointer[0] = (<cpp_enrichment.Enrichment *> self._inst).calc(conv.dict_to_map_int_dbl(input))
        elif isinstance(input, pyne.material._Material):
            in_mat = input
            output.mat_pointer[0] = (<cpp_enrichment.Enrichment *> self._inst).calc(<pyne.cpp_material.Material> in_mat.mat_pointer[0])

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
        return (<cpp_enrichment.Enrichment *> self._inst).PoverF(x_F, x_P, x_W)


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
        return (<cpp_enrichment.Enrichment *> self._inst).WoverF(x_F, x_P, x_W)

