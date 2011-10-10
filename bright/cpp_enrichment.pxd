"""Cython header for  library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

from pyne cimport std
from pyne cimport cpp_material
from pyne cimport material

#from bright cimport cpp_fccomp
#cimport bright.cpp_fccomp

from bright.cpp_fccomp cimport FCComp

#cdef extern from "fccomp.h" namespace "bright":
#    cdef cppclass FCComp:
#        pass

cdef extern from "enrichment.h" namespace "bright":

    cdef cppclass EnrichmentParameters:
        # Constructors
        EnrichmentParameters() except +

        # Attributes
        double alpha_0
        double Mstar_0

        int j
        int k

        double N0
        double M0

        double xP_j
        double xW_j

    EnrichmentParameters fillUraniumEnrichmentDefaults() except +

    cdef cppclass Enrichment(FCComp): 
        # Constructors
        Enrichment() except +
        Enrichment(std.string) except +
        Enrichment(EnrichmentParameters, std.string) except +

        # Attributes
        double alpha_0
        double Mstar_0
        double Mstar
        cpp_material.Material mat_tail

        int j
        int k
        double xP_j
        double xW_j

        double N
        double M
        double N0
        double M0

        double TotalPerFeed
        double SWUperFeed
        double SWUperProduct

        # Methods
        void initialize(EnrichmentParameters) except +
        void calc_params () except +
        cpp_material.Material calc () except +
        cpp_material.Material calc (map[int, double]) except +
        cpp_material.Material calc (cpp_material.Material) except +

        double PoverF (double, double, double) except +
        double WoverF (double, double, double) except +

        # The following are methods I am too lazy to expose to Python
        # FIXME
        #double get_alphastar_i (double)

        #double get_Ei (double)
        #double get_Si (double)
        #void FindNM()

        #double xP_i(int)
        #double xW_i(int)
        #void SolveNM()
        #void Comp2UnitySecant()
        #void Comp2UnityOther()
        #double deltaU_i_OverG(int)
        #void LoverF()
        #void MstarOptimize()



