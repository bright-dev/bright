"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
cimport std

import isoname
cimport cpp_mass_stream
cimport cpp_bright

cdef class FCComp:
    cdef cpp_bright.FCComp * fccomp_pointer


cdef class EnrichmentParameters:
    cdef cpp_bright.EnrichmentParameters * ep_pointer


cdef class Enrichment(FCComp):
    cdef cpp_bright.Enrichment * e_pointer


cdef class Reprocess(FCComp):
    cdef cpp_bright.Reprocess * r_pointer


cdef class Storage(FCComp):
    cdef cpp_bright.Storage * s_pointer


cdef class FluencePoint:
    cdef cpp_bright.FluencePoint * fp_pointer


cdef class ReactorParameters:
    cdef cpp_bright.ReactorParameters * rp_pointer


cdef class Reactor1G(FCComp):
    cdef cpp_bright.Reactor1G * r1g_pointer


cdef class LightWaterReactor1G(Reactor1G):
    cdef cpp_bright.LightWaterReactor1G * lwr1g_pointer


cdef class FastReactor1G(Reactor1G):
    cdef cpp_bright.FastReactor1G * fr1g_pointer


cdef class FuelFabrication(FCComp):
    cdef cpp_bright.FuelFabrication * ff_pointer


cdef class ReactorMG(FCComp):
    cdef cpp_bright.ReactorMG * rmg_pointer


