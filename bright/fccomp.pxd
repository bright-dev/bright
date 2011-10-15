"""Python header for  library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc


from pyne cimport std
from pyne cimport nucname
from pyne cimport cpp_material
from pyne cimport material

cimport cpp_fccomp

cdef class FCComp:
    cdef void * _inst
    cdef public bint _free_inst

    cdef public object _params_prior_calc
    cdef public object _params_after_calc
    cdef public object _track_params
