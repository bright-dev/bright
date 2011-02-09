"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
cimport std
#cimport converters

cimport isoname
cimport cpp_mass_stream
#cimport mass_stream_wrapper
#import mass_stream_wrapper
#cimport "../MassStream/mass_stream" as mass_stream
#cimport "../MassStream/mass_stream_wrapper" as mass_stream_wrapper
#cimport src.MassStream.mass_stream_wrapper as mass_stream_wrapper
cimport cpp_bright

cdef class FCComp:
    cdef cpp_bright.FCComp * fccomp_pointer
