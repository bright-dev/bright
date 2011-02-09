# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# Local imports
cimport std

cdef extern from "converters.pyx" void set_py_to_cpp(set, cpp_set[int])
