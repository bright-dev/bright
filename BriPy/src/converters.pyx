# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# Local imports
cimport std


#
# Set conversions
#

cdef void py_to_cpp_set(set pyset, cpp_set[int] cppset):
    cppset.clear()
    for item in pyset:
        cppset.insert(item)


cdef set cpp_to_py_set(cpp_set[int] cppset):
    pyset = set()
    cdef cpp_set[int].iterator setiter = cppset.begin()
    while setiter != cppset.end():
        pyset.add(deref(setiter))
        inc(setiter)
    return pyset
