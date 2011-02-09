# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# Local imports
cimport std


#
# Map conversions
#

# <string, int> conversions

cdef void dict_to_map(dict pydict, cpp_map[std.string, int] cppmap):
    #cdef cpp_map[std.string, int] cppdict = cpp_map[std.string, int]()
    #for key, value in pydict.items():
    #    cppdict[std.string(key)] = value
    #cppmap = cpppdict
    cppmap.clear()
    for key, value in pydict.items():
        cppmap[std.string(key)] = value


cdef dict map_to_dict(cpp_map[std.string, int] cppmap):
    pydict = {}
    cdef map[std.string, int].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        LLzz[deref(mapiter).first.c_str()] = deref(mapiter).second
        inc(mapiter)

    return pydict


# <int, string> conversions

cdef dict map_to_dict(cpp_map[int, std.string] cppmap):
    pydict = {}
    cdef map[int, std.string].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        LLzz[deref(mapiter).first] = deref(mapiter).second.c_str()
        inc(mapiter)

    return pydict



#
# Set conversions
#

# Template primitive type sets

cdef void py_to_cpp_set(set pyset, cpp_set[T] cppset):
    cppset.clear()
    for item in pyset:
        cppset.insert(item)


cdef set cpp_to_py_set(cpp_set[T] cppset):
    pyset = set()
    cdef cpp_set[T].iterator setiter = cppset.begin()
    while setiter != cppset.end():
        pyset.add(deref(setiter))
        inc(setiter)
    return pyset


# String sets

cdef void py_to_cpp_set(set pyset, cpp_set[std.string] cppset):
    cppset.clear()
    cdef std.string s

    for item in pyset:
        s = std.string(item)
        cppset.insert(s)


cdef set cpp_to_py_set(cpp_set[std.string] cppset):
    pyset = set()
    cdef cpp_set[std.string].iterator setiter = cppset.begin()

    while setiter != cppset.end():
        pyset.add(deref(setiter).c_str())
        inc(setiter)

    return pyset
