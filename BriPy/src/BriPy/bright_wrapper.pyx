"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
cimport std
cimport cpp_bright

import isoname


#
# FCComps Configuration namespace
#

class isos2track(object):
    def __get_value__(self):
        value = set()
        cdef cpp_set[int].iterator isos_iter = cpp_bright.isos2track.begin()
        while isos_iter != cpp_bright.isos2track.end():
            value.add(deref(isos_iter))
            inc(isos_iter)
        return value

    def __set_value__(self, value):
        cdef int iso_zz
        cpp_bright.isos2track.clear()
        for iso in value:
            iso_zz = isoname.mixed_2_zzaaam(iso)
            cpp_bright.isos2track.insert(iso_zz)

    value = property(__get_value__, __set_value__)

# Make isos2track a singleton
isos2track = isos2track().value

verbosity = cpp_bright.verbosity

