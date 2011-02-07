"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
cimport std
cimport cpp_mass_stream


cdef class MassStream:
    cdef cpp_mass_stream.MassStream * ms_pointer

    def __cinit__(self, isovec=None, float mass=-1.0, char * name=''):
        cdef cpp_map[int, double] cd

        if isinstance(isovec, dict):
            # Mass Stream from dict
            cd = cpp_map[int, double]()
            for key, value in isovec.items():
                cd[key] = value
            self.ms_pointer = new cpp_mass_stream.MassStream(cd, mass, std.string(name))

        elif isinstance(isovec, basestring):
            # Mass Stream from file
            self.ms_pointer = new cpp_mass_stream.MassStream(<char *> isovec, mass, std.string(name))

        elif isovec is None:
            # Empty mass stream
            self.ms_pointer = new cpp_mass_stream.MassStream()

        else:
            # Bad MassStream 
            raise TypeError("The mass stream isotopic vector must be a dict, str, or None.")


    def __dealloc__(self):
        del self.ms_pointer


    property mass:
        def __get__(self):
            return self.ms_pointer.mass

        def __set__(self, double mass):
            self.ms_pointer.mass = mass
