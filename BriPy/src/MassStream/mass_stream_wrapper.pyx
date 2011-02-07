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
        """MassStream C++ constuctor."""
        cdef cpp_map[int, double] comp_dict

        if isinstance(isovec, dict):
            # Mass Stream from dict
            comp_dict = cpp_map[int, double]()
            for key, value in isovec.items():
                comp_dict[key] = value
            self.ms_pointer = new cpp_mass_stream.MassStream(comp_dict, mass, std.string(name))

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
        """MassStream C++ destructor."""
        del self.ms_pointer


    #
    # Class Attributes
    #

    property comp:
        def __get__(self):
            comp_dict = {}
            cdef cpp_map[int, double].iterator comp_iter = self.ms_pointer.comp.begin()
            while comp_iter != self.ms_pointer.comp.end():
                comp_dict[deref(comp_iter).first] = deref(comp_iter).second
            return comp_dict

        def __set__(self, dict comp):
            cdef cpp_map[int, double] comp_dict = cpp_map[int, double]()
            for key, value in comp.items():
                comp_dict[key] = value
            self.ms_pointer.comp = comp_dict


    property mass:
        def __get__(self):
            return self.ms_pointer.mass

        def __set__(self, double mass):
            self.ms_pointer.mass = mass


    property name:
        def __get__(self):
            cdef std.string ms_name = self.ms_pointer.name
            return ms_name.c_str()

        def __set__(self, char * name):
            self.ms_pointer.name = std.string(name)
