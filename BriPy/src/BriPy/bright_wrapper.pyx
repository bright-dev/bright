"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
cimport std
cimport converters

cimport isoname
cimport cpp_mass_stream
cimport mass_stream_wrapper
cimport cpp_bright

#from MassStream import MassStream
#import MassStream

import os

#######################################
### FCComps Configuration namespace ###
#######################################

class isos2track(object):
    def __get_value__(self):
        value = set()
        cdef cpp_set[int].iterator isos_iter = cpp_bright.isos2track.begin()
        while isos_iter != cpp_bright.isos2track.end():
            value.add(deref(isos_iter))
            inc(isos_iter)
        return value

    def __set_value__(self, value):
        """cdef int iso_zz
        cpp_bright.isos2track.clear()
        for iso in value:
            iso_zz = isoname.mixed_2_zzaaam(iso)
            cpp_bright.isos2track.insert(iso_zz)
        """
        s = set([isoname.mixed_2_zzaaam(v) for v in value])
        converters.set_py_to_cpp(s, cpp_bright.isos2track)

    value = property(__get_value__, __set_value__)

# Make isos2track a singleton
isos2track = isos2track().value

# Load isos2track from file functions
def load_isos2track_hdf5(char * filename, char * datasetname="", bint clear=False):
    """This convience function tries to load the isos2track set from a dataset 
    in an HDF5 file.  The dataset *must* be of integer type.  String-based
    nuclide names are currently not supported. 

    Args:
        * filename (str): Path to the data library.
        * dataset (str):  Dataset name to grab nuclides from.
        * clear (bool):   Flag that if set removes the currrent entries
          from isos2track prior to loading in new values.

    If the dataset argument is not provided or empty, the function tries to 
    load from various default datasets in the following order::

        "/isos2track"  
        "/Isos2Track"
        "/isostrack"   
        "/IsosTrack"
        "/isotrack"   
        "/IsoTrack"    
        "/ToIso"
        "/ToIsos"
        "/ToIso_zz"
        "/ToIso_MCNP"
        "/FromIso"  
        "/FromIsos"  
        "/FromIso_zz" 
        "/FromIso_MCNP"
    """
    cpp_bright.load_isos2track_hdf5(std.string(filename), std.string(datasetname), clear)


def load_isos2track_text(char * filename, bint clear=False):
    """This convience function tries to load the isos2track set from a text
    file.  The nuclide names may use any naming convention.  Mixing different
    conventions in the same file is allowed.  Whitespace is required between
    isotopic names.

    Args:
        * filename (str): Path to the data library.
        * clear (bool):   Flag that if set removes the currrent entries
          from isos2track prior to loading in new values.
    """
    cpp_bright.load_isos2track_text(std.string(filename), clear)


# Simple settings
verbosity = cpp_bright.verbosity
write_hdf5 = cpp_bright.write_hdf5
write_text = cpp_bright.write_text

class output_filename(object):
    def __get_value__(self):
        cdef std.string value = cpp_bright.output_filename
        return value.c_str()

    def __set_value__(self, char * value):
        cpp_bright.output_filename = std.string(value)

    value = property(__get_value__, __set_value__)

# Make isos2track a singleton
output_filename = output_filename().value



####################
### FCComp Class ###
####################

cdef class FCComp:
    """Base Fuel Cycle Component Class.

    Args:
        * paramlist (sequence of str): A set of parameter names (str) that the component will track.
        * name (str): The name of the fuel cycle component instance.

    Note that this automatically calls the protected :meth:`initialize` C function.
    """

    cdef cpp_bright.FCComp * fccomp_pointer

    def __cinit__(self, params=None, char * name=""):
        cdef cpp_set[std.string] param_set

        if params is None:
            self.fccomp_pointer = new cpp_bright.FCComp(std.string(name))
        else:
            for p in params:
                param_set.insert(std.string(p))
            self.fccomp_pointer = new cpp_bright.FCComp(param_set, std.string(name))


    def __dealloc__(self):
        del self.fccomp_pointer


    #
    # Class Attributes
    #

    property name:
        def __get__(self):
            cdef std.string n = self.fccomp_pointer.name
            return n.c_str()

        def __set__(self, char * n):
            self.fccomp_pointer.name = std.string(n)


    property natural_name:
        def __get__(self):
            cdef std.string n = self.fccomp_pointer.natural_name
            return n.c_str()

        def __set__(self, char * n):
            self.fccomp_pointer.natural_name = std.string(n)


    property IsosIn:
        def __get__(self):
            cdef mass_stream_wrapper.MassStream py_ms = mass_stream_wrapper.MassStream()
            py_ms.ms_pointer[0] = self.fccomp_pointer.IsosIn
            return py_ms

        def __set__(self, mass_stream_wrapper.MassStream ms):
            self.fccomp_pointer.IsosIn = <cpp_mass_stream.MassStream> ms.ms_pointer[0]


    property IsosOut:
        def __get__(self):
            cdef mass_stream_wrapper.MassStream py_ms = mass_stream_wrapper.MassStream()
            py_ms.ms_pointer[0] = self.fccomp_pointer.IsosOut
            return py_ms

        def __set__(self, mass_stream_wrapper.MassStream ms):
            self.fccomp_pointer.IsosOut = <cpp_mass_stream.MassStream> ms.ms_pointer[0]


    property ParamsIn:
        def __get__(self):
            pi = {}
            cdef cpp_map[std.string, double].iterator pi_iter = self.fccomp_pointer.ParamsIn.begin()
            while pi_iter != self.fccomp_pointer.ParamsIn.end():
                pi[deref(pi_iter).first.c_str()] = deref(pi_iter).second
                inc(pi_iter)
            return pi

        def __set__(self, dict pi):
            """cdef cpp_map[std.string, double] cpp_pi = self.fccomps_pointer.ParamsIn
            cpp_pi.clear()
            for key, value in pi.items():
                cpp_pi[std.string(key)] = <double> value
            """
            pass

#            cdef std.string cpp_key
#            cdef cpp_map[std.string, double] cpp_pi = cpp_map[std.string, double]()
#            for key, value in pi.items():
#                #cpp_key = std.string(key)
#               #cpp_pi[cpp_key] = <double> value
#               cpp_pi[std.string(key)] = value
#            self.fccomps_pointer.ParamsIn = cpp_pi
