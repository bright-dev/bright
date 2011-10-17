"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free

cimport numpy as np
import numpy as np

# pyne imports 
from pyne cimport std
from pyne cimport stlconverters as conv
from pyne import stlconverters as conv

cimport pyne.nucname
import pyne.nucname

cimport cpp_bright

import os

import bright.__init__

_local_dir = os.path.split(bright.__init__.__file__)[0]

lib = os.path.join(_local_dir, 'lib')
includes = os.path.join(_local_dir, 'includes')
bright_data = os.path.join(_local_dir, 'data')


######################################
### bright Configuration namespace ###
######################################

# Expose the C-code start up routine
def bright_start():
    # Specifiy the BRIGHT_DATA directory
    if "BRIGHT_DATA" not in os.environ:
        os.environ['BRIGHT_DATA'] = bright_data

    # Call the C-version of bright_start
    cpp_bright.bright_start()


# Run the appropriate start-up routines
bright_start()


#######################################
### FCComps Configuration namespace ###
#######################################

cdef class BrightConf:

    cdef object _track_nucs

    def __cinit__(self):
        self._track_nucs = None

    # From bright namespace

    property BRIGHT_DATA:
        def __get__(self):
            cdef std.string value = cpp_bright.BRIGHT_DATA
            return value.c_str()

        def __set__(self, char * value):
            cpp_bright.BRIGHT_DATA = std.string(value)
        

    property track_nucs:
        def __get__(self):
            cdef conv._SetInt proxy

            if self._track_nucs is None:
                proxy = conv.SetInt(False, False)
                proxy.set_ptr = &cpp_bright.track_nucs
                self._track_nucs = proxy

            return self._track_nucs

        def __set__(self, value):
            cdef cpp_set[int] s

            if isinstance(value, conv._SetInt):
                cpp_bright.track_nucs = deref((<conv._SetInt> value).set_ptr)
            elif hasattr(value, '__len__'):
                s = cpp_set[int]()
                for nuc in value:
                    s.insert(pyne.nucname.zzaaam(nuc))
                cpp_bright.track_nucs = s
            else:
                raise TypeError('{0} cannot be converted to a C++ set.'.format(type(value)))

            self._track_nucs = None


    property track_nucs_order:
        def __get__(self):
            return conv.vector_to_array_1d_int(cpp_bright.track_nucs_order)

        def __set__(self, value):
            s = set([pyne.nucname.zzaaam(v) for v in value])
            a = np.array(s)
            a.sort()
            cpp_bright.track_nucs = conv.py_to_cpp_set_int(s)
            cpp_bright.track_nucs_order = conv.array_to_vector_1d_int(a)


    property verbosity:
        def __get__(self):
            return cpp_bright.verbosity

        def __set__(self, int value):
            cpp_bright.verbosity = value


    property write_hdf5:
        def __get__(self):
            return cpp_bright.write_hdf5

        def __set__(self, bint value):
            cpp_bright.write_hdf5 = value


    property write_text:
        def __get__(self):
            return cpp_bright.write_text

        def __set__(self, bint value):
            cpp_bright.write_text = value


    property output_filename:
        def __get__(self):
            cdef std.string value = cpp_bright.output_filename
            return value.c_str()

        def __set__(self, char * value):
            cpp_bright.output_filename = std.string(value)


# Make a singleton of the Bright config object
bright_conf = BrightConf()


# Load track_nucs from file functions
def load_track_nucs_hdf5(char * filename, char * datasetname="", bint clear=False):
    """This convience function tries to load the track_nucs set from a dataset 
    in an HDF5 file.  The dataset *must* be of integer type.  String-based
    nuclide names are currently not supported. 

    Args:
        * filename (str): Path to the data library.
        * dataset (str):  Dataset name to grab nuclides from.
        * clear (bool):   Flag that if set removes the currrent entries
          from track_nucs prior to loading in new values.

    If the dataset argument is not provided or empty, the function tries to 
    load from various default datasets in the following order::

        "/track_nucs"  
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
    cpp_bright.load_track_nucs_hdf5(std.string(filename), std.string(datasetname), clear)


def load_track_nucs_text(char * filename, bint clear=False):
    """This convience function tries to load the track_nucs set from a text
    file.  The nuclide names may use any naming convention.  Mixing different
    conventions in the same file is allowed.  Whitespace is required between
    isotopic names.

    Args:
        * filename (str): Path to the data library.
        * clear (bool):   Flag that if set removes the currrent entries
          from track_nucs prior to loading in new values.
    """
    cpp_bright.load_track_nucs_text(std.string(filename), clear)



def sort_track_nucs():
    """This function sorts the track_nucs and places the results in track_nucs_order."""
    cpp_bright.sort_track_nucs()



