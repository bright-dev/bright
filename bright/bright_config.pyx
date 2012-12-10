"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free
from libcpp.string cimport string as std_string

cimport numpy as np
import numpy as np

# pyne imports 
from pyne cimport stlconverters as conv
from pyne import stlconverters as conv

cimport pyne.nucname
import pyne.nucname

cimport cpp_bright

import os

import bright.__init__

prefix = os.path.split(bright.__init__.__file__)[0]

lib = os.path.join(prefix, 'lib')
includes = os.path.join(prefix, 'include')
bright_data = os.path.join(prefix, 'data')


######################################
### bright Configuration namespace ###
######################################

# Expose the C-code start up routine
def bright_start():
    """Handles bright initialization.  Reads in the 'BRIGHT_DATA' 
    environment variable."""
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
    """A class whose attributes expose C++ bright namespace variables.
    The bright_conf object is a singleton instance of this class."""

    cdef object _track_nucs

    def __cinit__(self):
        self._track_nucs = None

    # From bright namespace

    property BRIGHT_DATA:
        """Overide for directory path which is (by default) read in from an
        environmental variable of the same name."""
        def __get__(self):
            cdef std_string value = cpp_bright.BRIGHT_DATA
            return value.c_str()

        def __set__(self, char * value):
            cpp_bright.BRIGHT_DATA = std_string(value)
        

    property track_nucs:
        """Set of nuclides (in zzaaam form, see pyne) which components should track."""
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
        """Array nuclides which determines the order of track_nucs."""
        def __get__(self):
            return conv.vector_to_array_1d_int(cpp_bright.track_nucs_order)

        def __set__(self, value):
            s = set([pyne.nucname.zzaaam(v) for v in value])
            a = np.array(s)
            a.sort()
            cpp_bright.track_nucs = conv.py_to_cpp_set_int(s)
            cpp_bright.track_nucs_order = conv.array_to_vector_1d_int(a)


    property verbosity:
        """Determines the at which to print messages. Lower numbers mean fewer messages (default 0)."""
        def __get__(self):
            return cpp_bright.verbosity

        def __set__(self, int value):
            cpp_bright.verbosity = value


    property write_hdf5:
        """Boolean flag for whether to write binary HDF5 output."""
        def __get__(self):
            return cpp_bright.write_hdf5

        def __set__(self, bint value):
            cpp_bright.write_hdf5 = value


    property write_text:
        """Boolean flag for whether to write flat text file output."""
        def __get__(self):
            return cpp_bright.write_text

        def __set__(self, bint value):
            cpp_bright.write_text = value


    property output_filename:
        """Path to outputh file."""
        def __get__(self):
            cdef std_string value = cpp_bright.output_filename
            return value.c_str()

        def __set__(self, char * value):
            cpp_bright.output_filename = std_string(value)


# Make a singleton of the Bright config object
bright_conf = BrightConf()


# Load track_nucs from file functions
def load_track_nucs_hdf5(char * filename, char * datasetname="", bint clear=False):
    """This convience function tries to load the track_nucs set from a dataset 
    in an HDF5 file.  The dataset *must* be of integer type.  String-based
    nuclide names are currently not supported. 

    Parameters
    ----------
    filename : str
        Path to the data library.
    dataset : str, optional  
        Dataset name to grab nuclides from.
    clear : bool, optional
        Flag that if set removes the currrent entries from track_nucs prior to loading in new values.

    Notes
    -----
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
    cpp_bright.load_track_nucs_hdf5(std_string(filename), std_string(datasetname), clear)


def load_track_nucs_text(char * filename, bint clear=False):
    """This convience function tries to load the track_nucs set from a text
    file.  The nuclide names may use any naming convention.  Mixing different
    conventions in the same file is allowed.  Whitespace is required between
    isotopic names.

    Parameters
    ----------
    filename : str 
        Path to the data library.
    clear : bool
        Flag that if set removes the currrent entries from track_nucs prior to loading in new values.

    """
    cpp_bright.load_track_nucs_text(std_string(filename), clear)



def sort_track_nucs():
    """This function sorts the track_nucs and places the result in track_nucs_order."""
    cpp_bright.sort_track_nucs()




