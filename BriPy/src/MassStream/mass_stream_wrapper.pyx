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

    #
    # Class Methods
    #

    def norm_comp_dict(self):
        """Normalizes the comp(osition), preserving the mass of the isotopic vector as mass."""
        self.ms_pointer.norm_comp_dict()


    def load_from_hdf5(self, char * filename, char * groupname, int row=-1):
        """A MassStream object may be initialized from an HDF5 file.
        The HDF5 representation of a MassStream is a group that holds several 
        extendable array datasets.  One array is entitled "Mass" while the other datasets
        are nuclide names in LLAAAM form ("U235", "NP237", *etc*).  For example::

            File.h5 (file)
                |-- MassStream (group)
                    |-- Mass (array)
                    |-- H1 (array)
                    |-- O16 (array)
                    |-- U235 (array)
                    |-- PU239 (array)
                    |-- ...

        The arrays are all of length N, where each row typically represents a different 
        fuel cycle pass.  The sum of all of the nuclide arrays should sum to one, like 
        MassStream.comp. 

        Args:
            * filename  (str): Path to HDF5 file that contains the data to read in.    
            * groupname (str): Path to HDF5 group that represents the data. 
              In the above example, groupname = "/MassStream".    

        Keyword Args:
            * row (int): The index of the arrays from which to read the data.  This 
              ranges from 0 to N-1.  Defaults to the last element of the array.
              Negative indexing is allowed (row[-N] = row[0]).

        Usage:
            This function loads data into a pre-existing :class:`MassStream`.  
            Initialization is therefore a two-step process::

                ms = MassStream()
                ms.load_from_hdf5("afile.h5", "/foo/bar/ms", -3)
        """
        self.ms_pointer.load_from_hdf5(std.string(filename), std.string(groupname), row)


    def load_from_text(self, char * filename):
        """A MassStream object may be initialized from a simple text file.
        The text representation of MassStreams are nuclide identifiers in the 
        first column and mass or weight values in the second column.  For example, 
        for natural uranium::

            922340  0.000055
            U235    0.00720
            92238   0.992745

        Data in this file must be whitespace separated.  Any valid nuclide naming
        scheme may be used for any isotope.

        Args:
            * filename (str): Path to HDF5 file that contains the data to read in.    

        Usage:
            This function loads data into a pre-existing MassStream.  
            Initialization is therefore a two-step process::

            ms = MassStream()
            ms.load_from_text("natu.h5")

        This function is most often called implicitly by the MassStream constructor.
        """
        self.ms_pointer.load_from_text(filename)


    def Print(self):
        """This prints a string representation of the MassStream to stdout.  Print is 
        particularly useful in C++.  In Python, this method simply duplicates 
        the functionality you would get from the built-in str() function.
        """
        self.ms_pointer.Print()


    def Normalize(self):
        """This convenience function normalizes the mass stream by setting its mass = 1.0."""
        self.ms_pointer.Normalize()


    def multByMass(self):
        """This function multiplies comp by mass and returns the resultant isotopic vector.

        Returns:
            * isovec(dict): For a MassStream ms, 

              .. math:: \mbox{isovec[iso]} = \mbox{ms.comp[iso]} \times \mbox{ms.mass}
        """
        isovec = {}
        cdef cpp_map[int, double] cpp_isovec = self.ms_pointer.multByMass()
        cdef cpp_map[int, double].iterator isovec_iter = cpp_isovec.begin()
        while isovec_iter != cpp_isovec.end():
            isovec[deref(isovec_iter).first] = deref(isovec_iter).second
        return isovec


    def atomic_weight(self):
        """This method returns the atomic weight of the comp of this MassStream.  Note that this is 
        only a rough estimate since this function is not yet coupled with measured atomic weights.

        Returns:
            * atomic_weight (float): Atomic weight in [amu]."""
        return self.ms_pointer.atomic_weight()
