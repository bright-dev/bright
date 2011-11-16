"""Python wrapper for fluence point."""
cimport cpp_fluence_point
cimport fluence_point


##########################
### FluencePoint Class ###
##########################


cdef class FluencePoint:
    """This class holds three simple data points that represent a fluence point.
    """

    def __cinit__(self):
        self.fp_pointer = new cpp_fluence_point.FluencePoint()

    def __dealloc__(self):
        del self.fp_pointer


    #
    # Class Attributes
    #

    property f:
        """Index (int) of fluence immediately lower than the value of F."""
        def __get__(self):
            return self.fp_pointer.f

        def __set__(self, int value):
            self.fp_pointer.f = value


    property F:
        """Fluence value itself (float).  In units of [neutrons/kilobarn], abbr [n/kb]."""
        def __get__(self):
            return self.fp_pointer.F

        def __set__(self, double value):
            self.fp_pointer.F = value


    property m:
        """The slope (float) dBU/dF between points f and f+1.  Has the odd units of [MWd kb / kgIHM n]"""
        def __get__(self):
            return self.fp_pointer.m

        def __set__(self, double value):
            self.fp_pointer.m = value



