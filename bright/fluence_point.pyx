"""Python wrapper for fluence point."""
cimport cpp_fluence_point
cimport fluence_point


##########################
### FluencePoint Class ###
##########################


cdef class FluencePoint:
    """This class holds three simple data points that represent a fluence point.

    Attributes:
        * f (int): Index of Reactor1G.F immediately lower than the value of F (int).
        * F (float): Fluence value itself (float). In units of [n/kb] or [neutrons/kilobarn].
        * m (float): The slope dBU/dF between points f and f+1 (float). 
          Has the odd units of [MWd kb / kgIHM n].
    """

    def __cinit__(self):
        self.fp_pointer = new cpp_fluence_point.FluencePoint()

    def __dealloc__(self):
        del self.fp_pointer


    #
    # Class Attributes
    #

    property f:
        def __get__(self):
            return self.fp_pointer.f

        def __set__(self, int value):
            self.fp_pointer.f = value


    property F:
        def __get__(self):
            return self.fp_pointer.F

        def __set__(self, double value):
            self.fp_pointer.F = value


    property m:
        def __get__(self):
            return self.fp_pointer.m

        def __set__(self, double value):
            self.fp_pointer.m = value



