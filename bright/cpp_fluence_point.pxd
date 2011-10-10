"""Cython header for fluence point library."""

cdef extern from "fluence_point.h" namespace "bright":

    cdef cppclass FluencePoint:
        # Constructors        
        FluencePoint() except +

        # Attributes
        int f
        double F
        double m


