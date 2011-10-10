"""Python header for fluence point library."""
cimport cpp_fluence_point

cdef class FluencePoint:
    cdef cpp_fluence_point.FluencePoint * fp_pointer

