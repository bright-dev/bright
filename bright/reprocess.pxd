"""Python header for reprocess library."""
cimport fccomp

cdef class Reprocess(fccomp.FCComp):
    cdef public object _sepeff


