"""Python header for reprocess library."""
cimport fccomp

cdef int convert_sepeff_key(object key)

cdef class Reprocess(fccomp.FCComp):
    cdef public object _sepeff


