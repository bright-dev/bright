"""Python header for fuel fabrication library."""
cimport fccomp

cdef class FuelFabrication(fccomp.FCComp):
    cdef public object _materials
    cdef public object _mass_weights_in
    cdef public object _mass_weights_out
    cdef public object _deltaRs


