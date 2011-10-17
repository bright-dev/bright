"""Python header for reactor1g library."""
cimport fccomp

cdef class Reactor1G(fccomp.FCComp):
    cdef public object _fuel_chemical_form
    cdef public object _coolant_chemical_form
    cdef public object _niF
    cdef public object _niC
    cdef public object _miF
    cdef public object _miC
    cdef public object _NiF
    cdef public object _NiC
