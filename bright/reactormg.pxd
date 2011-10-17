"""Python header for reactormg library."""
cimport fccomp

cdef class ReactorMG(fccomp.FCComp):
    cdef public object _chemical_form_fuel
    cdef public object _chemical_form_clad
    cdef public object _chemical_form_cool

    cdef public object _I
    cdef public object _J
    cdef public object _K

    cdef public object _K_ind


