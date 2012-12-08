"""Python header for reactor parameters library."""
cimport cpp_reactor_parameters

cdef class ReactorParameters:
    cdef cpp_reactor_parameters.ReactorParameters * rp_pointer
    cdef public object _fuel_form
    cdef public object _cladding_form
    cdef public object _coolant_form
