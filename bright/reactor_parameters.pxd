"""Python header for reactor parameters library."""
cimport cpp_reactor_parameters

cdef class ReactorParameters:
    cdef cpp_reactor_parameters.ReactorParameters * rp_pointer
