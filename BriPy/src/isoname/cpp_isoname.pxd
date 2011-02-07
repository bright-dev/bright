"""C++ wrapper for isoname library."""
cimport std

cdef extern from "isoname.h" namespace "isoname":
    std.string CurrentForm(std.string)

    int zzaaam_2_MCNP(int)
