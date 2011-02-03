cdef extern from "string" namespace "std":
    cdef cppclass cpp_string:
        char* c_str()

cdef extern from "isoname.h" namespace "isoname":
    cpp_string CurrentForm(cpp_string)

    int zzaaam_2_MCNP(int)
