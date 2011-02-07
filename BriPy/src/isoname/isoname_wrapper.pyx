"""Python wrapper for isoname library."""
#cimport cpp_isoname

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char * c_str()


cdef extern from "isoname.h" namespace "isoname":
    string CurrentForm(string)

    int zzaaam_2_MCNP(int)


def o_CurrentForm(char* name):
#    cdef std_string q = std_string(name)
    cdef string cpp_CurrentForm = CurrentForm(string(name))
    return cpp_CurrentForm.c_str()


def zzaaam_2_MCNP(int iso):
    return zzaaam_2_MCNP(iso)

