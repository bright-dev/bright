#cdef extern from "string" namespace "std":
#    cdef cppclass cpp_string:
#        char* c_str()

#cdef extern from "string":
#    ctypedef struct std_string "std::string":
#        pass
#
#    cdef std_string charp_to_stdstring "std::string"(char*)
##    cdef char* stdstring_to_charp "std::string".c_str()


cdef extern from "<string>" namespace "std":
    cdef cppclass std_string:
        string()
        string(char *)
        char * c_str()


cdef extern from "isoname.h" namespace "isoname":
    std_string c_CurrentForm(std_string)

    int zzaaam_2_MCNP(int)


def CurrentForm(char* name):
    cdef std_string q = std_string(name)
    cdef std_string cpp_CurrentForm = c_CurrentForm(q)
    return cpp_CurrentForm.c_str()
