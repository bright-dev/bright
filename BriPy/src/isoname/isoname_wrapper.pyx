"""Python wrapper for isoname library."""
cimport std
cimport cpp_isoname

def CurrentForm(char* name):
    cdef std.string cpp_CurrentForm = cpp_isoname.CurrentForm(std.string(name))
    return cpp_CurrentForm.c_str()


def zzaaam_2_MCNP(int iso):
    return cpp_isoname.zzaaam_2_MCNP(iso)

