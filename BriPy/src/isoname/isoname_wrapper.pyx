"""Python wrapper for isoname library."""
cimport std
cimport cpp_isoname


#def CurrentForm(char* nuc):
#    cdef std.string cpp_CurrentForm = cpp_isoname.CurrentForm(std.string(nuc))
#    return cpp_CurrentForm.c_str()
    

#def CurrentForm(int nuc):
#    cdef std.string cpp_CurrentForm = cpp_isoname.CurrentForm(nuc)
#    return cpp_CurrentForm.c_str()


def CurrentForm(nuc):
    """Find the current form of a nuclide."""
    cdef std.string cpp_CurrentForm 

    if isinstance(nuc, basestring):
        cpp_CurrentForm = cpp_isoname.CurrentForm(std.string(nuc))
    elif isinstance(nuc, int):
        cpp_CurrentForm = cpp_isoname.CurrentForm(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return cpp_CurrentForm.c_str()
    

def zzaaam_2_MCNP(int nuc):
    return cpp_isoname.zzaaam_2_MCNP(nuc)

