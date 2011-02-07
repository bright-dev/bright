"""Python wrapper for isoname library."""
cimport std
cimport cpp_isoname

#
# Current Form Function
#

def CurrentForm(nuc):
    """Find the current form of a nuclide.

    Args:
        * nuc (int or str): Nuclide whose 

    """
    cdef std.string cpp_CurrentForm 

    if isinstance(nuc, basestring):
        cpp_CurrentForm = cpp_isoname.CurrentForm(std.string(nuc))
    elif isinstance(nuc, int):
        cpp_CurrentForm = cpp_isoname.CurrentForm(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return cpp_CurrentForm.c_str()


#
# LLAAAM_2_* Functions
#

def LLAAAM_2_zzaaam(char * nuc):
    """Converts a nuclide """
    return cpp_isoname.LLAAAM_2_zzaaam(std.string(nuc))


def LLAAAM_2_MCNP(char * nuc):
    return cpp_isoname.LLAAAM_2_MCNP(std.string(nuc))


def zzaaam_2_MCNP(int nuc):
    return cpp_isoname.zzaaam_2_MCNP(nuc)

