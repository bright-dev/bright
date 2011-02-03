"""Python wrapper for isoname library."""
cimport cpp_isoname

def zzaaam_2_MCNP(int iso):
    return cpp_isoname.zzaaam_2_MCNP(iso)
