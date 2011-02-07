"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map
from libcpp.utility cimport pair
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
cimport std
cimport cpp_isoname


#
# Conversion dictionaries
#
LLzz = {}
cdef map[std.string, int].iterator LLzz_iter = cpp_isoname.LLzz.begin()
while LLzz_iter != cpp_isoname.LLzz.end():
    LLzz[deref(LLzz_iter).first.c_str()] = deref(LLzz_iter).second
    inc(LLzz_iter)

zzLL = {}
cdef map[int, std.string].iterator zzLL_iter = cpp_isoname.zzLL.begin()
while zzLL_iter != cpp_isoname.zzLL.end():
    zzLL[deref(zzLL_iter).first] = deref(zzLL_iter).second.c_str()
    inc(zzLL_iter)


#
# Current Form Function
#

def CurrentForm(nuc):
    """Find the current form of a nuclide.

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * form_flag (str): The form identifier string from ["zzaaam", "LLAAAM", "MCNP"].
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
    """Converts a nuclide from LLAAAM form to its zzaaam form. 

    Args:
        * nuc (str): Input nuclide in LLAAAM form.

    Returns:
        * newnuc (int): Output nuclide in zzaaam form.
    """
    return cpp_isoname.LLAAAM_2_zzaaam(std.string(nuc))


def LLAAAM_2_MCNP(char * nuc):
    """Converts a nuclide from LLAAAM form to its MCNP form. 

    Args:
        * nuc (str): Input nuclide in LLAAAM form.

    Returns:
        * newnuc (int): Output nuclide in MCNP form.
    """
    return cpp_isoname.LLAAAM_2_MCNP(std.string(nuc))

#
# zzaaam_2_* Functions 
#

def zzaaam_2_LLAAAM(int nuc):
    """Converts a nuclide from zzaaam form to its LLAAAM form. 

    Args:
        * nuc (str): Input nuclide in zzaaam form.

    Returns:
        * newnuc (int): Output nuclide in LLAAAM form.
    """
    cdef std.string cpp_LLAAAM = cpp_isoname.zzaaam_2_LLAAAM(nuc)
    return cpp_LLAAAM.c_str()


def zzaaam_2_MCNP(int nuc):
    """Converts a nuclide from zzaaam form to its MCNP form. 

    Args:
        * nuc (str): Input nuclide in zzaaam form.

    Returns:
        * newnuc (int): Output nuclide in MCNP form.
    """
    return cpp_isoname.zzaaam_2_MCNP(nuc)


#
# MCNP_2_* Functions
#

def MCNP_2_zzaaam(int nuc):
    """Converts a nuclide from MCNP form to its zzaaam form. 

    Args:
        * nuc (str): Input nuclide in MCNP form.

    Returns:
        * newnuc (int): Output nuclide in zzaaam form.
    """
    return cpp_isoname.MCNP_2_zzaaam(nuc)


def MCNP_2_LLAAAM(int nuc):
    """Converts a nuclide from MCNP form to its LLAAAM form. 

    Args:
        * nuc (str): Input nuclide in MCNP form.

    Returns:
        * newnuc (int): Output nuclide in LLAAAM form.
    """
    cdef std.string cpp_LLAAAM = cpp_isoname.MCNP_2_LLAAAM(nuc)
    return cpp_LLAAAM.c_str()


#
# mixed_2_* Functions
#

def mixed_2_zzaaam(nuc):
    """Converts an arbitrary nuclide and its zzaaam form. 

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * newnuc (int): Output nuclide in zzaaam form.
    """

    if isinstance(nuc, basestring):
        newnuc = cpp_isoname.mixed_2_zzaaam(std.string(nuc))
    elif isinstance(nuc, int):
        newnuc = cpp_isoname.mixed_2_zzaaam(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return newnuc


def mixed_2_LLAAAM(nuc):
    """Converts an arbitrary nuclide and its LLAAAM form. 

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * newnuc (int): Output nuclide in LLAAAM form.
    """
    cdef std.string newnuc

    if isinstance(nuc, basestring):
        newnuc = cpp_isoname.mixed_2_LLAAAM(std.string(nuc))
    elif isinstance(nuc, int):
        newnuc = cpp_isoname.mixed_2_LLAAAM(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return newnuc.c_str()


def mixed_2_zzaaam(nuc):
    """Converts an arbitrary nuclide and its MCNP form. 

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * newnuc (int): Output nuclide in MCNP form.
    """

    if isinstance(nuc, basestring):
        newnuc = cpp_isoname.mixed_2_MCNP(std.string(nuc))
    elif isinstance(nuc, int):
        newnuc = cpp_isoname.mixed_2_MCNP(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return newnuc


#
# Helper Functions
#

def nuc_weight_zzaaam(int nuc):
    """Calculates the weight of a nuclide in [amu].

    Args:
        * nuc (int): Input nuclide.

    Returns:
        * weight (float): Atomic weight of this nuclide [amu].
    """
    return cpp_isoname.nuc_weight_zzaaam(nuc)


def nuc_weight(nuc):
    """Calculates the weight of a nuclide in [amu].

    Args:
        * nuc (int or str): Input nuclide.

    Returns:
        * weight (float): Atomic weight of this nuclide [amu].
    """
    if isinstance(nuc, basestring):
        weight = cpp_isoname.nuc_weight(std.string(nuc))
    elif isinstance(nuc, int):
        weight = cpp_isoname.nuc_weight(<int> nuc)
    else:
        raise TypeError("Nuclide not a string ot integer.")

    return weight
