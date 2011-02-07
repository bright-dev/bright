"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
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
# Elemental string sets
#

LAN = set()
cdef cpp_set[std.string].iterator LAN_iter = cpp_isoname.LAN.begin()
while LAN_iter != cpp_isoname.LAN.end():
    LAN.add(deref(LAN_iter).c_str())
    inc(LAN_iter)

ACT = set()
cdef cpp_set[std.string].iterator ACT_iter = cpp_isoname.ACT.begin()
while ACT_iter != cpp_isoname.ACT.end():
    ACT.add(deref(ACT_iter).c_str())
    inc(ACT_iter)

TRU = set()
cdef cpp_set[std.string].iterator TRU_iter = cpp_isoname.TRU.begin()
while TRU_iter != cpp_isoname.TRU.end():
    TRU.add(deref(TRU_iter).c_str())
    inc(TRU_iter)

MA = set()
cdef cpp_set[std.string].iterator MA_iter = cpp_isoname.MA.begin()
while MA_iter != cpp_isoname.MA.end():
    MA.add(deref(MA_iter).c_str())
    inc(MA_iter)

FP = set()
cdef cpp_set[std.string].iterator FP_iter = cpp_isoname.FP.begin()
while FP_iter != cpp_isoname.FP.end():
    FP.add(deref(FP_iter).c_str())
    inc(FP_iter)


#
# Elemental integer sets
#

lan = set()
cdef cpp_set[int].iterator lan_iter = cpp_isoname.lan.begin()
while lan_iter != cpp_isoname.lan.end():
    lan.add(deref(lan_iter))
    inc(lan_iter)

act = set()
cdef cpp_set[int].iterator act_iter = cpp_isoname.act.begin()
while act_iter != cpp_isoname.act.end():
    act.add(deref(act_iter))
    inc(act_iter)

tru = set()
cdef cpp_set[int].iterator tru_iter = cpp_isoname.tru.begin()
while tru_iter != cpp_isoname.tru.end():
    tru.add(deref(tru_iter))
    inc(tru_iter)

ma = set()
cdef cpp_set[int].iterator ma_iter = cpp_isoname.ma.begin()
while ma_iter != cpp_isoname.ma.end():
    ma.add(deref(ma_iter))
    inc(ma_iter)

fp = set()
cdef cpp_set[int].iterator fp_iter = cpp_isoname.fp.begin()
while fp_iter != cpp_isoname.fp.end():
    fp.add(deref(fp_iter))
    inc(fp_iter)


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


def mixed_2_MCNP(nuc):
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


#
# (a*)_2_(b*)_List Functions
#

def LLAAAM_2_zzaaam_List(nuclist):
    """Converts a list of LLAAAM form to a list of zzaaam form.

    Args:
        * `nuclist` (str list): List of LLAAAM nuclides.

    Returns:
        * `newnuclist` (int list): List of zzaaam nuclides.
    """
    return [LLAAAM_2_zzaaam(nuc) for nuc in nuclist]


def LLAAAM_2_MCNP_List(nuclist):
    """Converts a list of LLAAAM form to a list of MCNP form.

    Args:
        * `nuclist` (str list): List of LLAAAM nuclides.

    Returns:
        * `newnuclist` (int list): List of MCNP nuclides.
    """
    return [LLAAAM_2_MCNP(nuc) for nuc in nuclist]


def zzaaam_2_LLAAAM_List(nuclist):
    """Converts a list of zzaaam form to a list of LLAAAM form.

    Args:
        * `nuclist` (int list): List of zzaaam nuclides.

    Returns:
        * `newnuclist` (str list): List of LLAAAM nuclides.
    """
    return [zzaaam_2_LLAAAM(nuc) for nuc in nuclist]


def zzaaam_2_MCNP_List(nuclist):
    """Converts a list of zzaaam form to a list of MCNP form.

    Args:
        * `nuclist` (int list): List of zzaaam nuclides.

    Returns:
        * `newnuclist` (int list): List of MCNP nuclides.
    """
    return [zzaaam_2_MCNP(nuc) for nuc in nuclist]


def MCNP_2_LLAAAM_List(nuclist):
    """Converts a list of MCNP form to a list of LLAAAM form.

    Args:
        * `nuclist` (int list): List of MCNP nuclides.

    Returns:
        * `newnuclist` (str list): List of LLAAAM nuclides.
    """
    return [MCNP_2_LLAAAM(nuc) for nuc in nuclist]


def MCNP_2_zzaaam_List(nuclist):
    """Converts a list of MCNP form to a list of zzaaam form.

    Args:
        * `nuclist` (int list): List of MCNP nuclides.

    Returns:
        * `newnuclist` (int list): List of zzaaam nuclides.
    """
    return [MCNP_2_zzaaam(nuc) for nuc in nuclist]


#
# mixed_2_*_List Functions
#

def RearRemoveDuplicates(l):
    """Removes duplicate entries from list l, starting from the back. 
    Used internally in the [form]_2_[form]_List() functions.

    Args:
       * `l` (list): input list.

    Returns:
       * `l` (list): input with duplicates removed.
    """

    for n in range(len(l)-1, -1, -1):
        if 1 < l.count(l[n]):
            l.pop(n)
    return l
    

def mixed_2_zzaaam_List(nuclist):
    """Converts a list of mixed form to a list of zzaaam form.

    Args:
        * `nuclist` (str or int list): List of nuclides of mixed form.

    Returns:
        * `newnuclist` (int list): List of zzaaam nuclides.
    """
    return RearRemoveDuplicates( [mixed_2_zzaaam(nuc) for nuc in nuclist] )


def mixed_2_LLAAAM_List(nuclist):
    """Converts a list of mixed form to a list of LLAAAM form.

    Args:
        * `nuclist` (str or int list): List of nuclides of mixed form.

    Returns:
        * `newnuclist` (str list): List of LLAAAM nuclides.
    """
    return RearRemoveDuplicates( [mixed_2_LLAAAM(nuc) for nuc in nuclist] )


def mixed_2_MCNP_List(nuclist):
    """Converts a list of mixed form to a list of MCNP form.

    Args:
        * `nuclist` (str or int list): List of nuclides of mixed form.

    Returns:
        * `newnuclist` (int list): List of MCNP nuclides.
    """
    return RearRemoveDuplicates( [mixed_2_MCNP(nuc) for nuc in nuclist] )

#
# isovec_keys_2_* Functions
#

def isovec_keys_2_zzaaam(isovec):
    """Converts all keys of an isotopic vector dictionary to zzaaam form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of zzaaam (int) form.
    """
    newvec = {}

    for iso in isovec.keys():
        newvec[mixed_2_zzaaam(iso)] = isovec[iso]

    return newvec


def isovec_keys_2_LLAAAM(isovec):
    """Converts all keys of an isotopic vector dictionary to LLAAAM form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of LLAAAM (str) form.
    """
    newvec = {}

    for iso in isovec.keys():
        newvec[mixed_2_LLAAAM(iso)] = isovec[iso]

    return newvec


def isovec_keys_2_MCNP(isovec):
    """Converts all keys of an isotopic vector dictionary to MCNP form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of MCNP (int) form.
    """
    newvec = {}

    for iso in isovec.keys():
        newvec[mixed_2_MCNP(iso)] = isovec[iso]

    return newvec

