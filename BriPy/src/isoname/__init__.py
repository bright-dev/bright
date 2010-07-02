import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from isoname import *

##################################
### (a*)_2_(b*)_List Functions ###
##################################

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


################################
### mixed_2_*_List Functions ###
################################
def RearRemoveDuplicates(l):
    """
    Removes duplicate entries from list l, starting from the back. 
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


#################################
### isovec_keys_2_* Functions ###
#################################
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

