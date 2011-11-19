from pyne import nucname

def load_nuc_file(path):
    """Takes a file that contains whitespace separated nuclide names and 
    returns the zzaaam representation as a sorted list."""
    with open(path, 'r') as f:
        s = f.read()

    nuc_list = [nucname.zzaaam(nuc) for nuc in s.split()]
    nuc_list.sort()
    return nuc_list


def temperature_flag(T):
    """Converts a temperature into the proper continuous energy 
    flag used in ACE and MCNP files.

    Parameters
    ----------
    T : number
        Temperature, multiple of 300 K.

    Returns
    -------
    temp_flag : 3-character string
        E.g. '06c' for 600 K.
    """

    t = int(T)

    # Check temperature value validity
    if t%300 != 0:
        raise ValueError("The temperature value must be a multiple of 300 K!")
    elif t <= 0:
        raise ValueError("The temperature value must be positive!")
    elif 9999 < t:
        raise ValueError("The temperature value must less than 10000 K!")

    # Make the temperature flag
    temp_flag = "{0:02}c".format(t/100)

    return temp_flag
