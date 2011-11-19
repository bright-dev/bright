from pyne import nucname

def load_nuc_file(path):
    """Takes a file that contains whitespace separated nuclide names and 
    returns the zzaaam representation as a sorted list."""
    with open(path, 'r') as f:
        s = f.read()

    nuc_list = [nucname.zzaaam(nuc) for nuc in s.split()]
    nuc_list.sort()
    return nuc_list

