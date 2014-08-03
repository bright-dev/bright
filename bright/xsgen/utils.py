from __future__ import print_function

import os
import subprocess
from pprint import pformat

from pyne import nucname

USE_COLOR = (os.name is 'posix')
DEFAULT_RC_FILE = "defaultrc.py"
DEFAULT_PLUGINS = ("bright.xsgen.pre", "bright.xsgen.make_mc_input", 
                   "bright.xsgen.run_transport", "bright.xsgen.post", 
                   "bright.xsgen.test")

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


class RemoteConnection(object):
    def __init__(self, url='', user='', dir=''):
        self.url  = url
        self.user = user
        self.dir  = dir

    def run(self, cmd):
        callcmd = 'ssh {user}@{url} \"{remcmd}\"'.format(remcmd=cmd, **self.__dict__)
        return subprocess.call(callcmd, shell=True)

    def put(self, loc_file, rem_file):
        callcmd = "rsync -rh --partial --progress --rsh=ssh {lf} {user}@{url}:{rf}"
        callcmd = callcmd.format(lf=loc_file, rf=rem_file, **self.__dict__)
        return subprocess.call(callcmd, shell=True)

    def get(self, rem_file, loc_file):
        callcmd = "rsync -rh --partial --progress --rsh=ssh {user}@{url}:{rf} {lf}"
        callcmd = callcmd.format(lf=loc_file, rf=rem_file, **self.__dict__)
        return subprocess.call(callcmd, shell=True)

class NotSpecified():
    def __repr__(self):
        return "NotSpecified"

class RunControl(object):
    """A composable configuration class for xdress. Unlike argparse.Namespace,
    this keeps the object dictionary (__dict__) separate from the run control
    attributes dictionary (_dict)."""

    def __init__(self, **kwargs):
        """Parameters
        -------------
        kwargs : optional
        Items to place into run control.

        """
        self._dict = {}
        for k, v in kwargs.items():
            setattr(self, k, v)
        self._updaters = {}

    def __getattr__(self, key):
        if key in self._dict:
            return self._dict[key]
        elif key in self.__dict__:
            return self.__dict__[key]
        elif key in self.__class__.__dict__:
            return self.__class__.__dict__[key]
        else:
            msg = "RunControl object has no attribute {0!r}.".format(key)
            raise AttributeError(msg)

    def get(self, key, default=None):
        if key in self._dict:
            return self._dict[key]
        elif key in self.__dict__:
            return self.__dict__[key]
        elif key in self.__class__.__dict__:
            return self.__class__.__dict__[key]
        else:
            return default


    def __setattr__(self, key, value):
        if key.startswith('_'):
            self.__dict__[key] = value
        else:
            if value is NotSpecified and key in self:
                return
            self._dict[key] = value

    def __delattr__(self, key):
        if key in self._dict:
            del self._dict[key]
        elif key in self.__dict__:
            del self.__dict__[key]
        elif key in self.__class__.__dict__:
            del self.__class__.__dict__[key]
        else:
            msg = "RunControl object has no attribute {0!r}.".format(key)
            raise AttributeError(msg)

    def __iter__(self):
        return iter(self._dict)

    def __repr__(self):
        keys = sorted(self._dict.keys())
        s = ", ".join(["{0!s}={1!r}".format(k, self._dict[k]) for k in keys])
        return "{0}({1})".format(self.__class__.__name__, s)

    def _pformat(self):
        keys = sorted(self._dict.keys())
        f = lambda k: "{0!s}={1}".format(k, pformat(self._dict[k], indent=2))
        s = ",\n ".join(map(f, keys))
        return "{0}({1})".format(self.__class__.__name__, s)

    def __contains__(self, key):
        return key in self._dict or key in self.__dict__ or \
                                    key in self.__class__.__dict__

    def __eq__(self, other):
        if hasattr(other, '_dict'):
            return self._dict == other._dict
        elif isinstance(other, Mapping):
            return self._dict == other
        else:
            return NotImplemented

    def __ne__(self, other):
        if hasattr(other, '_dict'):
            return self._dict != other._dict
        elif isinstance(other, Mapping):
            return self._dict != other
        else:
            return NotImplemented

    def _update(self, other):
        """Updates the rc with values from another mapping. If this rc has
        if a key is in self, other, and self._updaters, then the updaters
        value is called to perform the update. This function should return
        a copy to be safe and not update in-place.
        """
        if hasattr(other, '_dict'):
            other = other._dict
        elif not hasattr(other, 'items'):
            other = dict(other)
        for k, v in other.items():
            if v is NotSpecified:
                pass
            elif k in self._updaters and k in self:
                v = self._updaters[k](getattr(self, k), v)
            setattr(self, k, v)

def exec_file(filename, glb=None, loc=None):
    """A function equivalent to the Python 2.x execfile statement."""
    import io
    with io.open(filename, 'r') as f:
        src = f.read()
    exec(compile(src, filename, "exec"), glb, loc)

nyansep = r'~\_/' * 17 + '~=[,,_,,]:3'

def indent(s, n=4, join=True):
    """Indents all lines in the string or list s by n spaces."""
    spaces = " " * n
    lines = s.splitlines() if isinstance(s, basestring) else s
    lines = lines or ()
    if join:
        return '\n'.join([spaces + l for l in lines if l is not None])
    else:
        return [spaces + l for l in lines if l is not None]

#
# filesystem helpers
#

def newoverwrite(s, filename, verbose=False):
    """Useful for not forcing re-compiles and thus playing nicely with the
    build system.  This is acomplished by not writing the file if the existsing
    contents are exactly the same as what would be written out.

    Parameters
    ----------
    s : str
        string contents of file to possible
    filename : str
        Path to file.
    vebose : bool, optional
        prints extra message

    """
    if os.path.isfile(filename):
        with io.open(filename, 'rb') as f:
            old = f.read()
        if s == old:
            return
    with io.open(filename, 'wb') as f:
        f.write(s.encode())
    if verbose:
        print("  wrote " + filename)

def newcopyover(f1, f2, verbose=False):
    """Useful for not forcing re-compiles and thus playing nicely with the
    build system.  This is acomplished by not writing the file if the existsing
    contents are exactly the same as what would be written out.

    Parameters
    ----------
    f1 : str
        Path to file to copy from
    f2 : str
        Path to file to copy over
    vebose : bool, optional
        prints extra message

    """
    if os.path.isfile(f1):
        with io.open(f1, 'r') as f:
            s = f.read()
        return newoverwrite(s, f2, verbose)

def writenewonly(s, filename, verbose=False):
    """Only writes the contents of the string to a file if the file does not exist.
    Useful for not tocuhing files.

    Parameters
    ----------
    s : str
        string contents of file to possible
    filename : str
        Path to file.
    vebose : bool, optional
        prints extra message

    """
    if os.path.isfile(filename):
        return
    with open(filename, 'w') as f:
        f.write(str(s))
    if verbose:
        print("  wrote " + filename)

def ensuredirs(f):
    """For a file path, ensure that its directory path exists."""
    d = os.path.split(f)[0]
    if not os.path.isdir(d):
        os.makedirs(d)

def touch(filename):
    """Opens a file and updates the mtime, like the posix command of the same name."""
    with io.open(filename, 'a') as f:
        os.utime(filename, None)
