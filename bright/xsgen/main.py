import argparse
import plugins
import os

DEFAULT_RC_FILE="rc"
DEFAULT_PLUGINS="base"

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

def main():
    preparser = argparse.ArgumentParser()
    preparser.add_argument("infile")
    preparser.add_argument("--plugins", default=NotSpecified, nargs="+")
    preparser.add_argument("--rc", default=NotSpecified, nargs="+")
    prens = preparser.parse_known_args()[0]
    prens.abspath = os.path.abspath(prens.infile)
    predefaultrc = RunControl(rc=DEFAULT_RC_FILE, plugins=DEFAULT_PLUGINS)
    prerc = RunControl()
    prerc._update(predefaultrc)
    prerc.rc = prens.rc
    rcdict = {}
    if os.path.isfile(prerc.rc):
        exec_file(prerc.rc, rcdict, rcdict)
        prerc.rc = rcdict['rc'] if 'rc' in rcdict else NotSpecified
        prerc.plugins = rcdict['plugins'] if 'plugins' in rcdict else NotSpecified
    prerc._update([(k, v) for k, v in prens.__dict__.items()])

    plugins = Plugins(prerc.plugins)
    parser = plugins.build_cli()
    # if argcomplete is not None and prerc.bash_completion:
        # argcomplete.autocomplete(parser)
    ns = parser.parse_args()
    rc = plugins.merge_rcs()
    rc._update(rcdict)
    rc._update([(k, v) for k, v in ns.__dict__.items()])
    plugins.setup()
    plugins.execute()
    plugins.teardown()


if __name__ == "__main__":
    main()
