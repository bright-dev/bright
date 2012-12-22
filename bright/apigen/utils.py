"""Helper functions for bright API generation."""

def indent(s, n):
    """Indents all lines in the string or list s by n spaces."""
    spaces = " " * n
    lines = s.splitlines() if isinstance(s, basestring) else s
    return '\n'.join([spaces + l for l in lines])


def expand_methods_default_args(methods):
    """This function takes a collection of method tuples and expands all of the default
    arguments, returning a set of all methods possible."""
    methitems = set()
    for mkey, mrtn in methods:
        mname, margs = mkey[0], mkey[1:]
        havedefaults = [3 == len(arg) for arg in margs]
        if any(havedefaults):
            # expand default arguments  
            n = havedefaults.index(True)
            mitems = [((mname, margs[:n]), mrtn)] + \
                     [((mname, margs[:i]), mrtn) for i in range(n+1, len(margs))]
            methitems.update(mitems)
        else:
            # no default args
            methitems.add((mkey, mrtn))
    return methitems
