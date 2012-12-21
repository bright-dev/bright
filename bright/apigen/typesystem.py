"""Implements a simple type system for API generation."""
from collections import Sequence

BASE_TYPES = set(['char', 'str', 'int32', 'int64', 'uint32', 'uint64', 'float32', 
                  'float64', 'complex128', 'void', 'bool'])

type_aliases = {
    'i': 'int32',
    'i4': 'int32',
    'i8': 'int64',
    'int': 'int32',
    'ui': 'uint32',
    'ui4': 'uint32',
    'ui8': 'uint64',
    'uint': 'uint32',
    'f': 'float64',
    'f4': 'float32',
    'f8': 'float64',
    'float': 'float64',
    'complex': 'complex128',
    'b': 'bool',
    'v': 'void',
    's': 'str',
    'string': 'str',
    # 'c' has char / complex ambiquity, not included
    }


# template types are types whose instantiations are based on meta-types
# this dict maps their names to meta-type names in order.
template_types = {
    'map': ('key_type', 'value_type'),
    'dict': ('key_type', 'value_type'),
    'set': ('value_type',),
    'vector': ('value_type',),
    }


# This is a mapping from refinement type names to the parent types.
# The parent types may either be base types, compound types, template 
# types, or other refined types!
refined_types = {
    'nucid': 'int32',
    'nucname': 'str',
    }

def _raise_type_error(t):
    raise TypeError("type of {0!r} could not be determined".format(t))

def _resolve_dependent_type(tname, tinst=None):
    depkey = [k for k in refined_types if k[0] == tname][0]
    depval = refined_types[depkey]
    istemplated = any([isinstance(x, basestring) for x in depkey[1:]])
    if tinst is None:
        return depkey
    elif istemplated:
        assert len(tinst) == len(depkey)
        typemap = {k: tinst[i] for i, k in enumerate(depkey[1:], 1) \
                                                    if isinstance(k, basestring)}
        for k in typemap:
            if k in type_aliases:
                raise TypeError('template type {0} already exists'.format(k))
        type_aliases.update(typemap)
        resotype = canon(depval), (tname,) + \
                        tuple([canon(k) for k in depkey[1:] if k in typemap]) + \
                        tuple([(k[0], canon(k[1]), instval) \
                            for k, instval in zip(depkey[1:], tinst[1:])
                            if k not in typemap])
        for k in typemap:
            del type_aliases[k]
        return resotype
    else:
        assert len(tinst) == len(depkey)
        return canon(depval), (tname,) + tuple([(kname, canon(ktype), instval) \
                        for (kname, ktype), instval in zip(depkey[1:], tinst[1:])])

def canon(t):
    """Turns the type into a canonical form."""
    if isinstance(t, basestring):
        if t in BASE_TYPES:
            return t
        elif t in type_aliases:
            return canon(type_aliases[t])
        elif t in refined_types:
            return (canon(refined_types[t]), t)
        elif t in set([k[0] for k in refined_types if not isinstance(k, basestring)]):
            return _resolve_dependent_type(t)
        else:
            _raise_type_error(t)
            # BELOW this would be for complicated string representations, 
            # such as 'char *' or 'map<nucid, double>'.  Would need to write
            # the parse_type() function and that might be a lot of work.
            #parse_type(t)  
    elif isinstance(t, Sequence):
        t0 = t[0]
        tlen = len(t)
        if 0 == tlen:
            _raise_type_error(t)
        last_val = 0 if tlen == 1 else t[-1]
        if isinstance(t0, basestring):
            if t0 in template_types:
                templen = len(template_types[t0])
                last_val = 0 if tlen == 1 + templen else t[-1]
                filledt = [t0] + [canon(tt) for tt in t[1:1+templen]] + [last_val]
                return tuple(filledt)
            elif t0 in set([k[0] for k in refined_types if not isinstance(k, basestring)]):
                return _resolve_dependent_type(t0, t)
            else:
                #if 2 < tlen:
                #    _raise_type_error(t)
                return (canon(t0), last_val)
        elif isinstance(t0, Sequence):
            t00 = t0[0]
            if isinstance(t00, basestring):
                # template or independent refinement type
                return (canon(t0), last_val)
            elif isinstance(t00, Sequence):
                # zOMG dependent type
                return _resolve_dependent_type(t00, t0)
            # BELOW is for possible compound types...
            #return (tuple([canon(subt) for subt in t[0]]), last_val)
        else:
            _raise_type_error(t)
    else:
        _raise_type_error(t)

#################### Type System Above This Line ##########################


_cython_c_base_types = {
    'char': 'char',
    'str': 'std_string',
    'int32': 'int',
    'uint32': 'extra_types.uint',  # 'unsigned int'
    'float32': 'float',
    'float64': 'double',
    'complex128': 'extra_types.complex_t',
    }


def cython_ctype(t):
    """Given a type t, returns the cooresponding Cython C type declaration."""
    t = canon(t)
    if isinstance(t, basestring):
        if  t in BASE_TYPES:
            return _cython_c_base_types[t]
    # must be tuple below this line
    if 2 == len(t):
        if 0 == t[1]:
            return _cython_c_base_types[t[0]]
