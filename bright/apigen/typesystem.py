"""Implements a simple type system for API generation."""
from collections import Sequence, Set

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


def isdependent(t):
    """Returns whether t is a dependent type or not."""
    deptypes = set([k[0] for k in refined_types if not isinstance(k, basestring)])
    if isinstance(t, basestring):
        return t in deptypes
    if isinstance(t, Sequence):
        return isdependent(t[0])
    return False


def isrefinement(t):
    """Returns whether t is a refined type."""
    if isinstance(t, basestring):
        return t in refined_types
    return isdependent(t)
        

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
        elif isdependent(t):
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
            elif isdependent(t0):
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

_cython_c_template_types = {
    'map': 'cpp_map',
    'dict': 'dict',
    'set': 'cpp_set',
    'vector': 'cpp_vector',
    }

def cython_ctype(t):
    """Given a type t, returns the cooresponding Cython C type declaration."""
    t = canon(t)
    if isinstance(t, basestring):
        if  t in BASE_TYPES:
            return _cython_c_base_types[t]
    # must be tuple below this line
    tlen = len(t)
    if 2 == tlen:
        if 0 == t[1]:
            return cython_ctype(t[0])
        elif isrefinement(t[1]):
            return cython_ctype(t[0])
        else:
            last = '[{0}]'.format(t[-1]) if isinstance(t[-1], int) else t[-1]
            return cython_ctype(t[0]) + ' {0}'.format(last)
    elif 3 <= tlen:
        assert t[0] in template_types
        assert len(t) == len(template_types[t[0]]) + 2
        template_name = _cython_c_template_types[t[0]]
        assert template_name is not NotImplemented
        template_filling = ', '.join([cython_ctype(x) for x in t[1:-1]])
        cyct = '{0}[{1}]'.format(template_name, template_filling)
        if 0 != t[-1]:
            last = '[{0}]'.format(t[-1]) if isinstance(t[-1], int) else t[-1]
            cyct += ' {0}'.format(last)
        return cyct


_cython_cimport_base_types = {
    'char': None,
    'str': ('libcpp.string', 'string', 'std_string'),
    'int32': None,
    'uint32': ('pyne', 'extra_types'),  # 'unsigned int'
    'float32': None,
    'float64': None,
    'complex128': ('pyne', 'extra_types'),
    }

_cython_cimport_template_types = {
    'map': ('libcpp.map', 'map', 'cpp_map'),
    'dict': None,
    'set': ('libcpp.set', 'set', 'cpp_set'),
    'vector': ('libcpp.vector', 'vector', 'cpp_vector'),
    }

def cython_cimport_tuples(t, seen=None):
    """Given a type t, and possibily previously seen cimport tuples, return 
    the set of all seen cimport tuples."""
    t = canon(t)
    if seen is None:
        seen = set()
    if isinstance(t, basestring):
        if  t in BASE_TYPES:
            seen.add(_cython_cimport_base_types[t])
            seen.discard(None)
            return seen
    # must be tuple below this line
    tlen = len(t)
    if 2 == tlen:
        return cython_cimport_tuples(t[0], seen)
    elif 3 <= tlen:
        assert t[0] in template_types
        seen.add(_cython_cimport_template_types[t[0]])
        for x in t[1:-1]:
            cython_cimport_tuples(x, seen)
        seen.discard(None)
        return seen

_cython_cimport_cases = {
    1: lambda tup: "cimport {0}".format(*tup),
    2: lambda tup: "from {0} cimport {1}".format(*tup),
    3: lambda tup: "from {0} cimport {1} as {2}".format(*tup),
    }

def cython_cimports(x):
    """Retuns the cimport lines associtated with a type or a set of seen tuples."""
    if not isinstance(x, Set):
        x = cython_cimport_tuples(x)
    return set([_cython_cimport_cases[len(tup)](tup) for tup in x])
