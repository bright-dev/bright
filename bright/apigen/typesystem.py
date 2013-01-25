"""Implements a simple type system for API generation."""
import functools
from collections import Sequence, Set

def _ishashable(x):
    try:
        hash(x)
        return True
    except TypeError:
        return False 

def _memoize(obj):
    # based off code from http://wiki.python.org/moin/PythonDecoratorLibrary
    cache = obj.cache = {}
    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = args + tuple(sorted(kwargs.items()))
        hashable = _ishashable(key)
        if hashable:
            if key not in cache:
                cache[key] = obj(*args, **kwargs)
            return cache[key]
        else:
            return obj(*args, **kwargs)
    return memoizer


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
    'pair': ('key_type', 'value_type'),
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

@_memoize
def isdependent(t):
    """Returns whether t is a dependent type or not."""
    deptypes = set([k[0] for k in refined_types if not isinstance(k, basestring)])
    if isinstance(t, basestring):
        return t in deptypes
    if isinstance(t, Sequence):
        return isdependent(t[0])
    return False


@_memoize
def isrefinement(t):
    """Returns whether t is a refined type."""
    if isinstance(t, basestring):
        return t in refined_types
    return isdependent(t)
        

def _raise_type_error(t):
    raise TypeError("type of {0!r} could not be determined".format(t))

@_memoize
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
            canon.cache.pop((k,), None)
        return resotype
    else:
        assert len(tinst) == len(depkey)
        return canon(depval), (tname,) + tuple([(kname, canon(ktype), instval) \
                        for (kname, ktype), instval in zip(depkey[1:], tinst[1:])])

@_memoize
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
    'void': 'void', 
    }

_cython_c_template_types = {
    'map': 'cpp_map',
    'dict': 'dict',
    'pair': 'cpp_pair',
    'set': 'cpp_set',
    'vector': 'cpp_vector',
    }

@_memoize
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
    'void': None, 
    }

_cython_cimport_template_types = {
    'map': ('libcpp.map', 'map', 'cpp_map'),
    'dict': None,
    'pair': ('libcpp.utility', 'pair', 'cpp_pair'),
    'set': ('libcpp.set', 'set', 'cpp_set'),
    'vector': ('libcpp.vector', 'vector', 'cpp_vector'),
    }

_cython_cyimport_base_types = {
    'char': None,
    'str': None,
    'int32': None,
    'uint32': None,
    'float32': None,
    'float64': None,
    'complex128': ('pyne', 'stlconverters', 'conv'),  # for py2c_complex()
    'void': None, 
    }

_cython_cyimport_template_types = {
    'map': ('pyne', 'stlconverters', 'conv'),
    'dict': None,
    'pair': ('pyne', 'stlconverters', 'conv'),
    'set': ('pyne', 'stlconverters', 'conv'),
    'vector': ('pyne', 'stlconverters', 'conv'),
    }

@_memoize
def cython_cimport_tuples(t, seen=None, inc=frozenset(['c', 'cy'])):
    """Given a type t, and possibily previously seen cimport tuples, return 
    the set of all seen cimport tuples."""
    t = canon(t)
    if seen is None:
        seen = set()
    if isinstance(t, basestring):
        if  t in BASE_TYPES:
            if 'c' in inc:
                seen.add(_cython_cimport_base_types[t])
            if 'cy' in inc:
                seen.add(_cython_cyimport_base_types[t])
            seen.discard(None)
            return seen
    # must be tuple below this line
    tlen = len(t)
    if 2 == tlen:
        return cython_cimport_tuples(t[0], seen, inc)
    elif 3 <= tlen:
        assert t[0] in template_types
        if 'c' in inc:
            seen.add(_cython_cimport_template_types[t[0]])
        if 'cy' in inc:
            seen.add(_cython_cyimport_template_types[t[0]])
        for x in t[1:-1]:
            cython_cimport_tuples(x, seen, inc)
        seen.discard(None)
        return seen

_cython_cimport_cases = {
    1: lambda tup: "cimport {0}".format(*tup),
    2: lambda tup: "from {0} cimport {1}".format(*tup),
    3: lambda tup: "from {0} cimport {1} as {2}".format(*tup),
    }

@_memoize
def cython_cimports(x, inc=frozenset(['c', 'cy'])):
    """Retuns the cimport lines associtated with a type or a set of seen tuples.
    """
    if not isinstance(x, Set):
        x = cython_cimport_tuples(x, inc=inc)
    return set([_cython_cimport_cases[len(tup)](tup) for tup in x])


_cython_cy_base_types = {
    'char': 'char',
    'str': 'char *',
    'int32': 'int',
    'uint32': 'long',  # 'unsigned int'
    'float32': 'float',
    'float64': 'float',
    'complex128': 'object',
    }


_cython_cy_template_types = {
    'map': 'conv._Map{key_type}{value_type}',
    'dict': 'dict',
    'pair': 'conv._Pair{value_type}',
    'set': 'conv._Set{value_type}',
    'vector': 'conv._Vector{value_type}',
    }


_cython_template_class_names = {
    # base types
    'char': 'Str',
    'str': 'Str',
    'int32': 'Int',
    'uint32': 'UInt',
    'float32': 'Float',
    'float64': 'Double',
    'complex128': 'Complex',
    # template types
    'map': 'Map{key_type}{value_type}',
    'dict': 'Dict',
    'pair': 'Pair{value_type}',
    'set': 'Set{value_type}',
    'vector': 'Vector{value_type}',    
    }


@_memoize
def _fill_cycyt(cycyt, t):
    """Helper for cython_cytype()."""
    d = {}
    for key, x in zip(template_types[t[0]], t[1:-1]):
        if isinstance(x, basestring):
            val = _cython_template_class_names[x]
        elif x[0] in BASE_TYPES:
            val = _cython_template_class_names[x[0]]
        else: 
            val, _ = _fill_cycyt(_cython_template_class_names[x[0]], x)
        d[key] = val
    return cycyt.format(**d), t
    

@_memoize
def cython_cytype(t):
    """Given a type t, returns the cooresponding Cython type."""
    t = canon(t)
    if isinstance(t, basestring):
        if  t in BASE_TYPES:
            return _cython_cy_base_types[t]
    # must be tuple below this line
    tlen = len(t)
    if 2 == tlen:
        if 0 == t[1]:
            return cython_cytype(t[0])
        elif isrefinement(t[1]):
            return cython_cytype(t[0])
        else:
            last = '[{0}]'.format(t[-1]) if isinstance(t[-1], int) else t[-1]
            return cython_cytype(t[0]) + ' {0}'.format(last)
    elif 3 <= tlen:
        assert t[0] in template_types
        assert len(t) == len(template_types[t[0]]) + 2
        template_name = _cython_cy_template_types[t[0]]
        assert template_name is not NotImplemented        
        cycyt = _cython_cy_template_types[t[0]]
        cycyt, t = _fill_cycyt(cycyt, t)
        if 0 != t[-1]:
            last = '[{0}]'.format(t[-1]) if isinstance(t[-1], int) else t[-1]
            cycyt += ' {0}'.format(last)
        return cycyt


_cython_py_base_types = {
    'char': 'str',
    'str': 'str',
    'int32': 'int',
    'uint32': 'int',  # 'unsigned int'
    'float32': 'float',
    'float64': 'float',
    'complex128': 'object',
    }


_cython_py_template_types = {
    'map': 'conv.Map{key_type}{value_type}',
    'dict': 'dict',
    'pair': 'conv.Pair{value_type}',
    'set': 'conv.Set{value_type}',
    'vector': 'conv.Vector{value_type}',
    }

@_memoize
def _fill_cypyt(cypyt, t):
    """Helper for cython_pytype()."""
    d = {}
    for key, x in zip(template_types[t[0]], t[1:-1]):
        if isinstance(x, basestring):
            val = _cython_template_class_names[x]
        elif x[0] in BASE_TYPES:
            val = _cython_template_class_names[x[0]]
        else: 
            val, _ = _fill_cypyt(_cython_template_class_names[x[0]], x)
        d[key] = val
    return cypyt.format(**d), t
    

@_memoize
def cython_pytype(t):
    """Given a type t, returns the cooresponding Python type."""
    t = canon(t)
    if isinstance(t, basestring):
        if  t in BASE_TYPES:
            return _cython_py_base_types[t]
    # must be tuple below this line
    tlen = len(t)
    if 2 == tlen:
        if 0 == t[1]:
            return cython_pytype(t[0])
        elif isrefinement(t[1]):
            return cython_pytype(t[0])
        else:
            # FIXME last is ignored for strings, but what about other types?
            #last = '[{0}]'.format(t[-1]) if isinstance(t[-1], int) else t[-1]
            #return cython_pytype(t[0]) + ' {0}'.format(last)
            return cython_pytype(t[0])
    elif 3 <= tlen:
        assert t[0] in template_types
        assert len(t) == len(template_types[t[0]]) + 2
        template_name = _cython_py_template_types[t[0]]
        assert template_name is not NotImplemented        
        cypyt = _cython_py_template_types[t[0]]
        cypyt, t = _fill_cypyt(cypyt, t)
        # FIXME last is ignored for strings, but what about other types?
        #if 0 != t[-1]:
        #    last = '[{0}]'.format(t[-1]) if isinstance(t[-1], int) else t[-1]
        #    cypyt += ' {0}'.format(last)
        return cypyt


_cython_c2py_conv = {
    # Has tuple form of (copy, [view, [cached_view]])
    # base types
    'char': ('str(<char *> {var})',),
    'str': ('str(<char *> {var}.c_str())',),
    'int32': ('int({var})',),
    'uint32': ('int({var})',),
    'float32': ('float({var})',),
    'float64': ('float({var})',),
    'complex128': ('complex(float({var}.re), float({var}.im))',),
    # template types
    'map': ('{pytype}({var})', 
           ('{proxy_name} = {pytype}(False, False)\n'
            '{proxy_name}.map_ptr = &{var}\n'),
           ('if {cache_name} is None:\n'
            '    {proxy_name} = {pytype}(False, False)\n'
            '    {proxy_name}.map_ptr = &{var}\n'
            '    {cache_name} = {proxy_name}\n'
            )),
    'dict': ('dict({var})',),
    'pair': ('{pytype}({var})',
            ('{proxy_name} = {pytype}(False, False)\n'
             '{proxy_name}.pair_ptr = &{var}\n'),
            ('if {cache_name} is None:\n'
             '    {proxy_name} = {pytype}(False, False)\n'
             '    {proxy_name}.pair_ptr = &{var}\n'
             '    {cache_name} = {proxy_name}\n'
             )),
    'set': ('{pytype}({var})',
           ('{proxy_name} = {pytype}(False, False)\n'
            '{proxy_name}.vector_ptr = &{var}\n'),
           ('if {cache_name} is None:\n'
            '    {proxy_name} = {pytype}(False, False)\n'
            '    {proxy_name}.set_ptr = &{var}\n'
            '    {cache_name} = {proxy_name}\n'
            )),
    'vector': ('{pytype}({var})',
              ('{proxy_name} = {pytype}(False, False)\n'
               '{proxy_name}.pair_ptr = &{var}\n'),
              ('if {cache_name} is None:\n'
               '    {proxy_name} = {pytype}(False, False)\n'
               '    {proxy_name}.vector_ptr = &{var}\n'
               '    {cache_name} = {proxy_name}\n'
               )),
    }

@_memoize
def cython_c2py(name, t, view=True, cached=True, inst_name=None, proxy_name=None, 
                cache_name=None):
    """Given a varibale name and type, returns cython code (declaration, body, 
    and return) to convert the variable from C/C++ to Python."""
    import pdb; pdb.set_trace()
    tkey = t
    while not isinstance(tkey, basestring):
        tkey = tkey[0]
    c2pyt = _cython_c2py_conv[tkey]
    ind = int(view) + int(cached)
    if cached and not view:
        raise ValueError('cached views require view=True.')
    if c2pyt is NotImplemented:
        raise NotImplementedError('conversion from C/C++ to Python for ' + \
                                  t + 'has not been implemented for when ' + \
                                  'view={0}, cached={1}'.format(view, cached))
    cyt = cython_cytype(t)
    pyt = cython_pytype(t)
    var = name if inst_name is None else "{0}.{1}".format(inst_name, name)
    cache_name = "self._{0}".format(name) if cache_name is None else cache_name
    proxy_name = "{0}_proxy".format(name) if proxy_name is None else proxy_name
    iscached = False
    if 1 == len(c2pyt) or ind == 0:
        decl = body = None
        rtn = c2pyt[0].format(var=var, pytype=pyt)
    elif ind == 1:
        decl = "cdef {0} {1}".format(cyt, proxy_name)
        body = c2pyt[1].format(var=var, pytype=pyt, proxy_name=proxy_name)
        rtn = proxy_name
    elif ind == 2:
        decl = "cdef {0} {1}".format(cyt, proxy_name)
        body = c2pyt[2].format(var=var, cache_name=cache_name, pytype=pyt,
                            proxy_name=proxy_name)
        rtn = cache_name
        iscached = True
    return decl, body, rtn, iscached


_cython_py2c_conv = {
    # Has tuple form of (body or return,  return or False)
    # base types
    'char': ('<char{last}> {var}', False),
    'str': ('std_string(<char *> {var})', False),
    'int32': ('{var}', False),
    'uint32': ('<extra_types.uint> long({var})', False),
    'float32': ('<float> {var}', False),
    'float64': ('<double> {var}', False),
    'complex128': ('conv.py2c_complex({var})', False),
    # template types
    'map': ('{proxy_name} = {pytype}({var}, not isinstance({var}, {cytype}))',
            '{proxy_name}.map_ptr[0]'),
    'dict': ('dict({var})', False),
    'pair': ('{proxy_name} = {pytype}({var}, not isinstance({var}, {cytype}))',
             '{proxy_name}.pair_ptr[0]'),
    'set': ('{proxy_name} = {pytype}({var}, not isinstance({var}, {cytype}))', 
            '{proxy_name}.set_ptr[0]'),
    'vector': ('{proxy_name} = {pytype}({var}, not isinstance({var}, {cytype}))', 
               '{proxy_name}.vector_ptr[0]'),
    # refinement types
    'nucid': ('nucname.zzaaam({var})', False),
    'nucname': ('nucname.name({var})', False),
    }

@_memoize
def cython_py2c(name, t, inst_name=None, proxy_name=None):
    """Given a varibale name and type, returns cython code (declaration, body, 
    and return) to convert the variable from Python to C/C++."""
    t = canon(t)
    if isinstance(t, basestring) or 0 == t[-1] or isrefinement(t[-1]):
        last = ''
    elif isinstance(t[-1], int):
        last = ' [{0}]'.format(t[-1])
    else:
        last = ' ' + t[-1]
    tkey = t
    tinst = None
    while not isinstance(tkey, basestring):
        tinst = tkey
        tkey = tkey[1] if (0 < len(tkey) and isrefinement(tkey[1])) else tkey[0]
    py2ct = _cython_py2c_conv[tkey]
    if py2ct is NotImplemented:
        raise NotImplementedError('conversion from C/C++ to Python for ' + \
                                  t + 'has not been implemented for ')
    body_template, rtn_template = py2ct
    cyt = cython_cytype(t)
    pyt = cython_pytype(t)
    var = name if inst_name is None else "{0}.{1}".format(inst_name, name)
    proxy_name = "{0}_proxy".format(name) if proxy_name is None else proxy_name
    template_kw = dict(var=var, proxy_name=proxy_name, pytype=pyt, cytype=cyt, 
                       last=last)
    nested = False
    if isdependent(tkey):
        tsig = [ts for ts in refined_types if ts[0] == tkey][0]
        for ts, ti in zip(tsig[1:], tinst[1:]):
            if isinstance(ts, basestring):
                template_kw[ts] = cython_ctype(ti)
            else:
                template_kw[ti[0]] = ti[2]
        vartype = refined_types[tsig]
        if vartype in tsig[1:]:
            vartype = tinst[tsig.index(vartype)][1]
        if isrefinement(vartype):
            nested = True
            vdecl, vbody, vrtn = cython_py2c(var, vartype)
            template_kw['var'] = vrtn
    body_filled = body_template.format(**template_kw)
    if rtn_template:
        decl = "cdef {0} {1}".format(cyt, proxy_name)
        body = body_filled
        rtn = rtn_template.format(**template_kw)
    else:
        decl = body = None
        rtn = body_filled
    if nested:
        decl = '' if decl is None else decl
        vdecl = '' if vdecl is None else vdecl
        decl = (vdecl + '\n' + decl).strip()
        decl = None if 0 == len(decl) else decl
        body = '' if body is None else body
        vbody = '' if vbody is None else vbody
        body = (vbody + '\n' + body).strip()
        body = None if 0 == len(body) else body
    return decl, body, rtn
 


######################  Some utility functions for the typesystem #############



def register_class(name, template_args=None, cython_c_type=None, 
                   cython_cimport=None, cython_cy_type=None, cython_py_type=None,
                   cython_template_class_name=None, cython_cyimport=None, 
                   cython_c2py=None, cython_py2c=None):
    """Classes are user specified types.  This function will add a class to 
    the type system so that it may be used normally with the rest of the 
    type system.

    """
    # register the class name
    isbase = True
    if template_args is None: 
        BASE_TYPES.add(name)  # normal class        
    elif isinstance(template_args, Sequence):
        if 0 == len(template_args):
            BASE_TYPES.add(name)  # normal class
        elif isinstance(template_args, basestring):
            _raise_type_error(name)
        else:
            template_types[name] = tuple(template_args)  # templated class...
            isbase = False

    # Register with Cython C/C++ types
    if (cython_c_type is not None) or (cython_cy_type is not None):
        if isinstance(cython_cimport, basestring):
            cython_cimport = (cython_cimport,)

        if isinstance(cython_cyimport, basestring):
            cython_cyimport = (cython_cyimport,)

        if isinstance(cython_c2py, basestring):
            cython_c2py = (cython_c2py,)
        cython_c2py = None if cython_c2py is None else tuple(cython_c2py)

        if isinstance(cython_py2c, basestring):
            cython_py2c = (cython_py2c, False)

        if isbase:
            _cython_c_base_types[name] = cython_c_type
            _cython_cy_base_types[name] = cython_cy_type
            _cython_py_base_types[name] = cython_py_type
            _cython_cimport_base_types[name] = cython_cimport
            _cython_cyimport_base_types[name] = cython_cyimport
        else:
            _cython_c_template_types[name] = cython_c_type
            _cython_cy_template_types[name] = cython_cy_type
            _cython_py_template_types[name] = cython_py_type
            _cython_cimport_template_types[name] = cython_cimport
            _cython_cyimport_template_types[name] = cython_cyimport

        _cython_c2py_conv[name] = cython_c2py
        _cython_py2c_conv[name] = cython_py2c
        _cython_template_class_names[name] = cython_template_class_name


def deregister_class(name):
    """This function will remove previously registered classes from the type
    system.
    """
    isbase = name in BASE_TYPES
    if not isbase and name not in template_types:
        _raise_type_error(name)

    if isbase:
        BASE_TYPES.remove(name)
        _cython_c_base_types.pop(name, None)
        _cython_cy_base_types.pop(name, None)
        _cython_py_base_types.pop(name, None)
        _cython_cimport_base_types.pop(name, None)
        _cython_cyimport_base_types.pop(name, None)
    else:
        template_types.pop(name, None)
        _cython_c_template_types.pop(name, None)
        _cython_cy_template_types.pop(name, None)
        _cython_py_template_types.pop(name, None)
        _cython_cimport_template_types.pop(name, None)
        _cython_cyimport_template_types.pop(name, None)

    _cython_c2py_conv.pop(name, None)
    _cython_py2c_conv.pop(name, None)
    _cython_template_class_names.pop(name, None)

    # clear all caches
    funcs = [isdependent, isrefinement, _resolve_dependent_type, canon, 
             cython_ctype, cython_cimport_tuples, cython_cimports, _fill_cycyt,
             cython_cytype, _fill_cypyt, cython_pytype, cython_c2py, cython_py2c]
    for f in funcs:
        f.cache.clear()
