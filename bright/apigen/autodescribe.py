"""This module creates Bright facility descriptions from C++ source code.
"""
import os
import linecache
import subprocess

try:
    from clang import cindex
except ImportError:
    from .clang.v3_1 import cindex

import pyne

def describe(filename, classname=None, parser='clang'):
    """Automatically describes a class in a file.

    Parameters
    ----------
    filename : str
        The path to the file.
    classname : str or None, optional
        The classname, a 'None' value will attempt to infer this from the 
        filename.
    parser : str
        The parser / AST to use to use for the C++ file.  Currently only
        'clang' is supported, though others (such as gccxml) may be 
        implemented in the future.

    Returns
    -------
    desc : dict
        A dictionary describing the class which may be used to generate
        API bindings.
    """
    if classname is None:
        classname = os.path.split(filename)[-1].rsplit('.', 1)[0].capitalize()
    describers = {'clang': clang_describe}
    describer = describers[parser]
    desc = describer(filename, classname)
    return desc


def clang_describe(filename, classname):
    """Use clang to describe the class."""
    index = cindex.Index.create()
    tu = index.parse(filename, args=['-cc1', '-I' + pyne.includes])
    #onlyin = set([filename, filename.replace('.cpp', '.h')])
    onlyin = set([filename.replace('.cpp', '.h')])
    describer = ClangClassDescriber(classname, onlyin=onlyin, verbose=True)
    describer.visit(tu.cursor)
    print describer.desc
    return describer.desc


def clang_is_loc_in_range(location, source_range):
    """Returns whether a given Clang location is part of a source file range."""
    if source_range is None or location is None:
        return False
    start = source_range.start
    stop = source_range.end
    file = location.file
    if file != start.file or file != stop.file:
        return False
    line = location.line
    if line < start.line or stop.line < line:
        return False
    return start.column <= location.column <= stop.column


def clang_range_str(source_range):
    """Get the text present on a source range."""
    start = source_range.start
    stop = source_range.end
    filename = start.file.name
    if filename != stop.file.name:
        msg = 'range spans multiple files: {0!r} & {1!r}'
        msg = msg.format(filename, stop.file.name)
        raise ValueError(msg)
    lines = [linecache.getline(filename, n) for n in range(start.line, stop.line+1)]
    lines[-1] = lines[-1][:stop.column-1]  # stop slice must come first for 
    lines[0] = lines[0][start.column-1:]   # len(lines) == 1
    s = "".join(lines)
    return s
    


class ClangClassDescriber(object):

    _funckinds = set(['function_decl', 'cxx_method', 'constructor', 'destructor'])

    def __init__(self, classname, onlyin=None, verbose=False):
        self.desc = {'name': classname, 'attrs': {}, 'methods': {}}
        self.classname = classname
        self.verbose = verbose
        onlyin = [onlyin] if isinstance(onlyin, basestring) else onlyin
        self.onlyin = set() if onlyin is None else set(onlyin)
        self._currfunc = []  # this must be a stack to handle nested functions
        self._currfuncsig = None
        self._currfuncarg = None
        self._currclass = []  # this must be a stack to handle nested classes  

    def __del__(self):
        linecache.clearcache()

    def _pprint(self, node, typename):
        if self.verbose:
            print("{0}: {1}".format(typename, node.displayname))

    def visit(self, root):
        for node in root.get_children():
            if not node.location.file or node.location.file.name not in self.onlyin:
                continue  # Ignore AST elements not from the desired source files
            kind = node.kind.name.lower()
            meth_name = 'visit_' + kind
            meth = getattr(self, meth_name, None)
            if meth is not None:
                meth(node)
            if hasattr(node, 'get_children'):
                self.visit(node)

            # reset the current function and class
            if kind in self._funckinds and node.spelling == self._currfunc[-1]:
                _key, _value = self._currfuncsig
                _key = (_key[0],) + tuple([tuple(k) for k in _key[1:]])
                self.desc['methods'][_key] = _value
                self._currfunc.pop()
                self._currfuncsig = None
            elif 'class_decl' == kind and node.spelling == self._currclass[-1]:
                self._currclass.pop()
            elif 'unexposed_expr' == kind and node.spelling == self._currfuncarg:
                self._currfuncarg = None

    def visit_class_decl(self, node):
        self._pprint(node, "Class")
        self._currclass.append(node.spelling)  # This could also be node.displayname

    def visit_function_decl(self, node):
        self._pprint(node, "Function")
        self._currfunc.append(node.spelling)  # This could also be node.displayname
        rtntype = node.type.get_result()
        rtnname = clang_canonize(rtntype)
        self._currfuncsig = ([node.spelling], rtnname)

    visit_cxx_method = visit_function_decl

    def visit_constructor(self, node):
        self._pprint(node, "Constructor")
        self._currfunc.append(node.spelling)  # This could also be node.displayname
        self._currfuncsig = ([node.spelling], None)

    def visit_destructor(self, node):
        self._pprint(node, "Destructor")
        self._currfunc.append(node.spelling)  # This could also be node.displayname
        self._currfuncsig = ([node.spelling], None)

    def visit_parm_decl(self, node):
        self._pprint(node, "Function Argument")
        name = node.spelling
        t = clang_canonize(node.type)
        self._currfuncsig[0].append([name, t])
        self._currfuncarg = name

    def visit_field_decl(self, node):
        self._pprint(node, "Field")

    def visit_var_decl(self, node):
        self._pprint(node, "Variable")

    def visit_unexposed_expr(self, node):
        self._pprint(node, "Default Parameter (Unexposed Expression)")
        # a little hacky reading from the file, 
        # but Clang doesn't expose this data...
        if self._currfuncsig is None:
            return
        currarg = self._currfuncsig[0][-1]
        assert currarg[0] == self._currfuncarg
        r = node.extent
        default_val = clang_range_str(r)
        if 2 == len(currarg):
            currarg.append(default_val)
        elif 3 == len(currarg):
            currarg[2] = default_val




def clang_find_class(node, classname, namespace=None):
    """Find the node for a given class underneath the current node.
    """
    if namespace is None:
        nsdecls = [node]
    else:
        nsdecls = [n for n in clang_find_declarations(node) if n.spelling == namespace]
    classnode = None
    for nsnode in nsdecls[::-1]:
        decls = [n for n in clang_find_declarations(nsnode) if n.spelling == classname]
        if 0 < len(decls):
            assert 1 == len(decls)
            classnode = decls[0]
            break
    if classnode is None:
        msg = "the class {0} could not be found in {1}".format(classname, filename)
        raise ValueError(msg)
    return classnode


def clang_find_declarations(node):
    """Finds declarations one level below the Clang node."""
    return [n for n in node.get_children() if n.kind.is_declaration()]

def clang_find_attributes(node):
    """Finds attributes one level below the Clang node."""
    return [n for n in node.get_children() if n.kind.is_attribute()]


# maps Clang TypeKinds to typesystem types
clang_base_typekinds = {
    cindex.TypeKind.VOID: 'void',
    cindex.TypeKind.BOOL: 'bool',
    cindex.TypeKind.CHAR_U: 'char',
    cindex.TypeKind.UCHAR: 'char',
    cindex.TypeKind.UINT: 'uint32',
    cindex.TypeKind.ULONG: 'uint64',
    cindex.TypeKind.INT: 'int32',
    cindex.TypeKind.LONG: 'int64',
    cindex.TypeKind.FLOAT: 'float32',
    cindex.TypeKind.DOUBLE: 'float64',
    cindex.TypeKind.COMPLEX: 'complex128',
    }

def clang_canonize(t):
    """For a Clang type t, return the cooresponding typesystem name.
    """
    kind = t.kind
    if kind in clang_base_typekinds:
        name = clang_base_typekinds[kind]
    elif kind == cindex.TypeKind.UNEXPOSED:
        name = t.get_declaration().spelling
    elif kind == cindex.TypeKind.TYPEDEF:
        name = "<fixme>"
    else:
        name = "<error:{0}>".format(kind)
    return name
