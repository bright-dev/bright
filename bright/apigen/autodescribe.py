"""This module creates Bright facility descriptions from C++ source code.
"""
import os
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


def clang_loc_in_range(location, source_range):
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

class ClangClassDescriber(object):

    def __init__(self, classname, onlyin=None, verbose=False):
        self.desc = {'name': classname, 'attrs': {}, 'methods': {}}
        self.classname = classname
        self.verbose = verbose
        onlyin = [onlyin] if isinstance(onlyin, basestring) else onlyin
        self.onlyin = set() if onlyin is None else set(onlyin)
        self._currfunc = None
        self._currfuncrange = None
        self._currclass = []  # this must be a stack to handle nested classes  

    def _pprint(self, node, typename):
        if self.verbose:
            print("{0}: {1}".format(typename, node.displayname) + "  CURRCLASS={0}".format(self._currclass))

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
            #print "  " + node.displayname
            if not clang_loc_in_range(node.location, self._currfuncrange):
                self._currfunc = None
                self._currfuncrange = None
            if 'class_decl' == kind and node.spelling == self._currclass[-1]:
                self._currclass.pop()

    def visit_class_decl(self, node):
        self._currclass.append(node.spelling)  # This could also be node.displayname
        self._pprint(node, "Class")

    def visit_function_decl(self, node):
        self._currfunc = node.spelling  # This could also be node.displayname
        self._currfuncrange = node.extent
        self._pprint(node, "Function")

    visit_cxx_method = visit_function_decl

    def visit_constructor(self, node):
        self._currfunc = node.spelling  # This could also be node.displayname
        self._currfuncrange = node.extent
        self._pprint(node, "Constructor")

    def visit_destructor(self, node):
        self._currfunc = node.spelling  # This could also be node.displayname
        self._currfuncrange = node.extent
        self._pprint(node, "Destructor")

    def visit_parm_decl(self, node):
        #print node.displayname, self._currfunc, self._currclass
        self._pprint(node, "Function Argument")

    def visit_field_decl(self, node):
        self._pprint(node, "Field")

    def visit_var_decl(self, node):
        self._pprint(node, "Variable")




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

