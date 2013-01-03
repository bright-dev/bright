"""This module creates Bright facility descriptions from C++ source code.
"""
import os
import subprocess

try:
    from clang import cindex
except ImportError:
    from .clang.v3_1 import cindex

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
    tu = index.parse(filename, args=['-cc1'])
    classnode = clang_find_class(tu.cursor, classname, namespace='bright')
    for node in clang_find_declarations(classnode):
        print node.spelling
    return tu


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

