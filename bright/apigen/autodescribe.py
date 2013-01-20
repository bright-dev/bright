"""This module creates Bright facility descriptions from C++ source code.
"""
import faulthandler
faulthandler.enable()
import os
import linecache
import subprocess
from pprint import pprint, pformat

# CLang conditional imports
try:
    from clang import cindex
except ImportError:
    try:
        from .clang.v3_1 import cindex
    except ImportError:
        pass

# GCC-XML conditional imports
try:
    from lxml import etree
except ImportError:
    try:
        # Python 2.5
        import xml.etree.cElementTree as etree
    except ImportError:
        try:
          # Python 2.5
          import xml.etree.ElementTree as etree
        except ImportError:
            try:
                # normal cElementTree install
                import cElementTree as etree
            except ImportError:
                try:
                  # normal ElementTree install
                  import elementtree.ElementTree as etree
                except ImportError:
                    pass
import tempfile

# Other imports
import pyne

def describe(filename, classname=None, parser='gccxml', verbose=False):
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
        'clang' and 'gccxml' are supported, though others may be 
        implemented in the future.

    Returns
    -------
    desc : dict
        A dictionary describing the class which may be used to generate
        API bindings.
    """
    if classname is None:
        classname = os.path.split(filename)[-1].rsplit('.', 1)[0].capitalize()
    describers = {'clang': clang_describe, 'gccxml': gccxml_describe}
    describer = describers[parser]
    desc = describer(filename, classname, verbose=verbose)
    return desc


#
# GCC-XML Describers
#


def gccxml_describe(filename, classname, verbose=False):
    """Use GCC-XML to describe the class."""
    f = tempfile.NamedTemporaryFile()
    cmd = ['gccxml', filename, '-fxml=' + f.name, '-I' + pyne.includes]
    subprocess.call(cmd)
    f.seek(0)
    root = etree.parse(f)
    onlyin = set([filename, filename.replace('.cpp', '.h')])
    #onlyin = set([filename.replace('.cpp', '.h')])
    describer = GccxmlClassDescriber(classname, root, onlyin=onlyin, verbose=verbose)
    describer.visit(root)
    f.close()
    print describer.desc
    return describer.desc


class GccxmlClassDescriber(object):

    def __init__(self, classname, root=None, onlyin=None, verbose=False):
        self.desc = {'name': classname, 'attrs': {}, 'methods': {}}
        self.classname = classname
        self.verbose = verbose
        self._root = root
        onlyin = [onlyin] if isinstance(onlyin, basestring) else onlyin
        onlyin = set() if onlyin is None else set(onlyin)
        self.onlyin = set([root.find("File[@name='{0}']".format(oi)).attrib['id'] \
                           for oi in onlyin])
        self._currfunc = []  # this must be a stack to handle nested functions
        self._currfuncsig = None
        self._currfuncarg = None
        self._currclass = []  # this must be a stack to handle nested classes  
        self._indent = -1

    def __str__(self):
        return pformat(self.desc)

    def __del__(self):
        linecache.clearcache()

    def _pprint(self, node):
        if self.verbose:
            print("{0}{1}: {2}".format(self._indent * "  ", node.tag, 
                                       node.attrib.get('name', None)))

    def visit(self, node):
        self._indent += 1
        for child in node.iterfind('*'):
            if child.attrib.get('file', None) in self.onlyin:
                tag = child.tag.lower()
                meth_name = 'visit_' + tag
                meth = getattr(self, meth_name, None)
                if meth is not None:
                    meth(child)
        self._indent -= 1

    def visit_class(self, node):
        self._pprint(node)
        self.visit(node)  # Walk farther down the tree

    def visit_constructor(self, node):
        self._pprint(node)
        self.visit(node)  # Walk farther down the tree

    def visit_destructor(self, node):
        self._pprint(node)
        self.visit(node)  # Walk farther down the tree

    def visit_method(self, node):
        self._pprint(node)
        self.visit(node)  # Walk farther down the tree

    def visit_argument(self, node):
        self._pprint(node)
        self.visit(node)  # Walk farther down the tree

    def visit_field(self, node):
        self._pprint(node)
        self.visit(node)  # Walk farther down the tree


#
# Clang Describers
#

def clang_describe(filename, classname, verbose=False):
    """Use clang to describe the class."""
    index = cindex.Index.create()
    tu = index.parse(filename, args=['-cc1', '-I' + pyne.includes])
    #onlyin = set([filename, filename.replace('.cpp', '.h')])
    onlyin = set([filename.replace('.cpp', '.h')])
    describer = ClangClassDescriber(classname, onlyin=onlyin, verbose=verbose)
    describer.visit(tu.cursor)
    from pprint import pprint; pprint(describer.desc)
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

    def __init__(self, classname, root=None, onlyin=None, verbose=False):
        self.desc = {'name': classname, 'attrs': {}, 'methods': {}}
        self.classname = classname
        self.verbose = verbose
        onlyin = [onlyin] if isinstance(onlyin, basestring) else onlyin
        self.onlyin = set() if onlyin is None else set(onlyin)
        self._currfunc = []  # this must be a stack to handle nested functions
        self._currfuncsig = None
        self._currfuncarg = None
        self._currclass = []  # this must be a stack to handle nested classes  

    def __str__(self):
        return pformat(self.desc)

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
        rtnname = ClangTypeVisitor(verbose=self.verbose).visit(rtntype)
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
        t = ClangTypeVisitor(verbose=self.verbose).visit(node)
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

    ##########

    def visit_type_ref(self, cur):
        self._pprint(cur, "type ref")

    def visit_template_ref(self, cur):
        self._pprint(cur, "template ref")

    def visit_template_type_parameter(self, cur):
        self._pprint(cur, "template type param")

    def visit_template_non_type_parameter(self, cur):
        self._pprint(cur, "template non-type param")

    def visit_template_template_parameter(self, cur):
        self._pprint(cur, "template template param")

    def visit_class_template(self, cur):
        self._pprint(cur, "class template")

    def visit_class_template_partial_specialization(self, cur):
        self._pprint(cur, "class template partial specialization")


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


class ClangTypeVisitor(object):
    """For a Clang type located at a root node, compute the cooresponding 
    typesystem type.
    """

    def __init__(self, verbose=False):
        self.type = []
        self.verbose = verbose
        self.namespace = []  # this must be a stack to handle nested namespaces
        self._atrootlevel = True
        self._currtype = []

    def _pprint(self, node, typename):
        if self.verbose:
            msg = "{0}: {1}"
            if isinstance(node, cindex.Type):
                msg = msg.format(typename, node.kind.spelling)
            elif isinstance(node, cindex.Cursor):
                msg = msg.format(typename, node.displayname)
            else:
                msg = msg.format(typename, node)
            print(msg)

    def visit(self, root):
        """Takes a root type."""
        atrootlevel = self._atrootlevel

        if isinstance(root, cindex.Type):
            typekind = root.kind.name.lower()
            methname = 'visit_' + typekind
            meth = getattr(self, methname, None)
            if meth is not None and root.kind != cindex.TypeKind.INVALID:
                meth(root)
        elif isinstance(root, cindex.Cursor):
            self.visit(root.type)
            for child in root.get_children():
                kindname = child.kind.name.lower()
                methname = 'visit_' + kindname
                meth = getattr(self, methname, None)
                if meth is not None:
                    meth(child)
                if hasattr(child, 'get_children'):
                    self._atrootlevel = False
                    self.visit(child)
                    self._atrootlevel = atrootlevel
                else:
                    self.visit(child.type)

        if self._atrootlevel:
            currtype = self._currtype
            currtype = currtype[0] if 1 == len(currtype) else tuple(currtype)
            self.type = [self.type, currtype] if isinstance(self.type, basestring) \
                        else list(self.type) + [currtype]
            self._currtype = []
            self.type = self.type[0] if 1 == len(self.type) else tuple(self.type)
            return self.type

    def visit_void(self, typ):
        self._pprint(typ, "void")
        self._currtype.append("void")

    def visit_bool(self, typ):
        self._pprint(typ, "boolean")
        self._currtype.append("bool")

    def visit_char_u(self, typ):
        self._pprint(typ, "character")
        self._currtype.append("char")

    visit_uchar = visit_char_u

    def visit_uint(self, typ):
        self._pprint(typ, "unsigned integer, 32-bit")
        self._currtype.append("uint32")

    def visit_ulong(self, typ):
        self._pprint(typ, "unsigned integer, 64-bit")
        self._currtype.append("uint64")

    def visit_int(self, typ):
        self._pprint(typ, "integer, 32-bit")
        self._currtype.append("int32")

    def visit_long(self, typ):
        self._pprint(typ, "integer, 64-bit")
        self._currtype.append("int64")

    def visit_float(self, typ):
        self._pprint(typ, "float, 32-bit")
        self._currtype.append("float32")

    def visit_double(self, typ):
        self._pprint(typ, "float, 64-bit")
        self._currtype.append("float64")

    def visit_complex(self, typ):
        self._pprint(typ, "complex, 128-bit")
        self._currtype.append("complex128")

    def visit_unexposed(self, typ):
        self._pprint(typ, "unexposed")
        #typ = typ.get_canonical()
        decl = typ.get_declaration()
        self._currtype.append(decl.spelling)
        print "   canon: ",  typ.get_canonical().get_declaration().displayname
        #import pdb; pdb.set_trace()        
        #self.visit(decl)
        #self.visit(typ.get_canonical().get_declaration())
        #self.visit(typ.get_canonical())

    def visit_typedef(self, typ):
        self._pprint(typ, "typedef")
        decl = typ.get_declaration()
        t = decl.underlying_typedef_type
        #self.visit(t.get_canonical())

    def visit_record(self, typ):
        self._pprint(typ, "record")
        self.visit(typ.get_declaration())

    def visit_invalid(self, typ):
        self._pprint(typ, "invalid")
        self.visit(typ.get_declaration())

    def visit_namespace_ref(self, cur):
        self._pprint(cur, "namespace")
        if self._atrootlevel:
            self.namespace.append(cur.displayname)

    def visit_type_ref(self, cur):
        self._pprint(cur, "type ref")
        self._currtype.append(cur.displayname)
#        print "    cur type kin =", cur.type.kind
        #self.visit(cur.type)
        #self.visit(cur)

    def visit_template_ref(self, cur):
        self._pprint(cur, "template ref")
        self._currtype.append(cur.displayname)
        #self.visit(cur)

        #import pdb; pdb.set_trace()
#        self.visit(cur)
        print "   canon: ",  cur.type.get_canonical().get_declaration().displayname

    def visit_template_type_parameter(self, cur):
        self._pprint(cur, "template type param")

    def visit_template_non_type_parameter(self, cur):
        self._pprint(cur, "template non-type param")

    def visit_template_template_parameter(self, cur):
        self._pprint(cur, "template template param")

    def visit_function_template(self, cur):
        self._pprint(cur, "function template")

#    def visit_class_template(self, cur):
#        self._pprint(cur, "class template")

    def visit_class_template_partial_specialization(self, cur):
        self._pprint(cur, "class template partial specialization")

#    def visit_var_decl(self, cur):
#        self._pprint(cur, "variable")



def clang_canonize(t):
    kind = t.kind
    if kind in clang_base_typekinds:
        name = clang_base_typekinds[kind]
    elif kind == cindex.TypeKind.UNEXPOSED:
        name = t.get_declaration().spelling
    elif kind == cindex.TypeKind.TYPEDEF:
        print [n.displayname for n in t.get_declaration().get_children()]
        print [n.kind.name for n in t.get_declaration().get_children()]
        name = "<fixme>"
    else:
        name = "<error:{0}>".format(kind)
    return name
