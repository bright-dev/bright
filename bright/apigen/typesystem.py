"""Implements a simple type system for API generation."""


def canon(t):


_cython_c_base_types = {
    'str': 'std_string',
    'int': 'int',
    'uint': 'extra_types.uint',  # 'unsigned int'
    'float': 'float',
    'double': 'double',
    'complex': 'extra_types.complex_t',
    }


def cython_ctype(t):
    """Given a type t, returns the cooresponding Cython C type declaration."""
    t = canon(t)
