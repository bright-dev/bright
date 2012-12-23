from bright.apigen import typesystem as ts
from bright.apigen import cygen as cg

from nose.tools import assert_equal

toaster_desc = {
    'name': 'Toaster',
    'header_filename': 'toaster.h',
    'cpppxd_filename': 'cpp_toaster.pxd',
    'namespace': 'bright',
    'library_docstring': "I am the Toaster lib! Hear me sizzle!", 
    'parents': ['FCComp'],
    'attrs': {
        'nslices': 'uint',
        'toastiness': 'str',
        'rate': 'float',
        },
    'methods': {
        ('Toaster',): None,
        ('~Toaster',): None, 
        ('make_toast', ('when', 'str'), ('nslices', 'uint', 1)): 'int',
        },
    }


exp_cpppxd = cg.AUTOGEN_WARNING + \
"""from libcpp.string cimport string as std_string
from pyne cimport extra_types

cdef extern from "toaster.h" namespace "bright":

    cdef cppclass Toaster(FCComp):
        # constructors
        Toaster() except +
        ~Toaster() except +

        # attributes
        extra_types.uint nslices
        double rate
        std_string toastiness

        # methods
        int make_toast(std_string) except +
        int make_toast(std_string, extra_types.uint) except +
"""

def test_gencpppxd():
    obs = cg.gencpppxd(toaster_desc).splitlines()
    exp = exp_cpppxd.splitlines()
    assert_equal(len(obs), len(exp))
    for o, e in zip(obs, exp):
        assert_equal(o, e)


exp_pxd = cg.AUTOGEN_WARNING + \
"""cimport cpp_toaster
cimport fccomp

cdef class Toaster(fccomp.FCComp):
    cdef cpp_toaster.Toaster * _inst
    cdef public bint _free_inst
"""

def test_genpxd():
    ts.register_class('FCComp', cython_c_type='fccomp.FCComp', 
                      cython_cimport='fccomp')
    obs = cg.genpxd(toaster_desc).splitlines()
    ts.deregister_class('FCComp')
    exp = exp_pxd.splitlines()
    assert_equal(len(obs), len(exp))
    for o, e in zip(obs, exp):
        assert_equal(o, e)
