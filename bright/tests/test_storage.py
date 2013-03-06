"""Storage Component tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, with_setup

import os
import warnings
import tables as tb
import numpy as np

import bright
import bright.storage
import pyne.material

Storage = bright.storage.Storage
Material = pyne.material.Material
bright_conf = bright.bright_conf

#
# Fixtures
#

def teardown_storage():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "s.h5"]:
            os.remove(f)

#
# Tests that the storage component constructors work.
#

@with_setup(None, teardown_storage)
def test_Storage_1():
    s = Storage()
    assert_equal(s.name, '')
    assert_equal(s.track_params, set(["Mass"]))

@with_setup(None, teardown_storage)
def test_Storage_2():
    s = Storage(n="s")
    assert_equal(s.name, 's')
    assert_equal(s.track_params, set(["Mass"]))


#
# Tests that the fuel cycle component attributes work.
#

@with_setup(None, teardown_storage)
def test_decay_time():
    s = Storage()
    s.decay_time = 0.0
    assert_equal(s.decay_time, 0.0)
    s.decay_time = 628        
    assert_equal(s.decay_time, 628.0)

@with_setup(None, teardown_storage)
def test_track_params():
    s = Storage()
    assert_equal(s.track_params, set(["Mass"]))
    s.track_params = set(["Om nom nom"])
    assert_equal(s.track_params, set(["Om nom nom"]))


#                        
# Tests that the fuel cycle component methods work.
#

@with_setup(None, teardown_storage)
def test_calc_1():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    s = Storage()
    s.decay_time = 0.0
    s.mat_feed = Material({942390: 1.0})
    s.calc()
    assert_equal(s.mat_prod.mass, 1.0)
    assert_almost_equal(s.mat_prod.comp[942390], 1.0) 

@with_setup(None, teardown_storage)
def test_calc_2():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    s = Storage()
    s.decay_time = 0.0
    s.calc(Material({942390: 1.0}))
    assert_equal(s.mat_prod.mass, 1.0)
    assert_equal(s.mat_prod.comp[942390], 1.0) 

@with_setup(None, teardown_storage)
def test_calc_3():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    s = Storage()
    s.mat_feed = Material({942390: 1.0})
    s.calc(t=24110*365.25*24*3600)
    assert(s.mat_prod.mass < 1.0)
    assert_almost_equal(s.mat_prod.comp[942390], 0.5, 3) 

@with_setup(None, teardown_storage)
def test_calc_4():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    s = Storage()
    s.calc(Material({942390: 1.0}), 24110*365.25*24*3600)
    assert(s.mat_prod.mass < 1.0)
    assert_almost_equal(s.mat_prod.comp[942390], 0.5, 3) 

@with_setup(None, teardown_storage)
def test_calc_params():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    s = Storage()
    s.calc(Material({942390: 1.0}), 24110*365.25*24*3600)
    s.calc_params()
    assert_equal(s.params_prior_calc["Mass"],  1.00)
    assert(0.5 < s.params_after_calc["Mass"] < 1.0)
        

if __name__ == "__main__":
    nose.main()
