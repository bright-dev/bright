"""Reprocessing Component tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, with_setup

import os
import warnings
import tables as tb
import numpy as np

import bright
import bright.reprocess
import pyne.material

bright_conf = bright.bright_conf
Reprocess = bright.reprocess.Reprocess
Material = pyne.material.Material


#
# Fixtures
#

def teardown_rep():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "r.h5"] or (".h5" in f):
            os.remove(f)

#
# Tests that the reprocessing component constructors work.
#

@with_setup(None, teardown_rep)
def test_Reprocess_1():
    r = Reprocess()
    assert_equal(r.name, '')
    for value in r.sepeff.values():
        assert_equal(value, 1.0)
    assert_equal(r.track_params, set(["Mass"]))

@with_setup(None, teardown_rep)
def test_Reprocess_2():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "94239": 0.99})
    assert_equal(r.name, '')
    assert_almost_equal(r.sepeff[922350], 0.9)
    assert_almost_equal(r.sepeff[922380], 0.999)
    assert_almost_equal(r.sepeff[942390], 0.99)
    assert_equal(r.track_params, set(["Mass"]))

@with_setup(None, teardown_rep)
def test_Reprocess_3():
    bright_conf.track_nucs = set([922350])
    r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "94239": 0.99})
    assert_equal(r.name, '')
    assert_equal(r.sepeff, {922350: 0.9})
    assert_equal(r.track_params, set(["Mass"]))

@with_setup(None, teardown_rep)
def test_Reprocess_4():
    bright_conf.track_nucs = set([922350])
    r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999}, name="r")
    assert_equal(r.name, 'r')
    assert_equal(r.sepeff, {922350: 0.9})
    assert_equal(r.track_params, set(["Mass"]))


#
# Tests that the fuel cycle component attributes work.
#


@with_setup(None, teardown_rep)
def test_sepeff():
    r = Reprocess()
    for value in r.sepeff.values():
        assert_equal(value, 1.0)
    r.sepeff = {922350: 0.9}
    assert_equal(r.sepeff, {922350: 0.9})


@with_setup(None, teardown_rep)
def test_track_params():
    r = Reprocess()
    assert_equal(r.track_params, set(["Mass"]))
    r.track_params = set(["Om nom nom"])
    assert_equal(r.track_params, set(["Om nom nom"]))
                        

#
# Tests that the fuel cycle component methods work.
#

@with_setup(None, teardown_rep)
def test_calc_1():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "94239": 0.99})
    r.mat_feed = Material({942390: 1.0})
    r.calc()
    assert_equal(r.mat_prod.mass, 0.99)
    assert_equal(r.mat_prod.comp[942390], 1.0) # Recall ms.comp is normalized


@with_setup(None, teardown_rep)
def test_calc_2():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "94239": 0.99})
    r.calc(Material({942390: 1.0}))
    assert_equal(r.mat_prod.mass, 0.99)
    assert_equal(r.mat_prod.comp[942390], 1.0) # Recall ms.comp is normalized

@with_setup(None, teardown_rep)
def test_initialize_1():
    bright_conf.track_nucs = set()
    r = Reprocess()
    bright_conf.track_nucs = set([922350, 922380, 942390])
    assert_equal(r.sepeff, {})        
    r.initialize({92: 0.99, 942390: 0.9})
    assert_almost_equal(r.sepeff[922350], 0.99)
    assert_almost_equal(r.sepeff[922380], 0.99)
    assert_almost_equal(r.sepeff[942390], 0.9)
        
@with_setup(None, teardown_rep)
def test_initialize_2():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999})
    r.initialize({92: 0.99, 942390: 0.9})
    assert_almost_equal(r.sepeff[922350], 0.99)
    assert_almost_equal(r.sepeff[922380], 0.99)
    assert_almost_equal(r.sepeff[942390], 0.9)
        
@with_setup(None, teardown_rep)
def test_initialize_3():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "PU2390": 0.99})
    r.initialize({92: 0.99})
    assert_almost_equal(r.sepeff[922350], 0.99)
    assert_almost_equal(r.sepeff[922380], 0.99)
    assert_almost_equal(r.sepeff[942390], 1.0)

@with_setup(None, teardown_rep)
def test_calc_params():
    bright_conf.track_nucs = set([922350, 922380, 942390])
    r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "94239": 0.99})
    r.calc(Material({942390: 1.0}))
    r.calc_params()
    assert_equal(r.params_prior_calc["Mass"],  1.00)
    assert_equal(r.params_after_calc["Mass"], 0.99)
        

if __name__ == "__main__":
    nose.main()
