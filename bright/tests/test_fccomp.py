"""FCComp tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, with_setup

import os
import warnings
import tables as tb
import numpy as np

import bright
import bright.fccomp

from pyne.material import Material

FCComp = bright.fccomp.FCComp
bright_conf = bright.bright_conf

#
# Fixtures
#

def setup_fccomp():
    pass
    
def teardown_fccomp():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "fcc.h5", "fuel_cycle.h5"]:
            os.remove(f)

#
# Tests
#


@with_setup(None, teardown_fccomp)
def test_FCComp_1():
    fcc = FCComp()
    assert_equal(fcc.name, '')
    assert_equal(fcc.track_params, set())


@with_setup(None, teardown_fccomp)
def test_FCComp_2():
    fcc = FCComp(set(["Spam", "Spam", "Eggs", "Spam"]))
    assert_equal(fcc.name, '')
    assert_equal(fcc.track_params, set(["Eggs", "Spam"]))

@with_setup(None, teardown_fccomp)
def test_FCComp_3():
    fcc = FCComp(set(), "Waldo")
    assert_equal(fcc.name, 'Waldo')
    assert_equal(fcc.track_params, set())


@with_setup(None, teardown_fccomp)
def test_mat_feed_empty():
    fcc = FCComp()
    assert_equal(fcc.mat_feed.comp, {})
    assert_equal(fcc.mat_feed.mass, -1.0)
    assert_equal(fcc.mat_feed.name, '')


@with_setup(None, teardown_fccomp)
def test_mat_feed_filled():
    fcc = FCComp()
    mat = Material({922350: 1.0}, -1, "mat")
    fcc.mat_feed = mat
    assert_equal(fcc.mat_feed.comp, {922350: 1.0})
    assert_equal(fcc.mat_feed.mass, 1.0)
    assert_equal(fcc.mat_feed.name, 'mat')

@with_setup(None, teardown_fccomp)
def test_mat_prod_empty():
    fcc = FCComp()
    assert_equal(fcc.mat_prod.comp, {})
    assert_equal(fcc.mat_prod.mass, -1.0)
    assert_equal(fcc.mat_prod.name, '')

@with_setup(None, teardown_fccomp)
def test_mat_prod_filled():
    fcc = FCComp()
    mat = Material({922350: 1.0}, -1, "mat")
    fcc.mat_prod = mat
    assert_equal(fcc.mat_prod.comp, {922350: 1.0})
    assert_equal(fcc.mat_prod.mass, 1.0)
    assert_equal(fcc.mat_prod.name, 'mat')

@with_setup(None, teardown_fccomp)
def test_params_prior_calc_empty():
    fcc = FCComp()
    assert_equal(fcc.params_prior_calc, {})

@with_setup(None, teardown_fccomp)
def test_params_prior_calc_filled():
    fcc = FCComp(set(["Mass"]))
    fcc.params_prior_calc = {"Mass": 1.0}
    assert_equal(fcc.params_prior_calc, {"Mass": 1.0})

@with_setup(None, teardown_fccomp)
def test_params_after_calc_empty():
    fcc = FCComp()
    assert_equal(fcc.params_after_calc, {})

@with_setup(None, teardown_fccomp)
def test_params_after_calc_filled():
    fcc = FCComp(set(["Mass"]))
    fcc.params_after_calc = {"Mass": 1.0}
    assert_equal(fcc.params_after_calc, {"Mass": 1.0})

@with_setup(None, teardown_fccomp)
def test_pass_num():
    fcc = FCComp()
    assert_equal(fcc.pass_num, 0)
    fcc.pass_num = 42
    assert_equal(fcc.pass_num, 42)

@with_setup(None, teardown_fccomp)
def test_name():
    fcc = FCComp()
    assert_equal(fcc.name, '')
    fcc.name = "fcc"
    assert_equal(fcc.name, 'fcc')

@with_setup(None, teardown_fccomp)
def test_natural_name1():
    fcc = FCComp(set(), "Word")
    assert_equal(fcc.natural_name, 'Word')

@with_setup(None, teardown_fccomp)
def test_natural_name2():
    fcc = FCComp(set(), "Word to your mother")
    assert_equal(fcc.natural_name, 'Word_to_your_mother')

@with_setup(None, teardown_fccomp)
def test_natural_name3():
    fcc = FCComp(set(), "1 isthe ")
    assert_equal(fcc.natural_name, '_1_isthe_')

@with_setup(None, teardown_fccomp)
def test_natural_name4():
    fcc = FCComp(set(), "\t Try\nMe...$")
    assert_equal(fcc.natural_name, '__Try_Me')

@with_setup(None, teardown_fccomp)
def test_track_params():
    fcc = FCComp()
    assert_equal(fcc.track_params, set())
    fcc.track_params = set(["Dave"])
    assert_equal(fcc.track_params, set(["Dave"]))


#
# Test Methods
#

@with_setup(None, teardown_fccomp)
def test_calc_params():
    fcc = FCComp(set(["Mass"]))
    fcc.calc_params()
    assert_equal(fcc.params_prior_calc["Mass"],  0.0)
    assert_equal(fcc.params_after_calc["Mass"], 0.0)

@with_setup(None, teardown_fccomp)
def test_write_mat_pass():
    bright_conf.track_nucs = set([922350])
    fcc = FCComp()
    fcc.mat_feed  = Material({922350: 1.0})
    fcc.mat_prod = Material({922350: 0.5})
    fcc.write_mat_pass()

@with_setup(None, teardown_fccomp)
def test_write_params_pass():
    fcc = FCComp(set(["Mass"]))
    fcc.calc_params()
    fcc.write_params_pass()

@with_setup(None, teardown_fccomp)
def test_write_text():
    bright_conf.track_nucs = set([922350])
    fcc = FCComp(set(["Mass"]))
    fcc.mat_feed  = Material({922350: 1.0})
    fcc.mat_prod = Material({922350: 0.5})
    fcc.calc_params()
    fcc.write_text()

@with_setup(None, teardown_fccomp)
def test_write_hdf5_1():
    bright_conf.track_nucs = set([922350])
    bright_conf.write_hdf5 = True
    fcc = FCComp(set(), 'fcc')
    fcc.mat_feed  = Material({922350: 1.0})
    fcc.mat_prod = Material({922350: 0.5})
    fcc.pass_num = 1
    fcc.write_hdf5()

@with_setup(None, teardown_fccomp)
def test_write_hdf5_2():
    bright_conf.track_nucs = set([922350])
    bright_conf.write_hdf5 = True
    fcc = FCComp(set(["Mass"]), 'fcc')
    fcc.mat_feed  = Material({922350: 1.0})
    fcc.mat_prod = Material({922350: 0.5})
    fcc.calc_params()
    fcc.pass_num = 1
    fcc.write_hdf5()

@with_setup(None, teardown_fccomp)
def test_write_1():
    "Text only."
    bright_conf.track_nucs = set([922350])
    fcc = FCComp(set(["Mass"]))
    fcc.mat_feed  = Material({922350: 1.0})
    fcc.mat_prod = Material({922350: 0.5})
    fcc.write()
    
@with_setup(None, teardown_fccomp)
def test_write_2():
    "HDF5 only."
    bright_conf.track_nucs = set([922350])
    bright_conf.write_hdf5 = True
    bright_conf.write_text = False
    fcc = FCComp(set(["Mass"]), 'fcc')
    fcc.mat_feed  = Material({922350: 1.0})
    fcc.mat_prod = Material({922350: 0.5})
    fcc.write()

@with_setup(None, teardown_fccomp)
def test_write_3():
    "HDF5 & Text output."
    bright_conf.track_nucs = set([922350])
    bright_conf.write_hdf5 = True
    bright_conf.write_text = True
    fcc = FCComp(set(["Mass"]), 'fcc')
    fcc.mat_feed  = Material({922350: 1.0})
    fcc.mat_prod = Material({922350: 0.5})
    fcc.write()


if __name__ == "__main__":
    nose.main()
