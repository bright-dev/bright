"""Storage Component tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false

import os
import warnings
import tables as tb
import numpy as np

import bright
import mass_stream

Storage = bright.Storage
MassStream = mass_stream.MassStream
bright_config = bright.bright_config

class TestStorageConstructors(TestCase):
    """Tests that the storage component constructors work."""

    @classmethod
    def teardown_class(cls):
        for f in os.listdir('.'):
            if "Isos.txt" in f:
                os.remove(f)
            elif "Params.txt" in f:
                os.remove(f)
            elif f in [".h5", "s.h5"]:
                os.remove(f)

    def test_Storage_1(self):
        s = Storage()
        assert_equal(s.name, '')
        assert_equal(s.track_params, set(["Mass"]))

    def test_Storage_2(self):
        s = Storage(name="s")
        assert_equal(s.name, 's')
        assert_equal(s.track_params, set(["Mass"]))


class TestStorageAttributes(TestCase):
    """Tests that the fuel cycle component attributes work."""

    @classmethod
    def teardown_class(cls):
        for f in os.listdir('.'):
            if "Isos.txt" in f:
                os.remove(f)
            elif "Params.txt" in f:
                os.remove(f)
            elif f in [".h5", "s.h5"]:
                os.remove(f)

    def test_decay_time(self):
        s = Storage()
        s.decay_time = 0.0
        assert_equal(s.decay_time, 0.0)
        s.decay_time = 628        
        assert_equal(s.decay_time, 628.0)

    def test_track_params(self):
        s = Storage()
        assert_equal(s.track_params, set(["Mass"]))
        s.track_params = set(["Om nom nom"])
        assert_equal(s.track_params, set(["Om nom nom"]))
                        

class TestStorageMethods(TestCase):
    """Tests that the fuel cycle component methods work."""

    @classmethod
    def teardown_class(cls):
        for f in os.listdir('.'):
            if "Isos.txt" in f:
                os.remove(f)
            elif "Params.txt" in f:
                os.remove(f)
            elif f in [".h5", "s.h5"]:
                os.remove(f)

    def test_calc_1(self):
        bright_config.track_nucs = set([922350, 922380, 942390])
        s = Storage()
        s.decay_time = 0.0
        s.mat_feed = MassStream({942390: 1.0})
        s.calc()
        assert_equal(s.mat_prod.mass, 1.0)
        assert_almost_equal(s.mat_prod.comp[942390], 1.0) 

    def test_calc_2(self):
        bright_config.track_nucs = set([922350, 922380, 942390])
        s = Storage()
        s.decay_time = 0.0
        s.calc(MassStream({942390: 1.0}))
        assert_equal(s.mat_prod.mass, 1.0)
        assert_equal(s.mat_prod.comp[942390], 1.0) 

    def test_calc_3(self):
        bright_config.track_nucs = set([922350, 922380, 942390])
        s = Storage()
        s.mat_feed = MassStream({942390: 1.0})
        s.calc(decay_time=24110*365.25*24*3600)
        assert(s.mat_prod.mass < 1.0)
        assert_almost_equal(s.mat_prod.comp[942390], 0.5, 3) 

    def test_calc_4(self):
        bright_config.track_nucs = set([922350, 922380, 942390])
        s = Storage()
        s.calc(MassStream({942390: 1.0}), 24110*365.25*24*3600)
        assert(s.mat_prod.mass < 1.0)
        assert_almost_equal(s.mat_prod.comp[942390], 0.5, 3) 

    def test_calc_params(self):
        bright_config.track_nucs = set([922350, 922380, 942390])
        s = Storage()
        s.calc(MassStream({942390: 1.0}), 24110*365.25*24*3600)
        s.calc_params()
        assert_equal(s.params_prior_calc["Mass"],  1.00)
        assert(0.5 < s.params_after_calc["Mass"] < 1.0)
        

if __name__ == "__main__":
    nose.main()
