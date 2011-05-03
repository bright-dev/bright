"""FCComp tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false

import os
import warnings
import tables as tb
import numpy as np

import bright
from mass_stream import MassStream

FCComp = bright.FCComp
bright_config = bright.bright_config

class TestFCCompConstructors(TestCase):
    """Tests that the fuel cycle component constructors work."""

    @classmethod
    def teardown_class(cls):
        for f in os.listdir('.'):
            if "Isos.txt" in f:
                os.remove(f)
            elif "Params.txt" in f:
                os.remove(f)
            elif f in [".h5", "fcc.h5", "fuel_cycle.h5"]:
                os.remove(f)

    def test_FCComp_1(self):
        fcc = FCComp()
        assert_equal(fcc.name, '')
        assert_equal(fcc.track_params, set())

    def test_FCComp_2(self):
        fcc = FCComp(["Spam", "Spam", "Eggs", "Spam"])
        assert_equal(fcc.name, '')
        assert_equal(fcc.track_params, ["Eggs", "Spam"])

    def test_FCComp_2(self):
        fcc = FCComp(set(), "Waldo")
        assert_equal(fcc.name, 'Waldo')
        assert_equal(fcc.track_params, set())


class TestFCCompAttributes(TestCase):
    """Tests that the fuel cycle component attributes work."""

    @classmethod
    def teardown_class(cls):
        for f in os.listdir('.'):
            if "Isos.txt" in f:
                os.remove(f)
            elif "Params.txt" in f:
                os.remove(f)
            elif f in [".h5", "fcc.h5", "fuel_cycle.h5"]:
                os.remove(f)

    def test_ms_feed_Empty(self):
        fcc = FCComp()
        assert_equal(fcc.ms_feed.comp, {})
        assert_equal(fcc.ms_feed.mass, -1.0)
        assert_equal(fcc.ms_feed.name, '')

    def test_ms_feed_Filled(self):
        fcc = FCComp()
        ms = MassStream({922350: 1.0}, -1, "ms")
        fcc.ms_feed = ms
        assert_equal(fcc.ms_feed.comp, {922350: 1.0})
        assert_equal(fcc.ms_feed.mass, 1.0)
        assert_equal(fcc.ms_feed.name, 'ms')

    def test_ms_prod_Empty(self):
        fcc = FCComp()
        assert_equal(fcc.ms_prod.comp, {})
        assert_equal(fcc.ms_prod.mass, -1.0)
        assert_equal(fcc.ms_prod.name, '')

    def test_ms_prod_Filled(self):
        fcc = FCComp()
        ms = MassStream({922350: 1.0}, -1, "ms")
        fcc.ms_prod = ms
        assert_equal(fcc.ms_prod.comp, {922350: 1.0})
        assert_equal(fcc.ms_prod.mass, 1.0)
        assert_equal(fcc.ms_prod.name, 'ms')

    def test_params_prior_calc_Empty(self):
        fcc = FCComp()
        assert_equal(fcc.params_prior_calc, {})

    def test_params_prior_calc_Filled(self):
        fcc = FCComp(set(["Mass"]))
        fcc.params_prior_calc = {"Mass": 1.0}
        assert_equal(fcc.params_prior_calc, {"Mass": 1.0})

    def test_params_after_calc_Empty(self):
        fcc = FCComp()
        assert_equal(fcc.params_after_calc, {})

    def test_params_after_calc_Filled(self):
        fcc = FCComp(set(["Mass"]))
        fcc.params_after_calc = {"Mass": 1.0}
        assert_equal(fcc.params_after_calc, {"Mass": 1.0})

    def test_pass_num(self):
        fcc = FCComp()
        assert_equal(fcc.pass_num, 0)
        fcc.pass_num = 42
        assert_equal(fcc.pass_num, 42)

    def test_name(self):
        fcc = FCComp()
        assert_equal(fcc.name, '')
        fcc.name = "fcc"
        assert_equal(fcc.name, 'fcc')

    def test_natural_name1(self):
        fcc = FCComp(set(), "Word")
        assert_equal(fcc.natural_name, 'Word')

    def test_natural_name2(self):
        fcc = FCComp(set(), "Word to your mother")
        assert_equal(fcc.natural_name, 'Word_to_your_mother')

    def test_natural_name3(self):
        fcc = FCComp(set(), "1 isthe ")
        assert_equal(fcc.natural_name, '_1_isthe_')

    def test_natural_name4(self):
        fcc = FCComp(set(), "\t Try\nMe...$")
        assert_equal(fcc.natural_name, '__Try_Me')

    def test_track_params(self):
        fcc = FCComp()
        assert_equal(fcc.track_params, set())
        fcc.track_params = set(["Dave"])
        assert_equal(fcc.track_params, set(["Dave"]))


class TestFCCompMethods(TestCase):
    """Tests that the fuel cycle component methods work."""

    @classmethod
    def teardown_class(cls):
        for f in os.listdir('.'):
            if "Isos.txt" in f:
                os.remove(f)
            elif "Params.txt" in f:
                os.remove(f)
            elif f in [".h5", "fcc.h5", "fuel_cycle.h5"]:
                os.remove(f)

    def test_calc_params(self):
        fcc = FCComp(set(["Mass"]))
        fcc.calc_params()
        assert_equal(fcc.params_prior_calc["Mass"],  0.0)
        assert_equal(fcc.params_after_calc["Mass"], 0.0)

    def test_write_ms_pass(self):
        bright_config.track_isos = set([922350])
        fcc = FCComp()
        fcc.ms_feed  = MassStream({922350: 1.0})
        fcc.ms_prod = MassStream({922350: 0.5})
        fcc.write_ms_pass()

    def test_write_params_pass(self):
        fcc = FCComp(set(["Mass"]))
        fcc.calc_params()
        fcc.write_params_pass()

    def test_write_text(self):
        bright_config.track_isos = set([922350])
        fcc = FCComp(set(["Mass"]))
        fcc.ms_feed  = MassStream({922350: 1.0})
        fcc.ms_prod = MassStream({922350: 0.5})
        fcc.calc_params()
        fcc.write_text()

    def test_write_hdf5_1(self):
        bright_config.track_isos = set([922350])
        bright_config.write_hdf5 = True
        fcc = FCComp(set(), 'fcc')
        fcc.ms_feed  = MassStream({922350: 1.0})
        fcc.ms_prod = MassStream({922350: 0.5})
        fcc.pass_num = 1
        fcc.write_hdf5()

    def test_write_hdf5_2(self):
        bright_config.track_isos = set([922350])
        bright_config.write_hdf5 = True
        fcc = FCComp(set(["Mass"]), 'fcc')
        fcc.ms_feed  = MassStream({922350: 1.0})
        fcc.ms_prod = MassStream({922350: 0.5})
        fcc.calc_params()
        fcc.pass_num = 1
        fcc.write_hdf5()

    def test_write_1(self):
        """Text only."""
        bright_config.track_isos = set([922350])
        fcc = FCComp(set(["Mass"]))
        fcc.ms_feed  = MassStream({922350: 1.0})
        fcc.ms_prod = MassStream({922350: 0.5})
        fcc.write()
    
    def test_write_2(self):
        """HDF5 only."""
        bright_config.track_isos = set([922350])
        bright_config.write_hdf5 = True
        bright_config.write_text = False
        fcc = FCComp(set(["Mass"]), 'fcc')
        fcc.ms_feed  = MassStream({922350: 1.0})
        fcc.ms_prod = MassStream({922350: 0.5})
        fcc.write()

    def test_write_3(self):
        """HDF5 & Text output."""
        bright_config.track_isos = set([922350])
        bright_config.write_hdf5 = True
        bright_config.write_text = True
        fcc = FCComp(set(["Mass"]), 'fcc')
        fcc.ms_feed  = MassStream({922350: 1.0})
        fcc.ms_prod = MassStream({922350: 0.5})
        fcc.write()


if __name__ == "__main__":
    nose.main()
