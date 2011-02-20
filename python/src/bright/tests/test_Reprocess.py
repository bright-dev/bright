"""Reprocessing Component tests"""

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

bright_config = bright.bright_config
Reprocess = bright.Reprocess
MassStream = mass_stream.MassStream

class TestReprocessConstructors(TestCase):
    """Tests that the reprocessing component constructors work."""

    @classmethod
    def teardown_class(cls):
        for f in os.listdir('.'):
            if "Isos.txt" in f:
                os.remove(f)
            elif "Params.txt" in f:
                os.remove(f)
            elif f in [".h5", "r.h5"] or (".h5" in f):
                os.remove(f)

    def test_Reprocess_1(self):
        r = Reprocess()
        assert_equal(r.name, '')
        for value in r.sepeff.values():
            assert_equal(value, 1.0)
        assert_equal(r.track_params, set(["Mass"]))

    def test_Reprocess_2(self):
        bright_config.track_isos = set([922350, 922380, 942390])
        r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "94239": 0.99})
        assert_equal(r.name, '')
        assert_almost_equal(r.sepeff[922350], 0.9)
        assert_almost_equal(r.sepeff[922380], 0.999)
        assert_almost_equal(r.sepeff[942390], 0.99)
        assert_equal(r.track_params, set(["Mass"]))

    def test_Reprocess_3(self):
        bright_config.track_isos = set([922350])
        r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "94239": 0.99})
        assert_equal(r.name, '')
        assert_equal(r.sepeff, {922350: 0.9})
        assert_equal(r.track_params, set(["Mass"]))

    def test_Reprocess_4(self):
        bright_config.track_isos = set([922350])
        r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999}, name="r")
        assert_equal(r.name, 'r')
        assert_equal(r.sepeff, {922350: 0.9})
        assert_equal(r.track_params, set(["Mass"]))


class TestReprocessAttributes(TestCase):
    """Tests that the fuel cycle component attributes work."""

    @classmethod
    def teardown_class(cls):
        for f in os.listdir('.'):
            if "Isos.txt" in f:
                os.remove(f)
            elif "Params.txt" in f:
                os.remove(f)
            elif f in [".h5", "r.h5"] or (".h5" in f):
                os.remove(f)

    def test_sepeff(self):
        r = Reprocess()
        for value in r.sepeff.values():
            assert_equal(value, 1.0)
        r.sepeff = {922350: 0.9}
        assert_equal(r.sepeff, {922350: 0.9})

    def test_track_params(self):
        r = Reprocess()
        assert_equal(r.track_params, set(["Mass"]))
        r.track_params = set(["Om nom nom"])
        assert_equal(r.track_params, set(["Om nom nom"]))
                        

class TestReprocessMethods(TestCase):
    """Tests that the fuel cycle component methods work."""

    @classmethod
    def teardown_class(cls):
        for f in os.listdir('.'):
            if "Isos.txt" in f:
                os.remove(f)
            elif "Params.txt" in f:
                os.remove(f)
            elif f in [".h5", "r.h5"]:
                os.remove(f)

    def test_doCalc_1(self):
        bright_config.track_isos = set([922350, 922380, 942390])
        r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "94239": 0.99})
        r.ms_feed = MassStream({942390: 1.0})
        r.doCalc()
        assert_equal(r.ms_prod.mass, 0.99)
        assert_equal(r.ms_prod.comp[942390], 1.0) # Recall ms.comp is normalized

    def test_doCalc_2(self):
        bright_config.track_isos = set([922350, 922380, 942390])
        r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "94239": 0.99})
        r.doCalc(MassStream({942390: 1.0}))
        assert_equal(r.ms_prod.mass, 0.99)
        assert_equal(r.ms_prod.comp[942390], 1.0) # Recall ms.comp is normalized

    def test_initialize_1(self):
        bright_config.track_isos = set()
        r = Reprocess()
        bright_config.track_isos = set([922350, 922380, 942390])
        assert_equal(r.sepeff, {})        
        r.initialize({92: 0.99, 942390: 0.9})
        assert_almost_equal(r.sepeff[922350], 0.99)
        assert_almost_equal(r.sepeff[922380], 0.99)
        assert_almost_equal(r.sepeff[942390], 0.9)
        
    def test_initialize_2(self):
        bright_config.track_isos = set([922350, 922380, 942390])
        r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999})
        r.initialize({92: 0.99, 942390: 0.9})
        assert_almost_equal(r.sepeff[922350], 0.99)
        assert_almost_equal(r.sepeff[922380], 0.99)
        assert_almost_equal(r.sepeff[942390], 0.9)
        
    def test_initialize_3(self):
        bright_config.track_isos = set([922350, 922380, 942390])
        r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "PU2390": 0.99})
        r.initialize({92: 0.99})
        assert_almost_equal(r.sepeff[922350], 0.99)
        assert_almost_equal(r.sepeff[922380], 0.99)
        assert_almost_equal(r.sepeff[942390], 1.0)

    def test_setParams(self):
        bright_config.track_isos = set([922350, 922380, 942390])
        r = Reprocess(sepeff={"U235": 0.9, "922380": 0.999, "94239": 0.99})
        r.doCalc(MassStream({942390: 1.0}))
        r.setParams()
        assert_equal(r.params_prior_calc["Mass"],  1.00)
        assert_equal(r.params_after_calc["Mass"], 0.99)
        

if __name__ == "__main__":
    nose.main()
