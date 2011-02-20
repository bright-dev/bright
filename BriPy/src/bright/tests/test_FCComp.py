"""FCComp tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false

import os
import warnings
import tables as tb
import numpy as np

import BriPy
from mass_stream import MassStream

FCComp = BriPy.FCComp
bright_config = BriPy.bright_config

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
        assert_equal(fcc.params2track, set())

    def test_FCComp_2(self):
        fcc = FCComp(["Spam", "Spam", "Eggs", "Spam"])
        assert_equal(fcc.name, '')
        assert_equal(fcc.params2track, ["Eggs", "Spam"])

    def test_FCComp_2(self):
        fcc = FCComp(set(), "Waldo")
        assert_equal(fcc.name, 'Waldo')
        assert_equal(fcc.params2track, set())


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

    def test_IsosIn_Empty(self):
        fcc = FCComp()
        assert_equal(fcc.IsosIn.comp, {})
        assert_equal(fcc.IsosIn.mass, -1.0)
        assert_equal(fcc.IsosIn.name, '')

    def test_IsosIn_Filled(self):
        fcc = FCComp()
        ms = MassStream({922350: 1.0}, -1, "ms")
        fcc.IsosIn = ms
        assert_equal(fcc.IsosIn.comp, {922350: 1.0})
        assert_equal(fcc.IsosIn.mass, 1.0)
        assert_equal(fcc.IsosIn.name, 'ms')

    def test_IsosOut_Empty(self):
        fcc = FCComp()
        assert_equal(fcc.IsosOut.comp, {})
        assert_equal(fcc.IsosOut.mass, -1.0)
        assert_equal(fcc.IsosOut.name, '')

    def test_IsosOut_Filled(self):
        fcc = FCComp()
        ms = MassStream({922350: 1.0}, -1, "ms")
        fcc.IsosOut = ms
        assert_equal(fcc.IsosOut.comp, {922350: 1.0})
        assert_equal(fcc.IsosOut.mass, 1.0)
        assert_equal(fcc.IsosOut.name, 'ms')

    def test_ParamsIn_Empty(self):
        fcc = FCComp()
        assert_equal(fcc.ParamsIn, {})

    def test_ParamsIn_Filled(self):
        fcc = FCComp(set(["Mass"]))
        fcc.ParamsIn = {"Mass": 1.0}
        assert_equal(fcc.ParamsIn, {"Mass": 1.0})

    def test_ParamsOut_Empty(self):
        fcc = FCComp()
        assert_equal(fcc.ParamsOut, {})

    def test_ParamsOut_Filled(self):
        fcc = FCComp(set(["Mass"]))
        fcc.ParamsOut = {"Mass": 1.0}
        assert_equal(fcc.ParamsOut, {"Mass": 1.0})

    def test_PassNum(self):
        fcc = FCComp()
        assert_equal(fcc.PassNum, 0)
        fcc.PassNum = 42
        assert_equal(fcc.PassNum, 42)

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

    def test_params2track(self):
        fcc = FCComp()
        assert_equal(fcc.params2track, set())
        fcc.params2track = set(["Dave"])
        assert_equal(fcc.params2track, set(["Dave"]))


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

    def test_setParams(self):
        fcc = FCComp(set(["Mass"]))
        fcc.setParams()
        assert_equal(fcc.ParamsIn["Mass"],  0.0)
        assert_equal(fcc.ParamsOut["Mass"], 0.0)

    def test_writeIsoPass(self):
        bright_config.isos2track = set([922350])
        fcc = FCComp()
        fcc.IsosIn  = MassStream({922350: 1.0})
        fcc.IsosOut = MassStream({922350: 0.5})
        fcc.writeIsoPass()

    def test_writeParamPass(self):
        fcc = FCComp(set(["Mass"]))
        fcc.setParams()
        fcc.writeParamPass()

    def test_writeText(self):
        bright_config.isos2track = set([922350])
        fcc = FCComp(set(["Mass"]))
        fcc.IsosIn  = MassStream({922350: 1.0})
        fcc.IsosOut = MassStream({922350: 0.5})
        fcc.setParams()
        fcc.writeText()

    def test_writeHDF5_1(self):
        bright_config.isos2track = set([922350])
        bright_config.write_hdf5 = True
        fcc = FCComp(set(), 'fcc')
        fcc.IsosIn  = MassStream({922350: 1.0})
        fcc.IsosOut = MassStream({922350: 0.5})
        fcc.PassNum = 1
        fcc.writeHDF5()

    def test_writeHDF5_2(self):
        bright_config.isos2track = set([922350])
        bright_config.write_hdf5 = True
        fcc = FCComp(set(["Mass"]), 'fcc')
        fcc.IsosIn  = MassStream({922350: 1.0})
        fcc.IsosOut = MassStream({922350: 0.5})
        fcc.setParams()
        fcc.PassNum = 1
        fcc.writeHDF5()

    def test_writeout_1(self):
        """Text only."""
        bright_config.isos2track = set([922350])
        fcc = FCComp(set(["Mass"]))
        fcc.IsosIn  = MassStream({922350: 1.0})
        fcc.IsosOut = MassStream({922350: 0.5})
        fcc.writeout()
    
    def test_writeout_2(self):
        """HDF5 only."""
        bright_config.isos2track = set([922350])
        bright_config.write_hdf5 = True
        bright_config.write_text = False
        fcc = FCComp(set(["Mass"]), 'fcc')
        fcc.IsosIn  = MassStream({922350: 1.0})
        fcc.IsosOut = MassStream({922350: 0.5})
        fcc.writeout()

    def test_writeout_3(self):
        """HDF5 & Text output."""
        bright_config.isos2track = set([922350])
        bright_config.write_hdf5 = True
        bright_config.write_text = True
        fcc = FCComp(set(["Mass"]), 'fcc')
        fcc.IsosIn  = MassStream({922350: 1.0})
        fcc.IsosOut = MassStream({922350: 0.5})
        fcc.writeout()


if __name__ == "__main__":
    nose.main()
