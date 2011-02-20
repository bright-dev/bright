"""LightWaterReactor1G Component tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false

import os
import warnings
import tables as tb
import numpy as np

import BriPy
import mass_stream


LightWaterReactor1G = BriPy.LightWaterReactor1G
MassStream = mass_stream.MassStream
bright_config = BriPy.bright_config

def general_teardown():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "lwr.h5"]:
            os.remove(f)

class TestLightWaterReactorConstructors(TestCase):
    """Tests that the storage component constructors work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_LWRDefaults(self):
        lwrd = BriPy.LWRDefaults()
        assert_equal(lwrd.batches, 3)
        assert_equal(lwrd.flux, 4.0*(10.0**14))
        assert_equal(lwrd.FuelForm["IHM"], 1.0)
        assert_equal(lwrd.FuelForm["O16"], 2.0)
        assert_equal(lwrd.CoolantForm["O16"], 1.0)
        assert_equal(lwrd.CoolantForm["H1"], 2.0)
        assert_equal(lwrd.CoolantForm["B10"], 0.199 * 550 * (10.0**-6))
        assert_equal(lwrd.CoolantForm["B11"], 0.801 * 550 * (10.0**-6))
        assert_equal(lwrd.FuelDensity, 10.7)
        assert_equal(lwrd.CoolantDensity, 0.73)
        assert_equal(lwrd.pnl, 0.98)
        assert_equal(lwrd.BUt, 0.0)
        assert_true(lwrd.useDisadvantage)
        assert_equal(lwrd.LatticeType, "Cylindrical")
        assert_true(lwrd.HydrogenRescale)
        assert_equal(lwrd.Radius, 0.412)
        assert_equal(lwrd.Length, 1.33)
        assert_equal(lwrd.open_slots, 25.0)
        assert_equal(lwrd.total_slots, 289.0)

    def test_LightWaterReactor1G_1(self):
        lwr = LightWaterReactor1G()
        assert_equal(lwr.name, '')
        assert_equal(lwr.params2track, set(["ACT", "BUd", "FP", "LAN", "TRU", "U"]))
        assert_equal(lwr.B, 3)
        assert_equal(lwr.phi, 4.0*(10.0**14))
        assert_equal(lwr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(lwr.FuelChemicalForm["O16"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["O16"], 1.0)
        assert_equal(lwr.CoolantChemicalForm["H1"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["B10"], 0.199 * 550 * (10.0**-6))
        assert_equal(lwr.CoolantChemicalForm["B11"], 0.801 * 550 * (10.0**-6))
        assert_equal(lwr.rhoF, 10.7)
        assert_equal(lwr.rhoC, 0.73)
        assert_equal(lwr.P_NL, 0.98)
        assert_equal(lwr.TargetBU, 0.0)
        assert_true(lwr.useZeta)
        assert_equal(lwr.Lattice, "Cylindrical")
        assert_true(lwr.H_XS_Rescale)
        assert_equal(lwr.r, 0.412)
        assert_equal(lwr.l, 1.33)
        assert_equal(lwr.S_O, 25.0)
        assert_equal(lwr.S_T, 289.0)

    def test_LightWaterReactor1G_2(self):
        lf = os.getenv("BRIGHT_DATA") + "/LWR.h5"
        lwr = LightWaterReactor1G(libfile=lf)
        assert_equal(lwr.libfile, lf)
        assert_equal(lwr.name, '')
        assert_equal(lwr.params2track, set(["ACT", "BUd", "FP", "LAN", "TRU", "U"]))
        assert_equal(lwr.B, 3)
        assert_equal(lwr.phi, 4.0*(10.0**14))
        assert_equal(lwr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(lwr.FuelChemicalForm["O16"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["O16"], 1.0)
        assert_equal(lwr.CoolantChemicalForm["H1"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["B10"], 0.199 * 550 * (10.0**-6))
        assert_equal(lwr.CoolantChemicalForm["B11"], 0.801 * 550 * (10.0**-6))
        assert_equal(lwr.rhoF, 10.7)
        assert_equal(lwr.rhoC, 0.73)
        assert_equal(lwr.P_NL, 0.98)
        assert_equal(lwr.TargetBU, 0.0)
        assert_true(lwr.useZeta)
        assert_equal(lwr.Lattice, "Cylindrical")
        assert_true(lwr.H_XS_Rescale)
        assert_equal(lwr.r, 0.412)
        assert_equal(lwr.l, 1.33)
        assert_equal(lwr.S_O, 25.0)
        assert_equal(lwr.S_T, 289.0)

    def test_LightWaterReactor1G_3(self):
        lf = os.getenv("BRIGHT_DATA") + "/LWR.h5"
        lwr = LightWaterReactor1G(libfile=lf, name="lwr")
        assert_equal(lwr.libfile, lf)
        assert_equal(lwr.name, 'lwr')
        assert_equal(lwr.params2track, set(["ACT", "BUd", "FP", "LAN", "TRU", "U"]))
        assert_equal(lwr.B, 3)
        assert_equal(lwr.phi, 4.0*(10.0**14))
        assert_equal(lwr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(lwr.FuelChemicalForm["O16"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["O16"], 1.0)
        assert_equal(lwr.CoolantChemicalForm["H1"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["B10"], 0.199 * 550 * (10.0**-6))
        assert_equal(lwr.CoolantChemicalForm["B11"], 0.801 * 550 * (10.0**-6))
        assert_equal(lwr.rhoF, 10.7)
        assert_equal(lwr.rhoC, 0.73)
        assert_equal(lwr.P_NL, 0.98)
        assert_equal(lwr.TargetBU, 0.0)
        assert_true(lwr.useZeta)
        assert_equal(lwr.Lattice, "Cylindrical")
        assert_true(lwr.H_XS_Rescale)
        assert_equal(lwr.r, 0.412)
        assert_equal(lwr.l, 1.33)
        assert_equal(lwr.S_O, 25.0)
        assert_equal(lwr.S_T, 289.0)

    def test_LightWaterReactor1G_4(self):
        rp = BriPy.LWRDefaults()
        rp.BUt = 50.0
        lwr = LightWaterReactor1G(reactor_parameters=rp)
        assert_equal(lwr.name, '')
        assert_equal(lwr.params2track, set(["ACT", "BUd", "FP", "LAN", "TRU", "U"]))
        assert_equal(lwr.B, 3)
        assert_equal(lwr.phi, 4.0*(10.0**14))
        assert_equal(lwr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(lwr.FuelChemicalForm["O16"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["O16"], 1.0)
        assert_equal(lwr.CoolantChemicalForm["H1"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["B10"], 0.199 * 550 * (10.0**-6))
        assert_equal(lwr.CoolantChemicalForm["B11"], 0.801 * 550 * (10.0**-6))
        assert_equal(lwr.rhoF, 10.7)
        assert_equal(lwr.rhoC, 0.73)
        assert_equal(lwr.P_NL, 0.98)
        assert_equal(lwr.TargetBU, 50.0)
        assert_true(lwr.useZeta)
        assert_equal(lwr.Lattice, "Cylindrical")
        assert_true(lwr.H_XS_Rescale)
        assert_equal(lwr.r, 0.412)
        assert_equal(lwr.l, 1.33)
        assert_equal(lwr.S_O, 25.0)
        assert_equal(lwr.S_T, 289.0)

    def test_LightWaterReactor1G_5(self):
        rp = BriPy.LWRDefaults()
        rp.BUt = 50.0
        lwr = LightWaterReactor1G(reactor_parameters=rp, name='lwr')
        assert_equal(lwr.name, 'lwr')
        assert_equal(lwr.params2track, set(["ACT", "BUd", "FP", "LAN", "TRU", "U"]))
        assert_equal(lwr.B, 3)
        assert_equal(lwr.phi, 4.0*(10.0**14))
        assert_equal(lwr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(lwr.FuelChemicalForm["O16"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["O16"], 1.0)
        assert_equal(lwr.CoolantChemicalForm["H1"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["B10"], 0.199 * 550 * (10.0**-6))
        assert_equal(lwr.CoolantChemicalForm["B11"], 0.801 * 550 * (10.0**-6))
        assert_equal(lwr.rhoF, 10.7)
        assert_equal(lwr.rhoC, 0.73)
        assert_equal(lwr.P_NL, 0.98)
        assert_equal(lwr.TargetBU, 50.0)
        assert_true(lwr.useZeta)
        assert_equal(lwr.Lattice, "Cylindrical")
        assert_true(lwr.H_XS_Rescale)
        assert_equal(lwr.r, 0.412)
        assert_equal(lwr.l, 1.33)
        assert_equal(lwr.S_O, 25.0)
        assert_equal(lwr.S_T, 289.0)

    def test_LightWaterReactor1G_6(self):
        lf = os.getenv("BRIGHT_DATA") + "/LWR.h5"
        rp = BriPy.LWRDefaults()
        rp.BUt = 50.0
        lwr = LightWaterReactor1G(libfile=lf, reactor_parameters=rp)
        assert_equal(lwr.libfile, lf)
        assert_equal(lwr.name, '')
        assert_equal(lwr.params2track, set(["ACT", "BUd", "FP", "LAN", "TRU", "U"]))
        assert_equal(lwr.B, 3)
        assert_equal(lwr.phi, 4.0*(10.0**14))
        assert_equal(lwr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(lwr.FuelChemicalForm["O16"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["O16"], 1.0)
        assert_equal(lwr.CoolantChemicalForm["H1"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["B10"], 0.199 * 550 * (10.0**-6))
        assert_equal(lwr.CoolantChemicalForm["B11"], 0.801 * 550 * (10.0**-6))
        assert_equal(lwr.rhoF, 10.7)
        assert_equal(lwr.rhoC, 0.73)
        assert_equal(lwr.P_NL, 0.98)
        assert_equal(lwr.TargetBU, 50.0)
        assert_true(lwr.useZeta)
        assert_equal(lwr.Lattice, "Cylindrical")
        assert_true(lwr.H_XS_Rescale)
        assert_equal(lwr.r, 0.412)
        assert_equal(lwr.l, 1.33)
        assert_equal(lwr.S_O, 25.0)
        assert_equal(lwr.S_T, 289.0)

    def test_LightWaterReactor1G_7(self):
        lf = os.getenv("BRIGHT_DATA") + "/LWR.h5"
        rp = BriPy.LWRDefaults()
        rp.BUt = 50.0
        lwr = LightWaterReactor1G(libfile=lf, reactor_parameters=rp, name='lwr')
        assert_equal(lwr.libfile, lf)
        assert_equal(lwr.name, 'lwr')
        assert_equal(lwr.params2track, set(["ACT", "BUd", "FP", "LAN", "TRU", "U"]))
        assert_equal(lwr.B, 3)
        assert_equal(lwr.phi, 4.0*(10.0**14))
        assert_equal(lwr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(lwr.FuelChemicalForm["O16"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["O16"], 1.0)
        assert_equal(lwr.CoolantChemicalForm["H1"], 2.0)
        assert_equal(lwr.CoolantChemicalForm["B10"], 0.199 * 550 * (10.0**-6))
        assert_equal(lwr.CoolantChemicalForm["B11"], 0.801 * 550 * (10.0**-6))
        assert_equal(lwr.rhoF, 10.7)
        assert_equal(lwr.rhoC, 0.73)
        assert_equal(lwr.P_NL, 0.98)
        assert_equal(lwr.TargetBU, 50.0)
        assert_true(lwr.useZeta)
        assert_equal(lwr.Lattice, "Cylindrical")
        assert_true(lwr.H_XS_Rescale)
        assert_equal(lwr.r, 0.412)
        assert_equal(lwr.l, 1.33)
        assert_equal(lwr.S_O, 25.0)
        assert_equal(lwr.S_T, 289.0)


class TestLightWaterReactor1GAttributes(TestCase):
    """Tests that the fuel cycle component attributes work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_params2track(self):
        lwr = LightWaterReactor1G()
        assert_equal(lwr.params2track, set(["ACT", "BUd", "FP", "LAN", "TRU", "U"]))
        lwr.params2track = set(["Mass"])
        assert_equal(lwr.params2track, set(["Mass"]))

class TestLightWaterReactor1GMethods(TestCase):
    """Tests that the fuel cycle component methods work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_setParams(self):
        lf = os.getenv("BRIGHT_DATA") + "/LWR.h5"
        BriPy.load_isos2track_hdf5(lf)
        rp = BriPy.LWRDefaults()
        rp.BUt = 50.0
        lwr = LightWaterReactor1G(libfile=lf, reactor_parameters=rp, name='lwr')
        lwr.doCalc(MassStream({922350: 0.05, 922380:0.95}))
        lwr.setParams()
        assert_equal(lwr.ParamsIn["BUd"],  0.0)
        assert_equal(lwr.ParamsOut["BUd"], lwr.BUd)
        assert_equal(lwr.ParamsIn["U"],  lwr.InU.mass)
        assert_equal(lwr.ParamsOut["U"], lwr.OutU.mass)
        assert_equal(lwr.ParamsIn["TRU"],  lwr.InTRU.mass)
        assert_equal(lwr.ParamsOut["TRU"], lwr.OutTRU.mass)
        assert_equal(lwr.ParamsIn["ACT"],  lwr.InACT.mass)
        assert_equal(lwr.ParamsOut["ACT"], lwr.OutACT.mass)
        assert_equal(lwr.ParamsIn["LAN"],  lwr.InLAN.mass)
        assert_equal(lwr.ParamsOut["LAN"], lwr.OutLAN.mass)
        assert_equal(lwr.ParamsIn["FP"],  1.0 - lwr.InACT.mass - lwr.InLAN.mass)
        assert_equal(lwr.ParamsOut["FP"], 1.0 - lwr.OutACT.mass - lwr.OutLAN.mass)
        

# Put Integral tests here, if desired.


if __name__ == "__main__":
    nose.main()
