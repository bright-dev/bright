"""FastReactor1G Component tests"""

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

FastReactor1G = bright.FastReactor1G
MassStream = mass_stream.MassStream
bright_config = bright.bright_config


def general_teardown():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "fr.h5", "fuel_cycle.h5"]:
            os.remove(f)

class TestFastReactorConstructors(TestCase):
    """Tests that the storage component constructors work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_FRDefaults(self):
        frd = bright.FRDefaults()
        assert_equal(frd.batches, 3)
        assert_equal(frd.flux, 2.0*(10.0**15))
        assert_equal(frd.FuelForm["IHM"], 1.0)
        assert_equal(frd.CoolantForm["NA23"], 1.0)
        assert_equal(frd.FuelDensity, 18.0)
        assert_equal(frd.CoolantDensity, 0.927)
        assert_equal(frd.pnl, 0.65)
        assert_equal(frd.BUt, 0.0)
        assert_false(frd.useDisadvantage)
        assert_equal(frd.LatticeType, "Cylindrical")
        assert_false(frd.HydrogenRescale)
        assert_equal(frd.Radius, 0.3115)
        assert_equal(frd.Length, 0.956)
        assert_equal(frd.open_slots, 19.0)
        assert_equal(frd.total_slots, 163.0)

    def test_FastReactor1G_1(self):
        fr = FastReactor1G()
        assert_equal(fr.name, '')
        assert_equal(fr.params2track, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
        assert_equal(fr.B, 3)
        assert_equal(fr.phi, 2.0*(10.0**15))
        assert_equal(fr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(fr.CoolantChemicalForm["NA23"], 1.0)
        assert_equal(fr.rhoF, 18.0)
        assert_equal(fr.rhoC, 0.927)
        assert_equal(fr.P_NL, 0.65)
        assert_equal(fr.TargetBU, 0.0)
        assert_false(fr.useZeta)
        assert_equal(fr.Lattice, "Cylindrical")
        assert_false(fr.H_XS_Rescale)
        assert_equal(fr.r, 0.3115)
        assert_equal(fr.l, 0.956)
        assert_equal(fr.S_O, 19.0)
        assert_equal(fr.S_T, 163.0)

    def test_FastReactor1G_2(self):
        lf = os.getenv("BRIGHT_DATA") + "/FR.h5"
        fr = FastReactor1G(libfile=lf)
        assert_equal(fr.libfile, lf)
        assert_equal(fr.name, '')
        assert_equal(fr.params2track, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
        assert_equal(fr.B, 3)
        assert_equal(fr.phi, 2.0*(10.0**15))
        assert_equal(fr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(fr.CoolantChemicalForm["NA23"], 1.0)
        assert_equal(fr.rhoF, 18.0)
        assert_equal(fr.rhoC, 0.927)
        assert_equal(fr.P_NL, 0.65)
        assert_equal(fr.TargetBU, 0.0)
        assert_false(fr.useZeta)
        assert_equal(fr.Lattice, "Cylindrical")
        assert_false(fr.H_XS_Rescale)
        assert_equal(fr.r, 0.3115)
        assert_equal(fr.l, 0.956)
        assert_equal(fr.S_O, 19.0)
        assert_equal(fr.S_T, 163.0)

    def test_FastReactor1G_3(self):
        lf = os.getenv("BRIGHT_DATA") + "/FR.h5"
        fr = FastReactor1G(libfile=lf, name="fr")
        assert_equal(fr.libfile, lf)
        assert_equal(fr.name, 'fr')
        assert_equal(fr.params2track, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
        assert_equal(fr.B, 3)
        assert_equal(fr.phi, 2.0*(10.0**15))
        assert_equal(fr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(fr.CoolantChemicalForm["NA23"], 1.0)
        assert_equal(fr.rhoF, 18.0)
        assert_equal(fr.rhoC, 0.927)
        assert_equal(fr.P_NL, 0.65)
        assert_equal(fr.TargetBU, 0.0)
        assert_false(fr.useZeta)
        assert_equal(fr.Lattice, "Cylindrical")
        assert_false(fr.H_XS_Rescale)
        assert_equal(fr.r, 0.3115)
        assert_equal(fr.l, 0.956)
        assert_equal(fr.S_O, 19.0)
        assert_equal(fr.S_T, 163.0)

    def test_FastReactor1G_4(self):
        rp = bright.FRDefaults()
        rp.BUt = 140.0
        fr = FastReactor1G(reactor_parameters=rp)
        assert_equal(fr.name, '')
        assert_equal(fr.params2track, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
        assert_equal(fr.B, 3)
        assert_equal(fr.phi, 2.0*(10.0**15))
        assert_equal(fr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(fr.CoolantChemicalForm["NA23"], 1.0)
        assert_equal(fr.rhoF, 18.0)
        assert_equal(fr.rhoC, 0.927)
        assert_equal(fr.P_NL, 0.65)
        assert_equal(fr.TargetBU, 140.0)
        assert_false(fr.useZeta)
        assert_equal(fr.Lattice, "Cylindrical")
        assert_false(fr.H_XS_Rescale)
        assert_equal(fr.r, 0.3115)
        assert_equal(fr.l, 0.956)
        assert_equal(fr.S_O, 19.0)
        assert_equal(fr.S_T, 163.0)

    def test_FastReactor1G_5(self):
        rp = bright.FRDefaults()
        rp.BUt = 140.0
        fr = FastReactor1G(reactor_parameters=rp, name='fr')
        assert_equal(fr.name, 'fr')
        assert_equal(fr.params2track, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
        assert_equal(fr.B, 3)
        assert_equal(fr.phi, 2.0*(10.0**15))
        assert_equal(fr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(fr.CoolantChemicalForm["NA23"], 1.0)
        assert_equal(fr.rhoF, 18.0)
        assert_equal(fr.rhoC, 0.927)
        assert_equal(fr.P_NL, 0.65)
        assert_equal(fr.TargetBU, 140.0)
        assert_false(fr.useZeta)
        assert_equal(fr.Lattice, "Cylindrical")
        assert_false(fr.H_XS_Rescale)
        assert_equal(fr.r, 0.3115)
        assert_equal(fr.l, 0.956)
        assert_equal(fr.S_O, 19.0)
        assert_equal(fr.S_T, 163.0)

    def test_FastReactor1G_6(self):
        lf = os.getenv("BRIGHT_DATA") + "/FR.h5"
        rp = bright.FRDefaults()
        rp.BUt = 140.0
        fr = FastReactor1G(libfile=lf, reactor_parameters=rp)
        assert_equal(fr.libfile, lf)
        assert_equal(fr.name, '')
        assert_equal(fr.params2track, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
        assert_equal(fr.B, 3)
        assert_equal(fr.phi, 2.0*(10.0**15))
        assert_equal(fr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(fr.CoolantChemicalForm["NA23"], 1.0)
        assert_equal(fr.rhoF, 18.0)
        assert_equal(fr.rhoC, 0.927)
        assert_equal(fr.P_NL, 0.65)
        assert_equal(fr.TargetBU, 140.0)
        assert_false(fr.useZeta)
        assert_equal(fr.Lattice, "Cylindrical")
        assert_false(fr.H_XS_Rescale)
        assert_equal(fr.r, 0.3115)
        assert_equal(fr.l, 0.956)
        assert_equal(fr.S_O, 19.0)
        assert_equal(fr.S_T, 163.0)

    def test_FastReactor1G_7(self):
        lf = os.getenv("BRIGHT_DATA") + "/FR.h5"
        rp = bright.FRDefaults()
        rp.BUt = 140.0
        fr = FastReactor1G(libfile=lf, reactor_parameters=rp, name='fr')
        assert_equal(fr.libfile, lf)
        assert_equal(fr.name, 'fr')
        assert_equal(fr.params2track, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
        assert_equal(fr.B, 3)
        assert_equal(fr.phi, 2.0*(10.0**15))
        assert_equal(fr.FuelChemicalForm["IHM"], 1.0)
        assert_equal(fr.CoolantChemicalForm["NA23"], 1.0)
        assert_equal(fr.rhoF, 18.0)
        assert_equal(fr.rhoC, 0.927)
        assert_equal(fr.P_NL, 0.65)
        assert_equal(fr.TargetBU, 140.0)
        assert_false(fr.useZeta)
        assert_equal(fr.Lattice, "Cylindrical")
        assert_false(fr.H_XS_Rescale)
        assert_equal(fr.r, 0.3115)
        assert_equal(fr.l, 0.956)
        assert_equal(fr.S_O, 19.0)
        assert_equal(fr.S_T, 163.0)


class TestFastReactor1GAttributes(TestCase):
    """Tests that the fuel cycle component attributes work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_params2track(self):
        fr = FastReactor1G()
        assert_equal(fr.params2track, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
        fr.params2track = set(["Mass"])
        assert_equal(fr.params2track, set(["Mass"]))

class TestFastReactor1GMethods(TestCase):
    """Tests that the fuel cycle component methods work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_setParams(self):
        lf = os.getenv("BRIGHT_DATA") + "/FR.h5"
        bright.load_track_isos_hdf5(lf)
        rp = bright.FRDefaults()
        rp.BUt = 140.0
        fr = FastReactor1G(libfile=lf, reactor_parameters=rp, name='fr')
        fr.doCalc(MassStream({922350: 0.30, 922380: 0.70}))
        fr.setParams()
        assert_equal(fr.params_prior_calc["BUd"],  0.0)
        assert_equal(fr.params_after_calc["BUd"], fr.BUd)
        assert_equal(fr.params_prior_calc["TRUCR"],  0.0)
        assert_equal(fr.params_after_calc["TRUCR"], fr.TruCR)
        assert_equal(fr.params_prior_calc["P_NL"],  0.0)
        assert_equal(fr.params_after_calc["P_NL"], fr.P_NL)
        assert_equal(fr.params_prior_calc["U"],  fr.InU.mass)
        assert_equal(fr.params_after_calc["U"], fr.OutU.mass)
        assert_equal(fr.params_prior_calc["TRU"],  fr.InTRU.mass)
        assert_equal(fr.params_after_calc["TRU"], fr.OutTRU.mass)
        assert_equal(fr.params_prior_calc["ACT"],  fr.InACT.mass)
        assert_equal(fr.params_after_calc["ACT"], fr.OutACT.mass)
        assert_equal(fr.params_prior_calc["LAN"],  fr.InLAN.mass)
        assert_equal(fr.params_after_calc["LAN"], fr.OutLAN.mass)
        assert_equal(fr.params_prior_calc["FP"],  1.0 - fr.InACT.mass - fr.InLAN.mass)
        assert_equal(fr.params_after_calc["FP"], 1.0 - fr.OutACT.mass - fr.OutLAN.mass)
        

# Put Integral tests here, if desired.


if __name__ == "__main__":
    nose.main()
