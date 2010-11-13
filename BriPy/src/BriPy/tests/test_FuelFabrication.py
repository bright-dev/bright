"""Fuel Fabrication Component and Helper Class tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, set_trace

import os
import warnings
import tables as tb
import numpy as np

import BriPy


Reactor1G = BriPy.Reactor1G
FluencePoint = BriPy.FluencePoint
ReactorParameters = BriPy.ReactorParameters
MassStream = BriPy.MassStream
FuelFabrication = BriPy.FuelFabrication

default_rp = BriPy.ReactorParameters()
default_rp.batches = 3
default_rp.flux = 2*(10**14)
default_rp.FuelForm = {"IHM": 1.0, "O16": 2.0}
default_rp.CoolantForm = {"H1": 2.0, "O16": 1.0}
default_rp.FuelDensity = 10.7
default_rp.CoolantDensity = 0.73
default_rp.pnl = 0.98
default_rp.BUt = 50.0
default_rp.useDisadvantage = True
default_rp.LatticeType = 'Cylindrical'
default_rp.HydrogenRescale = True
default_rp.Radius = 0.411
default_rp.Length = 0.7
default_rp.open_slots = 123
default_rp.total_slots = 180

def general_teardown():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "ff.h5"]:
            os.remove(f)


class TestFuelFabricationConstructors(TestCase):
    """Tests that the FuelFabrication component constructors work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_FuelFabrication_1(self):
        ff = FuelFabrication()
        assert_equal(ff.name, '')
        assert_equal(ff.params2track, [])

    def test_FuelFabrication_2(self):
        ff = FuelFabrication("ff")
        assert_equal(ff.name, 'ff')
        assert_equal(ff.params2track, [])

    def test_FuelFabrication_3(self):
        ff = FuelFabrication(["Mass"])
        assert_equal(ff.name, '')
        assert_equal(ff.params2track, ["Mass"])

    def test_FuelFabrication_4(self):
        ff = FuelFabrication(["Mass"], 'ff')
        assert_equal(ff.name, 'ff')
        assert_equal(ff.params2track, ["Mass"])

    def test_FuelFabrication_5(self):
        # Reactor to use
        rp = ReactorParameters()
        r1g = Reactor1G(rp, "r1g")

        # Mass streams to use
        u235 = MassStream({922350: 1.0}, 1.0, "U-235")
        u238 = MassStream({922380: 1.0}, 1.0, "U-238")
        mss = {"U235": u235, "U238": u238}

        # Mass weights to use
        mws = {"U235": -1.0, "U238": -1.0}

        # Fuel Fabrication Facility
        ff = FuelFabrication(mss, mws, r1g)

        keys = ["U235", "U238"]
        assert_equal(set(ff.mass_streams.keys()), set(keys))

        for iso in keys:
            assert_equal(ff.mass_streams[iso].mass, 1.0)
            assert_equal(ff.mass_streams[iso].comp[BriPy.LLAAAM_2_zzaaam(iso)], 1.0)

        assert_equal(ff.mass_weights_in, mws)

        assert_equal(set(ff.params2track), set(["Weight_U235", "deltaR_U235", "Weight_U238", "deltaR_U238"]))

    def test_FuelFabrication_6(self):
        # Reactor to use
        rp = ReactorParameters()
        r1g = Reactor1G(rp, "r1g")

        # Mass streams to use
        u235 = MassStream({922350: 1.0}, 1.0, "U-235")
        u238 = MassStream({922380: 1.0}, 1.0, "U-238")
        mss = {"U235": u235, "U238": u238}

        # Mass weights to use
        mws = {"U235": -1.0, "U238": -1.0}

        # Fuel Fabrication Facility
        ff = FuelFabrication(mss, mws, r1g, ["Mass"])

        keys = ["U235", "U238"]
        assert_equal(set(ff.mass_streams.keys()), set(keys))

        for iso in keys:
            assert_equal(ff.mass_streams[iso].mass, 1.0)
            assert_equal(ff.mass_streams[iso].comp[BriPy.LLAAAM_2_zzaaam(iso)], 1.0)

        assert_equal(ff.mass_weights_in, mws)

        assert_equal(set(ff.params2track), set(["Mass", "Weight_U235", "deltaR_U235", "Weight_U238", "deltaR_U238"]))

    def test_FuelFabrication_7(self):
        # Reactor to use
        rp = ReactorParameters()
        r1g = Reactor1G(rp, "r1g")

        # Mass streams to use
        u235 = MassStream({922350: 1.0}, 1.0, "U-235")
        u238 = MassStream({922380: 1.0}, 1.0, "U-238")
        mss = {"U235": u235, "U238": u238}

        # Mass weights to use
        mws = {"U235": -1.0, "U238": -1.0}

        # Fuel Fabrication Facility
        ff = FuelFabrication(mss, mws, r1g, ["Mass"], "ff")

        keys = ["U235", "U238"]
        assert_equal(set(ff.mass_streams.keys()), set(keys))

        for iso in keys:
            assert_equal(ff.mass_streams[iso].mass, 1.0)
            assert_equal(ff.mass_streams[iso].comp[BriPy.LLAAAM_2_zzaaam(iso)], 1.0)

        assert_equal(ff.mass_weights_in, mws)

        assert_equal(set(ff.params2track), set(["Mass", "Weight_U235", "deltaR_U235", "Weight_U238", "deltaR_U238"]))

        assert_equal(ff.name, "ff")

class TestFuelFabricationAttributes(TestCase):
    """Tests that the FuelFabrication basic data attributes work."""

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        BriPy.load_isos2track_hdf5(libfile)

        r1g = Reactor1G(default_rp)
        r1g.loadLib(libfile)
        cls.r1g = r1g

        u235 = MassStream({922350: 1.0}, 1.0, "U-235")
        u238 = MassStream({922380: 1.0}, 1.0, "U-238")
        mss = {"U235": u235, "U238": u238}
        cls.mss = mss

        mws = {"U235": -1.0, "U238": -1.0}
        cls.mws = mws

        cls.ff = FuelFabrication(mss, mws, r1g)

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_mass_weights_in(self):
        assert_equal(self.ff.mass_weights_in, self.mws)
        self.ff.mass_weights_in = {"U235": -1.0, "U238": -1.0, "O16": 0.125}
        assert_equal(self.ff.mass_weights_in, {"U235": -1.0, "U238": -1.0, "O16": 0.125})

    def test_mass_streams(self):
        keys = ["U235", "U238"]
        assert_equal(set(self.ff.mass_streams.keys()), set(keys))

        for iso in keys:
            assert_equal(self.ff.mass_streams[iso].mass, 1.0)
            assert_equal(self.ff.mass_streams[iso].comp[BriPy.LLAAAM_2_zzaaam(iso)], 1.0)

        u235 = MassStream({922350: 1.0}, 1.0, "U-235")
        u238 = MassStream({922380: 1.0}, 1.0, "U-238")
        o16  = MassStream({80160:  1.0}, 1.0, "O-16")
        mss = {"U235": u235, "U238": u238, "O16": o16}
        self.ff.mass_streams = mss

        keys = ["U235", "U238", "O16"]
        assert_equal(set(self.ff.mass_streams.keys()), set(keys))

        for iso in keys:
            assert_equal(self.ff.mass_streams[iso].mass, 1.0)
            assert_equal(self.ff.mass_streams[iso].comp[BriPy.LLAAAM_2_zzaaam(iso)], 1.0)



class TestFuelFabricationMethodss(TestCase):
    """Tests that the FuelFabrication methods work."""

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        BriPy.load_isos2track_hdf5(libfile)

        r1g = Reactor1G(default_rp)
        r1g.loadLib(libfile)
        cls.r1g = r1g

        u235 = MassStream({922350: 1.0}, 1.0, "U-235")
        u238 = MassStream({922380: 1.0}, 1.0, "U-238")
        mss = {"U235": u235, "U238": u238}
        cls.mss = mss

        mws = {"U235": -1.0, "U238": -1.0}
        cls.mws = mws

        cls.ff = FuelFabrication(mss, mws, cls.r1g)

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_calc_deltaRs(self):
        #self.ff.calc_deltaRs()

        keys = ["U235", "U238"]
        assert_equal(set(self.ff.deltaRs.keys()), set(keys))

        assert(self.ff.deltaRs["U238"] <= self.ff.deltaRs["U235"])

        print self.ff.deltaRs
        

if __name__ == "__main__":
    nose.main()
