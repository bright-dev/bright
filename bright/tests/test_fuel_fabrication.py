"""Fuel Fabrication Tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
assert_almost_equal, assert_true, assert_false, set_trace, with_setup

import os
import warnings
import tables as tb
import numpy as np


from bright import bright_conf, load_track_nucs_hdf5
from bright.reactor_parameters import ReactorParameters, lwr_defaults
from bright.reactor1g import Reactor1G
from bright.fuel_fabrication import FuelFabrication
from pyne.material import Material
from pyne import nucname

default_rp = ReactorParameters()
default_rp.batches = 3
default_rp.flux = 2*(10**14)
default_rp.fuel_form = {"IHM": 1.0, "O16": 2.0}
default_rp.coolant_form = {"H1": 2.0, "O16": 1.0}
default_rp.fuel_density = 10.7
default_rp.coolant_density = 0.73
default_rp.pnl = 0.98
default_rp.BUt = 50.0
default_rp.use_disadvantage_factor = True
default_rp.lattice_type = 'Cylindrical'
default_rp.rescale_hydrogen = True
default_rp.fuel_radius = 0.411
default_rp.unit_cell_pitch = 0.7
default_rp.open_slots = 123
default_rp.total_slots = 180

#
# Fixtures
#
mws = None
ff = None
u235 = None
u238 = None
r1g = None
mats = None

def setup_ff():
    global mws, ff, u235, u238, r1g, mats
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    load_track_nucs_hdf5(libfile)

    r1g = Reactor1G(rp=default_rp)
    r1g.loadlib(libfile)

    u235 = Material({922350: 1.0}, 1.0, name="U-235")
    u238 = Material({922380: 1.0}, 1.0, name="U-238")
    mats = {"U235": u235, "U238": u238}

    mws = {"U235": -1.0, "U238": -1.0}

    ff = FuelFabrication(mats=mats, mws_in=mws, r=r1g)


def teardown_ff():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "ff.h5"]:
            os.remove(f)

def teardown_ff_clear():
    global mws, ff, u235, u238, r1g, mats
    teardown_ff()
    mws = None
    ff = None
    u235 = None
    u238 = None
    r1g = None
    mats = None



#
# Tests that the FuelFabrication component constructors work.
#

@with_setup(None, teardown_ff)
def test_FuelFabrication_1():
    ff = FuelFabrication()
    assert_equal(ff.name, '')
    assert_equal(ff.track_params, set())

@with_setup(None, teardown_ff)
def test_FuelFabrication_2():
    ff = FuelFabrication(n="ff")
    assert_equal(ff.name, 'ff')
    assert_equal(ff.track_params, set())

@with_setup(None, teardown_ff)
def test_FuelFabrication_3():
    ff = FuelFabrication(paramtrack=set(["Mass"]))
    assert_equal(ff.name, '')
    assert_equal(ff.track_params, set(["Mass"]))

@with_setup(None, teardown_ff)
def test_FuelFabrication_4():
    ff = FuelFabrication(paramtrack=set(["Mass"]), n='ff')
    assert_equal(ff.name, 'ff')
    assert_equal(ff.track_params, set(["Mass"]))

@with_setup(None, teardown_ff)
def test_FuelFabrication_5():
    # Reactor to use
    rp = ReactorParameters()
    r1g = Reactor1G(rp=rp, n="r1g")

    # Mass streams to use
    u235 = Material({922350: 1.0}, 1.0, name="U-235")
    u238 = Material({922380: 1.0}, 1.0, name="U-238")
    mats = {"U235": u235, "U238": u238}

    # Mass weights to use
    mws = {"U235": -1.0, "U238": -1.0}

    # Fuel Fabrication Facility
    ff = FuelFabrication(mats=mats, mws_in=mws, r=r1g)

    keys = ["U235", "U238"]
    print ff.materials
    assert_equal(set(ff.materials.keys()), set(keys))

    for iso in keys:
        assert_equal(ff.materials[iso].mass, 1.0)
        assert_equal(ff.materials[iso].comp[nucname.zzaaam(iso)], 1.0)

    assert_equal(ff.mass_weights_in, mws)

    assert_equal(ff.track_params, set(["Weight_U235", "deltaR_U235", "Weight_U238", "deltaR_U238"]))

    assert_equal(ff.reactor.name, "r1g")
    r1g.name = "r1g name"
    ff.initialize(mats, mws, r1g)
    assert_equal(ff.reactor.name, "r1g name")

@with_setup(None, teardown_ff)
def test_FuelFabrication_6():
    # Reactor to use
    rp = ReactorParameters()
    r1g = Reactor1G(rp=rp, n="r1g")

    # Mass streams to use
    u235 = Material({922350: 1.0}, 1.0, name="U-235")
    u238 = Material({922380: 1.0}, 1.0, name="U-238")
    mats = {"U235": u235, "U238": u238}

    # Mass weights to use
    mws = {"U235": -1.0, "U238": -1.0}

    # Fuel Fabrication Facility
    ff = FuelFabrication(mats=mats, mws_in=mws, r=r1g, paramtrack=set(["Mass"]))

    keys = ["U235", "U238"]
    assert_equal(set(ff.materials.keys()), set(keys))

    for iso in keys:
        assert_equal(ff.materials[iso].mass, 1.0)
        assert_equal(ff.materials[iso].comp[nucname.zzaaam(iso)], 1.0)

    assert_equal(ff.mass_weights_in, mws)

    assert_equal(ff.track_params, set(["Mass", "Weight_U235", "deltaR_U235", "Weight_U238", "deltaR_U238"]))

    assert_equal(ff.reactor.name, "r1g")
    r1g.name = "r1g name"
    ff.initialize(mats, mws, r1g)
    assert_equal(ff.reactor.name, "r1g name")

@with_setup(None, teardown_ff)
def test_FuelFabrication_7():
    # Reactor to use
    rp = ReactorParameters()
    r1g = Reactor1G(rp=rp, n="r1g")

    # Mass streams to use
    u235 = Material({922350: 1.0}, 1.0, name="U-235")
    u238 = Material({922380: 1.0}, 1.0, name="U-238")
    mats = {"U235": u235, "U238": u238}

    # Mass weights to use
    mws = {"U235": -1.0, "U238": -1.0}

    # Fuel Fabrication Facility
    ff = FuelFabrication(mats=mats, mws_in=mws, r=r1g, paramtrack=set(["Mass"]), n="ff")

    keys = ["U235", "U238"]
    assert_equal(set(ff.materials.keys()), set(keys))

    for iso in keys:
        assert_equal(ff.materials[iso].mass, 1.0)
        assert_equal(ff.materials[iso].comp[nucname.zzaaam(iso)], 1.0)

    assert_equal(ff.mass_weights_in, mws)

    assert_equal(ff.track_params, set(["Mass", "Weight_U235", "deltaR_U235", "Weight_U238", "deltaR_U238"]))

    assert_equal(ff.name, "ff")

    assert_equal(ff.reactor.name, "r1g")
    r1g.name = "r1g name"
    ff.initialize(mats, mws, r1g)
    assert_equal(ff.reactor.name, "r1g name")


#
# Tests that the FuelFabrication basic data attributes work.
#

@with_setup(setup_ff, teardown_ff)
def test_mass_weights_in():
    assert_equal(ff.mass_weights_in, mws)
    ff.mass_weights_in = {"U235": -1.0, "U238": -1.0, "O16": 0.125}
    assert_equal(ff.mass_weights_in, {"U235": -1.0, "U238": -1.0, "O16": 0.125})


@with_setup(None, teardown_ff_clear)
def test_materials():
    keys = ["U235", "U238"]
    assert_equal(set(ff.materials.keys()), set(keys))

    for iso in keys:
        assert_equal(ff.materials[iso].mass, 1.0)
        assert_equal(ff.materials[iso].comp[nucname.zzaaam(iso)], 1.0)

    u235 = Material({922350: 1.0}, 1.0, name="U-235")
    u238 = Material({922380: 1.0}, 1.0, name="U-238")
    o16  = Material({80160:  1.0}, 1.0, name="O-16")
    mats = {"U235": u235, "U238": u238, "O16": o16}
    ff.materials = mats

    keys = ["U235", "U238", "O16"]
    assert_equal(set(ff.materials.keys()), set(keys))

    for iso in keys:
        assert_equal(ff.materials[iso].mass, 1.0)
        assert_equal(ff.materials[iso].comp[nucname.zzaaam(iso)], 1.0)


#
# Tests that the FuelFabrication methods work.
#


@with_setup(setup_ff, teardown_ff)
def test_calc_deltaRs():
    ff.calc_deltaRs()

    keys = ["U235", "U238"]
    assert_equal(set(ff.deltaRs.keys()), set(keys))

    assert(ff.deltaRs["U238"] <= ff.deltaRs["U235"])


@with_setup(None, teardown_ff)
def test_calc_mass_ratios():
    ff.calc_mass_ratios()

    keys = ["U235", "U238"]
    assert_equal(set(ff.mass_weights_out.keys()), set(keys))

    assert(ff.mass_weights_out["U235"] <= ff.mass_weights_out["U238"])

@with_setup(None, teardown_ff)
def test_calc_mass_ratios():
    ff.calc_mass_ratios()
    core_input = ff.calc_core_input()

    assert_equal(core_input.mass, 1.0)

    assert_almost_equal(core_input.comp[922350], ff.mass_weights_out["U235"])
    assert_almost_equal(core_input.comp[922380], ff.mass_weights_out["U238"])

@with_setup(None, teardown_ff)
def test_calc1():
    core_input = ff.calc()

    assert_equal(core_input.mass, 1.0)

    assert_almost_equal(core_input.comp[922350], ff.mass_weights_out["U235"])
    assert_almost_equal(core_input.comp[922380], ff.mass_weights_out["U238"])

@with_setup(None, teardown_ff)
def test_calc2():
    r1g.name = "r1g name"

    core_input = ff.calc(mats, mws, r1g)

    assert_equal(ff.reactor.name, "r1g name")

    assert_equal(core_input.mass, 1.0)

    assert_almost_equal(core_input.comp[922350], ff.mass_weights_out["U235"])
    assert_almost_equal(core_input.comp[922380], ff.mass_weights_out["U238"])

@with_setup(None, teardown_ff_clear)
def test_calc_params():
    core_input = ff.calc()
    ff.calc_params()

    assert_equal(ff.params_prior_calc["Weight_U235"], -1.0)
    assert_equal(ff.params_prior_calc["Weight_U238"], -1.0)

    assert_equal(ff.params_after_calc["Weight_U235"], ff.mass_weights_out["U235"])
    assert_equal(ff.params_after_calc["Weight_U238"], ff.mass_weights_out["U238"])

    assert_equal(ff.params_prior_calc["deltaR_U235"], ff.deltaRs["U235"])
    assert_equal(ff.params_prior_calc["deltaR_U238"], ff.deltaRs["U238"])

    assert_equal(ff.params_after_calc["deltaR_U235"], ff.deltaRs["U235"])
    assert_equal(ff.params_after_calc["deltaR_U238"], ff.deltaRs["U238"])


if __name__ == "__main__":
    nose.main()
