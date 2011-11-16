"""FastReactor1G Tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
assert_almost_equal, assert_true, assert_false, with_setup

import os
import warnings
import tables as tb
import numpy as np

from bright import bright_conf, load_track_nucs_hdf5
from bright.reactor_parameters import fr_defaults
from bright.fast_reactor1g import FastReactor1G
from pyne.material import Material

#
# Fixtures
#

def teardown_fr1g():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "fr.h5", "fuel_cycle.h5"]:
            os.remove(f)

#
# Tests that the storage component constructors work.
#

@with_setup(None, teardown_fr1g)
def test_fr_defaults():
    frd = fr_defaults()
    assert_equal(frd.batches, 3)
    assert_equal(frd.flux, 2.0*(10.0**15))
    assert_equal(frd.fuel_form["IHM"], 1.0)
    assert_equal(frd.coolant_form["NA23"], 1.0)
    assert_equal(frd.fuel_density, 18.0)
    assert_equal(frd.coolant_density, 0.927)
    assert_equal(frd.pnl, 0.65)
    assert_equal(frd.BUt, 0.0)
    assert_false(frd.use_disadvantage_factor)
    assert_equal(frd.lattice_type, "Cylindrical")
    assert_false(frd.rescale_hydrogen)
    assert_equal(frd.fuel_radius, 0.3115)
    assert_equal(frd.unit_cell_pitch, 0.956)
    assert_equal(frd.open_slots, 19.0)
    assert_equal(frd.total_slots, 163.0)

@with_setup(None, teardown_fr1g)
def test_FastReactor1G_1():
    fr = FastReactor1G()
    assert_equal(fr.name, '')
    assert_equal(fr.track_params, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
    assert_equal(fr.B, 3)
    assert_equal(fr.phi, 2.0*(10.0**15))
    assert_equal(fr.fuel_chemical_form["IHM"], 1.0)
    assert_equal(fr.coolant_chemical_form["NA23"], 1.0)
    assert_equal(fr.rhoF, 18.0)
    assert_equal(fr.rhoC, 0.927)
    assert_equal(fr.P_NL, 0.65)
    assert_equal(fr.target_BU, 0.0)
    assert_false(fr.use_zeta)
    assert_equal(fr.lattice_flag, "Cylindrical")
    assert_false(fr.rescale_hydrogen_xs)
    assert_equal(fr.r, 0.3115)
    assert_equal(fr.l, 0.956)
    assert_equal(fr.S_O, 19.0)
    assert_equal(fr.S_T, 163.0)

@with_setup(None, teardown_fr1g)
def test_FastReactor1G_2():
    lf = os.getenv("BRIGHT_DATA") + "/FR.h5"
    fr = FastReactor1G(libfile=lf)
    assert_equal(fr.libfile, lf)
    assert_equal(fr.name, '')
    assert_equal(fr.track_params, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
    assert_equal(fr.B, 3)
    assert_equal(fr.phi, 2.0*(10.0**15))
    assert_equal(fr.fuel_chemical_form["IHM"], 1.0)
    assert_equal(fr.coolant_chemical_form["NA23"], 1.0)
    assert_equal(fr.rhoF, 18.0)
    assert_equal(fr.rhoC, 0.927)
    assert_equal(fr.P_NL, 0.65)
    assert_equal(fr.target_BU, 0.0)
    assert_false(fr.use_zeta)
    assert_equal(fr.lattice_flag, "Cylindrical")
    assert_false(fr.rescale_hydrogen_xs)
    assert_equal(fr.r, 0.3115)
    assert_equal(fr.l, 0.956)
    assert_equal(fr.S_O, 19.0)
    assert_equal(fr.S_T, 163.0)

@with_setup(None, teardown_fr1g)
def test_FastReactor1G_3():
    lf = os.getenv("BRIGHT_DATA") + "/FR.h5"
    fr = FastReactor1G(libfile=lf, name="fr")
    assert_equal(fr.libfile, lf)
    assert_equal(fr.name, 'fr')
    assert_equal(fr.track_params, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
    assert_equal(fr.B, 3)
    assert_equal(fr.phi, 2.0*(10.0**15))
    assert_equal(fr.fuel_chemical_form["IHM"], 1.0)
    assert_equal(fr.coolant_chemical_form["NA23"], 1.0)
    assert_equal(fr.rhoF, 18.0)
    assert_equal(fr.rhoC, 0.927)
    assert_equal(fr.P_NL, 0.65)
    assert_equal(fr.target_BU, 0.0)
    assert_false(fr.use_zeta)
    assert_equal(fr.lattice_flag, "Cylindrical")
    assert_false(fr.rescale_hydrogen_xs)
    assert_equal(fr.r, 0.3115)
    assert_equal(fr.l, 0.956)
    assert_equal(fr.S_O, 19.0)
    assert_equal(fr.S_T, 163.0)

@with_setup(None, teardown_fr1g)
def test_FastReactor1G_4():
    rp = fr_defaults()
    rp.BUt = 140.0
    fr = FastReactor1G(reactor_parameters=rp)
    assert_equal(fr.name, '')
    assert_equal(fr.track_params, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
    assert_equal(fr.B, 3)
    assert_equal(fr.phi, 2.0*(10.0**15))
    assert_equal(fr.fuel_chemical_form["IHM"], 1.0)
    assert_equal(fr.coolant_chemical_form["NA23"], 1.0)
    assert_equal(fr.rhoF, 18.0)
    assert_equal(fr.rhoC, 0.927)
    assert_equal(fr.P_NL, 0.65)
    assert_equal(fr.target_BU, 140.0)
    assert_false(fr.use_zeta)
    assert_equal(fr.lattice_flag, "Cylindrical")
    assert_false(fr.rescale_hydrogen_xs)
    assert_equal(fr.r, 0.3115)
    assert_equal(fr.l, 0.956)
    assert_equal(fr.S_O, 19.0)
    assert_equal(fr.S_T, 163.0)

@with_setup(None, teardown_fr1g)
def test_FastReactor1G_5():
    rp = fr_defaults()
    rp.BUt = 140.0
    fr = FastReactor1G(reactor_parameters=rp, name='fr')
    assert_equal(fr.name, 'fr')
    assert_equal(fr.track_params, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
    assert_equal(fr.B, 3)
    assert_equal(fr.phi, 2.0*(10.0**15))
    assert_equal(fr.fuel_chemical_form["IHM"], 1.0)
    assert_equal(fr.coolant_chemical_form["NA23"], 1.0)
    assert_equal(fr.rhoF, 18.0)
    assert_equal(fr.rhoC, 0.927)
    assert_equal(fr.P_NL, 0.65)
    assert_equal(fr.target_BU, 140.0)
    assert_false(fr.use_zeta)
    assert_equal(fr.lattice_flag, "Cylindrical")
    assert_false(fr.rescale_hydrogen_xs)
    assert_equal(fr.r, 0.3115)
    assert_equal(fr.l, 0.956)
    assert_equal(fr.S_O, 19.0)
    assert_equal(fr.S_T, 163.0)

@with_setup(None, teardown_fr1g)
def test_FastReactor1G_6():
    lf = os.getenv("BRIGHT_DATA") + "/FR.h5"
    rp = fr_defaults()
    rp.BUt = 140.0
    fr = FastReactor1G(libfile=lf, reactor_parameters=rp)
    assert_equal(fr.libfile, lf)
    assert_equal(fr.name, '')
    assert_equal(fr.track_params, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
    assert_equal(fr.B, 3)
    assert_equal(fr.phi, 2.0*(10.0**15))
    assert_equal(fr.fuel_chemical_form["IHM"], 1.0)
    assert_equal(fr.coolant_chemical_form["NA23"], 1.0)
    assert_equal(fr.rhoF, 18.0)
    assert_equal(fr.rhoC, 0.927)
    assert_equal(fr.P_NL, 0.65)
    assert_equal(fr.target_BU, 140.0)
    assert_false(fr.use_zeta)
    assert_equal(fr.lattice_flag, "Cylindrical")
    assert_false(fr.rescale_hydrogen_xs)
    assert_equal(fr.r, 0.3115)
    assert_equal(fr.l, 0.956)
    assert_equal(fr.S_O, 19.0)
    assert_equal(fr.S_T, 163.0)

@with_setup(None, teardown_fr1g)
def test_FastReactor1G_7():
    lf = os.getenv("BRIGHT_DATA") + "/FR.h5"
    rp = fr_defaults()
    rp.BUt = 140.0
    fr = FastReactor1G(libfile=lf, reactor_parameters=rp, name='fr')
    assert_equal(fr.libfile, lf)
    assert_equal(fr.name, 'fr')
    assert_equal(fr.track_params, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
    assert_equal(fr.B, 3)
    assert_equal(fr.phi, 2.0*(10.0**15))
    assert_equal(fr.fuel_chemical_form["IHM"], 1.0)
    assert_equal(fr.coolant_chemical_form["NA23"], 1.0)
    assert_equal(fr.rhoF, 18.0)
    assert_equal(fr.rhoC, 0.927)
    assert_equal(fr.P_NL, 0.65)
    assert_equal(fr.target_BU, 140.0)
    assert_false(fr.use_zeta)
    assert_equal(fr.lattice_flag, "Cylindrical")
    assert_false(fr.rescale_hydrogen_xs)
    assert_equal(fr.r, 0.3115)
    assert_equal(fr.l, 0.956)
    assert_equal(fr.S_O, 19.0)
    assert_equal(fr.S_T, 163.0)

#
# Tests that the fuel cycle component attributes work.
#

@with_setup(None, teardown_fr1g)
def test_track_params():
    fr = FastReactor1G()
    assert_equal(fr.track_params, set(["ACT", "BUd", "FP", "LAN", "P_NL", "TRU", "TRUCR", "U"]))
    fr.track_params = set(["Mass"])
    assert_equal(fr.track_params, set(["Mass"]))

#
# Tests that the fuel cycle component methods work.
#

@with_setup(None, teardown_fr1g)
def test_calc_params():
    lf = os.getenv("BRIGHT_DATA") + "/FR.h5"
    load_track_nucs_hdf5(lf)
    rp = fr_defaults()
    rp.BUt = 140.0
    fr = FastReactor1G(libfile=lf, reactor_parameters=rp, name='fr')
    fr.calc(Material({922350: 0.30, 922380: 0.70}))
    fr.calc_params()
    assert_equal(fr.params_prior_calc["BUd"],  0.0)
    assert_equal(fr.params_after_calc["BUd"], fr.BUd)
    assert_equal(fr.params_prior_calc["TRUCR"],  0.0)
    assert_equal(fr.params_after_calc["TRUCR"], fr.tru_cr)
    assert_equal(fr.params_prior_calc["P_NL"],  0.0)
    assert_equal(fr.params_after_calc["P_NL"], fr.P_NL)
    assert_equal(fr.params_prior_calc["U"],  fr.mat_feed_u.mass)
    assert_equal(fr.params_after_calc["U"], fr.mat_prod_u.mass)
    assert_equal(fr.params_prior_calc["TRU"],  fr.mat_feed_tru.mass)
    assert_equal(fr.params_after_calc["TRU"], fr.mat_prod_tru.mass)
    assert_equal(fr.params_prior_calc["ACT"],  fr.mat_feed_act.mass)
    assert_equal(fr.params_after_calc["ACT"], fr.mat_prod_act.mass)
    assert_equal(fr.params_prior_calc["LAN"],  fr.mat_feed_lan.mass)
    assert_equal(fr.params_after_calc["LAN"], fr.mat_prod_lan.mass)
    assert_equal(fr.params_prior_calc["FP"],  1.0 - fr.mat_feed_act.mass - fr.mat_feed_lan.mass)
    assert_equal(fr.params_after_calc["FP"], 1.0 - fr.mat_prod_act.mass - fr.mat_prod_lan.mass)
    

# Put Integral tests here, if desired.


if __name__ == "__main__":
    nose.main()
