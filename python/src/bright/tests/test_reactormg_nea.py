import os

import numpy as np

from bright import *
import isoname
from mass_stream import MassStream

import serpent

from char.m2py import convert_res, convert_dep


def run_nea():
    mw_out = {'U234':  1.778273E-04,
              'U235':  8.005087E-03,
              'U236':  4.842712E-03,
              'U238':  9.337753E-01,
              'NP237': 6.146369E-04,
              'PU238': 2.258085E-04,
              'PU239': 5.993966E-03,
              'PU240': 2.389836E-03,
              'PU241': 1.636997E-03,
              'PU242': 6.026438E-04,
              'AM241': 4.672990E-05,
              'AM243': 1.481173E-04,
              'MO95':  1.977351E-03,
              'TC99':  2.227474E-03,
              'RU101': 2.191019E-03,
              'RH103': 1.163463E-03,
              'AG109': 1.723287E-04,
              'CS133': 2.403604E-03,
              'ND143': 1.552958E-03,
              'ND145': 1.328163E-03,
              'SM147': 1.395532E-04,
              'SM149': 4.860031E-06,
              'SM150': 6.024706E-04,
              'SM151': 2.881857E-05,
              'SM152': 2.344635E-04,
              'EU153': 2.069572E-04,
             }
    mw_out = {isoname.LLAAAM_2_zzaaam(key): value for key, value in mw_out.items()}

    ms = MassStream(mw_out)
    return ms

def run_reactormg():
    # Init bright
    libfile = os.getenv("BRIGHT_DATA") + 'lwr_mg.h5'
    load_track_isos_hdf5(libfile)
    bright_config.write_text = False
    bright_config.write_hdf5 = False
    bright_config.verbosity = 100

    # Init reactor paramters
    rp = lwr_defaults()
    rp.batches = 4
    rp.flux = 4E+14

    rp.fuel_form = {"IHM": 1.0, "O16": 2.0}
    rp.cladding_form = {"ZR93": 0.5, "ZR95": 0.5}
    rp.coolant_form = {"H1": 2.0, "O16": 1.0}

    rp.fuel_density = 10.7
    rp.cladding_density = 5.87
    rp.coolant_density = 0.7295

    rp.pnl = 0.98
    rp.BUt = 40.0
    rp.use_disadvantage_factor = True
    rp.lattice_type = 'Spherical'
    rp.lattice_type = 'Planar'
    rp.rescale_hydrogen = True
    rp.burnup_via_constant = 'flux'

    rp.fuel_radius = 0.409575
    rp.void_radius = 0.41783
    rp.clad_radius = 0.47483
    rp.unit_cell_pitch = 1.25984

    rp.open_slots = 25
    rp.total_slots = 289

    rp.burn_times = np.linspace(0.0, 500.0, 50)

    # Init mass stream
    leu = MassStream({922340: 0.00032, 
                      922350: 0.036, 
                      922360: 0.00016,
                      922380: 0.96352})

    # Init ReactorMG
    rmg = ReactorMG(reactor_parameters=rp, name="rmg")
    rmg.loadlib(libfile)

    # Run the reactor
    rmg.calc(leu)
    #rmg.ms_feed = leu
    #rmg.burnup_core()
    return rmg


def test_regression():
    ms_nea = run_nea()
    rmg = run_reactormg()
    return rmg, ms_nea

def calc_diff(r, n, name=""):
    print "Summary for {0}:".format(name)
    print "Reactor: "
    print repr(r)
    print "NEA: "
    print repr(n)
    print "Fractional Diff: "
    diff = 1.0 - r / n
    print repr(diff)
    print
    return r, n, diff


if __name__ == "__main__":
    rmg, ms_nea = test_regression()

    r_BU, n_BU, diff_BU = calc_diff(rmg.BUd, 40.0, "Burnup")

    rmg_comp = rmg.ms_prod.mult_by_mass()
    nea_comp = ms_nea.mult_by_mass()

    for key in nea_comp.keys():
        if key in rmg_comp:
            r_, n_, diff_ = calc_diff(rmg_comp[key], nea_comp[key], "Mass of " + isoname.zzaaam_2_LLAAAM(key))

    mss = [MassStream({i: T_it[i][t] for i in T_it.keys()}) for t in range(len(rmg.burn_times))]
    masses = r_mass, s_mass, diff_mass = calc_diff(np.array([ms.mass for ms in mss]), 1.0, "Mass")
