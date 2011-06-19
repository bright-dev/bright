import os

import numpy as np

from bright import *
import isoname
from mass_stream import MassStream

import serpent

from char.m2py import convert_res, convert_dep


def run_serpent():
    bu_file = 'serp_bu'
    xs_file = 'serp_xs'

    if (not os.path.exists(bu_file + '_res.py')) or (not os.path.exists(bu_file + '_dep.py')):
        serpent.main(bu_file)
        convert_res(bu_file + '_res.m')
        convert_dep(bu_file + '_dep.m')

    if (not os.path.exists(xs_file + '_res.py')):
        serpent.main(xs_file)
        convert_res(xs_file + '_res.m')

    res_bu = {}
    execfile(bu_file + '_res.py', {}, res_bu)

    dep_bu = {}
    execfile(bu_file + '_dep.py', {}, dep_bu)

    res_xs = {}
    execfile(xs_file + '_res.py', {}, res_xs)

    return res_bu, dep_bu, res_xs


def run_reactormg():
    # Init bright
    libfile = os.getenv("BRIGHT_DATA") + 'lwr_mg.h5'
    load_track_isos_hdf5(libfile)
    bright_config.write_text = False
    bright_config.write_hdf5 = False
    bright_config.verbosity = 100

    # Init reactor paramters
    rp = lwr_defaults()
    rp.batches = 3
    rp.flux = 4*(10**14)

    rp.fuel_form = {"IHM": 1.0, "O16": 2.0}
    rp.cladding_form = {"ZR93": 0.5, "ZR95": 0.5}
    rp.coolant_form = {"H1": 2.0, "O16": 1.0}

    rp.fuel_density = 10.7
    rp.cladding_density = 5.87
    rp.coolant_density = 0.73

    rp.pnl = 0.98
    rp.BUt = 50.0
    rp.use_disadvantage_factor = True
    rp.lattice_type = 'Spherical'
    rp.lattice_type = 'Planar'
    rp.rescale_hydrogen = True

    rp.fuel_radius = 0.412
    rp.void_radius = 0.4205
    rp.clad_radius = 0.475
    rp.unit_cell_pitch = 1.33

    rp.open_slots = 25
    rp.total_slots = 289

    rp.burn_times = np.linspace(0.0, 365.0, 10)

    # Init mass stream
    leu = MassStream({922350: 0.05, 922380: 0.95})

    # Init ReactorMG
    rmg = ReactorMG(reactor_parameters=rp, name="rmg")
    rmg.loadlib(libfile)

    # Run the reactor
    #rmg.calc(leu)
    rmg.ms_feed = leu
    rmg.burnup_core()
    return rmg


def test_regression():
    res_bu, dep_bu, res_xs = run_serpent()
    rmg = run_reactormg()
    return rmg, res_bu, dep_bu, res_xs

if __name__ == "__main__":
    rmg, res_bu, dep_bu, res_xs = test_regression()

    print "Reactor k: ", rmg.k_t
    print "Serpent k: ", res_bu['SIX_FF_KEFF'][:, 0]
    print "Fractional Diff: ", 1.0 - rmg.k_t / res_bu['SIX_FF_KEFF'][:, 0]
    print

    r_phi = rmg.phi_tg /rmg.phi_t
    s_phi = res_bu['FLUX'][:, 2::2] / res_bu['FLUX'][:, np.newaxis, 0]
    print "Normalized Reactor phi: ", r_phi
    print "Normalized Serpent phi: ", s_phi
    print "Fractional Diff: ", 1.0 - r_phi / s_phi
    print

