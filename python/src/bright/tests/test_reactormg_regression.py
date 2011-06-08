import os

import numpy as np

from bright import *
import isoname
from mass_stream import MassStream

import serpent

from char.m2py import convert_res, convert_dep


def run_serpent():
    input_file = 'serpent_reactormg_regression'

    if (not os.path.exists(input_file + '_res.py')) or (not os.path.exists(input_file + '_dep.py')):
        serpent.main(input_file)
        convert_res(input_file + '_res.m')
        convert_dep(input_file + '_dep.m')

    res = {}
    exec(input_file + '_res.py', {}, res)

    dep = {}
    exec(input_file + '_dep.py', {}, dep)

    return res, dep


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
    rmg.calc(leu)
    return rmg


def test_regression():
    res, dep = run_serpent()
    rmg = run_reactormg()


if __name__ == "__main__":
    test_regression()
