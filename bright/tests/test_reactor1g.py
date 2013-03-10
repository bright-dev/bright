"""Reactor1G Class Tests"""

import nose 
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
assert_almost_equal, assert_true, assert_false, set_trace, with_setup

from numpy.testing import assert_array_equal, assert_array_almost_equal

import os
import warnings
import tables as tb
import numpy as np

import bright
import bright.reactor_parameters
import bright.fluence_point
import bright.reactor1g
import pyne.material

#
# Shortcuts
#

bright_conf = bright.bright_conf
Material = pyne.material.Material
FluencePoint = bright.fluence_point.FluencePoint
ReactorParameters = bright.reactor_parameters.ReactorParameters
Reactor1G = bright.reactor1g.Reactor1G

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

r1g = None

def setup_r1g():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G()
    r1g.loadlib(libfile)


def setup_r1g_attrs():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G(rp=default_rp, n='r1g')
    r1g.loadlib(libfile)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.fold_mass_weights()

def setup_r1g_data_attrs():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G(rp=default_rp, n='r1g')
    r1g.loadlib(libfile)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.fold_mass_weights()


def setup_r1g_discharge():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G(rp=default_rp, n='r1g')
    r1g.loadlib(libfile)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.calc()

def setup_r1g_sub():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G(rp=default_rp, n='r1g')
    r1g.loadlib(libfile)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.calc()
    r1g.calc_sub_mats()


def setup_r1g_zeta():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G(rp=default_rp, n='r1g')
    r1g.loadlib(libfile)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.fold_mass_weights()


def setup_r1g_init():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G(rp=default_rp, n='r1g')
    r1g.loadlib(libfile)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.fold_mass_weights()


def setup_r1g_transmute():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G(rp=default_rp, n='r1g')
    r1g.loadlib(libfile)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.fold_mass_weights()


def setup_r1g_basic():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G(rp=default_rp, n='r1g')
    r1g.loadlib(libfile)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.fold_mass_weights()


def setup_r1g_burnup():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G(rp=default_rp, n='r1g')
    r1g.loadlib(libfile)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.calc()


def setup_r1g_lattice():
    global r1g
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    r1g = Reactor1G(rp=default_rp, n='r1g')
    r1g.loadlib(libfile)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.calc()


def resetup_r1g_lattice():
    global r1g
    r1g.initialize(default_rp)
    r1g.mat_feed = Material({922350: 0.5, 922380: 0.5})
    r1g.calc()


def teardown_r1g():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "r1g.h5"]:
            os.remove(f)

def teardown_r1g_clear():
    global r1g
    teardown_r1g()
    r1g = None

#
# Tests that the Reactor1G component constructors work.
#


@with_setup(None, teardown_r1g)
def test_Reactor1G_1():
    r1g = Reactor1G()
    assert_equal(r1g.name, '')
    assert_equal(r1g.track_params, set())

@with_setup(None, teardown_r1g)
def test_Reactor1G_2():
    r1g = Reactor1G(n="r1g")
    assert_equal(r1g.name, 'r1g')
    assert_equal(r1g.track_params, set())

@with_setup(None, teardown_r1g)
def test_Reactor1G_3():
    r1g = Reactor1G(paramtrack=set(["Mass"]))
    assert_equal(r1g.name, '')
    assert_equal(r1g.track_params, set(["Mass"]))

@with_setup(None, teardown_r1g)
def test_Reactor1G_4():
    r1g = Reactor1G(paramtrack=set(["Mass"]), n='r1g')
    assert_equal(r1g.name, 'r1g')
    assert_equal(r1g.track_params, set(["Mass"]))

@with_setup(None, teardown_r1g)
def test_Reactor1G_5():
    rp = ReactorParameters()
    r1g = Reactor1G(rp=rp, n="r1g")
    assert_equal(r1g.name, 'r1g')
    assert_equal(r1g.track_params, set())
    assert_equal(r1g.B, 0)
    assert_equal(r1g.phi, 0.0)
    assert_equal(r1g.fuel_chemical_form, {})
    assert_equal(r1g.coolant_chemical_form, {})
    assert_equal(r1g.rhoF, 0.0)
    assert_equal(r1g.rhoC, 0.0)
    assert_equal(r1g.P_NL, 0.0)
    assert_equal(r1g.target_BU, 0.0)
    assert_false(r1g.use_zeta)
    assert_equal(r1g.lattice_flag, '')
    assert_false(r1g.rescale_hydrogen_xs)
    assert_equal(r1g.r, 0.0)
    assert_equal(r1g.l, 0.0)
    assert_equal(r1g.S_O, 0.0)
    assert_equal(r1g.S_T, 0.0)

@with_setup(None, teardown_r1g)
def test_Reactor1G_6():
    rp = ReactorParameters()
    r1g = Reactor1G(rp=rp, paramtrack=set(["Mass"]), n="r1g")
    assert_equal(r1g.name, 'r1g')
    assert_equal(r1g.track_params, set(["Mass"]))
    assert_equal(r1g.B, 0)
    assert_almost_equal(r1g.phi, 0.0)
    assert_equal(r1g.fuel_chemical_form, {})
    assert_equal(r1g.coolant_chemical_form, {})
    assert_almost_equal(r1g.rhoF, 0.0)
    assert_almost_equal(r1g.rhoC, 0.0)
    assert_almost_equal(r1g.P_NL, 0.0)
    assert_almost_equal(r1g.target_BU, 0.0)
    assert_false(r1g.use_zeta)
    assert_equal(r1g.lattice_flag, '')
    assert_false(r1g.rescale_hydrogen_xs)
    assert_almost_equal(r1g.r, 0.0)
    assert_almost_equal(r1g.l, 0.0)
    assert_almost_equal(r1g.S_O, 0.0)
    assert_almost_equal(r1g.S_T, 0.0)


#
# Tests that the Reactor1G parameter attributes work.
#


def test_B():
    r1g = Reactor1G()
    r1g.B = 3
    assert_equal(r1g.B, 3)

@with_setup(None, teardown_r1g)
def test_phi():
    r1g = Reactor1G()
    r1g.phi = 2*(10**14)
    assert_equal(r1g.phi, 2*(10**14))

@with_setup(None, teardown_r1g)
def test_fuel_chemical_form():
    r1g = Reactor1G()
    r1g.fuel_chemical_form = {"IHM": 1.0, "O16": 2.0}
    assert_equal(r1g.fuel_chemical_form, {"IHM": 1.0, "O16": 2.0})

@with_setup(None, teardown_r1g)
def test_coolant_chemical_form():
    r1g = Reactor1G()
    r1g.coolant_chemical_form = {"H1": 2.0, "O16": 1.0,
                                 "B10": 0.199 * 550 * 10.0**-6, 
                                 "B11": 0.801 * 550 * 10.0**-6}
    assert_equal(r1g.coolant_chemical_form, {"H1": 2.0, "O16": 1.0,
                                             "B10": 0.199 * 550 * 10.0**-6,
                                             "B11": 0.801 * 550 * 10.0**-6})

@with_setup(None, teardown_r1g)
def test_rhoF():
    r1g = Reactor1G()
    r1g.rhoF = 10.7
    assert_equal(r1g.rhoF, 10.7)

@with_setup(None, teardown_r1g)
def test_rhoC():
    r1g = Reactor1G()
    r1g.rhoC = 0.73
    assert_equal(r1g.rhoC, 0.73)

@with_setup(None, teardown_r1g)
def test_P_NL():
    r1g = Reactor1G()
    r1g.P_NL = 0.98
    assert_equal(r1g.P_NL, 0.98)

@with_setup(None, teardown_r1g)
def test_target_BU():
    r1g = Reactor1G()
    r1g.target_BU = 50.0
    assert_equal(r1g.target_BU, 50.0)

@with_setup(None, teardown_r1g)
def test_use_zeta():
    r1g = Reactor1G()
    r1g.use_zeta = True
    assert_true(r1g.use_zeta)

@with_setup(None, teardown_r1g)
def test_Lattice():
    r1g = Reactor1G()
    r1g.lattice_flag = 'Spherical'
    assert_equal(r1g.lattice_flag, 'Spherical')

@with_setup(None, teardown_r1g)
def test_rescale_hydrogen_xs():
    r1g = Reactor1G()
    r1g.rescale_hydrogen_xs = True
    assert_true(r1g.rescale_hydrogen_xs)

@with_setup(None, teardown_r1g)
def test_r():
    r1g = Reactor1G()
    r1g.r = 0.411
    assert_equal(r1g.r, 0.411)

@with_setup(None, teardown_r1g)
def test_l():
    r1g = Reactor1G()
    r1g.l = 0.7
    assert_equal(r1g.l, 0.7)

@with_setup(None, teardown_r1g)
def test_S_O():
    r1g = Reactor1G()
    r1g.S_O = 123
    assert_equal(r1g.S_O, 123)

@with_setup(None, teardown_r1g)
def test_S_T():
    r1g = Reactor1G()
    r1g.S_T = 180
    assert_equal(r1g.S_T, 180)

@with_setup(None, teardown_r1g)
def test_VF():
    rp = ReactorParameters()
    rp.fuel_radius = 0.5
    rp.unit_cell_pitch = 1.0
    rp.open_slots = 0
    rp.total_slots = 1
    r1g = Reactor1G(rp=rp, n='r1g')
    assert_almost_equal(r1g.VF, 3.14159265*0.25) 

@with_setup(None, teardown_r1g)
def test_VC():
    rp = ReactorParameters()
    rp.fuel_radius = 0.5
    rp.unit_cell_pitch = 1.0
    rp.open_slots = 0
    rp.total_slots = 1
    r1g = Reactor1G(rp=rp, n='r1g')
    assert_almost_equal(r1g.VC, 1.0 - 3.14159265*0.25) 




#
# Tests that the Reactor1G basic data attributes work.
#

@with_setup(setup_r1g, teardown_r1g)
def test_libfile():
    r1g.libfile = "It's Ruth!"
    assert_equal(r1g.libfile, "It's Ruth!")

@with_setup(None, teardown_r1g)
def test_F():
    assert_equal(r1g.F[0], 0.0)
    assert(1 < len(r1g.F))
    for f in range(1, len(r1g.F)):        
        assert(r1g.F[f-1] < r1g.F[f])

    old_F = r1g.F

    r1g.F = np.arange(0.0, 10.0)        
    assert_array_equal(r1g.F, np.arange(0.0, 10.0))

    r1g.F = old_F

@with_setup(None, teardown_r1g)
def test_BUi_F_():
    BUi_F_ = r1g.BUi_F_
    for i in BUi_F_.keys():
        assert_equal(BUi_F_[i][0], 0.0)
        assert_equal(len(r1g.F), len(BUi_F_[i]))
        for f in range(1, len(r1g.F)):        
            assert(BUi_F_[i][f-1] <= BUi_F_[i][f])


    old_BU = r1g.BUi_F_
    r1g.BUi_F_ = {1: np.arange(0.0, 10.0)}
    assert_equal(r1g.BUi_F_.keys(), [1])
    assert_array_equal(r1g.BUi_F_[1], np.arange(0.0, 10.0))
    r1g.BUi_F_ = old_BU


@with_setup(None, teardown_r1g)
def test_pi_F_():
    pi_F_ = r1g.pi_F_
    for i in pi_F_.keys():
        assert_equal(len(r1g.F), len(pi_F_[i]))

@with_setup(None, teardown_r1g)
def test_di_F_():
    di_F_ = r1g.di_F_
    for i in di_F_.keys():
        assert_equal(len(r1g.F), len(di_F_[i]))

# FIXME once we have bound map[int, map[int, vector[double]]]
#@with_setup(None, teardown_r1g_clear)
#def test_Tij_F_():
#    Tij_F_ = r1g.Tij_F_
#    jsos   = bright_conf.track_nucs
#    for i in Tij_F_.keys():
#        for j in jsos:
#            assert_equal(len(r1g.F), len(Tij_F_[i][j]))
#
#    old_T = r1g.Tij_F_
#    r1g.Tij_F_ = {1: {2: np.arange(0.0, 10.0)}}
#    assert_equal(r1g.Tij_F_.keys(), [1])
#    assert_equal(r1g.Tij_F_[1].keys(), [2])
#    assert_array_equal(r1g.Tij_F_[1][2], np.arange(0.0, 10.0))
#    r1g.Tij_F_ = old_T


#
# Tests that the Reactor1G calculated weight attributes work.
#

@with_setup(setup_r1g_attrs, teardown_r1g)
def test_A_IHM():
    assert_almost_equal(r1g.A_IHM / 236.5, 1.0, 3)

@with_setup(None, teardown_r1g)
def test_MWF():
    assert_almost_equal(r1g.MWF / 268.5, 1.0, 3)

@with_setup(None, teardown_r1g)
def test_MWC():
    assert_almost_equal(r1g.MWC / 18.0, 1.0, 2)

@with_setup(None, teardown_r1g)
def test_niF():
    assert_equal(r1g.niF[80160],  2.0)
    assert_equal(r1g.niF[922350], 0.5)
    assert_equal(r1g.niF[922380], 0.5)

@with_setup(None, teardown_r1g)
def test_niC():
    assert_equal(r1g.niC[10010],  2.0)
    assert_equal(r1g.niC[80160],  1.0)

@with_setup(None, teardown_r1g)
def test_miF():
    assert_almost_equal(r1g.miF[80160],  2.0 * 16  / 236.5, 3)
    assert_almost_equal(r1g.miF[922350], 0.5 * 235 / 236.5, 3)
    assert_almost_equal(r1g.miF[922380], 0.5 * 238 / 236.5, 3)

@with_setup(None, teardown_r1g)
def test_miC():
    #First calculate the relative volume
    rel_Vol = (r1g.rhoC * 268.5 * r1g.VC) / (r1g.rhoF * 18.0 * r1g.VF)
    assert_almost_equal(r1g.miC[10010],  rel_Vol * 2.0 * 1   / 236.5, 3)
    assert_almost_equal(r1g.miC[80160],  rel_Vol * 1.0 * 16  / 236.5, 3)
    
@with_setup(None, teardown_r1g)
def test_NiF():
    assert_almost_equal(r1g.NiF[80160]  / (2.0 * r1g.rhoF * 6.022*(10**23) / 268.5), 1.0, 3)
    assert_almost_equal(r1g.NiF[922350] / (0.5 * r1g.rhoF * 6.022*(10**23) / 268.5), 1.0, 3)
    assert_almost_equal(r1g.NiF[922380] / (0.5 * r1g.rhoF * 6.022*(10**23) / 268.5), 1.0, 3)

@with_setup(None, teardown_r1g_clear)
def test_NiC():
    assert_almost_equal(r1g.NiC[10010] / (2.0 * r1g.rhoC * 6.022*(10**23) / 18.0), 1.0, 2)
    assert_almost_equal(r1g.NiC[80160] / (1.0 * r1g.rhoC * 6.022*(10**23) / 18.0), 1.0, 2)


#
# Tests that the Reactor1G calculated data attributes work.
#

@with_setup(setup_r1g_data_attrs, teardown_r1g)
def test_BU_F_():
    BU_F_ = r1g.BU_F_
    BUi_F_ = r1g.BUi_F_
    assert_equal(BU_F_[0], 0.0)
    for f in range(1, len(r1g.F)):
        tmp_BU = 0.0
        for i in r1g.miF.keys():
            tmp_BU = tmp_BU + (r1g.miF[i] * BUi_F_[i][f])
        assert_almost_equal(BU_F_[f] / tmp_BU, 1.0, 4)

@with_setup(None, teardown_r1g)
def test_P_F_():
    P_F_ = r1g.P_F_
    pi_F_ = r1g.pi_F_
    for f in range(len(r1g.F)):
        tmp_P = 0.0
        for i in r1g.miF.keys():
            tmp_P = tmp_P + (r1g.P_NL * r1g.miF[i] * pi_F_[i][f])
        assert_almost_equal(P_F_[f] / tmp_P, 1.0, 4)

@with_setup(None, teardown_r1g)
def test_dF_F_():
    dF_F_ = r1g.dF_F_
    di_F_ = r1g.di_F_
    for f in range(len(r1g.F)):
        tmp_dF = 0.0
        for i in r1g.miF.keys():
            tmp_dF = tmp_dF + (r1g.miF[i] * di_F_[i][f])
        assert_almost_equal(dF_F_[f] / tmp_dF, 1.0, 4)

@with_setup(None, teardown_r1g)
def test_dC_F_():
    dC_F_ = r1g.dC_F_
    di_F_ = r1g.di_F_
    zeta_F_ = r1g.zeta_F_
    for f in range(len(r1g.F)):
        tmp_dC = 0.0
        for i in r1g.miC.keys():
            if i == 10010 and r1g.rescale_hydrogen_xs:
                tmp_dC = tmp_dC + (r1g.miC[i] * di_F_[i][f] * (1.36927-(0.01119*r1g.BU_F_[f])) )
            else:
                tmp_dC = tmp_dC + (r1g.miC[i] * di_F_[i][f])
        if r1g.use_zeta:
            tmp_dC = tmp_dC * zeta_F_[f]
        assert_almost_equal(dC_F_[f] / tmp_dC, 1.0, 4)

@with_setup(None, teardown_r1g)
def test_D_F_():
    D_F_  = r1g.D_F_
    dF_F_ = r1g.dF_F_
    dC_F_ = r1g.dC_F_
    for f in range(len(r1g.F)):
        assert_almost_equal(D_F_[f] / (dF_F_[f] + dC_F_[f]), 1.0, 4)

@with_setup(None, teardown_r1g)
def test_k_F_():
    P_F_ = r1g.P_F_
    D_F_ = r1g.D_F_
    k_F_ = r1g.k_F_
    for f in range(len(r1g.F)):
        assert_almost_equal(k_F_[f] / (P_F_[f]/D_F_[f]), 1.0, 4)

# FIXME once we bind Tij_F_
#@with_setup(None, teardown_r1g)
#def test_Mj_F_():
#    Mj_F_  = r1g.Mj_F_
#    Tij_F_ = r1g.Tij_F_
#    for j in Mj_F_.keys():
#        for f in range(len(r1g.F)):
#            tmp_Mj = 0.0
#            for i in r1g.miF.keys():
#                tmp_Mj = tmp_Mj + (r1g.miF[i] * Tij_F_[i][j][f])
#            assert_almost_equal(Mj_F_[j][f] / tmp_Mj, 1.0, 4)

@with_setup(None, teardown_r1g_clear)
def test_zeta_F_():
    zeta_F_ = r1g.zeta_F_
    for f in range(len(r1g.F)):
        assert(1.0 <= zeta_F_[f])


#
# Tests that the Reactor1G discharge attributes are right.
#

@with_setup(setup_r1g_discharge, teardown_r1g)
def test_fd():
    assert(0 <= r1g.fd)
    assert(r1g.fd <= len(r1g.F))

@with_setup(None, teardown_r1g)
def test_Fd():
    assert(r1g.F[r1g.fd] <= r1g.Fd)
    assert(r1g.Fd <= r1g.F[r1g.fd+1])

@with_setup(None, teardown_r1g)
def test_BUd():
    assert(r1g.BU_F_[r1g.fd] <= r1g.BUd)
    assert(r1g.BUd <= r1g.BU_F_[r1g.fd+1])

@with_setup(None, teardown_r1g_clear)
def test_k():
    assert_almost_equal(r1g.k, 1.0, 6)


#
# Tests that the Reactor1G sub-stream and transuranic conversion ratio
# attributes are right.
#

@with_setup(setup_r1g_sub, teardown_r1g)
def test_mat_feed_u():
    assert_equal(r1g.mat_feed_u.mass, 1.0)

@with_setup(None, teardown_r1g)
def test_mat_feed_tru():
    assert_equal(r1g.mat_feed_tru.mass, 0.0)

@with_setup(None, teardown_r1g)
def test_mat_feed_lan():
    assert_equal(r1g.mat_feed_lan.mass, 0.0)

@with_setup(None, teardown_r1g)
def test_mat_feed_act():
    assert_equal(r1g.mat_feed_act.mass, 1.0)

@with_setup(None, teardown_r1g)
def test_mat_prod_u():
    assert(r1g.mat_prod_u.mass < 1.0)

@with_setup(None, teardown_r1g)
def test_mat_prod_tru():
    assert(0.0 < r1g.mat_prod_tru.mass)

@with_setup(None, teardown_r1g)
def test_mat_prod_lan():
    assert(0.0 < r1g.mat_prod_lan.mass)

@with_setup(None, teardown_r1g)
def test_mat_prod_act():
    assert(r1g.mat_prod_act.mass < 1.0)

@with_setup(None, teardown_r1g)
def test_tru_cr():
    r1g.calc_tru_cr()
    tmp_tru_cr = 1.0 - (r1g.mat_feed_tru.mass - r1g.mat_prod_tru.mass) / (r1g.BUd / 931.46)
    assert_almost_equal(r1g.tru_cr / tmp_tru_cr, 1.0)

@with_setup(None, teardown_r1g_clear)
def test_deltaR():
    r1g.calc_deltaR()
    tmp_deltaR = r1g.batch_average(r1g.target_BU, "p") - r1g.batch_average(r1g.target_BU, "d")
    assert_almost_equal(r1g.deltaR / tmp_deltaR, 1.0)


#
# Tests that the Reactor1G calculated data attributes work.
#

@with_setup(setup_r1g_zeta, teardown_r1g)
def test_SigmaFa_F_():
    "Can not do full test due to private sigma_a_therm data member."
    assert_equal(len(r1g.SigmaFa_F_), len( r1g.F))

@with_setup(None, teardown_r1g)
def test_SigmaFtr_F_():
    "Can not do full test due to private sigma_a_therm data member."
    assert_equal(len(r1g.SigmaFtr_F_), len( r1g.F))

@with_setup(None, teardown_r1g)
def test_kappaF_F_():
    kappaF_F_ = r1g.kappaF_F_
    SigmaFa_F_ = r1g.SigmaFa_F_
    SigmaFtr_F_ = r1g.SigmaFtr_F_
    for f in range(len(r1g.F)):
        tmp_kappaF = np.sqrt(3.0 * SigmaFtr_F_[f] * SigmaFa_F_[f])
        assert_almost_equal(kappaF_F_[f] / tmp_kappaF, 1.0)

@with_setup(None, teardown_r1g)
def test_SigmaCa_F_():
    "Can not do full test due to private sigma_a_therm data member."
    assert_equal(len(r1g.SigmaCa_F_), len( r1g.F))

@with_setup(None, teardown_r1g)
def test_SigmaCtr_F_():
    "Can not do full test due to private sigma_a_therm data member."
    assert_equal(len(r1g.SigmaCtr_F_), len( r1g.F))

@with_setup(None, teardown_r1g_clear)
def test_kappaC_F_():
    kappaC_F_ = r1g.kappaC_F_
    SigmaCa_F_ = r1g.SigmaCa_F_
    SigmaCtr_F_ = r1g.SigmaCtr_F_
    for f in range(len(r1g.F)):
        tmp_kappaC = np.sqrt(3.0 * SigmaCtr_F_[f] * SigmaCa_F_[f])
        assert_almost_equal(kappaC_F_[f] / tmp_kappaC, 1.0)

#Test lattice_E_F_ here

#Test lattice_F_F_ here


#
# Tests that the fuel cycle component initialization methods work.
#

@with_setup(setup_r1g_init, teardown_r1g)
def test_initialize():
    rp = ReactorParameters()
    r1g.initialize(rp)
    assert_equal(r1g.B, 0)
    assert_equal(r1g.phi, 0.0)
    assert_equal(r1g.fuel_chemical_form, {})
    assert_equal(r1g.coolant_chemical_form, {})
    assert_equal(r1g.rhoF, 0.0)
    assert_equal(r1g.rhoC, 0.0)
    assert_equal(r1g.P_NL, 0.0)
    assert_equal(r1g.target_BU, 0.0)
    assert_false(r1g.use_zeta)
    assert_equal(r1g.lattice_flag, '')
    assert_false(r1g.rescale_hydrogen_xs)
    assert_equal(r1g.r, 0.0)
    assert_equal(r1g.l, 0.0)
    assert_equal(r1g.S_O, 0.0)
    assert_equal(r1g.S_T, 0.0)

@with_setup(None, teardown_r1g)
def test_loadlib():
    r1g.loadlib(os.getenv("BRIGHT_DATA") + '/FR.h5')
    r1g.loadlib(os.getenv("BRIGHT_DATA") + '/LWR.h5')

@with_setup(None, teardown_r1g_clear)
def test_fold_mass_weights():
    prevkey = r1g.miF.keys()
    assert(922380 in r1g.miF.keys())
    r1g.mat_feed = Material({922350: 0.5})
    r1g.fold_mass_weights()
    assert(922380 not in r1g.miF.keys())

    

#
# Tests that the fuel cycle component transmutation matrix methods work.
#

# FIXME once Tij_F_ has bindings
#@with_setup(setup_r1g_transmute, teardown_r1g)
#def test_calc_Mj_F_():
#    r1g.calc_Mj_F_()
#    # Test below is the same as test_Mj_F_()
#    Mj_F_  = r1g.Mj_F_
#    Tij_F_ = r1g.Tij_F_
#    for j in Mj_F_.keys():
#        for f in range(len(r1g.F)):
#            tmp_Mj = 0.0
#            for i in r1g.mat_feed.comp.keys():
#                tmp_Mj = tmp_Mj + (r1g.miF[i] * Tij_F_[i][j][f])
#            if tmp_Mj == 0.0:
#                assert_almost_equal(Mj_F_[j][f], 0.0, 4)
#            else:
#                assert_almost_equal(Mj_F_[j][f] / tmp_Mj, 1.0, 4)

@with_setup(setup_r1g_transmute, teardown_r1g_clear)
def test_calc_Mj_Fd_():
    r1g.BUd_bisection_method()
    r1g.calc_Mj_F_()
    r1g.calc_Mj_Fd_()
    assert(0.0 < r1g.mat_prod.mass)
    assert(r1g.mat_prod.mass < 1.0)



#
# Tests that the Reactor1G basic calculation methods work.
#


@with_setup(setup_r1g_basic, teardown_r1g)
def test_calc_mat_prod():
    r1g.BUd_bisection_method()
    r1g.calc_mat_prod()
    assert(0.0 < r1g.mat_prod.mass)
    assert(r1g.mat_prod.mass < 1.0)


@with_setup(None, teardown_r1g)
def test_calc_sub_mats():
    r1g.calc()
    r1g.calc_sub_mats()
    assert_equal(r1g.mat_feed_u.mass, 1.0)
    assert_equal(r1g.mat_feed_tru.mass, 0.0)
    assert_equal(r1g.mat_feed_lan.mass, 0.0)
    assert_equal(r1g.mat_feed_act.mass, 1.0)
    assert(r1g.mat_prod_u.mass < 1.0)
    assert(0.0 < r1g.mat_prod_tru.mass)
    assert(0.0 < r1g.mat_prod_lan.mass)
    assert(r1g.mat_prod_act.mass < 1.0)

@with_setup(None, teardown_r1g)
def test_calc_tru_cr():
    r1g.calc()
    tmp_tru_cr = 1.0 - (r1g.mat_feed_tru.mass - r1g.mat_prod_tru.mass) / (r1g.BUd / 931.46)
    assert_almost_equal(r1g.calc_tru_cr() / tmp_tru_cr, 1.0)

@with_setup(None, teardown_r1g)
def test_deltaR1():
    r1g.calc_deltaR()
    tmp_deltaR = r1g.batch_average(r1g.target_BU, "p") - r1g.batch_average(r1g.target_BU, "d")
    assert_almost_equal(r1g.deltaR / tmp_deltaR, 1.0)

@with_setup(None, teardown_r1g)
def test_deltaR2():
    r1g.calc_deltaR({922350: 0.5, 922380: 0.5, 80160: 0.125})
    tmp_deltaR = r1g.batch_average(r1g.target_BU, "p") - r1g.batch_average(r1g.target_BU, "d")
    assert_almost_equal(r1g.deltaR / tmp_deltaR, 1.0)

@with_setup(None, teardown_r1g_clear)
def test_deltaR3():
    ms = Material({922350: 0.5, 922380: 0.5, 80160: 0.125})
    r1g.calc_deltaR(ms)
    tmp_deltaR = r1g.batch_average(r1g.target_BU, "p") - r1g.batch_average(r1g.target_BU, "d")
    assert_almost_equal(r1g.deltaR / tmp_deltaR, 1.0)

#
# Tests that the Reactor1G burnup methods work.
#

@with_setup(setup_r1g_burnup, teardown_r1g)
def test_fluence_at_BU():
    fp = r1g.fluence_at_BU(80.0)
    assert(0 <= fp.f)
    assert(fp.f <= len(r1g.F))
    assert(r1g.F[fp.f] <= fp.F)
    assert(fp.F <= r1g.F[fp.f+1])
    assert(r1g.BU_F_[fp.f] <= 80.0)
    assert(80.0 <= r1g.BU_F_[fp.f+1])
    tmp_m = (r1g.BU_F_[fp.f+1] - r1g.BU_F_[fp.f]) / (r1g.F[fp.f+1] - r1g.F[fp.f])
    assert_equal(fp.m / tmp_m, 1.0)


@with_setup(None, teardown_r1g)
def test_batch_average():
    BUd = r1g.BUd
    p = r1g.batch_average(BUd, "P")
    d = r1g.batch_average(BUd, "D")
    k = r1g.batch_average(BUd, "K")
    kk = r1g.batch_average(BUd)
    assert_equal(k, kk)
    #assert_equal(p/d, k) # Averaging messes this up.

@with_setup(None, teardown_r1g)
def test_batch_average_k():
    BUd = r1g.BUd
    assert_almost_equal(r1g.batch_average_k(BUd), 1.0, 6)

@with_setup(None, teardown_r1g)
def test_calc_1():
    r1g.calc()
    assert(r1g.mat_prod.mass < 1.0)
    assert(r1g.mat_prod.comp[922350] < 0.5) 

@with_setup(None, teardown_r1g)
def test_calc_2():
    r1g.calc(Material({942390: 0.05, 922380: 0.95}))
    assert(r1g.mat_prod.mass < 1.0)
    assert(r1g.mat_prod.comp[942390] < 1.0) 


@with_setup(None, teardown_r1g)
def test_BUd_bisection_method():
    assert_almost_equal(r1g.k, 1.0, 5)
    r1g.B = 1
    r1g.BUd_bisection_method()
    assert_almost_equal(r1g.k, 1.0, 5)
    r1g.B = 3


@with_setup(None, teardown_r1g)
def test_run_P_NL():
    # Convergence is not gaurenteed!
    r1g.run_P_NL(0.99)
    assert_equal(r1g.P_NL, 0.99)
    assert_almost_equal(r1g.k, 1.0, 1)
    r1g.run_P_NL(0.98)


@with_setup(None, teardown_r1g_clear)
def test_calibrate_P_NL_to_BUd():
    r1g.calibrate_P_NL_to_BUd()
    assert_not_equal(r1g.P_NL, 0.98)
    assert_almost_equal(r1g.BUd / r1g.target_BU, 1.0, 1)
    

#
# Tests that the Reactor1G Lattice methods work. These are not exposed to Python directly :(
#


@with_setup(setup_r1g_lattice, teardown_r1g)
def test_lattice_E_planar():
    prev = np.array(r1g.lattice_E_F_)
    r1g.lattice_flag = "Planar"
    r1g.r = 0.5
    r1g.l = 1.0
    r1g.fold_mass_weights()
    curr = r1g.lattice_E_F_
    for f in range(len(r1g.F)):
        assert_not_equal(prev[f], curr[f])

@with_setup(resetup_r1g_lattice, teardown_r1g)
def test_lattice_F_planar():
    prev = np.array(r1g.lattice_F_F_)
    r1g.lattice_flag = "Planar"
    r1g.r = 0.5
    r1g.l = 1.0
    r1g.fold_mass_weights()
    curr = r1g.lattice_F_F_
    for f in range(len(r1g.F)):
        assert_not_equal(prev[f], curr[f])

@with_setup(resetup_r1g_lattice, teardown_r1g)
def test_lattice_E_spherical():
    prev = np.array(r1g.lattice_E_F_)
    r1g.lattice_flag = "Spherical"
    r1g.r = 0.5
    r1g.l = 1.0
    r1g.fold_mass_weights()
    curr = r1g.lattice_E_F_
    for f in range(len(r1g.F)):
        assert_not_equal(prev[f], curr[f])

@with_setup(resetup_r1g_lattice, teardown_r1g)
def test_lattice_F_spherical():
    prev = np.array(r1g.lattice_F_F_)
    r1g.lattice_flag = "Spherical"
    r1g.r = 0.5
    r1g.l = 1.0
    r1g.fold_mass_weights()
    curr = r1g.lattice_F_F_
    for f in range(len(r1g.F)):
        assert_not_equal(prev[f], curr[f])

@with_setup(resetup_r1g_lattice, teardown_r1g)
def test_lattice_E_cylindrical():
    prev = np.array(r1g.lattice_E_F_)
    r1g.lattice_flag = "Cylindrical"
    r1g.r = 0.5
    r1g.l = 1.0
    r1g.fold_mass_weights()
    curr = r1g.lattice_E_F_
    for f in range(len(r1g.F)):
        assert_not_equal(prev[f], curr[f])


@with_setup(resetup_r1g_lattice, teardown_r1g_clear)
def test_lattice_F_cylindrical():
    prev = np.array(r1g.lattice_F_F_)
    r1g.lattice_flag = "Cylindrical"
    r1g.r = 0.5
    r1g.l = 1.0
    r1g.fold_mass_weights()
    curr = r1g.lattice_F_F_
    for f in range(len(r1g.F)):
        assert_not_equal(prev[f], curr[f])

# Since the above are not exposed directly, 
# They implicitly test the following Reactor1G functions:
#   calc_zeta()
#   calc_zeta_planar()
#   calc_zeta_spherical()
#   calc_zeta_cylindrical()


if __name__ == "__main__":
    nose.main()
