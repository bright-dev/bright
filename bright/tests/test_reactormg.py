"""ReactorMG Tests"""

#import faulthandler
#faulthandler.enable()

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
assert_almost_equal, assert_true, assert_false, set_trace, with_setup

from numpy.testing import dec, assert_array_equal, assert_array_almost_equal

import os
import warnings
import tables as tb
import numpy as np

from bright import bright_conf, load_track_nucs_hdf5
from bright.reactor_parameters import ReactorParameters, lwr_defaults
from bright.fluence_point import FluencePoint
from bright.reactormg import ReactorMG
from pyne.material import Material
from pyne import nucname

#
# Shortcuts
#

#default_rp = ReactorParameters()
default_rp = lwr_defaults()

default_rp.batches = 3
default_rp.flux = 2*(10**14)

default_rp.fuel_form = {"IHM": 1.0, "O16": 2.0}
default_rp.cladding_form = {"ZR93": 0.5, "ZR95": 0.5}
default_rp.coolant_form = {"H1": 2.0, "O16": 1.0}

default_rp.fuel_density = 10.7
default_rp.cladding_density = 5.87
default_rp.coolant_density = 0.73

default_rp.pnl = 0.98
default_rp.BUt = 50.0

default_rp.use_disadvantage_factor = True
default_rp.rescale_hydrogen = True

#default_rp.lattice_type = 'Cylindrical'
default_rp.lattice_type = 'Spherical'


#default_rp.burn_times = np.linspace(0.0, 4200.0, 5)
#default_rp.burn_times = np.linspace(0.0, 100.0, 5)
default_rp.burn_times = np.linspace(0.0, 2100.0, 5)

default_rp.fuel_radius = 0.412
default_rp.void_radius = 0.4205
default_rp.clad_radius = 0.475
default_rp.unit_cell_pitch = 1.33


#default_rp.unit_cell_pitch = 0.7
#default_rp.unit_cell_pitch = 0.4751 * 2.0
#default_rp.unit_cell_pitch = 100.0
#default_rp.unit_cell_pitch = 0.4121 * 2.0


default_rp.open_slots = 25
default_rp.total_slots = 289


# Only Fuel
#default_rp.void_radius = 0.412
#default_rp.clad_radius = 0.412
#default_rp.unit_cell_pitch = 0.412 * np.sqrt(np.pi)
#default_rp.open_slots = 0

# Only Cool
#default_rp.fuel_radius = 0.001
#default_rp.void_radius = 0.0
#default_rp.clad_radius = 0.0
#default_rp.unit_cell_pitch = 1.33
#default_rp.open_slots = 0


#
# Fixtures
#

rmg = None

def setup_rmg_attr():
    global rmg
    libfile = os.getenv("BRIGHT_DATA") + '/lwr_base.h5'
    load_track_nucs_hdf5(libfile)
    rmg = ReactorMG()
    rmg.loadlib(libfile)


def teardown_rmg():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "rmg.h5"]:
            os.remove(f)

def teardown_rmg_clear():
    global rmg
    teardown_rmg()
    rmg = None

#
# Tests that the ReactorMG component constructors work.

@with_setup(None, teardown_rmg)
def test_ReactorMG_1():
    rmg = ReactorMG()
    assert_equal(rmg.name, '')
    assert_equal(rmg.track_params, set())

@with_setup(None, teardown_rmg)
def test_ReactorMG_2():
    rmg = ReactorMG(n="rmg")
    assert_equal(rmg.name, 'rmg')
    assert_equal(rmg.track_params, set())

@with_setup(None, teardown_rmg)
def test_ReactorMG_3():
    rmg = ReactorMG(paramtrack=set(["Mass"]))
    assert_equal(rmg.name, '')
    assert_equal(rmg.track_params, set(["Mass"]))

@with_setup(None, teardown_rmg)
def test_ReactorMG_4():
    rmg = ReactorMG(paramtrack=set(["Mass"]), n='rmg')
    assert_equal(rmg.name, 'rmg')
    assert_equal(rmg.track_params, set(["Mass"]))

@with_setup(None, teardown_rmg)
def test_ReactorMG_5():
    rp = ReactorParameters()
    rmg = ReactorMG(rp=rp, n="rmg")
    assert_equal(rmg.name, 'rmg')
    assert_equal(rmg.track_params, set())
    assert_equal(rmg.B, 0)
    assert_equal(rmg.flux, 0.0)
    assert_equal(rmg.chemical_form_fuel, {})
    assert_equal(rmg.chemical_form_clad, {})
    assert_equal(rmg.chemical_form_cool, {})
    assert_equal(rmg.rho_fuel, 0.0)
    assert_equal(rmg.rho_cool, 0.0)
    assert_equal(rmg.P_NL, 0.0)
    assert_equal(rmg.target_BU, 0.0)
    assert_false(rmg.use_zeta)
    assert_equal(rmg.lattice_flag, '')
    assert_false(rmg.rescale_hydrogen_xs)
    assert_equal(rmg.r_fuel, 0.0)
    assert_equal(rmg.pitch, 0.0)
    assert_equal(rmg.S_O, 0.0)
    assert_equal(rmg.S_T, 0.0)

@with_setup(None, teardown_rmg)
def test_ReactorMG_6():
    rp = ReactorParameters()
    rmg = ReactorMG(rp=rp, paramtrack=set(["Mass"]), n="rmg")
    assert_equal(rmg.name, 'rmg')
    assert_equal(rmg.track_params, set(["Mass"]))
    assert_equal(rmg.B, 0)
    assert_almost_equal(rmg.flux, 0.0)
    assert_equal(rmg.chemical_form_fuel, {})
    assert_equal(rmg.chemical_form_clad, {})
    assert_equal(rmg.chemical_form_cool, {})
    assert_almost_equal(rmg.rho_fuel, 0.0)
    assert_almost_equal(rmg.rho_cool, 0.0)
    assert_almost_equal(rmg.P_NL, 0.0)
    assert_almost_equal(rmg.target_BU, 0.0)
    assert_false(rmg.use_zeta)
    assert_equal(rmg.lattice_flag, '')
    assert_false(rmg.rescale_hydrogen_xs)
    assert_almost_equal(rmg.r_fuel, 0.0)
    assert_almost_equal(rmg.pitch, 0.0)
    assert_almost_equal(rmg.S_O, 0.0)
    assert_almost_equal(rmg.S_T, 0.0)



#
# Tests that the ReactorMG parameter attributes work.
#

@with_setup(None, teardown_rmg)
def test_B():
    rmg = ReactorMG()
    rmg.B = 3
    assert_equal(rmg.B, 3)

@with_setup(None, teardown_rmg)
def test_flux():
    rmg = ReactorMG()
    rmg.flux = 2*(10**14)
    assert_equal(rmg.flux, 2*(10**14))

@with_setup(None, teardown_rmg)
def test_chemical_form_fuel():
    rmg = ReactorMG()
    rmg.chemical_form_fuel = {"IHM": 1.0, "O16": 2.0}
    assert_equal(rmg.chemical_form_fuel, {"IHM": 1.0, "O16": 2.0})

@with_setup(None, teardown_rmg)
def test_chemical_form_cool():
    rmg = ReactorMG()
    rmg.chemical_form_cool = {"H1": 2.0, "O16": 1.0,
        "B10": 0.199 * 550 * 10.0**-6, 
        "B11": 0.801 * 550 * 10.0**-6}
    assert_equal(rmg.chemical_form_cool, {"H1": 2.0, "O16": 1.0,
        "B10": 0.199 * 550 * 10.0**-6,
        "B11": 0.801 * 550 * 10.0**-6})

@with_setup(None, teardown_rmg)
def test_rho_fuel():
    rmg = ReactorMG()
    rmg.rho_fuel = 10.7
    assert_equal(rmg.rho_fuel, 10.7)

@with_setup(None, teardown_rmg)
def test_rho_cool():
    rmg = ReactorMG()
    rmg.rho_cool = 0.73
    assert_equal(rmg.rho_cool, 0.73)

@with_setup(None, teardown_rmg)
def test_P_NL():
    rmg = ReactorMG()
    rmg.P_NL = 0.98
    assert_equal(rmg.P_NL, 0.98)

@with_setup(None, teardown_rmg)
def test_target_BU():
    rmg = ReactorMG()
    rmg.target_BU = 50.0
    assert_equal(rmg.target_BU, 50.0)

@with_setup(None, teardown_rmg)
def test_use_zeta():
    rmg = ReactorMG()
    rmg.use_zeta = True
    assert_true(rmg.use_zeta)

@with_setup(None, teardown_rmg)
def test_Lattice():
    rmg = ReactorMG()
    rmg.lattice_flag = 'Spherical'
    assert_equal(rmg.lattice_flag, 'Spherical')

@with_setup(None, teardown_rmg)
def test_rescale_hydrogen_xs():
    rmg = ReactorMG()
    rmg.rescale_hydrogen_xs = True
    assert_true(rmg.rescale_hydrogen_xs)

@with_setup(None, teardown_rmg)
def test_r():
    rmg = ReactorMG()
    rmg.r_fuel = 0.411
    assert_equal(rmg.r_fuel, 0.411)

@with_setup(None, teardown_rmg)
def test_l():
    rmg = ReactorMG()
    rmg.pitch = 0.7
    assert_equal(rmg.pitch, 0.7)

@with_setup(None, teardown_rmg)
def test_S_O():
    rmg = ReactorMG()
    rmg.S_O = 123
    assert_equal(rmg.S_O, 123)

@with_setup(None, teardown_rmg)
def test_S_T():
    rmg = ReactorMG()
    rmg.S_T = 180
    assert_equal(rmg.S_T, 180)

@with_setup(None, teardown_rmg)
def test_V_fuel():
    rp = ReactorParameters()
    rp.fuel_radius = 0.5
    rp.unit_cell_pitch = 1.0
    rp.open_slots = 0
    rp.total_slots = 1
    rmg = ReactorMG(rp=rp, n='rmg')
    assert_almost_equal(rmg.V_fuel, 3.14159265*0.25) 

@with_setup(None, teardown_rmg)
def test_V_cool():
    rp = ReactorParameters()
    rp.fuel_radius = 0.5
    rp.void_radius = 0.5
    rp.clad_radius = 0.5
    rp.unit_cell_pitch = 1.0
    rp.open_slots = 0
    rp.total_slots = 1
    rmg = ReactorMG(rp=rp, n='rmg')
    assert_almost_equal(rmg.V_cool, 1.0 - 3.14159265*0.25) 


#
# Tests that the ReactorMG basic data attributes work.
#

@with_setup(setup_rmg_attr, teardown_rmg)
def test_libfile():
    rmg.libfile = "It's Ruth!"
    assert_equal(rmg.libfile, "It's Ruth!")


@with_setup(None, teardown_rmg)
def test_npertubations():
    assert(0 < rmg.nperturbations)


@with_setup(None, teardown_rmg)
def test_perturbed_fields():
    pf = rmg.perturbed_fields
    for key, value in pf.items():
        assert_equal(len(value), 3)
        assert_equal(value[1] - value[0], value[2])


@with_setup(None, teardown_rmg)
def test_I():
    assert_not_equal(len(rmg.I), 0)
    for iso in rmg.I:
        assert_equal(nucname.current_form(iso), 'zzaaam')


@with_setup(None, teardown_rmg)
def test_J():
    assert_not_equal(len(rmg.J), 0)
    for iso in rmg.J:
        assert_equal(nucname.current_form(iso), 'zzaaam')


#@with_setup(None, teardown_rmg)
#def test_IJ():
#    # Bad dataset
#    print rmg.J
#    assert(rmg.I <= rmg.J)


@with_setup(None, teardown_rmg)
def test_G():
    assert(0 < rmg.G)


@with_setup(None, teardown_rmg)
def test_E_g():
    assert_equal(len(rmg.E_g), rmg.G+1)
    assert((rmg.E_g[1:] < rmg.E_g[:-1]).all())


@dec.skipif(True)
@with_setup(None, teardown_rmg)
def test_phi_g():
    phi_g = rmg.phi_g
    assert_equal(len(phi_g), rmg.nperturbations)
    assert_equal(len(phi_g[0]), rmg.G)


@dec.skipif(True)
@with_setup(None, teardown_rmg)
def test_phi():
    assert_equal(len(rmg.phi), rmg.nperturbations)
    assert_array_almost_equal(rmg.phi_g.sum(axis=1) / rmg.phi, np.ones(rmg.nperturbations), 5)


@with_setup(None, teardown_rmg)
def test_Phi():
    assert_equal(len(rmg.Phi), rmg.nperturbations)


@with_setup(None, teardown_rmg)
def test_time0():
    assert_equal(len(rmg.time0), rmg.nperturbations)


@with_setup(None, teardown_rmg)
def test_BU0():
    assert_equal(len(rmg.BU0), rmg.nperturbations)


@with_setup(None, teardown_rmg)
def test_Ti0():
    J = rmg.J
    Ti0 = rmg.Ti0
    nperturbations = rmg.nperturbations
    for iso in Ti0.keys():
        assert(iso in J)
        assert_equal(len(Ti0[iso]), nperturbations)


@dec.skipif(True)
@with_setup(None, teardown_rmg)
def test_sigma_f_pg():
    J = rmg.J
    G = rmg.G
    sigma_f_pg = rmg.sigma_f_pg
    nperturbations = rmg.nperturbations
    for iso in sigma_f_pg.keys():
        assert_equal(sigma_f_pg[iso].shape, (nperturbations, G))


@dec.skipif(True)
@with_setup(None, teardown_rmg)
def test_nubar_sigma_f_pg():
    J = rmg.J
    G = rmg.G
    nubar_sigma_f_pg = rmg.nubar_sigma_f_pg
    nperturbations = rmg.nperturbations
    for iso in nubar_sigma_f_pg.keys():
        assert_equal(nubar_sigma_f_pg[iso].shape, (nperturbations, G))



@dec.skipif(True)
@with_setup(None, teardown_rmg)
def test_sigma_s_pgh():
    J = rmg.J
    G = rmg.G
    sigma_s_pgh = rmg.sigma_s_pgh
    nperturbations = rmg.nperturbations
    for iso in sigma_s_pgh.keys():
        assert_equal(sigma_s_pgh[iso].shape, (nperturbations, G, G))



#
# Cannot implement the following tests without a woring native burnup method.
#
"""\
class TestReactorMGMutliGroupMethods(TestCase):
"Tests that the ReactorMG basic data attributes work."

@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + 'lwr_mg.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.05, 922380: 0.95})
#    rmg.mat_feed = MassStream({922350: 0.03, 922380: 0.97})
    rmg.burnup_core()
    rmg.burn_time = 0.0
    rmg.bt_s = 0

@classmethod
def teardown_class():
    general_teardown()


def test_calc_nearest_neighbors():
    assert_equal(len(rmg.nearest_neighbors), rmg.nperturbations)


def test_interpolate_cross_sections():
    # Grab data from C++
    G = rmg.G
    J = rmg.J
    S = rmg.S
    s_f_itg = rmg.sigma_f_itg
    ns_f_itg = rmg.nubar_sigma_f_itg
    s_s_itgh = rmg.sigma_s_itgh

    # Test Data
    assert_equal(len(s_f_itg), len(J))        
    assert_equal(len(ns_f_itg), len(J))
    assert_equal(len(s_s_itgh), len(J))

    for j in J:
        # Assert shapes
        assert_equal(s_f_itg[j].shape, (S, G))
        assert_equal(ns_f_itg[j].shape, (S, G))
        assert_equal(s_s_itgh[j].shape, (S, G, G))

        # Assert Values
        assert (0.0 <= s_f_itg[j]).all()
        assert (0.0 <= ns_f_itg[j]).all()
        assert (0.0 <= s_s_itgh[j]).all()


def test_calc_mass_weights():
    # FIXME
    
    # Grab data from C++
    G = rmg.G
    J = rmg.J
    S = rmg.S

    assert (0 <= rmg.A_HM_t).all()
    assert (0 <= rmg.MW_fuel_t).all()
    assert (18 == rmg.MW_cool_t).all()
    assert (94 <= rmg.MW_clad_t).all()

    for j in J:
        #for s in range(S):
        for s in range(1):
            # Test Fuel atomic density
            assert 0 <= rmg.n_fuel_it[j][s] 
            assert rmg.n_fuel_it[j][s] <= 3

            # Test cool and clad atomic denisty
            if j == 10010:
                assert_equal(2.0, rmg.n_cool_it[j][s])
            elif j == 80160:
                assert_equal(1.0, rmg.n_cool_it[j][s])
            elif j == 400930:
                assert_equal(0.5, rmg.n_clad_it[j][s])
            elif j == 400950:
                assert_equal(0.5, rmg.n_clad_it[j][s])
            else:
                assert_equal(0.0, rmg.n_cool_it[j][s])
                assert_equal(0.0, rmg.n_clad_it[j][s])

            # Test number densities
            if rmg.N_fuel_it[j][s] != 0.0: 
                assert_almost_equal(1.0, rmg.n_fuel_it[j][s] * rmg.rho_fuel * 6.0221415E+23 / rmg.MW_fuel_t[s] / rmg.N_fuel_it[j][s]) 

            if rmg.N_clad_it[j][s] != 0.0: 
                assert_almost_equal(1.0, rmg.n_clad_it[j][s] * rmg.rho_clad * 6.0221415E+23 / rmg.MW_clad_t[s] / rmg.N_clad_it[j][s])

            if rmg.N_cool_it[j][s] != 0.0: 
                assert_almost_equal(1.0, rmg.n_cool_it[j][s] * rmg.rho_cool * 6.0221415E+23 / rmg.MW_cool_t[s] / rmg.N_cool_it[j][s])



def test_fold_mass_weights():
    # FIXME

    # Grab data from C++
    G = rmg.G
    J = rmg.J
    S = rmg.S

    #for s in range(S):
    for s in range(1):
        assert (0 <= rmg.Sigma_t_tg[s]).all()
        assert (rmg.Sigma_f_tg[s] <= rmg.Sigma_t_tg[s]).all()

        assert (1.0 < rmg.nubar_Sigma_f_tg[s] / rmg.Sigma_f_tg[s]).all()
        assert (rmg.nubar_Sigma_f_tg[s] / rmg.Sigma_f_tg[s] <= 4.0).all()

        assert (0.0 <= rmg.chi_tg[s]).all()
        assert (rmg.chi_tg[s] <= 1.0).all()
        assert_almost_equal(1.0, np.sum(rmg.chi_tg[s]))

        assert (rmg.Sigma_gamma_tg[s] <= rmg.Sigma_t_tg[s]).all()
        assert (rmg.Sigma_2n_tg[s] <= rmg.Sigma_t_tg[s]).all()
        assert (rmg.Sigma_3n_tg[s] <= rmg.Sigma_t_tg[s]).all()
        assert (rmg.Sigma_alpha_tg[s] <= rmg.Sigma_t_tg[s]).all()
        assert (rmg.Sigma_proton_tg[s] <= rmg.Sigma_t_tg[s]).all()
        assert (rmg.Sigma_gamma_x_tg[s] <= rmg.Sigma_t_tg[s]).all()
        assert (rmg.Sigma_2n_x_tg[s] <= rmg.Sigma_t_tg[s]).all()


def test_assemble_multigroup_matrices():
    # Grab data from C++
    G = rmg.G
    J = rmg.J
    S = rmg.S

    dia = np.identity(G, dtype=bool)
    ndia = True - dia

    #for s in range(S):
    for s in range(1):
        # Check A diagonal is positive, and off-diagonal elements
        # are negative or zero.  This is a test of the physicality of 
        # The algorithm
        assert (0 < rmg.A_tgh[s][dia]).all()
        assert (rmg.A_tgh[s][ndia] <= 0).all()

        assert_array_almost_equal(rmg.F_tgh[s], 
            (rmg.nubar_Sigma_f_tg[s][:, np.newaxis] *  rmg.chi_tg[s]) )

        # Test that A_inv is truly the inverso of A
        assert_array_almost_equal(np.dot(rmg.A_inv_tgh[s], rmg.A_tgh[s]), np.identity(G))

        assert_array_almost_equal(rmg.A_inv_F_tgh[s], 
            np.dot(rmg.A_inv_tgh[s], rmg.F_tgh[s]) )

        assert not np.isnan(rmg.T_int_tij[s]).any()
        assert not np.isnan(rmg.M_tij[s]).any()



def test_calc_criticality():
    #raise nose.SkipTest
    # Grab data from C++
    G = rmg.G
    J = rmg.J
    S = rmg.S

    #print rmg.phi_tg[0]
    #print rmg.phi_t
    #print rmg.k_t
    #print rmg.P_NL

#    print rmg.A_inv_F_tgh[0]
#    print 
    s = 0
    k0 = 1.0
    phi0 = np.ones(G, dtype=float)
    phi0[2] = 2.0

    "-"\
    print "phi0 = ", phi0
    print
#    for i in range(100):
    for i in range(3):
        phi1 = (1.0/(0.98 * k0)) * np.dot(rmg.A_inv_F_tgh[s], phi0)
        k1 = k0 * sum(rmg.nubar_Sigma_f_tg[s] * phi1) / sum(rmg.nubar_Sigma_f_tg[s] * phi0)

        print "k1 = ", k1
        print "phi1 = ", phi1
        print 

        k0 = k1
        phi0 = phi1
    "-"

#    print rmg.phi_tg
#    print 
    print "s = ", s
    print
    print "Sigma_a_fuel_tg[s] = ", rmg.Sigma_a_fuel_tg[s]
    print
    print "Sigma_a_cool_tg[s] = ", rmg.Sigma_a_cool_tg[s]
    print
    print "zeta_tg[s] = ", rmg.zeta_tg[s]
    print
    print "k_t = ", rmg.k_t
    print
    print "phi_tg[s] = ", rmg.phi_tg[s]
    print
    print rmg.V_fuel
    print rmg.V_clad
    print rmg.V_cool

    raise TypeError 
"""

################################ Unset below



"""\
class TestReactorMGCalculatedWeightAttributes(TestCase):
    "Tests that the ReactorMG calculated weight attributes work."

@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.fold_mass_weights()

@classmethod
def teardown_class():
    general_teardown()

def test_A_IHM():
    assert_almost_equal(rmg.A_IHM / 236.5, 1.0, 4)

def test_MWF():
    assert_almost_equal(rmg.MWF / 268.5, 1.0, 4)

def test_MWC():
    assert_almost_equal(rmg.MWC / 18.0, 1.0, 4)

def test_niF():
    assert_equal(rmg.niF[80160],  2.0)
    assert_equal(rmg.niF[922350], 0.5)
    assert_equal(rmg.niF[922380], 0.5)

def test_niC():
    assert_equal(rmg.niC[10010],  2.0)
    assert_equal(rmg.niC[80160],  1.0)

def test_miF():
    assert_almost_equal(rmg.miF[80160],  2.0 * 16  / 236.5, 4)
    assert_almost_equal(rmg.miF[922350], 0.5 * 235 / 236.5, 4)
    assert_almost_equal(rmg.miF[922380], 0.5 * 238 / 236.5, 4)

def test_miC():
    #First calculate the relative volume
    rel_Vol = (rmg.rho_cool * 268.5 * rmg.V_cool) / (rmg.rho_fuel * 18.0 * rmg.V_fuel)
    assert_almost_equal(rmg.miC[10010],  rel_Vol * 2.0 * 1   / 236.5, 4)
    assert_almost_equal(rmg.miC[80160],  rel_Vol * 1.0 * 16  / 236.5, 4)
    
def test_NiF():
    assert_almost_equal(rmg.NiF[80160]  / (2.0 * rmg.rho_fuel * 6.022*(10**23) / 268.5), 1.0, 3)
    assert_almost_equal(rmg.NiF[922350] / (0.5 * rmg.rho_fuel * 6.022*(10**23) / 268.5), 1.0, 3)
    assert_almost_equal(rmg.NiF[922380] / (0.5 * rmg.rho_fuel * 6.022*(10**23) / 268.5), 1.0, 3)

def test_NiC():
    assert_almost_equal(rmg.NiC[10010] / (2.0 * rmg.rho_cool * 6.022*(10**23) / 18.0), 1.0, 3)
    assert_almost_equal(rmg.NiC[80160] / (1.0 * rmg.rho_cool * 6.022*(10**23) / 18.0), 1.0, 3)


class TestReactorMGCalculatedDataAttributes(TestCase):
"Tests that the ReactorMG calculated data attributes work."

@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.fold_mass_weights()

@classmethod
def teardown_class():
    general_teardown()

def test_BU_F_():
    BU_F_ = rmg.BU_F_
    BUi_F_ = rmg.BUi_F_
    assert_equal(BU_F_[0], 0.0)
    for f in range(1, len(rmg.F)):
        tmp_BU = 0.0
        for i in rmg.miF.keys():
            tmp_BU = tmp_BU + (rmg.miF[i] * BUi_F_[i][f])
        assert_almost_equal(BU_F_[f] / tmp_BU, 1.0, 4)

def test_P_F_():
    P_F_ = rmg.P_F_
    pi_F_ = rmg.pi_F_
    for f in range(len(rmg.F)):
        tmp_P = 0.0
        for i in rmg.miF.keys():
            tmp_P = tmp_P + (rmg.P_NL * rmg.miF[i] * pi_F_[i][f])
        assert_almost_equal(P_F_[f] / tmp_P, 1.0, 4)

def test_dF_F_():
    dF_F_ = rmg.dF_F_
    di_F_ = rmg.di_F_
    for f in range(len(rmg.F)):
        tmp_dF = 0.0
        for i in rmg.miF.keys():
            tmp_dF = tmp_dF + (rmg.miF[i] * di_F_[i][f])
        assert_almost_equal(dF_F_[f] / tmp_dF, 1.0, 4)

def test_dC_F_():
    dC_F_ = rmg.dC_F_
    di_F_ = rmg.di_F_
    zeta_F_ = rmg.zeta_F_
    for f in range(len(rmg.F)):
        tmp_dC = 0.0
        for i in rmg.miC.keys():
            if i == 10010 and rmg.rescale_hydrogen_xs:
                tmp_dC = tmp_dC + (rmg.miC[i] * di_F_[i][f] * (1.36927-(0.01119*rmg.BU_F_[f])) )
            else:
                tmp_dC = tmp_dC + (rmg.miC[i] * di_F_[i][f])
        if rmg.use_zeta:
            tmp_dC = tmp_dC * zeta_F_[f]
        assert_almost_equal(dC_F_[f] / tmp_dC, 1.0, 4)

def test_D_F_():
    D_F_  = rmg.D_F_
    dF_F_ = rmg.dF_F_
    dC_F_ = rmg.dC_F_
    for f in range(len(rmg.F)):
        assert_almost_equal(D_F_[f] / (dF_F_[f] + dC_F_[f]), 1.0, 4)

def test_k_F_():
    P_F_ = rmg.P_F_
    D_F_ = rmg.D_F_
    k_F_ = rmg.k_F_
    for f in range(len(rmg.F)):
        assert_almost_equal(k_F_[f] / (P_F_[f]/D_F_[f]), 1.0, 4)

def test_Mj_F_():
    Mj_F_  = rmg.Mj_F_
    Tij_F_ = rmg.Tij_F_
    for j in Mj_F_.keys():
        for f in range(len(rmg.F)):
            tmp_Mj = 0.0
            for i in rmg.miF.keys():
                tmp_Mj = tmp_Mj + (rmg.miF[i] * Tij_F_[i][j][f])
            assert_almost_equal(Mj_F_[j][f] / tmp_Mj, 1.0, 4)

def test_zeta_F_():
    zeta_F_ = rmg.zeta_F_
    for f in range(len(rmg.F)):
        assert(1.0 <= zeta_F_[f])


class TestReactorMGDischargeAttributes(TestCase):
"Tests that the ReactorMG discharge attributes are right."

@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.calc()

@classmethod
def teardown_class():
    general_teardown()

def test_fd():
    assert(0 <= rmg.fd)
    assert(rmg.fd <= len(rmg.F))

def test_Fd():
    assert(rmg.F[rmg.fd] <= rmg.Fd)
    assert(rmg.Fd <= rmg.F[rmg.fd+1])

def test_BUd():
    assert(rmg.BU_F_[rmg.fd] <= rmg.BUd)
    assert(rmg.BUd <= rmg.BU_F_[rmg.fd+1])

def test_k():
    assert_almost_equal(rmg.k, 1.0)


class TestReactorMGSubStreamAndtru_crAttributes(TestCase):
"Tests that the ReactorMG sub-stream and transuranic conversion ratio attributes are right."

@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.calc()
    rmg.calc_sub_mats()

@classmethod
def teardown_class():
    general_teardown()

def test_mat_feed_u():
    assert_equal(rmg.mat_feed_u.mass, 1.0)

def test_mat_feed_tru():
    assert_equal(rmg.mat_feed_tru.mass, 0.0)

def test_mat_feed_lan():
    assert_equal(rmg.mat_feed_lan.mass, 0.0)

def test_mat_feed_act():
    assert_equal(rmg.mat_feed_act.mass, 1.0)

def test_mat_prod_u():
    assert(rmg.mat_prod_u.mass < 1.0)

def test_mat_prod_tru():
    assert(0.0 < rmg.mat_prod_tru.mass)

def test_mat_prod_lan():
    assert(0.0 < rmg.mat_prod_lan.mass)

def test_mat_prod_act():
    assert(rmg.mat_prod_act.mass < 1.0)

def test_tru_cr():
    rmg.calc_tru_cr()
    tmp_tru_cr = 1.0 - (rmg.mat_feed_tru.mass - rmg.mat_prod_tru.mass) / (rmg.BUd / 931.46)
    assert_almost_equal(rmg.tru_cr / tmp_tru_cr, 1.0)

def test_deltaR():
    rmg.calc_deltaR()
    tmp_deltaR = rmg.batch_average(rmg.target_BU, "p") - rmg.batch_average(rmg.target_BU, "d")
    assert_almost_equal(rmg.deltaR / tmp_deltaR, 1.0)


class TestReactorMGThermalDisadvantageFactorAttributes(TestCase):
"Tests that the ReactorMG calculated data attributes work."

@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.fold_mass_weights()

@classmethod
def teardown_class():
    general_teardown()


def test_SigmaFa_F_():
    "Can not do full test do to private sigma_a_therm data member."
    assert_equal(len(rmg.SigmaFa_F_), len( rmg.F))

def test_SigmaFtr_F_():
    "Can not do full test do to private sigma_a_therm data member."
    assert_equal(len(rmg.SigmaFtr_F_), len( rmg.F))

def test_kappaF_F_():
    kappaF_F_ = rmg.kappaF_F_
    SigmaFa_F_ = rmg.SigmaFa_F_
    SigmaFtr_F_ = rmg.SigmaFtr_F_
    for f in range(len(rmg.F)):
        tmp_kappaF = np.sqrt(3.0 * SigmaFtr_F_[f] * SigmaFa_F_[f])
        assert_almost_equal(kappaF_F_[f] / tmp_kappaF, 1.0)

def test_SigmaCa_F_():
    "Can not do full test do to private sigma_a_therm data member."
    assert_equal(len(rmg.SigmaCa_F_), len( rmg.F))

def test_SigmaCtr_F_():
    "Can not do full test do to private sigma_a_therm data member."
    assert_equal(len(rmg.SigmaCtr_F_), len( rmg.F))

def test_kappaC_F_():
    kappaC_F_ = rmg.kappaC_F_
    SigmaCa_F_ = rmg.SigmaCa_F_
    SigmaCtr_F_ = rmg.SigmaCtr_F_
    for f in range(len(rmg.F)):
        tmp_kappaC = np.sqrt(3.0 * SigmaCtr_F_[f] * SigmaCa_F_[f])
        assert_almost_equal(kappaC_F_[f] / tmp_kappaC, 1.0)

#Test lattice_E_F_ here

#Test lattice_F_F_ here


class TestReactorMGInitializationMethods(TestCase):
"Tests that the fuel cycle component initialization methods work."

@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.fold_mass_weights()

@classmethod
def teardown_class():
    general_teardown()

def test_initialize():
    rp = ReactorParameters()
    rmg.initialize(rp)
    assert_equal(rmg.B, 0)
    assert_equal(rmg.flux, 0.0)
    assert_equal(rmg.chemical_form_fuel, {})
    assert_equal(rmg.chemical_form_cool, {})
    assert_equal(rmg.rho_fuel, 0.0)
    assert_equal(rmg.rho_cool, 0.0)
    assert_equal(rmg.P_NL, 0.0)
    assert_equal(rmg.target_BU, 0.0)
    assert_false(rmg.use_zeta)
    assert_equal(rmg.lattice_flag, '')
    assert_false(rmg.rescale_hydrogen_xs)
    assert_equal(rmg.r_fuel, 0.0)
    assert_equal(rmg.pitch, 0.0)
    assert_equal(rmg.S_O, 0.0)
    assert_equal(rmg.S_T, 0.0)

def test_loadlib():
    rmg.loadlib(os.getenv("BRIGHT_DATA") + '/FR.h5')
    rmg.loadlib(os.getenv("BRIGHT_DATA") + '/LWR.h5')

def test_fold_mass_weights():
    prevkey = rmg.miF.keys()
    assert(922380 in rmg.miF.keys())
    rmg.mat_feed = MassStream({922350: 0.5})
    rmg.fold_mass_weights()
    assert(922380 not in rmg.miF.keys())

    

class TestReactorMGTransmutationMatrixMethods(TestCase):
"Tests that the fuel cycle component transmutation matrix methods work."


@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.fold_mass_weights()

@classmethod
def teardown_class():
    general_teardown()

def test_calc_Mj_F_():
    rmg.calc_Mj_F_()
    # Test below is the same as test_Mj_F_()
    Mj_F_  = rmg.Mj_F_
    Tij_F_ = rmg.Tij_F_
    for j in Mj_F_.keys():
        for f in range(len(rmg.F)):
            tmp_Mj = 0.0
            for i in rmg.mat_feed.comp.keys():
                tmp_Mj = tmp_Mj + (rmg.miF[i] * Tij_F_[i][j][f])
            if tmp_Mj == 0.0:
                assert_almost_equal(Mj_F_[j][f], 0.0, 4)
            else:
                assert_almost_equal(Mj_F_[j][f] / tmp_Mj, 1.0, 4)

def test_calc_Mj_Fd_():
    rmg.BUd_bisection_method()
    rmg.calc_Mj_F_()
    rmg.calc_Mj_Fd_()
    assert(0.0 < rmg.mat_prod.mass)
    assert(rmg.mat_prod.mass < 1.0)



class TestReactorMGBasicCalculationMethods(TestCase):
"Tests that the ReactorMG basic calculation methods work."


@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.fold_mass_weights()

@classmethod
def teardown_class():
    general_teardown()

def test_calc_mat_prod():
    rmg.BUd_bisection_method()
    rmg.calc_mat_prod()
    assert(0.0 < rmg.mat_prod.mass)
    assert(rmg.mat_prod.mass < 1.0)


def test_calc_sub_mats():
    rmg.calc()
    rmg.calc_sub_mats()
    assert_equal(rmg.mat_feed_u.mass, 1.0)
    assert_equal(rmg.mat_feed_tru.mass, 0.0)
    assert_equal(rmg.mat_feed_lan.mass, 0.0)
    assert_equal(rmg.mat_feed_act.mass, 1.0)
    assert(rmg.mat_prod_u.mass < 1.0)
    assert(0.0 < rmg.mat_prod_tru.mass)
    assert(0.0 < rmg.mat_prod_lan.mass)
    assert(rmg.mat_prod_act.mass < 1.0)

def test_calc_tru_cr():
    rmg.calc()
    tmp_tru_cr = 1.0 - (rmg.mat_feed_tru.mass - rmg.mat_prod_tru.mass) / (rmg.BUd / 931.46)
    assert_almost_equal(rmg.calc_tru_cr() / tmp_tru_cr, 1.0)

def test_deltaR1():
    rmg.calc_deltaR()
    tmp_deltaR = rmg.batch_average(rmg.target_BU, "p") - rmg.batch_average(rmg.target_BU, "d")
    assert_almost_equal(rmg.deltaR / tmp_deltaR, 1.0)

def test_deltaR2():
    rmg.calc_deltaR({922350: 0.5, 922380: 0.5, 80160: 0.125})
    tmp_deltaR = rmg.batch_average(rmg.target_BU, "p") - rmg.batch_average(rmg.target_BU, "d")
    assert_almost_equal(rmg.deltaR / tmp_deltaR, 1.0)

def test_deltaR3():
    ms = MassStream({922350: 0.5, 922380: 0.5, 80160: 0.125})
    rmg.calc_deltaR(ms)
    tmp_deltaR = rmg.batch_average(rmg.target_BU, "p") - rmg.batch_average(rmg.target_BU, "d")
    assert_almost_equal(rmg.deltaR / tmp_deltaR, 1.0)



class TestReactorMGBurnupMethods(TestCase):
"Tests that the ReactorMG burnup methods work."


@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.calc()

@classmethod
def teardown_class():
    general_teardown()

def test_fluence_at_BU():
    fp = rmg.fluence_at_BU(80.0)
    assert(0 <= fp.f)
    assert(fp.f <= len(rmg.F))
    assert(rmg.F[fp.f] <= fp.F)
    assert(fp.F <= rmg.F[fp.f+1])
    assert(rmg.BU_F_[fp.f] <= 80.0)
    assert(80.0 <= rmg.BU_F_[fp.f+1])
    tmp_m = (rmg.BU_F_[fp.f+1] - rmg.BU_F_[fp.f]) / (rmg.F[fp.f+1] - rmg.F[fp.f])
    assert_equal(fp.m / tmp_m, 1.0)

def test_batch_average():
    BUd = rmg.BUd
    p = rmg.batch_average(BUd, "P")
    d = rmg.batch_average(BUd, "D")
    k = rmg.batch_average(BUd, "K")
    kk = rmg.batch_average(BUd)
    assert_equal(k, kk)
    #assert_equal(p/d, k) # Averaging messes this up.

def test_batch_average_k():
    BUd = rmg.BUd
    assert_almost_equal(rmg.batch_average_k(BUd), 1.0)

def test_calc_1():
    rmg.calc()
    assert(rmg.mat_prod.mass < 1.0)
    assert(rmg.mat_prod.comp[922350] < 0.5) 

def test_calc_2():
    rmg.calc(MassStream({942390: 0.05, 922380: 0.95}))
    assert(rmg.mat_prod.mass < 1.0)
    assert(rmg.mat_prod.comp[942390] < 1.0) 



class TestReactorMGBurnupMethods2(TestCase):
"Tests that the ReactorMG burnup methods work."


@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.calc()

@classmethod
def teardown_class():
    general_teardown()

def test_BUd_bisection_method():
    assert_almost_equal(rmg.k, 1.0, 5)
    rmg.B = 1
    rmg.BUd_bisection_method()
    assert_almost_equal(rmg.k, 1.0, 5)
    rmg.B = 3


class TestReactorMGBurnupMethods3(TestCase):
"Tests that the ReactorMG burnup methods work."


@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.calc()

@classmethod
def teardown_class():
    general_teardown()

def test_run_P_NL():
    # Convergence is not gaurenteed!
    rmg.r_fuelun_P_NL(0.99)
    assert_equal(rmg.P_NL, 0.99)
    assert_almost_equal(rmg.k, 1.0, 1)
    rmg.r_fuelun_P_NL(0.98)



class TestReactorMGBurnupMethods4(TestCase):
"Tests that the ReactorMG burnup methods work.


@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.calc()

@classmethod
def teardown_class():
    general_teardown()

def test_calibrate_P_NL_to_BUd():
    rmg.calibrate_P_NL_to_BUd()
    assert_not_equal(rmg.P_NL, 0.98)
    assert_almost_equal(rmg.BUd / rmg.target_BU, 1.0, 5)
    


class TestReactorMGLatticeMethods(TestCase):
"Tests that the ReactorMG burnup methods work. hese are not exposed to Python directly =("


@classmethod
def setup_class():
    libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
    bright.load_track_nucs_hdf5(libfile)
    rmg = ReactorMG(rp=default_rp, n='rmg')
    rmg.loadlib(libfile)
    rmg.mat_feed = MassStream({922350: 0.5, 922380: 0.5})
    rmg.calc()

@classmethod
def teardown_class():
    general_teardown()

def test_lattice_E_planar():
    prev = rmg.lattice_E_F_
    rmg.lattice_flag = "Planar"
    rmg.r_fuel = 0.5
    rmg.pitch = 1.0
    rmg.fold_mass_weights()
    curr = rmg.lattice_E_F_
    for f in range(len(rmg.F)):
        assert_not_equal(prev[f], curr[f])

def test_lattice_F_planar():
    prev = rmg.lattice_F_F_
    rmg.lattice_flag = "Planar"
    rmg.r_fuel = 0.5
    rmg.pitch = 1.0
    rmg.fold_mass_weights()
    curr = rmg.lattice_F_F_
    for f in range(len(rmg.F)):
        assert_not_equal(prev[f], curr[f])

def test_lattice_E_spherical():
    prev = rmg.lattice_E_F_
    rmg.lattice_flag = "Spherical"
    rmg.r_fuel = 0.5
    rmg.pitch = 1.0
    rmg.fold_mass_weights()
    curr = rmg.lattice_E_F_
    for f in range(len(rmg.F)):
        assert_not_equal(prev[f], curr[f])

def test_lattice_F_spherical():
    prev = rmg.lattice_F_F_
    rmg.lattice_flag = "Spherical"
    rmg.r_fuel = 0.5
    rmg.pitch = 1.0
    rmg.fold_mass_weights()
    curr = rmg.lattice_F_F_
    for f in range(len(rmg.F)):
        assert_not_equal(prev[f], curr[f])

def test_lattice_E_cylindrical():
    prev = rmg.lattice_E_F_
    rmg.lattice_flag = "Cylindrical"
    rmg.r_fuel = 0.5
    rmg.pitch = 1.0
    rmg.fold_mass_weights()
    curr = rmg.lattice_E_F_
    for f in range(len(rmg.F)):
        assert_not_equal(prev[f], curr[f])

def test_lattice_F_cylindrical():
    prev = rmg.lattice_F_F_
    rmg.lattice_flag = "Cylindrical"
    rmg.r_fuel = 0.5
    rmg.pitch = 1.0
    rmg.fold_mass_weights()
    curr = rmg.lattice_F_F_
    for f in range(len(rmg.F)):
        assert_not_equal(prev[f], curr[f])

# Since the above are not exposed directly, 
# They implicitly test the following ReactorMG functions:
#   calc_zeta()
#   calc_zeta_planar()
#   calc_zeta_spherical()
#   calc_zeta_cylindrical()

"""

if __name__ == "__main__":
    nose.main()
