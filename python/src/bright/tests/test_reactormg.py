"""ReactorMG Component and Helper Class tests"""

#import faulthandler
#faulthandler.enable()

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, set_trace

from numpy.testing import assert_array_equal, assert_array_almost_equal

import os
import warnings
import tables as tb
import numpy as np

import bright
import mass_stream
import isoname

bright_config = bright.bright_config
MassStream = mass_stream.MassStream
FluencePoint = bright.FluencePoint
ReactorParameters = bright.ReactorParameters
ReactorMG = bright.ReactorMG

#default_rp = bright.ReactorParameters()
default_rp = bright.lwr_defaults()

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
default_rp.lattice_type = 'Cylindrical'
default_rp.rescale_hydrogen = True

default_rp.burn_times = np.linspace(0.0, 4200.0, 5)
#default_rp.burn_times = np.linspace(0.0, 100.0, 5)

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
default_rp.void_radius = 0.412
default_rp.clad_radius = 0.412
default_rp.unit_cell_pitch = 0.412 * np.sqrt(np.pi)
default_rp.open_slots = 0

# Only Cool
default_rp.fuel_radius = 0.001
default_rp.void_radius = 0.0
default_rp.clad_radius = 0.0
default_rp.unit_cell_pitch = 1.33
default_rp.open_slots = 0



def general_teardown():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "rmg.h5"]:
            os.remove(f)



class TestReactorMGConstructors(TestCase):
    "Tests that the ReactorMG component constructors work."

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_ReactorMG_1(self):
        rmg = ReactorMG()
        assert_equal(rmg.name, '')
        assert_equal(rmg.track_params, set())

    def test_ReactorMG_2(self):
        rmg = ReactorMG(name="rmg")
        assert_equal(rmg.name, 'rmg')
        assert_equal(rmg.track_params, set())

    def test_ReactorMG_3(self):
        rmg = ReactorMG(track_params=set(["Mass"]))
        assert_equal(rmg.name, '')
        assert_equal(rmg.track_params, set(["Mass"]))

    def test_ReactorMG_4(self):
        rmg = ReactorMG(track_params=set(["Mass"]), name='rmg')
        assert_equal(rmg.name, 'rmg')
        assert_equal(rmg.track_params, set(["Mass"]))

    def test_ReactorMG_5(self):
        rp = ReactorParameters()
        rmg = ReactorMG(reactor_parameters=rp, name="rmg")
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

    def test_ReactorMG_6(self):
        rp = ReactorParameters()
        rmg = ReactorMG(reactor_parameters=rp, track_params=set(["Mass"]), name="rmg")
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
    



class TestReactorMGParameterAttributes(TestCase):
    "Tests that the ReactorMG parameter attributes work."

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_B(self):
        rmg = ReactorMG()
        rmg.B = 3
        assert_equal(rmg.B, 3)

    def test_flux(self):
        rmg = ReactorMG()
        rmg.flux = 2*(10**14)
        assert_equal(rmg.flux, 2*(10**14))

    def test_chemical_form_fuel(self):
        rmg = ReactorMG()
        rmg.chemical_form_fuel = {"IHM": 1.0, "O16": 2.0}
        assert_equal(rmg.chemical_form_fuel, {"IHM": 1.0, "O16": 2.0})

    def test_chemical_form_cool(self):
        rmg = ReactorMG()
        rmg.chemical_form_cool = {"H1": 2.0, "O16": 1.0,
            "B10": 0.199 * 550 * 10.0**-6, 
            "B11": 0.801 * 550 * 10.0**-6}
        assert_equal(rmg.chemical_form_cool, {"H1": 2.0, "O16": 1.0,
            "B10": 0.199 * 550 * 10.0**-6,
            "B11": 0.801 * 550 * 10.0**-6})

    def test_rho_fuel(self):
        rmg = ReactorMG()
        rmg.rho_fuel = 10.7
        assert_equal(rmg.rho_fuel, 10.7)

    def test_rho_cool(self):
        rmg = ReactorMG()
        rmg.rho_cool = 0.73
        assert_equal(rmg.rho_cool, 0.73)

    def test_P_NL(self):
        rmg = ReactorMG()
        rmg.P_NL = 0.98
        assert_equal(rmg.P_NL, 0.98)

    def test_target_BU(self):
        rmg = ReactorMG()
        rmg.target_BU = 50.0
        assert_equal(rmg.target_BU, 50.0)

    def test_use_zeta(self):
        rmg = ReactorMG()
        rmg.use_zeta = True
        assert_true(rmg.use_zeta)

    def test_Lattice(self):
        rmg = ReactorMG()
        rmg.lattice_flag = 'Spherical'
        assert_equal(rmg.lattice_flag, 'Spherical')

    def test_rescale_hydrogen_xs(self):
        rmg = ReactorMG()
        rmg.rescale_hydrogen_xs = True
        assert_true(rmg.rescale_hydrogen_xs)

    def test_r(self):
        rmg = ReactorMG()
        rmg.r_fuel = 0.411
        assert_equal(rmg.r_fuel, 0.411)

    def test_l(self):
        rmg = ReactorMG()
        rmg.pitch = 0.7
        assert_equal(rmg.pitch, 0.7)

    def test_S_O(self):
        rmg = ReactorMG()
        rmg.S_O = 123
        assert_equal(rmg.S_O, 123)

    def test_S_T(self):
        rmg = ReactorMG()
        rmg.S_T = 180
        assert_equal(rmg.S_T, 180)

    def test_V_fuel(self):
        rp = ReactorParameters()
        rp.fuel_radius = 0.5
        rp.unit_cell_pitch = 1.0
        rp.open_slots = 0
        rp.total_slots = 1
        rmg = ReactorMG(reactor_parameters=rp, name='rmg')
        assert_almost_equal(rmg.V_fuel, 3.14159265*0.25) 

    def test_V_cool(self):
        rp = ReactorParameters()
        rp.fuel_radius = 0.5
        rp.void_radius = 0.5
        rp.clad_radius = 0.5
        rp.unit_cell_pitch = 1.0
        rp.open_slots = 0
        rp.total_slots = 1
        rmg = ReactorMG(reactor_parameters=rp, name='rmg')
        assert_almost_equal(rmg.V_cool, 1.0 - 3.14159265*0.25) 



class TestReactorMGBasicDataAttributes(TestCase):
    "Tests that the ReactorMG basic data attributes work."

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + 'lwr_mg.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG()
        cls.rmg.loadlib(libfile)

    @classmethod
    def teardown_class(cls):
        general_teardown()


    def test_libfile(self):
        self.rmg.libfile = "It's Ruth!"
        assert_equal(self.rmg.libfile, "It's Ruth!")


    def test_npertubations(self):
        assert(0 < self.rmg.nperturbations)


    def test_perturbed_fields(self):
        pf = self.rmg.perturbed_fields
        for key, value in pf.items():
            assert_equal(len(value), 3)
            assert_equal(value[1] - value[0], value[2])


    def test_I(self):
        assert_not_equal(len(self.rmg.I), 0)
        for iso in self.rmg.I:
            assert_equal(isoname.CurrentForm(iso), 'zzaaam')


    def test_J(self):
        assert_not_equal(len(self.rmg.J), 0)
        for iso in self.rmg.J:
            assert_equal(isoname.CurrentForm(iso), 'zzaaam')


#    def test_IJ(self):
#        # Bad dataset
#        print self.rmg.J
#        assert(self.rmg.I <= self.rmg.J)


    def test_G(self):
        assert(0 < self.rmg.G)


    def test_E_g(self):
        assert_equal(len(self.rmg.E_g), self.rmg.G+1)
        assert((self.rmg.E_g[1:] < self.rmg.E_g[:-1]).all())


    def test_phi_g(self):
        phi_g = self.rmg.phi_g
        assert_equal(len(phi_g), self.rmg.nperturbations)
        assert_equal(len(phi_g[0]), self.rmg.G)


    def test_phi(self):
        assert_equal(len(self.rmg.phi), self.rmg.nperturbations)
        assert_array_almost_equal(self.rmg.phi_g.sum(axis=1) / self.rmg.phi, np.ones(self.rmg.nperturbations), 5)


    def test_Phi(self):
        assert_equal(len(self.rmg.Phi), self.rmg.nperturbations)


    def test_time0(self):
        assert_equal(len(self.rmg.time0), self.rmg.nperturbations)


    def test_BU0(self):
        assert_equal(len(self.rmg.BU0), self.rmg.nperturbations)


    def test_Ti0(self):
        J = self.rmg.J
        Ti0 = self.rmg.Ti0
        nperturbations = self.rmg.nperturbations
        for iso in Ti0.keys():
            assert(iso in J)
            assert_equal(len(Ti0[iso]), nperturbations)


    def test_sigma_f_pg(self):
        J = self.rmg.J
        G = self.rmg.G
        sigma_f_pg = self.rmg.sigma_f_pg
        nperturbations = self.rmg.nperturbations
        for iso in sigma_f_pg.keys():
            assert(iso in J)
            assert_equal(sigma_f_pg[iso].shape, (nperturbations, G))


    def test_nubar_sigma_f_pg(self):
        J = self.rmg.J
        G = self.rmg.G
        nubar_sigma_f_pg = self.rmg.nubar_sigma_f_pg
        nperturbations = self.rmg.nperturbations
        for iso in nubar_sigma_f_pg.keys():
            assert(iso in J)
            assert_equal(nubar_sigma_f_pg[iso].shape, (nperturbations, G))



    def test_sigma_s_pgh(self):
        J = self.rmg.J
        G = self.rmg.G
        sigma_s_pgh = self.rmg.sigma_s_pgh
        nperturbations = self.rmg.nperturbations
        for iso in sigma_s_pgh.keys():
            assert(iso in J)
            assert_equal(sigma_s_pgh[iso].shape, (nperturbations, G, G))


    def test_decay_matrix(self):
        assert not np.isnan(self.rmg.decay_matrix).any()


    def test_thermal_yield_matrix(self):
        assert not np.isnan(self.rmg.thermal_yield_matrix).any()


    def test_fast_yield_matrix(self):
        assert not np.isnan(self.rmg.fast_yield_matrix).any()


    def test_fission_product_yield_matrix(self):
        assert not np.isnan(self.rmg.fission_product_yield_matrix).any()


    """\

    def test_F(self):
        assert_equal(self.rmg.F[0], 0.0)
        assert(1 < len(self.rmg.F))
        for f in range(1, len(self.rmg.F)):        
            assert(self.rmg.F[f-1] < self.rmg.F[f])

        old_F = self.rmg.F

        self.rmg.F = np.arange(0.0, 10.0)        
        assert_array_equal(self.rmg.F, np.arange(0.0, 10.0))

        self.rmg.F = old_F

    def test_BUi_F_(self):
        BUi_F_ = self.rmg.BUi_F_
        for i in BUi_F_.keys():
            assert_equal(BUi_F_[i][0], 0.0)
            assert_equal(len(self.rmg.F), len(BUi_F_[i]))
            for f in range(1, len(self.rmg.F)):        
                assert(BUi_F_[i][f-1] <= BUi_F_[i][f])


        old_BU = self.rmg.BUi_F_
        self.rmg.BUi_F_ = {1: np.arange(0.0, 10.0)}
        assert_equal(self.rmg.BUi_F_.keys(), [1])
        assert_array_equal(self.rmg.BUi_F_[1], np.arange(0.0, 10.0))
        self.rmg.BUi_F_ = old_BU


    def test_Tij_F_(self):
        Tij_F_ = self.rmg.Tij_F_
        jsos   = bright_config.track_isos
        for i in Tij_F_.keys():
            for j in jsos:
                assert_equal(len(self.rmg.F), len(Tij_F_[i][j]))

        old_T = self.rmg.Tij_F_
        self.rmg.Tij_F_ = {1: {2: np.arange(0.0, 10.0)}}
        assert_equal(self.rmg.Tij_F_.keys(), [1])
        assert_equal(self.rmg.Tij_F_[1].keys(), [2])
        assert_array_equal(self.rmg.Tij_F_[1][2], np.arange(0.0, 10.0))
        self.rmg.Tij_F_ = old_T
    """





class TestReactorMGMutliGroupMethods(TestCase):
    "Tests that the ReactorMG basic data attributes work."

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + 'lwr_mg.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.05, 922380: 0.95})
        cls.rmg.burnup_core()
        cls.rmg.burn_time = 0.0
        cls.rmg.bt_s = 0

    @classmethod
    def teardown_class(cls):
        general_teardown()


    def test_calc_nearest_neighbors(self):
        assert_equal(len(self.rmg.nearest_neighbors), self.rmg.nperturbations)


    def test_interpolate_cross_sections(self):
        # Grab data from C++
        G = self.rmg.G
        J = self.rmg.J
        S = self.rmg.S
        s_f_itg = self.rmg.sigma_f_itg
        ns_f_itg = self.rmg.nubar_sigma_f_itg
        s_s_itgh = self.rmg.sigma_s_itgh

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


    def test_calc_mass_weights(self):
        # FIXME
        
        # Grab data from C++
        G = self.rmg.G
        J = self.rmg.J
        S = self.rmg.S

        assert (0 <= self.rmg.A_HM_t).all()
        assert (0 <= self.rmg.MW_fuel_t).all()
        assert (18 == self.rmg.MW_cool_t).all()
        assert (94 <= self.rmg.MW_clad_t).all()

        for j in J:
            #for s in range(S):
            for s in range(1):
                # Test Fuel atomic density
                assert 0 <= self.rmg.n_fuel_it[j][s] 
                assert self.rmg.n_fuel_it[j][s] <= 3

                # Test cool and clad atomic denisty
                if j == 10010:
                    assert_equal(2.0, self.rmg.n_cool_it[j][s])
                elif j == 80160:
                    assert_equal(1.0, self.rmg.n_cool_it[j][s])
                elif j == 400930:
                    assert_equal(0.5, self.rmg.n_clad_it[j][s])
                elif j == 400950:
                    assert_equal(0.5, self.rmg.n_clad_it[j][s])
                else:
                    assert_equal(0.0, self.rmg.n_cool_it[j][s])
                    assert_equal(0.0, self.rmg.n_clad_it[j][s])

                # Test number densities
                if self.rmg.N_fuel_it[j][s] != 0.0: 
                    assert_almost_equal(1.0, self.rmg.n_fuel_it[j][s] * self.rmg.rho_fuel * 6.0221415E+23 / self.rmg.MW_fuel_t[s] / self.rmg.N_fuel_it[j][s]) 

                if self.rmg.N_clad_it[j][s] != 0.0: 
                    assert_almost_equal(1.0, self.rmg.n_clad_it[j][s] * self.rmg.rho_clad * 6.0221415E+23 / self.rmg.MW_clad_t[s] / self.rmg.N_clad_it[j][s])

                if self.rmg.N_cool_it[j][s] != 0.0: 
                    assert_almost_equal(1.0, self.rmg.n_cool_it[j][s] * self.rmg.rho_cool * 6.0221415E+23 / self.rmg.MW_cool_t[s] / self.rmg.N_cool_it[j][s])



    def test_fold_mass_weights(self):
        # FIXME

        # Grab data from C++
        G = self.rmg.G
        J = self.rmg.J
        S = self.rmg.S
    
        #for s in range(S):
        for s in range(1):
            assert (0 <= self.rmg.Sigma_t_tg[s]).all()
            assert (self.rmg.Sigma_f_tg[s] <= self.rmg.Sigma_t_tg[s]).all()

            assert (1.0 < self.rmg.nubar_Sigma_f_tg[s] / self.rmg.Sigma_f_tg[s]).all()
            assert (self.rmg.nubar_Sigma_f_tg[s] / self.rmg.Sigma_f_tg[s] <= 4.0).all()

            assert (0.0 <= self.rmg.chi_tg[s]).all()
            assert (self.rmg.chi_tg[s] <= 1.0).all()
            assert_almost_equal(1.0, np.sum(self.rmg.chi_tg[s]))

            assert (self.rmg.Sigma_gamma_tg[s] <= self.rmg.Sigma_t_tg[s]).all()
            assert (self.rmg.Sigma_2n_tg[s] <= self.rmg.Sigma_t_tg[s]).all()
            assert (self.rmg.Sigma_3n_tg[s] <= self.rmg.Sigma_t_tg[s]).all()
            assert (self.rmg.Sigma_alpha_tg[s] <= self.rmg.Sigma_t_tg[s]).all()
            assert (self.rmg.Sigma_proton_tg[s] <= self.rmg.Sigma_t_tg[s]).all()
            assert (self.rmg.Sigma_gamma_x_tg[s] <= self.rmg.Sigma_t_tg[s]).all()
            assert (self.rmg.Sigma_2n_x_tg[s] <= self.rmg.Sigma_t_tg[s]).all()


    def test_assemble_multigroup_matrices(self):
        # Grab data from C++
        G = self.rmg.G
        J = self.rmg.J
        S = self.rmg.S

        dia = np.identity(G, dtype=bool)
        ndia = True - dia

        #for s in range(S):
        for s in range(1):
            # Check A diagonal is positive, and off-diagonal elements
            # are negative or zero.  This is a test of the physicality of 
            # The algorithm
            assert (0 < self.rmg.A_tgh[s][dia]).all()
            assert (self.rmg.A_tgh[s][ndia] <= 0).all()

            assert_array_almost_equal(self.rmg.F_tgh[s], 
                (self.rmg.nubar_Sigma_f_tg[s][:, np.newaxis] *  self.rmg.chi_tg[s]) )

            # Test that A_inv is truly the inverso of A
            assert_array_almost_equal(np.dot(self.rmg.A_inv_tgh[s], self.rmg.A_tgh[s]), np.identity(G))

            assert_array_almost_equal(self.rmg.A_inv_F_tgh[s], 
                np.dot(self.rmg.A_inv_tgh[s], self.rmg.F_tgh[s]) )

            assert not np.isnan(self.rmg.T_int_tij[s]).any()
            assert not np.isnan(self.rmg.M_tij[s]).any()



    def test_calc_criticality(self):
        #raise nose.SkipTest
        # Grab data from C++
        G = self.rmg.G
        J = self.rmg.J
        S = self.rmg.S

        #print self.rmg.phi_tg[0]
        #print self.rmg.phi_t
        #print self.rmg.k_t
        #print self.rmg.P_NL

#        print self.rmg.A_inv_F_tgh[0]
#        print 
        s = 0
        k0 = 1.0
        phi0 = np.ones(G, dtype=float)
        phi0[2] = 2.0

        """\
        print "phi0 = ", phi0
        print
#        for i in range(100):
        for i in range(3):
            phi1 = (1.0/(0.98 * k0)) * np.dot(self.rmg.A_inv_F_tgh[s], phi0)
            k1 = k0 * sum(self.rmg.nubar_Sigma_f_tg[s] * phi1) / sum(self.rmg.nubar_Sigma_f_tg[s] * phi0)

            print "k1 = ", k1
            print "phi1 = ", phi1
            print 

            k0 = k1
            phi0 = phi1
        """

#        print self.rmg.phi_tg
#        print 
        print "k_t = ", self.rmg.k_t
        print
        print "phi_tg[s] = ", self.rmg.phi_tg[s]
        print
        print self.rmg.V_fuel
        print self.rmg.V_clad
        print self.rmg.V_cool

        raise TypeError 
        

"""\

class TestReactorMGCalculatedWeightAttributes(TestCase):
    "Tests that the ReactorMG calculated weight attributes work."

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_A_IHM(self):
        assert_almost_equal(self.rmg.A_IHM / 236.5, 1.0, 4)

    def test_MWF(self):
        assert_almost_equal(self.rmg.MWF / 268.5, 1.0, 4)

    def test_MWC(self):
        assert_almost_equal(self.rmg.MWC / 18.0, 1.0, 4)

    def test_niF(self):
        assert_equal(self.rmg.niF[80160],  2.0)
        assert_equal(self.rmg.niF[922350], 0.5)
        assert_equal(self.rmg.niF[922380], 0.5)

    def test_niC(self):
        assert_equal(self.rmg.niC[10010],  2.0)
        assert_equal(self.rmg.niC[80160],  1.0)

    def test_miF(self):
        assert_almost_equal(self.rmg.miF[80160],  2.0 * 16  / 236.5, 4)
        assert_almost_equal(self.rmg.miF[922350], 0.5 * 235 / 236.5, 4)
        assert_almost_equal(self.rmg.miF[922380], 0.5 * 238 / 236.5, 4)

    def test_miC(self):
        #First calculate the relative volume
        rel_Vol = (self.rmg.rho_cool * 268.5 * self.rmg.V_cool) / (self.rmg.rho_fuel * 18.0 * self.rmg.V_fuel)
        assert_almost_equal(self.rmg.miC[10010],  rel_Vol * 2.0 * 1   / 236.5, 4)
        assert_almost_equal(self.rmg.miC[80160],  rel_Vol * 1.0 * 16  / 236.5, 4)
        
    def test_NiF(self):
        assert_almost_equal(self.rmg.NiF[80160]  / (2.0 * self.rmg.rho_fuel * 6.022*(10**23) / 268.5), 1.0, 3)
        assert_almost_equal(self.rmg.NiF[922350] / (0.5 * self.rmg.rho_fuel * 6.022*(10**23) / 268.5), 1.0, 3)
        assert_almost_equal(self.rmg.NiF[922380] / (0.5 * self.rmg.rho_fuel * 6.022*(10**23) / 268.5), 1.0, 3)

    def test_NiC(self):
        assert_almost_equal(self.rmg.NiC[10010] / (2.0 * self.rmg.rho_cool * 6.022*(10**23) / 18.0), 1.0, 3)
        assert_almost_equal(self.rmg.NiC[80160] / (1.0 * self.rmg.rho_cool * 6.022*(10**23) / 18.0), 1.0, 3)


class TestReactorMGCalculatedDataAttributes(TestCase):
    "Tests that the ReactorMG calculated data attributes work."

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_BU_F_(self):
        BU_F_ = self.rmg.BU_F_
        BUi_F_ = self.rmg.BUi_F_
        assert_equal(BU_F_[0], 0.0)
        for f in range(1, len(self.rmg.F)):
            tmp_BU = 0.0
            for i in self.rmg.miF.keys():
                tmp_BU = tmp_BU + (self.rmg.miF[i] * BUi_F_[i][f])
            assert_almost_equal(BU_F_[f] / tmp_BU, 1.0, 4)

    def test_P_F_(self):
        P_F_ = self.rmg.P_F_
        pi_F_ = self.rmg.pi_F_
        for f in range(len(self.rmg.F)):
            tmp_P = 0.0
            for i in self.rmg.miF.keys():
                tmp_P = tmp_P + (self.rmg.P_NL * self.rmg.miF[i] * pi_F_[i][f])
            assert_almost_equal(P_F_[f] / tmp_P, 1.0, 4)

    def test_dF_F_(self):
        dF_F_ = self.rmg.dF_F_
        di_F_ = self.rmg.di_F_
        for f in range(len(self.rmg.F)):
            tmp_dF = 0.0
            for i in self.rmg.miF.keys():
                tmp_dF = tmp_dF + (self.rmg.miF[i] * di_F_[i][f])
            assert_almost_equal(dF_F_[f] / tmp_dF, 1.0, 4)

    def test_dC_F_(self):
        dC_F_ = self.rmg.dC_F_
        di_F_ = self.rmg.di_F_
        zeta_F_ = self.rmg.zeta_F_
        for f in range(len(self.rmg.F)):
            tmp_dC = 0.0
            for i in self.rmg.miC.keys():
                if i == 10010 and self.rmg.rescale_hydrogen_xs:
                    tmp_dC = tmp_dC + (self.rmg.miC[i] * di_F_[i][f] * (1.36927-(0.01119*self.rmg.BU_F_[f])) )
                else:
                    tmp_dC = tmp_dC + (self.rmg.miC[i] * di_F_[i][f])
            if self.rmg.use_zeta:
                tmp_dC = tmp_dC * zeta_F_[f]
            assert_almost_equal(dC_F_[f] / tmp_dC, 1.0, 4)

    def test_D_F_(self):
        D_F_  = self.rmg.D_F_
        dF_F_ = self.rmg.dF_F_
        dC_F_ = self.rmg.dC_F_
        for f in range(len(self.rmg.F)):
            assert_almost_equal(D_F_[f] / (dF_F_[f] + dC_F_[f]), 1.0, 4)

    def test_k_F_(self):
        P_F_ = self.rmg.P_F_
        D_F_ = self.rmg.D_F_
        k_F_ = self.rmg.k_F_
        for f in range(len(self.rmg.F)):
            assert_almost_equal(k_F_[f] / (P_F_[f]/D_F_[f]), 1.0, 4)

    def test_Mj_F_(self):
        Mj_F_  = self.rmg.Mj_F_
        Tij_F_ = self.rmg.Tij_F_
        for j in Mj_F_.keys():
            for f in range(len(self.rmg.F)):
                tmp_Mj = 0.0
                for i in self.rmg.miF.keys():
                    tmp_Mj = tmp_Mj + (self.rmg.miF[i] * Tij_F_[i][j][f])
                assert_almost_equal(Mj_F_[j][f] / tmp_Mj, 1.0, 4)

    def test_zeta_F_(self):
        zeta_F_ = self.rmg.zeta_F_
        for f in range(len(self.rmg.F)):
            assert(1.0 <= zeta_F_[f])


class TestReactorMGDischargeAttributes(TestCase):
    "Tests that the ReactorMG discharge attributes are right."

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_fd(self):
        assert(0 <= self.rmg.fd)
        assert(self.rmg.fd <= len(self.rmg.F))

    def test_Fd(self):
        assert(self.rmg.F[self.rmg.fd] <= self.rmg.Fd)
        assert(self.rmg.Fd <= self.rmg.F[self.rmg.fd+1])

    def test_BUd(self):
        assert(self.rmg.BU_F_[self.rmg.fd] <= self.rmg.BUd)
        assert(self.rmg.BUd <= self.rmg.BU_F_[self.rmg.fd+1])

    def test_k(self):
        assert_almost_equal(self.rmg.k, 1.0)


class TestReactorMGSubStreamAndtru_crAttributes(TestCase):
    "Tests that the ReactorMG sub-stream and transuranic conversion ratio attributes are right."

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.calc()
        cls.rmg.calcSubStreams()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_ms_feed_u(self):
        assert_equal(self.rmg.ms_feed_u.mass, 1.0)

    def test_ms_feed_tru(self):
        assert_equal(self.rmg.ms_feed_tru.mass, 0.0)

    def test_ms_feed_lan(self):
        assert_equal(self.rmg.ms_feed_lan.mass, 0.0)

    def test_ms_feed_act(self):
        assert_equal(self.rmg.ms_feed_act.mass, 1.0)

    def test_ms_prod_u(self):
        assert(self.rmg.ms_prod_u.mass < 1.0)

    def test_ms_prod_tru(self):
        assert(0.0 < self.rmg.ms_prod_tru.mass)

    def test_ms_prod_lan(self):
        assert(0.0 < self.rmg.ms_prod_lan.mass)

    def test_ms_prod_act(self):
        assert(self.rmg.ms_prod_act.mass < 1.0)

    def test_tru_cr(self):
        self.rmg.calc_tru_cr()
        tmp_tru_cr = 1.0 - (self.rmg.ms_feed_tru.mass - self.rmg.ms_prod_tru.mass) / (self.rmg.BUd / 931.46)
        assert_almost_equal(self.rmg.tru_cr / tmp_tru_cr, 1.0)

    def test_deltaR(self):
        self.rmg.calc_deltaR()
        tmp_deltaR = self.rmg.batch_average(self.rmg.target_BU, "p") - self.rmg.batch_average(self.rmg.target_BU, "d")
        assert_almost_equal(self.rmg.deltaR / tmp_deltaR, 1.0)


class TestReactorMGThermalDisadvantageFactorAttributes(TestCase):
    "Tests that the ReactorMG calculated data attributes work."

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()


    def test_SigmaFa_F_(self):
        "Can not do full test do to private sigma_a_therm data member."
        assert_equal(len(self.rmg.SigmaFa_F_), len( self.rmg.F))

    def test_SigmaFtr_F_(self):
        "Can not do full test do to private sigma_a_therm data member."
        assert_equal(len(self.rmg.SigmaFtr_F_), len( self.rmg.F))

    def test_kappaF_F_(self):
        kappaF_F_ = self.rmg.kappaF_F_
        SigmaFa_F_ = self.rmg.SigmaFa_F_
        SigmaFtr_F_ = self.rmg.SigmaFtr_F_
        for f in range(len(self.rmg.F)):
            tmp_kappaF = np.sqrt(3.0 * SigmaFtr_F_[f] * SigmaFa_F_[f])
            assert_almost_equal(kappaF_F_[f] / tmp_kappaF, 1.0)

    def test_SigmaCa_F_(self):
        "Can not do full test do to private sigma_a_therm data member."
        assert_equal(len(self.rmg.SigmaCa_F_), len( self.rmg.F))

    def test_SigmaCtr_F_(self):
        "Can not do full test do to private sigma_a_therm data member."
        assert_equal(len(self.rmg.SigmaCtr_F_), len( self.rmg.F))

    def test_kappaC_F_(self):
        kappaC_F_ = self.rmg.kappaC_F_
        SigmaCa_F_ = self.rmg.SigmaCa_F_
        SigmaCtr_F_ = self.rmg.SigmaCtr_F_
        for f in range(len(self.rmg.F)):
            tmp_kappaC = np.sqrt(3.0 * SigmaCtr_F_[f] * SigmaCa_F_[f])
            assert_almost_equal(kappaC_F_[f] / tmp_kappaC, 1.0)

    #Test lattice_E_F_ here

    #Test lattice_F_F_ here


class TestReactorMGInitializationMethods(TestCase):
    "Tests that the fuel cycle component initialization methods work."

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_initialize(self):
        rp = ReactorParameters()
        self.rmg.initialize(rp)
        assert_equal(self.rmg.B, 0)
        assert_equal(self.rmg.flux, 0.0)
        assert_equal(self.rmg.chemical_form_fuel, {})
        assert_equal(self.rmg.chemical_form_cool, {})
        assert_equal(self.rmg.rho_fuel, 0.0)
        assert_equal(self.rmg.rho_cool, 0.0)
        assert_equal(self.rmg.P_NL, 0.0)
        assert_equal(self.rmg.target_BU, 0.0)
        assert_false(self.rmg.use_zeta)
        assert_equal(self.rmg.lattice_flag, '')
        assert_false(self.rmg.rescale_hydrogen_xs)
        assert_equal(self.rmg.r_fuel, 0.0)
        assert_equal(self.rmg.pitch, 0.0)
        assert_equal(self.rmg.S_O, 0.0)
        assert_equal(self.rmg.S_T, 0.0)

    def test_loadlib(self):
        self.rmg.loadlib(os.getenv("BRIGHT_DATA") + '/FR.h5')
        self.rmg.loadlib(os.getenv("BRIGHT_DATA") + '/LWR.h5')

    def test_fold_mass_weights(self):
        prevkey = self.rmg.miF.keys()
        assert(922380 in self.rmg.miF.keys())
        self.rmg.ms_feed = MassStream({922350: 0.5})
        self.rmg.fold_mass_weights()
        assert(922380 not in self.rmg.miF.keys())

        

class TestReactorMGTransmutationMatrixMethods(TestCase):
    "Tests that the fuel cycle component transmutation matrix methods work."


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_calc_Mj_F_(self):
        self.rmg.calc_Mj_F_()
        # Test below is the same as test_Mj_F_()
        Mj_F_  = self.rmg.Mj_F_
        Tij_F_ = self.rmg.Tij_F_
        for j in Mj_F_.keys():
            for f in range(len(self.rmg.F)):
                tmp_Mj = 0.0
                for i in self.rmg.ms_feed.comp.keys():
                    tmp_Mj = tmp_Mj + (self.rmg.miF[i] * Tij_F_[i][j][f])
                if tmp_Mj == 0.0:
                    assert_almost_equal(Mj_F_[j][f], 0.0, 4)
                else:
                    assert_almost_equal(Mj_F_[j][f] / tmp_Mj, 1.0, 4)

    def test_calc_Mj_Fd_(self):
        self.rmg.BUd_bisection_method()
        self.rmg.calc_Mj_F_()
        self.rmg.calc_Mj_Fd_()
        assert(0.0 < self.rmg.ms_prod.mass)
        assert(self.rmg.ms_prod.mass < 1.0)



class TestReactorMGBasicCalculationMethods(TestCase):
    "Tests that the ReactorMG basic calculation methods work."


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_calc_ms_prod(self):
        self.rmg.BUd_bisection_method()
        self.rmg.calc_ms_prod()
        assert(0.0 < self.rmg.ms_prod.mass)
        assert(self.rmg.ms_prod.mass < 1.0)


    def test_calcSubStreams(self):
        self.rmg.calc()
        self.rmg.calcSubStreams()
        assert_equal(self.rmg.ms_feed_u.mass, 1.0)
        assert_equal(self.rmg.ms_feed_tru.mass, 0.0)
        assert_equal(self.rmg.ms_feed_lan.mass, 0.0)
        assert_equal(self.rmg.ms_feed_act.mass, 1.0)
        assert(self.rmg.ms_prod_u.mass < 1.0)
        assert(0.0 < self.rmg.ms_prod_tru.mass)
        assert(0.0 < self.rmg.ms_prod_lan.mass)
        assert(self.rmg.ms_prod_act.mass < 1.0)

    def test_calc_tru_cr(self):
        self.rmg.calc()
        tmp_tru_cr = 1.0 - (self.rmg.ms_feed_tru.mass - self.rmg.ms_prod_tru.mass) / (self.rmg.BUd / 931.46)
        assert_almost_equal(self.rmg.calc_tru_cr() / tmp_tru_cr, 1.0)

    def test_deltaR1(self):
        self.rmg.calc_deltaR()
        tmp_deltaR = self.rmg.batch_average(self.rmg.target_BU, "p") - self.rmg.batch_average(self.rmg.target_BU, "d")
        assert_almost_equal(self.rmg.deltaR / tmp_deltaR, 1.0)

    def test_deltaR2(self):
        self.rmg.calc_deltaR({922350: 0.5, 922380: 0.5, 80160: 0.125})
        tmp_deltaR = self.rmg.batch_average(self.rmg.target_BU, "p") - self.rmg.batch_average(self.rmg.target_BU, "d")
        assert_almost_equal(self.rmg.deltaR / tmp_deltaR, 1.0)

    def test_deltaR3(self):
        ms = MassStream({922350: 0.5, 922380: 0.5, 80160: 0.125})
        self.rmg.calc_deltaR(ms)
        tmp_deltaR = self.rmg.batch_average(self.rmg.target_BU, "p") - self.rmg.batch_average(self.rmg.target_BU, "d")
        assert_almost_equal(self.rmg.deltaR / tmp_deltaR, 1.0)



class TestReactorMGBurnupMethods(TestCase):
    "Tests that the ReactorMG burnup methods work."


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_fluence_at_BU(self):
        fp = self.rmg.fluence_at_BU(80.0)
        assert(0 <= fp.f)
        assert(fp.f <= len(self.rmg.F))
        assert(self.rmg.F[fp.f] <= fp.F)
        assert(fp.F <= self.rmg.F[fp.f+1])
        assert(self.rmg.BU_F_[fp.f] <= 80.0)
        assert(80.0 <= self.rmg.BU_F_[fp.f+1])
        tmp_m = (self.rmg.BU_F_[fp.f+1] - self.rmg.BU_F_[fp.f]) / (self.rmg.F[fp.f+1] - self.rmg.F[fp.f])
        assert_equal(fp.m / tmp_m, 1.0)

    def test_batch_average(self):
        BUd = self.rmg.BUd
        p = self.rmg.batch_average(BUd, "P")
        d = self.rmg.batch_average(BUd, "D")
        k = self.rmg.batch_average(BUd, "K")
        kk = self.rmg.batch_average(BUd)
        assert_equal(k, kk)
        #assert_equal(p/d, k) # Averaging messes this up.

    def test_batch_average_k(self):
        BUd = self.rmg.BUd
        assert_almost_equal(self.rmg.batch_average_k(BUd), 1.0)

    def test_calc_1(self):
        self.rmg.calc()
        assert(self.rmg.ms_prod.mass < 1.0)
        assert(self.rmg.ms_prod.comp[922350] < 0.5) 

    def test_calc_2(self):
        self.rmg.calc(MassStream({942390: 0.05, 922380: 0.95}))
        assert(self.rmg.ms_prod.mass < 1.0)
        assert(self.rmg.ms_prod.comp[942390] < 1.0) 



class TestReactorMGBurnupMethods2(TestCase):
    "Tests that the ReactorMG burnup methods work."


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_BUd_bisection_method(self):
        assert_almost_equal(self.rmg.k, 1.0, 5)
        self.rmg.B = 1
        self.rmg.BUd_bisection_method()
        assert_almost_equal(self.rmg.k, 1.0, 5)
        self.rmg.B = 3


class TestReactorMGBurnupMethods3(TestCase):
    "Tests that the ReactorMG burnup methods work."


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_run_P_NL(self):
        # Convergence is not gaurenteed!
        self.rmg.r_fuelun_P_NL(0.99)
        assert_equal(self.rmg.P_NL, 0.99)
        assert_almost_equal(self.rmg.k, 1.0, 1)
        self.rmg.r_fuelun_P_NL(0.98)



class TestReactorMGBurnupMethods4(TestCase):
    "Tests that the ReactorMG burnup methods work.


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_calibrate_P_NL_to_BUd(self):
        self.rmg.calibrate_P_NL_to_BUd()
        assert_not_equal(self.rmg.P_NL, 0.98)
        assert_almost_equal(self.rmg.BUd / self.rmg.target_BU, 1.0, 5)
        


class TestReactorMGLatticeMethods(TestCase):
    "Tests that the ReactorMG burnup methods work. hese are not exposed to Python directly =("


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_isos_hdf5(libfile)
        cls.rmg = ReactorMG(reactor_parameters=default_rp, name='rmg')
        cls.rmg.loadlib(libfile)
        cls.rmg.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.rmg.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_lattice_E_planar(self):
        prev = self.rmg.lattice_E_F_
        self.rmg.lattice_flag = "Planar"
        self.rmg.r_fuel = 0.5
        self.rmg.pitch = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_E_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_F_planar(self):
        prev = self.rmg.lattice_F_F_
        self.rmg.lattice_flag = "Planar"
        self.rmg.r_fuel = 0.5
        self.rmg.pitch = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_F_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_E_spherical(self):
        prev = self.rmg.lattice_E_F_
        self.rmg.lattice_flag = "Spherical"
        self.rmg.r_fuel = 0.5
        self.rmg.pitch = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_E_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_F_spherical(self):
        prev = self.rmg.lattice_F_F_
        self.rmg.lattice_flag = "Spherical"
        self.rmg.r_fuel = 0.5
        self.rmg.pitch = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_F_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_E_cylindrical(self):
        prev = self.rmg.lattice_E_F_
        self.rmg.lattice_flag = "Cylindrical"
        self.rmg.r_fuel = 0.5
        self.rmg.pitch = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_E_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_F_cylindrical(self):
        prev = self.rmg.lattice_F_F_
        self.rmg.lattice_flag = "Cylindrical"
        self.rmg.r_fuel = 0.5
        self.rmg.pitch = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_F_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    # Since the above are not exposed directly, 
    # They implicitly test the following ReactorMG functions:
    #   calc_zeta()
    #   calc_zeta_planar()
    #   calc_zeta_spherical()
    #   calc_zeta_cylindrical()

"""\

if __name__ == "__main__":
    nose.main()
