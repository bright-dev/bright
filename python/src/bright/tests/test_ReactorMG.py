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

default_rp = bright.ReactorParameters()
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
default_rp.radius = 0.411
default_rp.pitch = 0.7
default_rp.open_slots = 123
default_rp.total_slots = 180


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
        assert_equal(rmg.fuel_chemical_form, {})
        assert_equal(rmg.coolant_chemical_form, {})
        assert_equal(rmg.rho_fuel, 0.0)
        assert_equal(rmg.rho_cool, 0.0)
        assert_equal(rmg.P_NL, 0.0)
        assert_equal(rmg.target_BU, 0.0)
        assert_false(rmg.use_zeta)
        assert_equal(rmg.lattice_flag, '')
        assert_false(rmg.rescale_hydrogen_xs)
        assert_equal(rmg.r, 0.0)
        assert_equal(rmg.l, 0.0)
        assert_equal(rmg.S_O, 0.0)
        assert_equal(rmg.S_T, 0.0)

    def test_ReactorMG_6(self):
        rp = ReactorParameters()
        rmg = ReactorMG(reactor_parameters=rp, track_params=set(["Mass"]), name="rmg")
        assert_equal(rmg.name, 'rmg')
        assert_equal(rmg.track_params, set(["Mass"]))
        assert_equal(rmg.B, 0)
        assert_almost_equal(rmg.flux, 0.0)
        assert_equal(rmg.fuel_chemical_form, {})
        assert_equal(rmg.coolant_chemical_form, {})
        assert_almost_equal(rmg.rho_fuel, 0.0)
        assert_almost_equal(rmg.rho_cool, 0.0)
        assert_almost_equal(rmg.P_NL, 0.0)
        assert_almost_equal(rmg.target_BU, 0.0)
        assert_false(rmg.use_zeta)
        assert_equal(rmg.lattice_flag, '')
        assert_false(rmg.rescale_hydrogen_xs)
        assert_almost_equal(rmg.r, 0.0)
        assert_almost_equal(rmg.l, 0.0)
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

    def test_fuel_chemical_form(self):
        rmg = ReactorMG()
        rmg.fuel_chemical_form = {"IHM": 1.0, "O16": 2.0}
        assert_equal(rmg.fuel_chemical_form, {"IHM": 1.0, "O16": 2.0})

    def test_coolant_chemical_form(self):
        rmg = ReactorMG()
        rmg.coolant_chemical_form = {"H1": 2.0, "O16": 1.0,
            "B10": 0.199 * 550 * 10.0**-6, 
            "B11": 0.801 * 550 * 10.0**-6}
        assert_equal(rmg.coolant_chemical_form, {"H1": 2.0, "O16": 1.0,
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
        rmg.r = 0.411
        assert_equal(rmg.r, 0.411)

    def test_l(self):
        rmg = ReactorMG()
        rmg.l = 0.7
        assert_equal(rmg.l, 0.7)

    def test_S_O(self):
        rmg = ReactorMG()
        rmg.S_O = 123
        assert_equal(rmg.S_O, 123)

    def test_S_T(self):
        rmg = ReactorMG()
        rmg.S_T = 180
        assert_equal(rmg.S_T, 180)

    def test_VF(self):
        rp = ReactorParameters()
        rp.radius = 0.5
        rp.pitch = 1.0
        rp.open_slots = 0
        rp.total_slots = 1
        rmg = ReactorMG(reactor_parameters=rp, name='rmg')
        assert_almost_equal(rmg.VF, 3.14159265*0.25) 

    def test_VC(self):
        rp = ReactorParameters()
        rp.radius = 0.5
        rp.pitch = 1.0
        rp.open_slots = 0
        rp.total_slots = 1
        rmg = ReactorMG(reactor_parameters=rp, name='rmg')
        assert_almost_equal(rmg.VC, 1.0 - 3.14159265*0.25) 



class TestReactorMGBasicDataAttributes(TestCase):
    "Tests that the ReactorMG basic data attributes work."

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/lwr_mg.h5'
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


    def test_sigma_a(self):
        J = self.rmg.J
        G = self.rmg.G
        sigma_a = self.rmg.sigma_a
        nperturbations = self.rmg.nperturbations
        for iso in sigma_a.keys():
            assert(iso in J)
            assert_equal(sigma_a[iso].shape, (nperturbations, G))


    def test_sigma_s(self):
        J = self.rmg.J
        G = self.rmg.G
        sigma_s = self.rmg.sigma_s
        nperturbations = self.rmg.nperturbations
        for iso in sigma_s.keys():
            assert(iso in J)
            assert_equal(sigma_s[iso].shape, (nperturbations, G))


    def test_sigma_f(self):
        J = self.rmg.J
        G = self.rmg.G
        sigma_f = self.rmg.sigma_f
        nperturbations = self.rmg.nperturbations
        for iso in sigma_f.keys():
            assert(iso in J)
            assert_equal(sigma_f[iso].shape, (nperturbations, G))


    def test_nubar_sigma_f(self):
        J = self.rmg.J
        G = self.rmg.G
        nubar_sigma_f = self.rmg.nubar_sigma_f
        nperturbations = self.rmg.nperturbations
        for iso in nubar_sigma_f.keys():
            assert(iso in J)
            assert_equal(nubar_sigma_f[iso].shape, (nperturbations, G))


    def test_nubar(self):
        J = self.rmg.J
        G = self.rmg.G
        nubar = self.rmg.nubar
        nperturbations = self.rmg.nperturbations
        for iso in nubar.keys():
            assert(iso in J)
            assert_equal(nubar[iso].shape, (nperturbations, G))


    def test_sigma_s_gh(self):
        J = self.rmg.J
        G = self.rmg.G
        sigma_s_gh = self.rmg.sigma_s_gh
        nperturbations = self.rmg.nperturbations
        for iso in sigma_s_gh.keys():
            assert(iso in J)
            assert_equal(sigma_s_gh[iso].shape, (nperturbations, G, G))

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


    def test_pi_F_(self):
        pi_F_ = self.rmg.pi_F_
        for i in pi_F_.keys():
            assert_equal(len(self.rmg.F), len(pi_F_[i]))

    def test_di_F_(self):
        di_F_ = self.rmg.di_F_
        for i in di_F_.keys():
            assert_equal(len(self.rmg.F), len(di_F_[i]))

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
    """\



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
        rel_Vol = (self.rmg.rho_cool * 268.5 * self.rmg.VC) / (self.rmg.rho_fuel * 18.0 * self.rmg.VF)
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
        assert_equal(self.rmg.fuel_chemical_form, {})
        assert_equal(self.rmg.coolant_chemical_form, {})
        assert_equal(self.rmg.rho_fuel, 0.0)
        assert_equal(self.rmg.rho_cool, 0.0)
        assert_equal(self.rmg.P_NL, 0.0)
        assert_equal(self.rmg.target_BU, 0.0)
        assert_false(self.rmg.use_zeta)
        assert_equal(self.rmg.lattice_flag, '')
        assert_false(self.rmg.rescale_hydrogen_xs)
        assert_equal(self.rmg.r, 0.0)
        assert_equal(self.rmg.l, 0.0)
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
        self.rmg.run_P_NL(0.99)
        assert_equal(self.rmg.P_NL, 0.99)
        assert_almost_equal(self.rmg.k, 1.0, 1)
        self.rmg.run_P_NL(0.98)



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
        self.rmg.r = 0.5
        self.rmg.l = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_E_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_F_planar(self):
        prev = self.rmg.lattice_F_F_
        self.rmg.lattice_flag = "Planar"
        self.rmg.r = 0.5
        self.rmg.l = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_F_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_E_spherical(self):
        prev = self.rmg.lattice_E_F_
        self.rmg.lattice_flag = "Spherical"
        self.rmg.r = 0.5
        self.rmg.l = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_E_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_F_spherical(self):
        prev = self.rmg.lattice_F_F_
        self.rmg.lattice_flag = "Spherical"
        self.rmg.r = 0.5
        self.rmg.l = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_F_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_E_cylindrical(self):
        prev = self.rmg.lattice_E_F_
        self.rmg.lattice_flag = "Cylindrical"
        self.rmg.r = 0.5
        self.rmg.l = 1.0
        self.rmg.fold_mass_weights()
        curr = self.rmg.lattice_E_F_
        for f in range(len(self.rmg.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_F_cylindrical(self):
        prev = self.rmg.lattice_F_F_
        self.rmg.lattice_flag = "Cylindrical"
        self.rmg.r = 0.5
        self.rmg.l = 1.0
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
