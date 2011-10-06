"""Reactor1G Component and Helper Class tests"""

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

bright_config = bright.bright_config
MassStream = mass_stream.MassStream
FluencePoint = bright.FluencePoint
ReactorParameters = bright.ReactorParameters
Reactor1G = bright.Reactor1G

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
default_rp.fuel_radius = 0.411
default_rp.unit_cell_pitch = 0.7
default_rp.open_slots = 123
default_rp.total_slots = 180


def general_teardown():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "r1g.h5"]:
            os.remove(f)


class TestFluencePoint(TestCase):
    """Tests the fluence point Reactor1G helper clasr1g."""

    def test_constructor(self):
        fp = FluencePoint()
        assert_equal(fp.f, 0)
        assert_equal(fp.F, 0.0)
        assert_equal(fp.m, 0.0)

    def test_f(self):
        fp = FluencePoint()
        fp.f = 10
        assert_equal(fp.f, 10)

    def test_F(self):
        fp = FluencePoint()
        fp.F = 10.0
        assert_equal(fp.F, 10.0)

    def test_m(self):
        fp = FluencePoint()
        fp.m = 10.0
        assert_equal(fp.m, 10.0)


class TestReactorParameters(TestCase):
    """Tests the reactor parameters Reactor1G helper clasr1g."""

    def test_constructor(self):
        rp = ReactorParameters()
        assert_equal(rp.batches, 0)
        assert_almost_equal(rp.flux, 0.0)
        assert_equal(rp.fuel_form, {})
        assert_equal(rp.coolant_form, {})
        assert_almost_equal(rp.fuel_density, 0.0)
        assert_almost_equal(rp.coolant_density, 0.0)
        assert_almost_equal(rp.pnl, 0.0)
        assert_almost_equal(rp.BUt, 0.0)
        assert_false(rp.use_disadvantage_factor)
        assert_equal(rp.lattice_type, '')
        assert_false(rp.rescale_hydrogen)
        assert_almost_equal(rp.fuel_radius, 0.0)
        assert_almost_equal(rp.unit_cell_pitch, 0.0)
        assert_almost_equal(rp.open_slots, 0.0)
        assert_almost_equal(rp.total_slots, 0.0)

    def test_batches(self):
        rp = ReactorParameters()
        rp.batches = 3
        assert_equal(rp.batches, 3)

    def test_flux(self):
        rp = ReactorParameters()
        rp.flux = 2*(10**14)
        assert_equal(rp.flux, 2*(10**14))

    def test_fuel_form(self):
        rp = ReactorParameters()
        rp.fuel_form = {"IHM": 1.0, "O16": 2.0}
        assert_equal(rp.fuel_form, {"IHM": 1.0, "O16": 2.0})

    def test_coolant_form(self):
        rp = ReactorParameters()
        rp.coolant_form = {"H1": 2.0, "O16": 1.0,
            "B10": 0.199 * 550 * 10.0**-6, 
            "B11": 0.801 * 550 * 10.0**-6}
        assert_equal(rp.coolant_form, {"H1": 2.0, "O16": 1.0,
            "B10": 0.199 * 550 * 10.0**-6,
            "B11": 0.801 * 550 * 10.0**-6})

    def test_fuel_density(self):
        rp = ReactorParameters()
        rp.fuel_density = 10.7
        assert_equal(rp.fuel_density, 10.7)

    def test_coolant_density(self):
        rp = ReactorParameters()
        rp.coolant_density = 0.73
        assert_equal(rp.coolant_density, 0.73)

    def test_pnl(self):
        rp = ReactorParameters()
        rp.pnl = 0.98
        assert_equal(rp.pnl, 0.98)

    def test_BUt(self):
        rp = ReactorParameters()
        rp.BUt = 50.0
        assert_equal(rp.BUt, 50.0)

    def test_use_disadvantage_factor(self):
        rp = ReactorParameters()
        rp.use_disadvantage_factor = True
        assert_true(rp.use_disadvantage_factor)

    def test_lattice_type(self):
        rp = ReactorParameters()
        rp.lattice_type = 'Spherical'
        assert_equal(rp.lattice_type, 'Spherical')

    def test_rescale_hydrogen(self):
        rp = ReactorParameters()
        rp.rescale_hydrogen = True
        assert_true(rp.rescale_hydrogen)

    def test_fuel_radius(self):
        rp = ReactorParameters()
        rp.fuel_radius = 0.411
        assert_equal(rp.fuel_radius, 0.411)

    def test_unit_cell_pitch(self):
        rp = ReactorParameters()
        rp.unit_cell_pitch = 0.7
        assert_equal(rp.unit_cell_pitch, 0.7)

    def test_open_slots(self):
        rp = ReactorParameters()
        rp.open_slots = 123
        assert_equal(rp.open_slots, 123)

    def test_total_slots(self):
        rp = ReactorParameters()
        rp.total_slots = 180
        assert_equal(rp.total_slots, 180)


class TestReactor1GConstructors(TestCase):
    """Tests that the Reactor1G component constructors work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_Reactor1G_1(self):
        r1g = Reactor1G()
        assert_equal(r1g.name, '')
        assert_equal(r1g.track_params, set())

    def test_Reactor1G_2(self):
        r1g = Reactor1G(name="r1g")
        assert_equal(r1g.name, 'r1g')
        assert_equal(r1g.track_params, set())

    def test_Reactor1G_3(self):
        r1g = Reactor1G(track_params=set(["Mass"]))
        assert_equal(r1g.name, '')
        assert_equal(r1g.track_params, set(["Mass"]))

    def test_Reactor1G_4(self):
        r1g = Reactor1G(track_params=set(["Mass"]), name='r1g')
        assert_equal(r1g.name, 'r1g')
        assert_equal(r1g.track_params, set(["Mass"]))

    def test_Reactor1G_5(self):
        rp = ReactorParameters()
        r1g = Reactor1G(reactor_parameters=rp, name="r1g")
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

    def test_Reactor1G_6(self):
        rp = ReactorParameters()
        r1g = Reactor1G(reactor_parameters=rp, track_params=set(["Mass"]), name="r1g")
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
    


class TestReactor1GParameterAttributes(TestCase):
    """Tests that the Reactor1G parameter attributes work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_B(self):
        r1g = Reactor1G()
        r1g.B = 3
        assert_equal(r1g.B, 3)

    def test_phi(self):
        r1g = Reactor1G()
        r1g.phi = 2*(10**14)
        assert_equal(r1g.phi, 2*(10**14))

    def test_fuel_chemical_form(self):
        r1g = Reactor1G()
        r1g.fuel_chemical_form = {"IHM": 1.0, "O16": 2.0}
        assert_equal(r1g.fuel_chemical_form, {"IHM": 1.0, "O16": 2.0})

    def test_coolant_chemical_form(self):
        r1g = Reactor1G()
        r1g.coolant_chemical_form = {"H1": 2.0, "O16": 1.0,
            "B10": 0.199 * 550 * 10.0**-6, 
            "B11": 0.801 * 550 * 10.0**-6}
        assert_equal(r1g.coolant_chemical_form, {"H1": 2.0, "O16": 1.0,
            "B10": 0.199 * 550 * 10.0**-6,
            "B11": 0.801 * 550 * 10.0**-6})

    def test_rhoF(self):
        r1g = Reactor1G()
        r1g.rhoF = 10.7
        assert_equal(r1g.rhoF, 10.7)

    def test_rhoC(self):
        r1g = Reactor1G()
        r1g.rhoC = 0.73
        assert_equal(r1g.rhoC, 0.73)

    def test_P_NL(self):
        r1g = Reactor1G()
        r1g.P_NL = 0.98
        assert_equal(r1g.P_NL, 0.98)

    def test_target_BU(self):
        r1g = Reactor1G()
        r1g.target_BU = 50.0
        assert_equal(r1g.target_BU, 50.0)

    def test_use_zeta(self):
        r1g = Reactor1G()
        r1g.use_zeta = True
        assert_true(r1g.use_zeta)

    def test_Lattice(self):
        r1g = Reactor1G()
        r1g.lattice_flag = 'Spherical'
        assert_equal(r1g.lattice_flag, 'Spherical')

    def test_rescale_hydrogen_xs(self):
        r1g = Reactor1G()
        r1g.rescale_hydrogen_xs = True
        assert_true(r1g.rescale_hydrogen_xs)

    def test_r(self):
        r1g = Reactor1G()
        r1g.r = 0.411
        assert_equal(r1g.r, 0.411)

    def test_l(self):
        r1g = Reactor1G()
        r1g.l = 0.7
        assert_equal(r1g.l, 0.7)

    def test_S_O(self):
        r1g = Reactor1G()
        r1g.S_O = 123
        assert_equal(r1g.S_O, 123)

    def test_S_T(self):
        r1g = Reactor1G()
        r1g.S_T = 180
        assert_equal(r1g.S_T, 180)

    def test_VF(self):
        rp = ReactorParameters()
        rp.fuel_radius = 0.5
        rp.unit_cell_pitch = 1.0
        rp.open_slots = 0
        rp.total_slots = 1
        r1g = Reactor1G(reactor_parameters=rp, name='r1g')
        assert_almost_equal(r1g.VF, 3.14159265*0.25) 

    def test_VC(self):
        rp = ReactorParameters()
        rp.fuel_radius = 0.5
        rp.unit_cell_pitch = 1.0
        rp.open_slots = 0
        rp.total_slots = 1
        r1g = Reactor1G(reactor_parameters=rp, name='r1g')
        assert_almost_equal(r1g.VC, 1.0 - 3.14159265*0.25) 


class TestReactor1GBasicDataAttributes(TestCase):
    """Tests that the Reactor1G basic data attributes work."""

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G()
        cls.r1g.loadlib(libfile)

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_libfile(self):
        self.r1g.libfile = "It's Ruth!"
        assert_equal(self.r1g.libfile, "It's Ruth!")

    def test_F(self):
        assert_equal(self.r1g.F[0], 0.0)
        assert(1 < len(self.r1g.F))
        for f in range(1, len(self.r1g.F)):        
            assert(self.r1g.F[f-1] < self.r1g.F[f])

        old_F = self.r1g.F

        self.r1g.F = np.arange(0.0, 10.0)        
        assert_array_equal(self.r1g.F, np.arange(0.0, 10.0))

        self.r1g.F = old_F

    def test_BUi_F_(self):
        BUi_F_ = self.r1g.BUi_F_
        for i in BUi_F_.keys():
            assert_equal(BUi_F_[i][0], 0.0)
            assert_equal(len(self.r1g.F), len(BUi_F_[i]))
            for f in range(1, len(self.r1g.F)):        
                assert(BUi_F_[i][f-1] <= BUi_F_[i][f])


        old_BU = self.r1g.BUi_F_
        self.r1g.BUi_F_ = {1: np.arange(0.0, 10.0)}
        assert_equal(self.r1g.BUi_F_.keys(), [1])
        assert_array_equal(self.r1g.BUi_F_[1], np.arange(0.0, 10.0))
        self.r1g.BUi_F_ = old_BU


    def test_pi_F_(self):
        pi_F_ = self.r1g.pi_F_
        for i in pi_F_.keys():
            assert_equal(len(self.r1g.F), len(pi_F_[i]))

    def test_di_F_(self):
        di_F_ = self.r1g.di_F_
        for i in di_F_.keys():
            assert_equal(len(self.r1g.F), len(di_F_[i]))

    def test_Tij_F_(self):
        Tij_F_ = self.r1g.Tij_F_
        jsos   = bright_config.track_nucs
        for i in Tij_F_.keys():
            for j in jsos:
                assert_equal(len(self.r1g.F), len(Tij_F_[i][j]))

        old_T = self.r1g.Tij_F_
        self.r1g.Tij_F_ = {1: {2: np.arange(0.0, 10.0)}}
        assert_equal(self.r1g.Tij_F_.keys(), [1])
        assert_equal(self.r1g.Tij_F_[1].keys(), [2])
        assert_array_equal(self.r1g.Tij_F_[1][2], np.arange(0.0, 10.0))
        self.r1g.Tij_F_ = old_T

class TestReactor1GCalculatedWeightAttributes(TestCase):
    """Tests that the Reactor1G calculated weight attributes work."""

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_A_IHM(self):
        assert_almost_equal(self.r1g.A_IHM / 236.5, 1.0, 4)

    def test_MWF(self):
        assert_almost_equal(self.r1g.MWF / 268.5, 1.0, 4)

    def test_MWC(self):
        assert_almost_equal(self.r1g.MWC / 18.0, 1.0, 4)

    def test_niF(self):
        assert_equal(self.r1g.niF[80160],  2.0)
        assert_equal(self.r1g.niF[922350], 0.5)
        assert_equal(self.r1g.niF[922380], 0.5)

    def test_niC(self):
        assert_equal(self.r1g.niC[10010],  2.0)
        assert_equal(self.r1g.niC[80160],  1.0)

    def test_miF(self):
        assert_almost_equal(self.r1g.miF[80160],  2.0 * 16  / 236.5, 4)
        assert_almost_equal(self.r1g.miF[922350], 0.5 * 235 / 236.5, 4)
        assert_almost_equal(self.r1g.miF[922380], 0.5 * 238 / 236.5, 4)

    def test_miC(self):
        #First calculate the relative volume
        rel_Vol = (self.r1g.rhoC * 268.5 * self.r1g.VC) / (self.r1g.rhoF * 18.0 * self.r1g.VF)
        assert_almost_equal(self.r1g.miC[10010],  rel_Vol * 2.0 * 1   / 236.5, 4)
        assert_almost_equal(self.r1g.miC[80160],  rel_Vol * 1.0 * 16  / 236.5, 4)
        
    def test_NiF(self):
        assert_almost_equal(self.r1g.NiF[80160]  / (2.0 * self.r1g.rhoF * 6.022*(10**23) / 268.5), 1.0, 3)
        assert_almost_equal(self.r1g.NiF[922350] / (0.5 * self.r1g.rhoF * 6.022*(10**23) / 268.5), 1.0, 3)
        assert_almost_equal(self.r1g.NiF[922380] / (0.5 * self.r1g.rhoF * 6.022*(10**23) / 268.5), 1.0, 3)

    def test_NiC(self):
        assert_almost_equal(self.r1g.NiC[10010] / (2.0 * self.r1g.rhoC * 6.022*(10**23) / 18.0), 1.0, 3)
        assert_almost_equal(self.r1g.NiC[80160] / (1.0 * self.r1g.rhoC * 6.022*(10**23) / 18.0), 1.0, 3)


class TestReactor1GCalculatedDataAttributes(TestCase):
    """Tests that the Reactor1G calculated data attributes work."""

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_BU_F_(self):
        BU_F_ = self.r1g.BU_F_
        BUi_F_ = self.r1g.BUi_F_
        assert_equal(BU_F_[0], 0.0)
        for f in range(1, len(self.r1g.F)):
            tmp_BU = 0.0
            for i in self.r1g.miF.keys():
                tmp_BU = tmp_BU + (self.r1g.miF[i] * BUi_F_[i][f])
            assert_almost_equal(BU_F_[f] / tmp_BU, 1.0, 4)

    def test_P_F_(self):
        P_F_ = self.r1g.P_F_
        pi_F_ = self.r1g.pi_F_
        for f in range(len(self.r1g.F)):
            tmp_P = 0.0
            for i in self.r1g.miF.keys():
                tmp_P = tmp_P + (self.r1g.P_NL * self.r1g.miF[i] * pi_F_[i][f])
            assert_almost_equal(P_F_[f] / tmp_P, 1.0, 4)

    def test_dF_F_(self):
        dF_F_ = self.r1g.dF_F_
        di_F_ = self.r1g.di_F_
        for f in range(len(self.r1g.F)):
            tmp_dF = 0.0
            for i in self.r1g.miF.keys():
                tmp_dF = tmp_dF + (self.r1g.miF[i] * di_F_[i][f])
            assert_almost_equal(dF_F_[f] / tmp_dF, 1.0, 4)

    def test_dC_F_(self):
        dC_F_ = self.r1g.dC_F_
        di_F_ = self.r1g.di_F_
        zeta_F_ = self.r1g.zeta_F_
        for f in range(len(self.r1g.F)):
            tmp_dC = 0.0
            for i in self.r1g.miC.keys():
                if i == 10010 and self.r1g.rescale_hydrogen_xs:
                    tmp_dC = tmp_dC + (self.r1g.miC[i] * di_F_[i][f] * (1.36927-(0.01119*self.r1g.BU_F_[f])) )
                else:
                    tmp_dC = tmp_dC + (self.r1g.miC[i] * di_F_[i][f])
            if self.r1g.use_zeta:
                tmp_dC = tmp_dC * zeta_F_[f]
            assert_almost_equal(dC_F_[f] / tmp_dC, 1.0, 4)

    def test_D_F_(self):
        D_F_  = self.r1g.D_F_
        dF_F_ = self.r1g.dF_F_
        dC_F_ = self.r1g.dC_F_
        for f in range(len(self.r1g.F)):
            assert_almost_equal(D_F_[f] / (dF_F_[f] + dC_F_[f]), 1.0, 4)

    def test_k_F_(self):
        P_F_ = self.r1g.P_F_
        D_F_ = self.r1g.D_F_
        k_F_ = self.r1g.k_F_
        for f in range(len(self.r1g.F)):
            assert_almost_equal(k_F_[f] / (P_F_[f]/D_F_[f]), 1.0, 4)

    def test_Mj_F_(self):
        Mj_F_  = self.r1g.Mj_F_
        Tij_F_ = self.r1g.Tij_F_
        for j in Mj_F_.keys():
            for f in range(len(self.r1g.F)):
                tmp_Mj = 0.0
                for i in self.r1g.miF.keys():
                    tmp_Mj = tmp_Mj + (self.r1g.miF[i] * Tij_F_[i][j][f])
                assert_almost_equal(Mj_F_[j][f] / tmp_Mj, 1.0, 4)

    def test_zeta_F_(self):
        zeta_F_ = self.r1g.zeta_F_
        for f in range(len(self.r1g.F)):
            assert(1.0 <= zeta_F_[f])


class TestReactor1GDischargeAttributes(TestCase):
    """Tests that the Reactor1G discharge attributes are right."""

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_fd(self):
        assert(0 <= self.r1g.fd)
        assert(self.r1g.fd <= len(self.r1g.F))

    def test_Fd(self):
        assert(self.r1g.F[self.r1g.fd] <= self.r1g.Fd)
        assert(self.r1g.Fd <= self.r1g.F[self.r1g.fd+1])

    def test_BUd(self):
        assert(self.r1g.BU_F_[self.r1g.fd] <= self.r1g.BUd)
        assert(self.r1g.BUd <= self.r1g.BU_F_[self.r1g.fd+1])

    def test_k(self):
        assert_almost_equal(self.r1g.k, 1.0)


class TestReactor1GSubStreamAndtru_crAttributes(TestCase):
    """Tests that the Reactor1G sub-stream and transuranic conversion ratio
    attributes are right."""

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.calc()
        cls.r1g.calcSubStreams()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_ms_feed_u(self):
        assert_equal(self.r1g.ms_feed_u.mass, 1.0)

    def test_ms_feed_tru(self):
        assert_equal(self.r1g.ms_feed_tru.mass, 0.0)

    def test_ms_feed_lan(self):
        assert_equal(self.r1g.ms_feed_lan.mass, 0.0)

    def test_ms_feed_act(self):
        assert_equal(self.r1g.ms_feed_act.mass, 1.0)

    def test_ms_prod_u(self):
        assert(self.r1g.ms_prod_u.mass < 1.0)

    def test_ms_prod_tru(self):
        assert(0.0 < self.r1g.ms_prod_tru.mass)

    def test_ms_prod_lan(self):
        assert(0.0 < self.r1g.ms_prod_lan.mass)

    def test_ms_prod_act(self):
        assert(self.r1g.ms_prod_act.mass < 1.0)

    def test_tru_cr(self):
        self.r1g.calc_tru_cr()
        tmp_tru_cr = 1.0 - (self.r1g.ms_feed_tru.mass - self.r1g.ms_prod_tru.mass) / (self.r1g.BUd / 931.46)
        assert_almost_equal(self.r1g.tru_cr / tmp_tru_cr, 1.0)

    def test_deltaR(self):
        self.r1g.calc_deltaR()
        tmp_deltaR = self.r1g.batch_average(self.r1g.target_BU, "p") - self.r1g.batch_average(self.r1g.target_BU, "d")
        assert_almost_equal(self.r1g.deltaR / tmp_deltaR, 1.0)


class TestReactor1GThermalDisadvantageFactorAttributes(TestCase):
    """Tests that the Reactor1G calculated data attributes work."""

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()


    def test_SigmaFa_F_(self):
        """Can not do full test do to private sigma_a_therm data member."""
        assert_equal(len(self.r1g.SigmaFa_F_), len( self.r1g.F))

    def test_SigmaFtr_F_(self):
        """Can not do full test do to private sigma_a_therm data member."""
        assert_equal(len(self.r1g.SigmaFtr_F_), len( self.r1g.F))

    def test_kappaF_F_(self):
        kappaF_F_ = self.r1g.kappaF_F_
        SigmaFa_F_ = self.r1g.SigmaFa_F_
        SigmaFtr_F_ = self.r1g.SigmaFtr_F_
        for f in range(len(self.r1g.F)):
            tmp_kappaF = np.sqrt(3.0 * SigmaFtr_F_[f] * SigmaFa_F_[f])
            assert_almost_equal(kappaF_F_[f] / tmp_kappaF, 1.0)

    def test_SigmaCa_F_(self):
        """Can not do full test do to private sigma_a_therm data member."""
        assert_equal(len(self.r1g.SigmaCa_F_), len( self.r1g.F))

    def test_SigmaCtr_F_(self):
        """Can not do full test do to private sigma_a_therm data member."""
        assert_equal(len(self.r1g.SigmaCtr_F_), len( self.r1g.F))

    def test_kappaC_F_(self):
        kappaC_F_ = self.r1g.kappaC_F_
        SigmaCa_F_ = self.r1g.SigmaCa_F_
        SigmaCtr_F_ = self.r1g.SigmaCtr_F_
        for f in range(len(self.r1g.F)):
            tmp_kappaC = np.sqrt(3.0 * SigmaCtr_F_[f] * SigmaCa_F_[f])
            assert_almost_equal(kappaC_F_[f] / tmp_kappaC, 1.0)

    #Test lattice_E_F_ here

    #Test lattice_F_F_ here


class TestReactor1GInitializationMethods(TestCase):
    """Tests that the fuel cycle component initialization methods work."""

    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_initialize(self):
        rp = ReactorParameters()
        self.r1g.initialize(rp)
        assert_equal(self.r1g.B, 0)
        assert_equal(self.r1g.phi, 0.0)
        assert_equal(self.r1g.fuel_chemical_form, {})
        assert_equal(self.r1g.coolant_chemical_form, {})
        assert_equal(self.r1g.rhoF, 0.0)
        assert_equal(self.r1g.rhoC, 0.0)
        assert_equal(self.r1g.P_NL, 0.0)
        assert_equal(self.r1g.target_BU, 0.0)
        assert_false(self.r1g.use_zeta)
        assert_equal(self.r1g.lattice_flag, '')
        assert_false(self.r1g.rescale_hydrogen_xs)
        assert_equal(self.r1g.r, 0.0)
        assert_equal(self.r1g.l, 0.0)
        assert_equal(self.r1g.S_O, 0.0)
        assert_equal(self.r1g.S_T, 0.0)

    def test_loadlib(self):
        self.r1g.loadlib(os.getenv("BRIGHT_DATA") + '/FR.h5')
        self.r1g.loadlib(os.getenv("BRIGHT_DATA") + '/LWR.h5')

    def test_fold_mass_weights(self):
        prevkey = self.r1g.miF.keys()
        assert(922380 in self.r1g.miF.keys())
        self.r1g.ms_feed = MassStream({922350: 0.5})
        self.r1g.fold_mass_weights()
        assert(922380 not in self.r1g.miF.keys())

        

class TestReactor1GTransmutationMatrixMethods(TestCase):
    """Tests that the fuel cycle component transmutation matrix methods work."""


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_calc_Mj_F_(self):
        self.r1g.calc_Mj_F_()
        # Test below is the same as test_Mj_F_()
        Mj_F_  = self.r1g.Mj_F_
        Tij_F_ = self.r1g.Tij_F_
        for j in Mj_F_.keys():
            for f in range(len(self.r1g.F)):
                tmp_Mj = 0.0
                for i in self.r1g.ms_feed.comp.keys():
                    tmp_Mj = tmp_Mj + (self.r1g.miF[i] * Tij_F_[i][j][f])
                if tmp_Mj == 0.0:
                    assert_almost_equal(Mj_F_[j][f], 0.0, 4)
                else:
                    assert_almost_equal(Mj_F_[j][f] / tmp_Mj, 1.0, 4)

    def test_calc_Mj_Fd_(self):
        self.r1g.BUd_bisection_method()
        self.r1g.calc_Mj_F_()
        self.r1g.calc_Mj_Fd_()
        assert(0.0 < self.r1g.ms_prod.mass)
        assert(self.r1g.ms_prod.mass < 1.0)



class TestReactor1GBasicCalculationMethods(TestCase):
    """Tests that the Reactor1G basic calculation methods work."""


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.fold_mass_weights()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_calc_ms_prod(self):
        self.r1g.BUd_bisection_method()
        self.r1g.calc_ms_prod()
        assert(0.0 < self.r1g.ms_prod.mass)
        assert(self.r1g.ms_prod.mass < 1.0)


    def test_calcSubStreams(self):
        self.r1g.calc()
        self.r1g.calcSubStreams()
        assert_equal(self.r1g.ms_feed_u.mass, 1.0)
        assert_equal(self.r1g.ms_feed_tru.mass, 0.0)
        assert_equal(self.r1g.ms_feed_lan.mass, 0.0)
        assert_equal(self.r1g.ms_feed_act.mass, 1.0)
        assert(self.r1g.ms_prod_u.mass < 1.0)
        assert(0.0 < self.r1g.ms_prod_tru.mass)
        assert(0.0 < self.r1g.ms_prod_lan.mass)
        assert(self.r1g.ms_prod_act.mass < 1.0)

    def test_calc_tru_cr(self):
        self.r1g.calc()
        tmp_tru_cr = 1.0 - (self.r1g.ms_feed_tru.mass - self.r1g.ms_prod_tru.mass) / (self.r1g.BUd / 931.46)
        assert_almost_equal(self.r1g.calc_tru_cr() / tmp_tru_cr, 1.0)

    def test_deltaR1(self):
        self.r1g.calc_deltaR()
        tmp_deltaR = self.r1g.batch_average(self.r1g.target_BU, "p") - self.r1g.batch_average(self.r1g.target_BU, "d")
        assert_almost_equal(self.r1g.deltaR / tmp_deltaR, 1.0)

    def test_deltaR2(self):
        self.r1g.calc_deltaR({922350: 0.5, 922380: 0.5, 80160: 0.125})
        tmp_deltaR = self.r1g.batch_average(self.r1g.target_BU, "p") - self.r1g.batch_average(self.r1g.target_BU, "d")
        assert_almost_equal(self.r1g.deltaR / tmp_deltaR, 1.0)

    def test_deltaR3(self):
        ms = MassStream({922350: 0.5, 922380: 0.5, 80160: 0.125})
        self.r1g.calc_deltaR(ms)
        tmp_deltaR = self.r1g.batch_average(self.r1g.target_BU, "p") - self.r1g.batch_average(self.r1g.target_BU, "d")
        assert_almost_equal(self.r1g.deltaR / tmp_deltaR, 1.0)



class TestReactor1GBurnupMethods(TestCase):
    """Tests that the Reactor1G burnup methods work."""


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_fluence_at_BU(self):
        fp = self.r1g.fluence_at_BU(80.0)
        assert(0 <= fp.f)
        assert(fp.f <= len(self.r1g.F))
        assert(self.r1g.F[fp.f] <= fp.F)
        assert(fp.F <= self.r1g.F[fp.f+1])
        assert(self.r1g.BU_F_[fp.f] <= 80.0)
        assert(80.0 <= self.r1g.BU_F_[fp.f+1])
        tmp_m = (self.r1g.BU_F_[fp.f+1] - self.r1g.BU_F_[fp.f]) / (self.r1g.F[fp.f+1] - self.r1g.F[fp.f])
        assert_equal(fp.m / tmp_m, 1.0)

    def test_batch_average(self):
        BUd = self.r1g.BUd
        p = self.r1g.batch_average(BUd, "P")
        d = self.r1g.batch_average(BUd, "D")
        k = self.r1g.batch_average(BUd, "K")
        kk = self.r1g.batch_average(BUd)
        assert_equal(k, kk)
        #assert_equal(p/d, k) # Averaging messes this up.

    def test_batch_average_k(self):
        BUd = self.r1g.BUd
        assert_almost_equal(self.r1g.batch_average_k(BUd), 1.0)

    def test_calc_1(self):
        self.r1g.calc()
        assert(self.r1g.ms_prod.mass < 1.0)
        assert(self.r1g.ms_prod.comp[922350] < 0.5) 

    def test_calc_2(self):
        self.r1g.calc(MassStream({942390: 0.05, 922380: 0.95}))
        assert(self.r1g.ms_prod.mass < 1.0)
        assert(self.r1g.ms_prod.comp[942390] < 1.0) 



class TestReactor1GBurnupMethods2(TestCase):
    """Tests that the Reactor1G burnup methods work."""


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_BUd_bisection_method(self):
        assert_almost_equal(self.r1g.k, 1.0, 5)
        self.r1g.B = 1
        self.r1g.BUd_bisection_method()
        assert_almost_equal(self.r1g.k, 1.0, 5)
        self.r1g.B = 3


class TestReactor1GBurnupMethods3(TestCase):
    """Tests that the Reactor1G burnup methods work."""


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_run_P_NL(self):
        # Convergence is not gaurenteed!
        self.r1g.run_P_NL(0.99)
        assert_equal(self.r1g.P_NL, 0.99)
        assert_almost_equal(self.r1g.k, 1.0, 1)
        self.r1g.run_P_NL(0.98)



class TestReactor1GBurnupMethods4(TestCase):
    """Tests that the Reactor1G burnup methods work."""


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_calibrate_P_NL_to_BUd(self):
        self.r1g.calibrate_P_NL_to_BUd()
        assert_not_equal(self.r1g.P_NL, 0.98)
        assert_almost_equal(self.r1g.BUd / self.r1g.target_BU, 1.0, 5)
        


class TestReactor1GLatticeMethods(TestCase):
    """Tests that the Reactor1G burnup methods work.
    These are not exposed to Python directly =("""


    @classmethod
    def setup_class(cls):
        libfile = os.getenv("BRIGHT_DATA") + '/LWR.h5'
        bright.load_track_nucs_hdf5(libfile)
        cls.r1g = Reactor1G(reactor_parameters=default_rp, name='r1g')
        cls.r1g.loadlib(libfile)
        cls.r1g.ms_feed = MassStream({922350: 0.5, 922380: 0.5})
        cls.r1g.calc()

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_lattice_E_planar(self):
        prev = self.r1g.lattice_E_F_
        self.r1g.lattice_flag = "Planar"
        self.r1g.r = 0.5
        self.r1g.l = 1.0
        self.r1g.fold_mass_weights()
        curr = self.r1g.lattice_E_F_
        for f in range(len(self.r1g.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_F_planar(self):
        prev = self.r1g.lattice_F_F_
        self.r1g.lattice_flag = "Planar"
        self.r1g.r = 0.5
        self.r1g.l = 1.0
        self.r1g.fold_mass_weights()
        curr = self.r1g.lattice_F_F_
        for f in range(len(self.r1g.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_E_spherical(self):
        prev = self.r1g.lattice_E_F_
        self.r1g.lattice_flag = "Spherical"
        self.r1g.r = 0.5
        self.r1g.l = 1.0
        self.r1g.fold_mass_weights()
        curr = self.r1g.lattice_E_F_
        for f in range(len(self.r1g.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_F_spherical(self):
        prev = self.r1g.lattice_F_F_
        self.r1g.lattice_flag = "Spherical"
        self.r1g.r = 0.5
        self.r1g.l = 1.0
        self.r1g.fold_mass_weights()
        curr = self.r1g.lattice_F_F_
        for f in range(len(self.r1g.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_E_cylindrical(self):
        prev = self.r1g.lattice_E_F_
        self.r1g.lattice_flag = "Cylindrical"
        self.r1g.r = 0.5
        self.r1g.l = 1.0
        self.r1g.fold_mass_weights()
        curr = self.r1g.lattice_E_F_
        for f in range(len(self.r1g.F)):
            assert_not_equal(prev[f], curr[f])

    def test_lattice_F_cylindrical(self):
        prev = self.r1g.lattice_F_F_
        self.r1g.lattice_flag = "Cylindrical"
        self.r1g.r = 0.5
        self.r1g.l = 1.0
        self.r1g.fold_mass_weights()
        curr = self.r1g.lattice_F_F_
        for f in range(len(self.r1g.F)):
            assert_not_equal(prev[f], curr[f])

    # Since the above are not exposed directly, 
    # They implicitly test the following Reactor1G functions:
    #   calc_zeta()
    #   calc_zeta_planar()
    #   calc_zeta_spherical()
    #   calc_zeta_cylindrical()


if __name__ == "__main__":
    nose.main()
