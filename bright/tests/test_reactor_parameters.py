"""Reactor Paramters Helper Class Tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, set_trace

from numpy.testing import assert_array_equal, assert_array_almost_equal

import os
import warnings
import tables as tb
import numpy as np

import bright.reactor_parameters
#from pyne import material

ReactorParameters = bright.reactor_parameters.ReactorParameters

#
# Tests the reactor parameters.
#

def test_constructor():
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

def test_batches():
    rp = ReactorParameters()
    rp.batches = 3
    assert_equal(rp.batches, 3)

def test_flux():
    rp = ReactorParameters()
    rp.flux = 2*(10**14)
    assert_equal(rp.flux, 2*(10**14))

def test_fuel_form():
    rp = ReactorParameters()
    rp.fuel_form = {"IHM": 1.0, "O16": 2.0}
    assert_equal(rp.fuel_form, {"IHM": 1.0, "O16": 2.0})

def test_coolant_form():
    rp = ReactorParameters()
    rp.coolant_form = {"H1": 2.0, "O16": 1.0,
                       "B10": 0.199 * 550 * 10.0**-6, 
                       "B11": 0.801 * 550 * 10.0**-6}
    assert_equal(rp.coolant_form, {"H1": 2.0, "O16": 1.0,
                                   "B10": 0.199 * 550 * 10.0**-6,
                                   "B11": 0.801 * 550 * 10.0**-6})

def test_fuel_density():
    rp = ReactorParameters()
    rp.fuel_density = 10.7
    assert_equal(rp.fuel_density, 10.7)

def test_coolant_density():
    rp = ReactorParameters()
    rp.coolant_density = 0.73
    assert_equal(rp.coolant_density, 0.73)

def test_pnl():
    rp = ReactorParameters()
    rp.pnl = 0.98
    assert_equal(rp.pnl, 0.98)

def test_BUt():
    rp = ReactorParameters()
    rp.BUt = 50.0
    assert_equal(rp.BUt, 50.0)

def test_use_disadvantage_factor():
    rp = ReactorParameters()
    rp.use_disadvantage_factor = True
    assert_true(rp.use_disadvantage_factor)

def test_lattice_type():
    rp = ReactorParameters()
    rp.lattice_type = 'Spherical'
    assert_equal(rp.lattice_type, 'Spherical')

def test_rescale_hydrogen():
    rp = ReactorParameters()
    rp.rescale_hydrogen = True
    assert_true(rp.rescale_hydrogen)

def test_fuel_radius():
    rp = ReactorParameters()
    rp.fuel_radius = 0.411
    assert_equal(rp.fuel_radius, 0.411)

def test_unit_cell_pitch():
    rp = ReactorParameters()
    rp.unit_cell_pitch = 0.7
    assert_equal(rp.unit_cell_pitch, 0.7)

def test_open_slots():
    rp = ReactorParameters()
    rp.open_slots = 123
    assert_equal(rp.open_slots, 123)

def test_total_slots():
    rp = ReactorParameters()
    rp.total_slots = 180
    assert_equal(rp.total_slots, 180)


