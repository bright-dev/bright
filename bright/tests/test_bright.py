"""Bright general tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, with_setup

import os
import warnings
import tables as tb
import numpy as np

from pyne import nucname
import bright

bright_conf = bright.bright_conf


#
# Fixtures
#

def setup_h5():
    if 'isos.h5' in os.listdir('.'):
        return
    f = tb.openFile('isos.h5', 'w')
    f.createArray(f.root, "ToIsos", np.array([92235, 922380, 10010]), "ToIsos")
    f.createArray(f.root, "NotIsos", np.array([92235, 922380, 10010]), "NotIsos")
    f.close()

def teardown_h5():
    os.remove('isos.h5')


def setup_txt():
    with open('isos.txt', 'w') as f:
        f.write('U-235, 922380\n10010}')

def teardown_txt():
    os.remove('isos.txt')


#
# Tests
#

def test_bright_start():
    current = os.getenv("BRIGHT_DATA")
    os.environ["BRIGHT_DATA"] = "/foo/bar"
    new = os.getenv("BRIGHT_DATA")
    bright.bright_start()
    assert_equal(new, "/foo/bar")
    os.environ["BRIGHT_DATA"] = current


def test_track_nucs():
    old_isolist = bright_conf.track_nucs
    new_isolist = [922350, 10010]
    bright_conf.track_nucs = set(new_isolist)
    assert_equal(bright_conf.track_nucs, set([10010, 922350]))
    bright_conf.track_nucs = old_isolist

def test_verbosity():
    old_verbosity = bright_conf.verbosity
    bright_conf.verbosity = 100
    assert_equal(bright_conf.verbosity, 100)
    bright.verbosity = old_verbosity

def test_write_hdf5():
    old_write = bright_conf.write_hdf5
    bright_conf.write_hdf5 = False
    assert_false(bright_conf.write_hdf5)
    bright_conf.write_hdf5 = 1
    assert_true(bright_conf.write_hdf5)
    bright_conf.write_hdf5 = old_write

def test_write_text():
    old_write = bright_conf.write_text
    bright_conf.write_text = False
    assert_false(bright_conf.write_text)
    bright_conf.write_text = 1
    assert_true(bright_conf.write_text)
    bright_conf.write_text = old_write
        
def test_output_filename():
    assert_equal( bright_conf.output_filename, 'fuel_cycle.h5')
    bright_conf.output_filename = 'new_name.h5'
    assert_equal( bright_conf.output_filename, 'new_name.h5')
        

@with_setup(setup_h5)
def test_load_track_nucs_hdf5_1():
    old_isos = bright_conf.track_nucs
    bright_conf.track_nucs = set([80160])
    bright.load_track_nucs_hdf5('isos.h5')
    assert_equal(bright_conf.track_nucs, set([10010, 80160, 922350, 922380]))
    bright_conf.track_nucs = old_isos


@with_setup(setup_h5)
def test_load_track_nucs_hdf5_2():
    old_isos = bright_conf.track_nucs
    bright_conf.track_nucs = set([80160])
    bright.load_track_nucs_hdf5('isos.h5', '/NotIsos')
    assert_equal(bright_conf.track_nucs, set([10010, 80160, 922350, 922380]))
    bright_conf.track_nucs = old_isos

@with_setup(setup_h5)
def test_load_track_nucs_hdf5_3():
    old_isos = bright_conf.track_nucs
    bright_conf.track_nucs = set([80160])
    bright.load_track_nucs_hdf5('isos.h5', '', True)
    assert_equal(bright_conf.track_nucs, set([10010, 922350, 922380]))
    bright_conf.track_nucs = old_isos
    
@with_setup(setup_h5, teardown_h5)
def test_load_track_nucs_hdf5_4():
    old_isos = bright_conf.track_nucs
    bright_conf.track_nucs = set([80160])
    bright.load_track_nucs_hdf5('isos.h5', '/NotIsos', True)
    assert_equal(bright_conf.track_nucs, set([10010, 922350, 922380]))
    bright_conf.track_nucs = old_isos


@with_setup(setup_txt)
def test_load_track_nucs_text_1():
    old_isos = bright_conf.track_nucs
    bright_conf.track_nucs = set([80160])
    bright.load_track_nucs_text('isos.txt')
    assert_equal(bright_conf.track_nucs, set([10010, 80160, 922350, 922380]))
    bright_conf.track_nucs = old_isos

@with_setup(setup_txt, teardown_txt)
def test_load_track_nucs_text_2():
    old_isos = bright_conf.track_nucs
    bright_conf.track_nucs = set([80160])
    bright.load_track_nucs_text('isos.txt', True)
    assert_equal(bright_conf.track_nucs, set([10010, 922350, 922380]))
    bright_conf.track_nucs = old_isos


if __name__ == "__main__":
    nose.main()
