"""Bright general tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false

import os
import warnings
import tables as tb
import numpy as np
import BriPy


class TestBright(TestCase):
    """Tests that the Bright general functions work."""

    def test_BrightStart(self):
        current = os.getenv("BRIGHT_DATA")
        os.environ["BRIGHT_DATA"] = "/foo/bar"
        new = os.getenv("BRIGHT_DATA")
        BriPy.BrightStart()
        assert_equal(new, "/foo/bar")
        os.environ["BRIGHT_DATA"] = current

    def test_isos2track(self):
        old_isolist = BriPy.isos2track()
        new_isolist = BriPy.mixed_2_zzaaam_List([92235, "H1"])
        BriPy.isos2track(new_isolist)
        assert_equal(BriPy.isos2track(), [10010, 922350])
        BriPy.isos2track(old_isolist)

    def test_verbosity(self):
        old_verbosity = BriPy.verbosity()
        BriPy.verbosity(100)
        assert_equal(BriPy.verbosity(), 100)
        BriPy.verbosity(old_verbosity)

    def test_write_hdf5(self):
        old_write = BriPy.write_hdf5()
        BriPy.write_hdf5(False)
        assert_false(BriPy.write_hdf5())
        BriPy.write_hdf5(1)
        assert_true(BriPy.write_hdf5())
        BriPy.write_hdf5(old_write)

    def test_write_text(self):
        old_write = BriPy.write_text()
        BriPy.write_text(False)
        assert_false(BriPy.write_text())
        BriPy.write_text(1)
        assert_true(BriPy.write_text())
        BriPy.write_text(old_write)
        
    def test_output_filename(self):
        assert_equal( BriPy.output_filename(), 'fuel_cycle.h5')
        BriPy.output_filename('new_name.h5')
        assert_equal( BriPy.output_filename(), 'new_name.h5')
        

class TestLoadFromHDF5(TestCase):
    """Tests isos2track can be loaded from an HDF5 file."""

    @classmethod
    def setup_class(cls):
        f = tb.openFile('isos.h5', 'w')
        f.createArray(f.root, "ToIsos", np.array([92235, 922380, 10010]), "ToIsos")
        f.createArray(f.root, "NotIsos", np.array([92235, 922380, 10010]), "NotIsos")
        f.close()

    @classmethod
    def teardown_class(cls):
        os.remove('isos.h5')

    def test_load_isos2track_hdf5_1(self):
        old_isos = BriPy.isos2track()
        BriPy.isos2track([80160])
        BriPy.load_isos2track_hdf5('isos.h5')
        assert_equal(BriPy.isos2track(), [10010, 80160, 922350, 922380])
        BriPy.isos2track(old_isos)

    def test_load_isos2track_hdf5_2(self):
        old_isos = BriPy.isos2track()
        BriPy.isos2track([80160])
        BriPy.load_isos2track_hdf5('isos.h5', '/NotIsos')
        assert_equal(BriPy.isos2track(), [10010, 80160, 922350, 922380])
        BriPy.isos2track(old_isos)

    def test_load_isos2track_hdf5_3(self):
        old_isos = BriPy.isos2track()
        BriPy.isos2track([80160])
        BriPy.load_isos2track_hdf5('isos.h5', '', True)
        assert_equal(BriPy.isos2track(), [10010, 922350, 922380])
        BriPy.isos2track(old_isos)

    def test_load_isos2track_hdf5_4(self):
        old_isos = BriPy.isos2track()
        BriPy.isos2track([80160])
        BriPy.load_isos2track_hdf5('isos.h5', '/NotIsos', True)
        assert_equal(BriPy.isos2track(), [10010, 922350, 922380])
        BriPy.isos2track(old_isos)

class TestLoadFromText(TestCase):
    """Tests isos2track can be loaded from a text file."""

    @classmethod
    def setup_class(cls):
        with open('isos.txt', 'w') as f:
            f.write('U-235, 922380\n10010}')

    @classmethod
    def teardown_class(cls):
        os.remove('isos.txt')

    def test_load_isos2track_text_1(self):
        old_isos = BriPy.isos2track()
        BriPy.isos2track([80160])
        BriPy.load_isos2track_text('isos.txt')
        assert_equal(BriPy.isos2track(), [10010, 80160, 922350, 922380])
        BriPy.isos2track(old_isos)

    def test_load_isos2track_text_2(self):
        old_isos = BriPy.isos2track()
        BriPy.isos2track([80160])
        BriPy.load_isos2track_text('isos.txt', True)
        assert_equal(BriPy.isos2track(), [10010, 922350, 922380])
        BriPy.isos2track(old_isos)

if __name__ == "__main__":
    nose.main()
