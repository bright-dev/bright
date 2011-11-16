"""Fluence Point Class Tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, set_trace

from numpy.testing import assert_array_equal, assert_array_almost_equal

import os
import warnings
import tables as tb
import numpy as np

import bright.fluence_point

FluencePoint = bright.fluence_point.FluencePoint

#
# Tests the fluence point reactor helper class.
#
def test_constructor():
    fp = FluencePoint()
    assert_equal(fp.f, 0)
    assert_equal(fp.F, 0.0)
    assert_equal(fp.m, 0.0)

def test_f():
    fp = FluencePoint()
    fp.f = 10
    assert_equal(fp.f, 10)

def test_F():
    fp = FluencePoint()
    fp.F = 10.0
    assert_equal(fp.F, 10.0)

def test_m():
    fp = FluencePoint()
    fp.m = 10.0
    assert_equal(fp.m, 10.0)


