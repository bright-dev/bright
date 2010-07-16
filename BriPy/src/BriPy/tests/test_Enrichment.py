"""Enrichment Component tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false

import os
import warnings
import tables as tb
import numpy as np

import BriPy
Enrichment = BriPy.Enrichment
MassStream = BriPy.MassStream
EnrichmentParameters = BriPy.EnrichmentParameters

def general_teardown():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "e.h5"]:
            os.remove(f)

class TestEnrichmentParameters(TestCase):
    """Tests the enrichment parameters helper class."""

    def test_constructor(self):
        ep = EnrichmentParameters()
        assert_equal(ep.alpha_0, 0.0)
        assert_equal(ep.Mstar_0, 0.0)
        assert_equal(ep.j, 0)
        assert_equal(ep.k, 0)
        assert_equal(ep.N0, 0.0)
        assert_equal(ep.M0, 0.0)
        assert_equal(ep.xP_j, 0.0)
        assert_equal(ep.xW_j, 0.0)

    def test_alpha_0(self):
        ep = EnrichmentParameters()
        ep.alpha_0 = 1.05
        assert_equal(ep.alpha_0, 1.05)

    def test_Mstar_0(self):
        ep = EnrichmentParameters()
        ep.Mstar_0 = 236.5
        assert_equal(ep.Mstar_0, 236.5)

    def test_j(self):
        ep = EnrichmentParameters()
        ep.j = 922350
        assert_equal(ep.j, 922350)

    def test_k(self):
        ep = EnrichmentParameters()
        ep.k = 922380
        assert_equal(ep.k, 922380)

    def test_N0(self):
        ep = EnrichmentParameters()
        ep.N0 = 30.0
        assert_equal(ep.N0, 30.0)

    def test_M0(self):
        ep = EnrichmentParameters()
        ep.M0 = 10.0
        assert_equal(ep.M0, 10.0)

    def test_xP_j(self):
        ep = EnrichmentParameters()
        ep.xP_j = 0.05
        assert_equal(ep.xP_j, 0.05)

    def test_xW_j(self):
        ep = EnrichmentParameters()
        ep.xW_j = 0.0025
        assert_equal(ep.xW_j, 0.0025)

    def test_UraniumEnrichmentDefualts(self):
        ep = BriPy.UraniumEnrichmentDefaults()
        assert_equal(ep.alpha_0, 1.05)
        assert_equal(ep.Mstar_0, 236.5)
        assert_equal(ep.j, 922350)
        assert_equal(ep.k, 922380)
        assert_equal(ep.N0, 30.0)
        assert_equal(ep.M0, 10.0)
        assert_equal(ep.xP_j, 0.05)
        assert_equal(ep.xW_j, 0.0025)


class TestEnrichmentConstructors(TestCase):
    """Tests that the Enrichment component constructors work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_Enrichment_1(self):
        e = Enrichment()
        assert_equal(e.name, '')
        assert_equal(e.params2track, ["M",  "MassFeed", "MassProduct", "MassTails", 
            "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"])

    def test_Enrichment_2(self):
        e = Enrichment("e")
        assert_equal(e.name, 'e')
        assert_equal(e.params2track, ["M",  "MassFeed", "MassProduct", "MassTails", 
            "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"])

    def test_Enrichment_3(self):
        ep = BriPy.UraniumEnrichmentDefaults()
        ep.xP_j = 0.1
        e = Enrichment(ep)
        assert_equal(e.name, '')
        assert_equal(e.params2track, ["M",  "MassFeed", "MassProduct", "MassTails", 
            "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"])
        assert_equal(e.alpha_0, 1.05)
        assert_equal(e.Mstar_0, 236.5)
        assert_equal(e.j, 922350)
        assert_equal(e.k, 922380)
        assert_equal(e.N0, 30.0)
        assert_equal(e.M0, 10.0)
        assert_equal(e.xP_j, 0.1)
        assert_equal(e.xW_j, 0.0025)

    def test_Enrichment_4(self):
        ep = BriPy.UraniumEnrichmentDefaults()
        ep.j = 922360
        e = Enrichment(ep, 'e')
        assert_equal(e.name, 'e')
        assert_equal(e.params2track, ["M",  "MassFeed", "MassProduct", "MassTails", 
            "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"])
        assert_equal(e.alpha_0, 1.05)
        assert_equal(e.Mstar_0, 236.5)
        assert_equal(e.j, 922360)
        assert_equal(e.k, 922380)
        assert_equal(e.N0, 30.0)
        assert_equal(e.M0, 10.0)
        assert_equal(e.xP_j, 0.05)
        assert_equal(e.xW_j, 0.0025)


class TestEnrichmentAttributes(TestCase):
    """Tests that enrichment the fuel cycle component attributes work."""

    @classmethod
    def teardown_class(cls):
        general_teardown()

    def test_params2track(self):
        e = Enrichment()
        assert_equal(e.params2track, ["M",  "MassFeed", "MassProduct", "MassTails", 
            "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"])
        e.params2track = ["Om nom nom"]
        assert_equal(e.params2track, ["Om nom nom"])
                        
    def test_alpha_0(self):
        e = Enrichment()
        e.alpha_0 = 1.06
        assert_equal(e.alpha_0, 1.06)

    def test_Mstar_0(self):
        e = Enrichment()
        e.Mstar_0 = 235.5
        assert_equal(e.Mstar_0, 235.5)

    def test_j(self):
        e = Enrichment()
        e.j = 922360
        assert_equal(e.j, 922360)

    def test_k(self):
        e = Enrichment()
        e.k = 922390
        assert_equal(e.k, 922390)

    def test_N0(self):
        e = Enrichment()
        e.N0 = 35.0
        assert_equal(e.N0, 35.0)

    def test_M0(self):
        e = Enrichment()
        e.M0 = 15.0
        assert_equal(e.M0, 15.0)

    def test_xP_j(self):
        e = Enrichment()
        e.xP_j = 0.1
        assert_equal(e.xP_j, 0.1)

    def test_xW_j(self):
        e = Enrichment()
        e.xW_j = 0.005
        assert_equal(e.xW_j, 0.005)

    def test_Mstar(self):
        e = Enrichment()
        e.Mstar = 235.5
        assert_equal(e.Mstar, 235.5)

    def test_IsosIn(self):
        e = Enrichment()
        ms = MassStream({922380: 0.5})        
        e.IsosIn = ms
        assert_equal(e.IsosIn.mass, 0.5)
        assert_equal(e.IsosIn.comp[922380], 1.0)

    def test_IsosOut(self):
        e = Enrichment()
        ms = MassStream({922380: 0.5})        
        e.IsosOut = ms
        assert_equal(e.IsosOut.mass, 0.5)
        assert_equal(e.IsosOut.comp[922380], 1.0)

    def test_IsosTail(self):
        e = Enrichment()
        ms = MassStream({922380: 0.5})        
        e.IsosTail = ms
        assert_equal(e.IsosTail.mass, 0.5)
        assert_equal(e.IsosTail.comp[922380], 1.0)    

    def test_N(self):
        e = Enrichment()
        e.N = 35.0
        assert_equal(e.N, 35.0)

    def test_M0(self):
        e = Enrichment()
        e.M0 = 15.0
        assert_equal(e.M0, 15.0)

    def test_TotalPerFeed(self):
        e = Enrichment()
        e.TotalPerFeed = 12.0
        assert_equal(e.TotalPerFeed, 12.0)

    def test_SWUperFeed(self):
        e = Enrichment()
        e.SWUperFeed = 101.0
        assert_equal(e.SWUperFeed, 101.0)

    def test_SWUperProduct(self):
        e = Enrichment()
        e.SWUperProduct = 1010.0
        assert_equal(e.SWUperProduct, 1010.0)


class TestEnrichmentMethods(TestCase):
    """Tests that the fuel cycle component methods work."""

    @classmethod
    def setup_class(cls):
        BriPy.isos2track([922350, 922360, 922380])
        BriPy.verbosity(0)

    @classmethod
    def teardown_class(cls):
        general_teardown()
        BriPy.verbosity(0)

    def test_initialize(self):
        e = Enrichment()
        assert_equal(e.xP_j, 0.05)
        ep = EnrichmentParameters()
        ep.xP_j = 0.1
        e.initialize(ep)
        assert_equal(e.alpha_0, 0.0)
        assert_equal(e.Mstar_0, 0.0)
        assert_equal(e.j, 0)
        assert_equal(e.k, 0)
        assert_equal(e.N0, 0.0)
        assert_equal(e.M0, 0.0)
        assert_equal(e.xP_j, 0.1)
        assert_equal(e.xW_j, 0.0)

    def test_doCalc_1(self):
        e = Enrichment()
        e.doCalc({922350: 0.01, 922380: 0.985, 922360: 0.005})
#        assert_equal(e.IsosOut.mass, 0.99)
#        assert_equal(e.IsosOut.comp[942390], 1.0) # Recall ms.comp is normalized

"""
    def test_doCalc_2(self):
        BriPy.isos2track([922350, 922380, 942390])
        e = Enrichment({"U235": 0.9, "922380": 0.999, "94239": 0.99})
        e.doCalc(MassStream({942390: 1.0}))
        assert_equal(e.IsosOut.mass, 0.99)
        assert_equal(e.IsosOut.comp[942390], 1.0) # Recall ms.comp is normalized

"""        

if __name__ == "__main__":
    nose.main()
