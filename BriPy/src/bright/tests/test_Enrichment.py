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
import mass_stream

Enrichment = BriPy.Enrichment
MassStream = mass_stream.MassStream
EnrichmentParameters = BriPy.EnrichmentParameters
bright_config = BriPy.bright_config

def general_teardown():
    for f in os.listdir('.'):
        if "Isos.txt" in f:
            os.remove(f)
        elif "Params.txt" in f:
            os.remove(f)
        elif f in [".h5", "e.h5", "fuel_cycle.h5"]:
            os.remove(f)

class TestEnrichmentParameters(TestCase):
    """Tests the enrichment parameters helper class."""

    def test_constructor(self):
        ep = EnrichmentParameters()
        assert_almost_equal(ep.alpha_0, 0.0)
        assert_almost_equal(ep.Mstar_0, 0.0)
        assert_equal(ep.j, 0)
        assert_equal(ep.k, 0)
        assert_almost_equal(ep.N0, 0.0)
        assert_almost_equal(ep.M0, 0.0)
        assert_almost_equal(ep.xP_j, 0.0)
        assert_almost_equal(ep.xW_j, 0.0)

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
        assert_equal(e.params2track, set(["M",  "MassFeed", "MassProduct", "MassTails", 
            "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"]))

    def test_Enrichment_2(self):
        e = Enrichment(name="e")
        assert_equal(e.name, 'e')
        assert_equal(e.params2track, set(["M",  "MassFeed", "MassProduct", "MassTails", 
            "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"]))

    def test_Enrichment_3(self):
        ep = BriPy.UraniumEnrichmentDefaults()
        ep.xP_j = 0.1
        e = Enrichment(enrich_params=ep)
        assert_equal(e.name, '')
        assert_equal(e.params2track, set(["M",  "MassFeed", "MassProduct", "MassTails", 
            "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"]))
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
        e = Enrichment(enrich_params=ep, name='e')
        assert_equal(e.name, 'e')
        assert_equal(e.params2track, set(["M",  "MassFeed", "MassProduct", "MassTails", 
            "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"]))
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
        assert_equal(e.params2track, set(["M",  "MassFeed", "MassProduct", "MassTails", 
            "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"]))
        e.params2track = set(["Om nom nom"])
        assert_equal(e.params2track, set(["Om nom nom"]))
                        
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
        bright_config.isos2track = set([922350, 922360, 922380])
        bright_config.verbosity = 0

    @classmethod
    def teardown_class(cls):
        general_teardown()
        bright_config.verbosity = 0

    def test_initialize(self):
        e = Enrichment()
        assert_equal(e.xP_j, 0.05)
        ep = EnrichmentParameters()
        ep.xP_j = 0.1
        e.initialize(ep)
        assert_almost_equal(e.alpha_0, 0.0)
        assert_almost_equal(e.Mstar_0, 0.0)
        assert_equal(e.j, 0)
        assert_equal(e.k, 0)
        assert_almost_equal(e.N0, 0.0)
        assert_almost_equal(e.M0, 0.0)
        assert_equal(e.xP_j, 0.1)
        assert_almost_equal(e.xW_j, 0.0)


    def test_doCalc_1(self):
        e = Enrichment()
        e.IsosIn = MassStream({922350: 0.01, 922380: 0.985, 922360: 0.005})
        e.doCalc()
        assert_almost_equal(e.IsosOut.comp[922350],  0.05,   6) 
        assert_almost_equal(e.IsosTail.comp[922350], 0.0025, 6)

    def test_doCalc_2(self):
        e = Enrichment()
        e.doCalc({922350: 0.01, 922380: 0.985, 922360: 0.005})
        assert_almost_equal(e.IsosOut.comp[922350],  0.05,   6) 
        assert_almost_equal(e.IsosTail.comp[922350], 0.0025, 6)

    def test_doCalc_3(self):
        e = Enrichment()
        ms = MassStream({922350: 0.01, 922380: 0.985, 922360: 0.005})
        e.doCalc(ms)
        assert_almost_equal(e.IsosOut.comp[922350],  0.05,   6) 
        assert_almost_equal(e.IsosTail.comp[922350], 0.0025, 6)

    def test_setParams(self):
        e = Enrichment()
        ms = MassStream({922350: 0.01, 922380: 0.985, 922360: 0.005})
        e.doCalc(ms)
        e.setParams()

        assert_equal(e.ParamsIn["MassFeed"],  e.IsosIn.mass)
        assert_equal(e.ParamsOut["MassFeed"], 0.0)

        assert_equal(e.ParamsIn["MassProduct"],  0.0)
        assert_equal(e.ParamsOut["MassProduct"], e.IsosOut.mass)

        assert_equal(e.ParamsIn["MassTails"],  0.0)
        assert_equal(e.ParamsOut["MassTails"], e.IsosTail.mass)

        assert_equal(e.ParamsIn["N"],  e.N)
        assert_equal(e.ParamsOut["N"], e.N)

        assert_equal(e.ParamsIn["M"],  e.M)
        assert_equal(e.ParamsOut["M"], e.M)

        assert_equal(e.ParamsIn["Mstar"],  e.Mstar)
        assert_equal(e.ParamsOut["Mstar"], e.Mstar)

        assert_equal(e.ParamsIn["TotalPerFeed"],  e.TotalPerFeed)
        assert_equal(e.ParamsOut["TotalPerFeed"], e.TotalPerFeed)

        assert_equal(e.ParamsIn["SWUperFeed"],  e.SWUperFeed)
        assert_equal(e.ParamsOut["SWUperFeed"], 0.0)

        assert_equal(e.ParamsIn["SWUperProduct"],  0.0)
        assert_equal(e.ParamsOut["SWUperProduct"], e.SWUperProduct)


class TestEnrichmentBenchmarks(TestCase):
    """Tests that the fuel cycle component methods work."""

    @classmethod
    def setup_class(cls):
        bright_config.isos2track = set([922320, 922340, 922350, 922360, 922380])
        bright_config.verbosity = 0

    @classmethod
    def teardown_class(cls):
        general_teardown()
        bright_config.verbosity = 0

    def test_SampleFeed(self):
        ep = BriPy.UraniumEnrichmentDefaults()
        ep.xP_j = 0.06
        e = Enrichment(enrich_params=ep, name='e')
        ms = MassStream({
            922320: 1.1 * (10.0**-9),
            922340: 0.00021,
            922350: 0.0092,
            922360: 0.0042,
            922380: 0.9863899989,
            })
        e.doCalc(ms)

        assert_almost_equal(e.IsosOut.comp[922350],  0.06,   5) 
        assert_almost_equal(e.IsosTail.comp[922350], 0.0025, 5)

        assert_almost_equal(e.IsosIn.mass   / 1.0,                 1.0)
        assert_almost_equal(e.IsosOut.mass  / 0.11652173913043479, 1.0)
        assert_almost_equal(e.IsosTail.mass / 0.88347826086956527, 1.0)

        assert_almost_equal(e.N / 26.92681830762287,  1.0, 4)
        assert_almost_equal(e.M / 16.677607227634965, 1.0, 4)

        assert_almost_equal(e.Mstar / 236.52972999999989, 1.0, 5)

        assert_almost_equal(e.TotalPerFeed  / 359.02324257153055,  1.0, 5)
        assert_almost_equal(e.SWUperFeed    / 0.93228086408760513, 1.0, 5)
        assert_almost_equal(e.SWUperProduct / 8.0009178634384011,  1.0, 5)

    def test_NU(self):
        ep = BriPy.UraniumEnrichmentDefaults()
        ep.xP_j = 0.05
        e = Enrichment(enrich_params=ep, name='e')
        ms = MassStream({
            922340: 0.000055,
            922350: 0.00720,
            922380: 0.992745,
            })
        e.doCalc(ms)

        assert_almost_equal(e.IsosOut.comp[922350],  0.05,   5) 
        assert_almost_equal(e.IsosTail.comp[922350], 0.0025, 5)

        assert_almost_equal(e.IsosIn.mass   / 1.0,             1.0)
        assert_almost_equal(e.IsosOut.mass  / 0.0989473684211, 1.0)
        assert_almost_equal(e.IsosTail.mass / 0.901052631579,  1.0)

        assert_almost_equal(e.N / 27.245874769544688, 1.0, 4)
        assert_almost_equal(e.M / 13.420267151618368, 1.0, 4)

        assert_almost_equal(e.Mstar / 236.51480600000014, 1.0, 5)

        assert_almost_equal(e.TotalPerFeed  / 289.94789485489417,  1.0, 5)
        assert_almost_equal(e.SWUperFeed    / 0.76126365178905142, 1.0, 5)
        assert_almost_equal(e.SWUperProduct / 7.6936220127616908,  1.0, 5)

    def test_VISION(self):
        ep = BriPy.UraniumEnrichmentDefaults()
        ep.xP_j = 0.055
        e = Enrichment(enrich_params=ep, name='e')
        ms = MassStream({
            922340: 0.000183963025893197,
            922350: 0.00818576605617839,
            922360: 0.00610641667100979,
            922380: 0.985523854246919,
            })
        e.doCalc(ms)

        assert_almost_equal(e.IsosOut.comp[922350],  0.055,   5) 
        assert_almost_equal(e.IsosTail.comp[922350], 0.0025, 5)

        assert_almost_equal(e.IsosIn.mass   / 1.0,                 1.0)
        assert_almost_equal(e.IsosOut.mass  / 0.10830030583196934, 1.0)
        assert_almost_equal(e.IsosTail.mass / 0.89169969416803063, 1.0)

        assert_almost_equal(e.N / 27.38162850698868, 1.0, 2)
        assert_almost_equal(e.M / 15.09646512546496, 1.0, 2)

        assert_almost_equal(e.Mstar / 236.53026, 1.0, 4)

        assert_almost_equal(e.TotalPerFeed  / 328.39281411626348,  1.0, 4)
        assert_almost_equal(e.SWUperFeed    / 0.8510218268267431,  1.0, 4)
        assert_almost_equal(e.SWUperProduct / 7.8579817507360037,  1.0, 4)


    def test_Tungsten(self):
        """This test comes from 'Multicomponent Isotope Separation in Matched
        Abundance Ratio Cascades Composed of Stages with Large Separation Factors' 
        by E. von Halle, 1987."""

        ep = EnrichmentParameters()
        ep.alpha_0 = 1.16306
        ep.Mstar_0 = 181.3
        ep.j = 741800
        ep.k = 741860
        ep.N0 = 30.0
        ep.M0 = 10.0
        ep.xP_j = 0.5109
        ep.xW_j = 0.00014

        e = Enrichment(enrich_params=ep, name='e')
        ms = MassStream({
            741800: 0.0014, 
            741820: 0.26416, 
            741830: 0.14409, 
            741840: 0.30618, 
            741860: 0.28417,
            })
        e.doCalc(ms)

        assert_almost_equal(e.IsosOut.comp[741800],  0.5109,  5) 
        assert_almost_equal(e.IsosTail.comp[741800], 0.00014, 5)

        assert_almost_equal(e.IsosIn.mass   / 1.0,                   1.0)
        assert_almost_equal(e.IsosOut.mass  / 0.0024669120526274574, 1.0)
        assert_almost_equal(e.IsosTail.mass / 0.99753308794737272,   1.0)

        assert_almost_equal(e.N / 43.557515688533513, 1.0, 2)
        assert_almost_equal(e.M / 11.49556481009056,  1.0, 2)

        assert_almost_equal(e.Mstar / 181.22106000000002, 1.0, 4)

        assert_almost_equal(e.TotalPerFeed  / 96.970110083845768, 1.0, 3)
        assert_almost_equal(e.SWUperFeed    / 2.2218643574439323, 1.0, 3)
        assert_almost_equal(e.SWUperProduct / 900.66622159370058, 1.0, 3)


if __name__ == "__main__":
    nose.main()
