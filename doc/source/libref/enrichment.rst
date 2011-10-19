.. _bright_enrichment:

**********
Enrichment
**********
Bright handles enrichment cascades in a very general manner.  The Enrichment 
object supports normal two isotope component mixture enrichment as well as 
multi-component cascades.  

Multi-component enrichment poses a special modeling issue because the problem 
is underdetermined by an order equal to the number of nuclides in the cascade 
minus two.  (Thus, two-component enrichment is *exactly* determined.) 
Therefore, an N-2 dimensional surface of possible solutions exists. 
The 'correct' solution is found by optimizing :math:`M^*`, the mass at which
the overall stage separation factor is unity.  This model is largely based off 
of the work of A. de la Garza, E. von Halle, and H. Wood.

All functionality may be found in the ``enrichment`` module::

    import bright.enrichment


.. currentmodule:: bright.enrichment

====================
The Enrichment Class
====================

.. autoclass:: Enrichment(enrich_params=None, name="")


    .. autoattribute:: mat_tail
    .. attribute:: track_params

        For Enrichment, the parameters that are automatically tracked are members of the following list: 
        ``["MassFeed", "MassProduct", "MassTails", "N", "M", "Mstar", "TotalPerFeed", "SWUperFeed", "SWUperProduct"]``.

    .. autoattribute:: alpha_0
    .. autoattribute:: Mstar_0
    .. autoattribute:: Mstar
    .. autoattribute:: j
    .. autoattribute:: k
    .. autoattribute:: xP_j
    .. autoattribute:: xW_j
    .. autoattribute:: N
    .. autoattribute:: M
    .. autoattribute:: N0
    .. autoattribute:: M0
    .. autoattribute:: TotalPerFeed
    .. autoattribute:: SWUperFeed
    .. autoattribute:: SWUperProduct


    .. automethod:: initialize(enrich_params)
    .. automethod:: calc(input=None)
    .. automethod:: calc_params()
    .. automethod:: PoverF(x_F, x_P, x_W)
    .. automethod:: WoverF(x_F, x_P, x_W)


==============================
The EnrichmentParameters Class
==============================
Enrichment has one major helper class, EnrichmentParameters, which 
aids in specifying the inputs necessary for the multi-component cascade 
model to run.  

.. autoclass:: EnrichmentParameters()

    .. autoattribute:: alpha_0
    .. autoattribute:: Mstar_0
    .. autoattribute:: j
    .. autoattribute:: k
    .. autoattribute:: xP_j
    .. autoattribute:: xW_j
    .. autoattribute:: N0
    .. autoattribute:: M0


===========================
Enrichment Helper Functions
===========================

.. autofunction:: uranium_enrichment_defaults()
