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

.. autoclass:: Enrichment(ep=None, n="")
    :members:

    .. attribute:: track_params

        For Enrichment, the parameters that are automatically tracked are members of the following list: 
        ``["MassFeed", "MassProduct", "MassTails", "N", "M", "Mstar", "TotalPerFeed", "SWUperFeed", "SWUperProduct"]``.

