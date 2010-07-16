****************
Enrichment Class
****************
Bright handles enrichment cascades in a very general maaner.  The Enrichment 
object supports normal two isotope component mixture enrichment as well as 
multi-component cascades.  

Multi-component enrichment poses a special modeling issue because the problem 
is underdetermined by an order equal to the number of nuclides in the cascade 
minus two.  (Thus, two-component enrichment is *exactly* determined.) 
Therefore, an N-2 dimensional surface of possible solutions exists. 
The 'correct' solution is found by optimizing :math:`M^*`, the mass at which
the overall stage separation factor is unity.
This is largely based off of the work of A. de la Garza, E. von Halle, 
and H. Wood.

.. currentmodule:: BriPy
    
.. class:: Enrichment([enrich_params[, name]])

    Enrichment Fuel Cycle Component Class.  Daughter of :class:`BriPy.FCComp` class.

    Args:
        * `enrich_params` (EnrichmentParameters): This specifies how the enrichment 
          cascade should be set up.  It is a :class:`BriPy.EnrichmentParameters`
          helper object.  If ``enrich_params`` is not specified, then the cascade 
          is initialized with ``UraniumEnrichmentDefaults``.
        * `name` (str): The name of the reprocessing fuel cycle component instance.

    Note that this automatically calls the public :meth:`initialize` C function.

.. _Enrichment_Attributes:

=====================
Enrichment Attributes
=====================
As a daughter class of :class:`BriPy.FCComp`, :class:`Enrichment` inherits all of 
the attributes of its parent.  The following is a listing of the additional 
attributes specific to this class.

.. attribute:: Enrichment.alpha_0

    The :math:`\alpha_0` attribute specifies the overall stage separation factor 
    for the cascade.  This should be set on initialization.  Values should be 
    greater than one.  Values less than one represent de-enrichment.

.. attribute:: Enrichment.Mstar_0

    The :math:`M^*_0` represents a first guess at what the :attr:`Mstar` is.  
    The value of :attr:`Mstar_0` on initialization should be in the ballpark 
    of the optimized result of :attr:`Mstar`.  However, :math:`M^*_0` must 
    always have a value between the wieghts of the :attr:`j` and :attr:`k`
    key components.

.. attribute:: Enrichment.Mstar

    The :math:`M^*` attribute represents the mass for which the adjusted 
    stage separation factor, :math:`\alpha^*_i`, is equal to one.  It is this
    value that is varied to achieve an optimized enrichment cascade.
     
.. attribute:: Enrichment.j

    This is an integer in zzaaam-form that represents the jth key component.
    This is the nuclide is prefferentially enriched in the product stream.
    For standard uranium cascades :attr:`j` is ``922350``, or U-235.

.. attribute:: Enrichment.k

    This is an integer in zzaaam-form that represents the kth key component.
    This is the nuclide is prefferentially enriched in the waste stream.
    For standard uranium cascades :attr:`k` is ``922380``, or U-238.

.. attribute:: Enrichment.IsosTail

    In addition to the :attr:`IsosIn <BriPy.FCComp.IsosIn>` and 
    :attr:`IsosOut <BriPy.FCComp.IsosOut>` mass streams, :class:`Enrichment`
    also has a tails or waste stream that is represented by this attribute.
    The mass of this stream and the :attr:`IsosOut <BriPy.FCComp.IsosOut>`
    product stream should always add up to the mass of the 
    :attr:`IsosIn <BriPy.FCComp.IsosIn>` feed stream.

.. attribute:: Enrichment.params2track

    For :class:`Enrichment`, the parameters that are automatically tracked are 
    members of the following list: 
    ``["MassFeed", "MassProduct", "MassTails", "N", "M", "Mstar", "TotalPerFeed", "SWUperFeed", "SWUperProduct"]``.

.. _Enrichment_Methods:

==================
Enrichment Methods
==================

.. method:: Enrichment.doCalc([input])

    This method performs an optimization calculation on :math:`M^*` and solves for 
    appropriate values for all :class:`Enrichment` attributes.  This includes the 
    product and waste streams flowing out of the the cascade as well.

    Args:
        * `input` (dict or MassStream): If input is present, it is set as the component's 
          :attr:`IsosIn <BriPy.FCComp.IsosIn>`.  If input is a isotopic dictionary 
          (zzaaam keys, float values), this dictionary is first converted into a MassStream 
          before being set as :attr:`IsosIn <BriPy.FCComp.IsosIn>`.

    Returns:
        * `output` (MassStream): :attr:`IsosOut <BriPy.FCComp.IsosOut>`.


.. method:: Enrichment.initialize(enrich_params)

    The :meth:`initialize` function takes an enrichment parameter object and sets
    the cooresponding :class:`Enrichment` attributes to the same value.

    Args:
        * `enrich_params` (EnrichmentParameters): A class containing the values to
          (re-)initialize an :class:`Enrichment` cascade with.


.. method:: Enrichment.setParams()

    Here the parameters for :class:`Enrichment` are set::

        self.ParamsIn["MassFeed"]  = self.IsosIn.mass
        self.ParamsOut["MassFeed"] = 0.0

        self.ParamsIn["MassProduct"]  = 0.0
        self.ParamsOut["MassProduct"] = self.IsosOut.mass

        self.ParamsIn["MassTails"]  = 0.0
        self.ParamsOut["MassTails"] = self.IsosTail.mass

        self.ParamsIn["N"]  = self.N
        self.ParamsOut["N"] = self.N

        self.ParamsIn["M"]  = self.M
        self.ParamsOut["M"] = self.M

        self.ParamsIn["Mstar"]  = self.Mstar
        self.ParamsOut["Mstar"] = self.Mstar

        self.ParamsIn["TotalPerFeed"]  = self.TotalPerFeed
        self.ParamsOut["TotalPerFeed"] = self.TotalPerFeed

        self.ParamsIn["SWUperFeed"]  = self.SWUperFeed
        self.ParamsOut["SWUperFeed"] = 0.0

        self.ParamsIn["SWUperProduct"]  = 0.0
        self.ParamsOut["SWUperProduct"] = self.SWUperProduct

