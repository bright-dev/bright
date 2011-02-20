****************
Enrichment Class
****************
Bright handles enrichment cascades in a very general manner.  The Enrichment 
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

.. currentmodule:: bright
    
.. class:: Enrichment([enrich_params[, name]])

    Enrichment Fuel Cycle Component Class.  Daughter of :class:`bright.FCComp` class.

    Args:
        * `enrich_params` (EnrichmentParameters): This specifies how the enrichment 
          cascade should be set up.  It is a :class:`bright.EnrichmentParameters`
          helper object.  If ``enrich_params`` is not specified, then the cascade 
          is initialized with ``UraniumEnrichmentDefaults``.
        * `name` (str): The name of the enrichment fuel cycle component instance.

    Note that this automatically calls the public :meth:`initialize` C function.

.. _Enrichment_Attributes:

=====================
Enrichment Attributes
=====================
As a daughter class of :class:`bright.FCComp`, :class:`Enrichment` inherits all of 
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
    always have a value between the weights of the :attr:`j` and :attr:`k`
    key components.

.. attribute:: Enrichment.Mstar

    The :math:`M^*` attribute represents the mass for which the adjusted 
    stage separation factor, :math:`\alpha^*_i`, is equal to one.  It is this
    value that is varied to achieve an optimized enrichment cascade.
     
.. attribute:: Enrichment.j

    This is an integer in zzaaam-form that represents the jth key component.
    This nuclide is preferentially enriched in the product stream.
    For standard uranium cascades :attr:`j` is ``922350``, or U-235.

.. attribute:: Enrichment.k

    This is an integer in zzaaam-form that represents the kth key component.
    This nuclide is preferentially enriched in the waste stream.
    For standard uranium cascades :attr:`k` is ``922380``, or U-238.

.. attribute:: Enrichment.xP_j

    This is the target enrichment of the :attr:`j`th isotope in the 
    product stream :attr:`ms_prod <bright.FCComp.ms_prod>`.  The
    :math:`x^P_j` value is set by the user at initialization or
    run-time.  For typical uranium vectors, this value is about 
    U-235 = 0.05.

.. attribute:: Enrichment.xW_j

    This is the target enrichment of the :attr:`j`th isotope in the 
    waste stream :attr:`ms_tail`.  The :math:`x^W_j` value is set 
    by the user at initialization or runtime.  For typical uranium vectors, 
    this value is about U-235 = 0.0025.

.. attribute:: Enrichment.N

    This is the number of enriching stages present in an ideal cascade.
    Along with :attr:`Mstar` and :attr:`M`, this number is optimized to 
    ensure that a product enrichment of :attr:`xP_j` is attained.  

.. attribute:: Enrichment.M

    This is the number of stripping stages present in an ideal cascade.
    Along with :attr:`Mstar` and :attr:`N`, this number is optimized to 
    ensure that a waste enrichment of :attr:`xW_j` is attained.  

.. attribute:: Enrichment.N0

    This is the number of enriching stages that is initially guessed
    by the user.

.. attribute:: Enrichment.M0

    This is the number of enriching stages that is initially guessed
    by the user.

.. attribute:: Enrichment.TotalPerFeed

    This represents the total flow rate of the cascade divided by the feed
    flow rate.  As such, it shows the mass of material needed in the 
    cascade to enrich an additional kilogram of feed.  Symbolically, 
    the total flow rate is given as :math:`L` while the feed rate is 
    :math:`F`.  Therefore, this quantity is sometimes seen as ``L-over-F``
    or as ``L/F``.  :attr:`TotalPerFeed` is the value that is minimized to
    form an optimized cascade.

.. attribute:: Enrichment.SWUperFeed

    This value denotes the number of separative work units (SWU) required
    per kg of feed for the specified cascade.   

.. attribute:: Enrichment.SWUperProduct

    This value is the number of separative work units (SWU) required
    to produce 1 [kg] of product in the specified cascade.   

.. attribute:: Enrichment.ms_tail

    In addition to the :attr:`ms_feed <bright.FCComp.ms_feed>` and 
    :attr:`ms_prod <bright.FCComp.ms_prod>` mass streams, :class:`Enrichment`
    also has a tails or waste stream that is represented by this attribute.
    The mass of this stream and the :attr:`ms_prod <bright.FCComp.ms_prod>`
    product stream should always add up to the mass of the 
    :attr:`ms_feed <bright.FCComp.ms_feed>` feed stream.

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
          :attr:`ms_feed <bright.FCComp.ms_feed>`.  If input is a isotopic dictionary 
          (zzaaam keys, float values), this dictionary is first converted into a MassStream 
          before being set as :attr:`ms_feed <bright.FCComp.ms_feed>`.

    Returns:
        * `output` (MassStream): :attr:`ms_prod <bright.FCComp.ms_prod>`.


.. method:: Enrichment.initialize(enrich_params)

    The :meth:`initialize` function takes an enrichment parameter object and sets
    the corresponding :class:`Enrichment` attributes to the same value.

    Args:
        * `enrich_params` (EnrichmentParameters): A class containing the values to
          (re-)initialize an :class:`Enrichment` cascade with.


.. method:: Enrichment.setParams()

    Here the parameters for :class:`Enrichment` are set::

        self.ParamsIn["MassFeed"]  = self.ms_feed.mass
        self.ParamsOut["MassFeed"] = 0.0

        self.ParamsIn["MassProduct"]  = 0.0
        self.ParamsOut["MassProduct"] = self.ms_prod.mass

        self.ParamsIn["MassTails"]  = 0.0
        self.ParamsOut["MassTails"] = self.ms_tail.mass

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



=========================
Enrichment Helper Classes
=========================

Enrichment has one major helper class, :class:`EnrichmentParameters`, that 
aids in specifying the inputs necessary for the multi-component cascade 
model to run.  Additonally, the :func:`UraniumEnrichmentDefaults` is
a sub-class of :class:`EnrichmentParameters` that provides stock values 
for uranium enrichment.


.. class:: EnrichmentParameters()

    This class is a collection of values that mirror the attributes in 
    :class:`Enrichment` that are required for the cascade model to run.
    In C-code this a simple ``struct``.  Like 
    :class:`ReactorParameters <bright.ReactorParameters>`, this class 
    takes no arguments on initialization.  An empty :class:`ErichmentParameters`
    instance has all values (weakly) set to zero. 

.. attribute:: EnrichmentParameters.alpha_0

    The :math:`\alpha_0` attribute specifies the overall stage separation factor 
    for the cascade.  This should be set on initialization.  Values should be 
    greater than one.  Values less than one represent de-enrichment.

.. attribute:: EnrichmentParameters.Mstar_0

    The :math:`M^*_0` represents a first guess at what the :attr:`Mstar` is.  
    The value of :attr:`Mstar_0` on initialization should be in the ballpark 
    of the optimized result of :attr:`Mstar`.  However, :math:`M^*_0` must 
    always have a value between the weights of the :attr:`j` and :attr:`k`
    key components.

.. attribute:: EnrichmentParameters.j

    This is an integer in zzaaam-form that represents the jth key component.
    This nuclide is preferentially enriched in the product stream.
    For standard uranium cascades :attr:`j` is ``922350``, or U-235.

.. attribute:: EnrichmentParameters.k

    This is an integer in zzaaam-form that represents the kth key component.
    This nuclide is preferentially enriched in the waste stream.
    For standard uranium cascades :attr:`k` is ``922380``, or U-238.

.. attribute:: EnrichmentParameters.xP_j

    This is the target enrichment of the :attr:`j`th isotope in the 
    product stream :attr:`ms_prod <bright.FCComp.ms_prod>`.  The
    :math:`x^P_j` value is set by the user at initialization or
    run-time.  For typical uranium vectors, this value is about 
    U-235 = 0.05.

.. attribute:: EnrichmentParameters.xW_j

    This is the target enrichment of the :attr:`j`th isotope in the 
    waste stream :attr:`ms_tail`.  The :math:`x^W_j` value is set 
    by the user at initialization or runtime.  For typical uranium vectors, 
    this value is about U-235 = 0.0025.

.. attribute:: EnrichmentParameters.N0

    This is the number of enriching stages that is initially guessed
    by the user.

.. attribute:: EnrichmentParameters.M0

    This is the number of enriching stages that is initially guessed
    by the user.


.. function:: UraniumEnrichmentDefaults()

    This function returns a deep copy of a UraniumEnrichmentDefaults
    class.  This houses values for the initial parameters that make
    it easy to specify a urnaium enrichment cascade.

    The values of this sub-class of :class:`EnrichmentParameters` are as
    follows::

        ued = bright.UraniumEnrichmentDefaults()

        ued.alpha_0 = 1.05
        ued.Mstar_0 = 236.5

        ued.j = 922350
        ued.k = 922380

        ued.xP_j = 0.05
        ued.xW_j = 0.0025

        ued.N0 = 30.0
        ued.M0 = 10.0

