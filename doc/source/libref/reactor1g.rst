.. _bright_reactor1g:

***************
Reactor1G Class
***************
Reactors are the most computationally complex of the nuclear fuel cycle components currently implemented.  
Bright handles nuclear reactors in two distinct object classes: a one neutron energy group methodology 
(implemented here) and a multi-group algorithm. Moreover, there are several different types of reactors.  
Each type has its own characteristic data library associated with it and takes on type-specific base case 
values.  The reactor types that have been fully implemented thus far are a light water reactor (LWR) and a 
fast reactor (FR).  You may read more about these on their own pages.

All one-group (1G) reactors share a common methodological backbone.  This page describes what is fundamentally 
the same about one group reactor objects via the Reactor1G class.  This is a subclass of FCComp.  Moreover, all 
one-group reactor types have :class:`bright.Reactor1G` as their parent.  The type-specific reactor objects turn 
out to be relatively simple since most of the computational effort is in Reactor1G.

All one energy group reactors are based on a algorithm published by the author in Nuclear Engineering & Design, 
"`A new method for rapid computation of transient fuel cycle material balances <http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6V4D-4W0SSGY-1&_user=10&_coverDate=10/31/2009&_rdoc=1&_fmt=high&_orig=search&_sort=d&_docanchor=&view=c&_acct=C000050221&_version=1&_urlVersion=0&_userid=10&md5=e52ececcd93400f84cc630ba20b01994>`_".

**Reactor1G Helper & Child Classes**

.. toctree::
    :maxdepth: 1

    fluence_point
    reactor_parameters
    light_water_reactor1g
    fast_reactor1g

.. currentmodule:: bright.reactor1g
    
.. autoclass:: Reactor1G(reactor_parameters=None, track_params=None, name="")

    .. autoattribute:: B
    .. autoattribute:: phi
    .. autoattribute:: fuel_chemical_form
    .. autoattribute:: coolant_chemical_form
    .. autoattribute:: rhoF
    .. autoattribute:: rhoC
    .. autoattribute:: P_NL
    .. autoattribute:: target_BU
    .. autoattribute:: use_zeta
    .. autoattribute:: lattice_flag
    .. autoattribute:: rescale_hydrogen_xs
    .. autoattribute:: r
    .. autoattribute:: l
    .. autoattribute:: S_O
    .. autoattribute:: S_T
    .. autoattribute:: VF
    .. autoattribute:: VC
    .. autoattribute:: libfile
    .. autoattribute:: F
    .. autoattribute:: BUi_F_
    .. autoattribute:: pi_F_
    .. autoattribute:: di_F_
    .. autoattribute:: Tij_F_
    .. autoattribute:: A_IHM
    .. autoattribute:: MWF
    .. autoattribute:: MWC
    .. autoattribute:: niF
    .. autoattribute:: niC
    .. autoattribute:: miF
    .. autoattribute:: miC
    .. autoattribute:: NiF
    .. autoattribute:: NiC

====================
Reactor1G Attributes
====================
    
--------------------------
Calculated Data Attributes
--------------------------
The following represents reactor that are calculated from the initial isotopics :attr:`ms_feed <bright.FCComp.ms_feed>`.
These attributes are assigned appropriate values when :meth:`Reactor1G.fold_mass_weights` is called.  Almost all of these
are C-vectors of doubles or floats.

.. attribute:: Reactor1G.dF_F_

    The neutron destruction rate of the fuel as a function of fluence.  This has units of [n/s] and is the linear
    combination of :attr:`di_F_` using :attr:`miF` as weights.

.. attribute:: Reactor1G.dC_F_

    The neutron destruction rate of the coolant as a function of fluence.  This has units of [n/s] and is the linear
    combination of :attr:`di_F_` using :attr:`miC` as weights.  If the disadvantage factor is used, then
    :attr:`zeta_F_` is multiplied by the linear combination before being assigned to :attr:`dC_F_`.
    

.. attribute:: Reactor1G.BU_F_

    The reactor burnup as a function of fluence.  This has units of [MWd/kgIHM] and is the linear
    combination of :attr:`BUi_F_` using :attr:`miF` as weights.

.. attribute:: Reactor1G.P_F_

    The full-core neutron production rate as a function of fluence.  This has units of [n/s] and is the linear
    combination of :attr:`pi_F_` using :attr:`miF` as weights. (Note: the coolant does not have a production rate).
    The linear combination is subsequently multiplied by the non-leakage probability, :attr:`P_NL`, before being 
    assigned to :attr:`P_F_`.

.. attribute:: Reactor1G.D_F_

    The full-core neutron destruction rate a function of fluence.  This has units of [n/s] and is the 
    sum of :attr:`dF_F_` and :attr:`dC_F_`.

.. attribute:: Reactor1G.k_F_

    The multiplication factor of the core.  Calculated from :attr:`P_F_` divided by :attr:`D_F_`.  This 
    attribute is unitless and not often used.

.. attribute:: Reactor1G.Mj_F_

    The transmutation matrix of the fuel (specifically, :attr:`ms_feed <bright.FCComp.ms_feed>`) into 
    the jth nuclide as a function of fluence.  Used with the discharge fluence :attr:`Fd` to calculate 
    :attr:`ms_prod <bright.FCComp.ms_prod>`.  This object is therefore a dictionary from zzaaam-integers to the 
    vectors of floats.

.. attribute:: Reactor1G.zeta_F_

    The thermal disadvantage factor as a function of fluence.  This attribute is unitless and is set when 
    :meth:`calc_zeta` is called.


--------------------
Discharge Attributes
--------------------
These attributes are ones that are solved for during a burnup calculation.  They are solved for in the 
process and relate important information about the state of the reactor at discharge.

.. attribute:: Reactor1G.fd

    The lower index of the discharge fluence (int).  The value of the discharge fluence is a float 
    that likely lies between two fluence points in :attr:`F`.  The :attr:`fd` attribute 
    gives the index of the lower of these fluence points.

.. attribute:: Reactor1G.Fd

    The discharge fluence (float).  May be used to calculate the amount of time that the 
    fuel was irradiated for.  Has units of [n/kb].

.. attribute:: Reactor1G.BUd

    The discharge burnup (float).  This has units of [MWd/kgIHM] and, unless something went very wrong, 
    should be rather close in value to :attr:`target_BU`.

.. attribute:: Reactor1G.k

    This is the multiplication factor of the reactor at discharge.  This should be very close in value to ``1.0``.

--------------------
SubStream Attributes
--------------------
Several parameters are dependent on knowing the mass or composition of specific substreams.  The ``In`` streams are derived from
:attr:`ms_feed <bright.FCComp.ms_feed>` while the ``Out`` SubStreams are derived from :attr:`ms_prod <BriPy.FCComp.ms_prod>`.
The easiest way to set these attributes is through the :meth:`Reactor1G.calcSubStreams` method.

.. attribute:: Reactor1G.ms_feed_u

    The input uranium mass stream, ``Reactor1G.ms_feed.get_u()``.

.. attribute:: Reactor1G.ms_feed_tru

    The input transuranic mass stream, ``Reactor1G.ms_feed.get_tru()``.

.. attribute:: Reactor1G.ms_feed_lan

    The input lanthanide mass stream, ``Reactor1G.ms_feed.get_lan()``.

.. attribute:: Reactor1G.ms_feed_act

    The input actinide mass stream, ``Reactor1G.ms_feed.get_act()``.

.. attribute:: Reactor1G.ms_prod_u

    The output uranium mass stream, ``Reactor1G.ms_prod.get_u()``.

.. attribute:: Reactor1G.ms_prod_tru

    The output transuranic mass stream, ``Reactor1G.ms_prod.get_tru()``.

.. attribute:: Reactor1G.ms_prod_lan

    The output lanthanide mass stream, ``Reactor1G.ms_prod.get_lan()``.

.. attribute:: Reactor1G.ms_prod_act

    The output actinide mass stream, ``Reactor1G.ms_prod.get_act()``.


----------------
Other Attributes
----------------

.. attribute:: Reactor1G.deltaR

    The :math:`\delta R` value of the core with ``ms_feed``.  This is equal to the 
    production rate minus the destruction rate at the target burnup::

        deltaR = batch_average(target_BU, "P") - batch_average(target_BU, "D")

    This is computed via the ``Reactor1G.calc_deltaR()`` method.

.. attribute:: Reactor1G.tru_cr

    The transuranic conversion ratio of the reactor (float).  This is set via :meth:`calc_tru_cr`.  


--------------------------------------
Thermal Disadvantage Factor Attributes
--------------------------------------
Many parameters go into the calculation of :attr:`zeta_F_` that are otherwise not needed.  Like
the disadvantage factor itself, these are functions of fluence and thus stored as vectors of 
float data.

.. attribute:: Reactor1G.SigmaFa_F_

    The fuel macroscopic absorption cross section.  In units of [1/cm].

.. attribute:: Reactor1G.SigmaFtr_F_

    The fuel macroscopic transport cross section.  In units of [1/cm].

.. attribute:: Reactor1G.kappaF_F_

    One over the thermal diffusion length of the fuel.  In units of [1/cm].


.. attribute:: Reactor1G.SigmaCa_F_

    The coolant macroscopic absorption cross section.  In units of [1/cm].

.. attribute:: Reactor1G.SigmaCtr_F_

    The coolant macroscopic transport cross section.  In units of [1/cm].

.. attribute:: Reactor1G.kappaC_F_

    One over the thermal diffusion length of the coolant.  In units of [1/cm].


.. attribute:: Reactor1G.lattice_E_F_

    The lattice function E(:attr:`F`).

.. attribute:: Reactor1G.lattice_F_F_

    The lattice function F(:attr:`F`).



.. _Reactor1G_Methods:

=================
Reactor1G Methods
=================
Along with all of the additional object attributes, :class:`Reactor1G` contains many more 
methods than other fuel cycle components.  All of the member functions are public.

------------------------------
Initialization Related Methods
------------------------------
Unlike other fuel cycle objects, the reactor component contains more than one setup related method.
Moreover, some of these will be often be called outside of a strict object instantiation context.
For instance, every time :attr:`ms_feed <bright.FCComp.ms_feed>` is changed, :meth:`fold_mass_weights`
should be called as well.

.. method:: Reactor1G.initialize(reactor_parameters)

    The :meth:`initialize` method for reactors copies all of the reactor specific parameters to this instance.
    Additionally, it calculates and sets the volumes :attr:`VF` and :attr:`VC`.

    Args:
        * `reactor_parameters` (:class:`ReactorParameters`): A special data structure that contains information
          on how to setup and run the reactor.

.. method:: Reactor1G.loadlib([libfile="Reactor.h5"])

    This method finds the HDF5 library for this reactor and extracts the necessary information from it.
    This method is typically called by the constructor of the child reactor type object.  It must be 
    called before attempting to do any real computation.

    Args: 
        * `libfile` (string): Path to the reactor library.


.. method:: Reactor1G.fold_mass_weights()

    This method performs the all-important task of doing the isotopically-weighted linear combination of raw data. 
    In a very real sense this is what makes this reactor *this specific reactor*.  The weights are taken 
    as the values of :attr:`ms_feed <bright.FCComp.ms_feed>`.  The raw data must have previously been 
    read in from :meth:`loadlib`.  

    .. warning::

        Anytime any reactor parameter whatsoever (:attr:`ms_feed <bright.FCComp.ms_feed>`, :attr:`P_NL`, *etc*) is 
        altered in any way, the :meth:`fold_mass_weights` function must be called to reset all of the resultant data.
        If you are unsure, please call this function anyway to be safe.  There is little harm in calling it twice by accident
        as it is computationally cheap to do so.  


----------------------------
Transmutation Matrix Methods
----------------------------
Computing the transmutation matrix is one of the more expensive :class:`Reactor1G` operations.  This is 
because transmutation is implicitly a function of three parameters: the input isotope ``i``, the output isotope
``j``, and the fluence :attr:`F`.  Therefore, for speedy execution times, the following functions should
only be called as needed (*ie* for generating :attr:`ms_prod <bright.FCComp.ms_prod>`).


.. method:: Reactor1G.calc_Mj_F_()

    This function calculates and sets the :attr:`Mj_F_` attribute from :attr:`ms_feed <bright.FCComp.ms_feed>` and the 
    raw reactor data :attr:`Tij_F_`.

.. method:: Reactor1G.calc_Mj_Fd_()

    This function evaluates :attr:`Mj_F_` calculated from :meth:`calc_Mj_F_` at the discharge fluence :attr:`Fd`.
    The resultant isotopic dictionary is then converted into the :attr:`ms_prod <bright.FCComp.ms_prod>` mass stream
    for this pass through the reactor.  Thus if ever you need to calculate :attr:`ms_prod <bright.FCComp.ms_prod>`
    without going through :meth:`calc`, use this function.


-------------------------
Basic Calculation Methods
-------------------------
The following functions represent basic calculations common to most reactor types.

.. method:: Reactor1G.calc_ms_prod()

    This is a convenience function that wraps the transmutation matrix methods.  It is equivalent to::

        #Wrapper to calculate discharge isotopics.
        calc_Mj_F_()
        calc_Mj_Fd_()

.. method:: Reactor1G.calcSubStreams()

    This sets possibly relevant reactor input and output substreams.  Specifically, it calculates the 
    attributes:

        * :attr:`ms_feed_u`
        * :attr:`ms_feed_tru`
        * :attr:`ms_feed_lan`
        * :attr:`ms_feed_act`
        * :attr:`ms_prod_u`
        * :attr:`ms_prod_tru`
        * :attr:`ms_prod_lan`
        * :attr:`ms_prod_act`

.. method:: Reactor1G.calc_tru_cr()

    This calculates and sets the transuranic conversion ratio :attr:`tru_cr` through the equation:

    .. math:: 

        \mbox{tru\_cr} = 1.0 - \frac{\mbox{ms\_feed\_tru.mass} - \mbox{ms\_prod\_tru.mass}}{\frac{\mbox{BUd}}{931.46}}

    Returns:
        * `tru_cr` (float): The value of the transuranic conversion ratio just calculated.

.. method:: Reactor1G.calc_deltaR([input])

    Calculates and sets the :math:`\\delta R` (:attr:`deltaR`) value of the reactor.  
    This is equal to the production rate minus the destruction rate at the target burnup::

        deltaR = batch_average(target_BU, "P") - batch_average(target_BU, "D")

    Args:
        * `input` (dict or MassStream): If input is present, it set as the component's 
          :attr:`ms_feed <bright.FCComp.ms_feed>`.  If input is a isotopic dictionary (zzaaam keys, float values), this
          dictionary is first converted into a MassStream before being set as :attr:`ms_feed <bright.FCComp.ms_feed>`.

    Returns:
        * `deltaR` (float): :attr:`deltaR`.


--------------
Burnup Methods    
--------------
These functions all relate to the calculation of burnup for the reactor from a given :attr:`ms_feed <bright.FCComp.ms_feed>`
stream.  The burnup functionality is separated out into so many functions so that the user has a fine 
degree of control over the burnup operation, if desired.  Higher level functions are also provided 
should this control be unnecessary for simple calculations.

.. method:: Reactor1G.fluence_at_BU(burnup)

    This function takes a burnup value  and returns a special fluence point object.  
    The fluence point is an amalgamation of data where the at which the burnup occurs.
    This object instance ``FP`` contains three pieces of information::
    
        FP.f    #Index immediately lower than where BU achieved (int)
        FP.F    #Fluence value itself (float)
        FP.m    #Slope dBU/dF between points f and f+1 (double)

    Args:
        * `burnup` (float): Burnup [MWd/kgIHM] at which to calculate the corresponding fluence.

    Returns:
        * `FP` (:class:`FluencePoint`): A class containing fluence information.


.. method:: Reactor1G.batch_average(BUd[, PDk_flag = "K"])

    Finds the batch-averaged ``P(F)``, ``D(F)``, or ``k(F)`` when at discharge burnup BUd.
    This function is typically iterated over until a BUd is found such that ``k(F) = 1.0 + err``.

    Args:
        * `BUd` (float): The discharge burnup [MWd/kgIHM] to obtain a batch-averaged value for.
        * `PDk_flag` (string): Flag that determines whether the neutron production rate "P" [n/s], 
          the neutron destruction rate "D" [n/s], or the multiplication factor "K" is reported in the output.

    Returns:
        * `PDk` (float): the batch averaged neutron production rate,
          neutron destruction rate, or the multiplication factor as determined by the input.

.. method:: Reactor1G.batch_average_k(BUd)

    Convenience function that calls :meth:`batch_average(BUd, "K") <batch_average>`.

    Args:
        * `BUd` (float): The discharge burnup [MWd/kgIHM] to obtain a batch-averaged value for.

    Returns:
        * `k` (float): the batch averaged multiplication factor.


.. method:: Reactor1G.BUd_bisection_method()

    Calculates the maximum discharge burnup via the Bisection Method for a given :attr:`ms_feed <bright.FCComp.ms_feed>`
    in this reactor.  This iterates over values of ``BUd`` to find a batch averaged multiplication factor 
    that is closest to ``1.0``.

    Other root finding methods for determining maximum discharge burnup are certainly possible.  
    However, with Bright's piecewise reactor data, the bisection method was found to return the most reliable results.


.. method:: Reactor1G.run_P_NL(pnl)

    Performs a reactor run for a specific non-leakage probability value.
    This requires that :attr:`ms_feed <bright.FCComp.ms_feed>` be (meaningfully) set and is
    for use with :meth:`calibrate_P_NL_to_BUd`.

    This function amounts to the following code::

        self.P_NL = pnl
        self.fold_mass_weights()
        self.BUd_bisection_method()

    Args:
        * `pnl` (float): The new non-leakage probability for the reactor.


.. method:: Reactor1G.calibrate_P_NL_to_BUd()

    Often times the non-leakage probability of a reactor is not known, though the input isotopics 
    and the target discharge burnup are.  This function handles that situation by
    calibrating the non-leakage probability of this reactor :attr:`P_NL` to hit its target burnup :attr:`target_BU`.
    Such a calibration proceeds by bisection method as well.  This function is extremely useful for 
    benchmarking calculations.


.. method:: Reactor1G.calc([input])

    Since many other methods provide the computational heavy-lifting of reactor calculations, 
    the :meth:`calc` method is relatively simple::

        self.ms_feed = input
        self.fold_mass_weights()
        self.BUd_bisection_method()
        self.calc_ms_prod()
        return self.ms_prod

    As you can see, all this function does is set burn an input stream to its maximum discharge burnup and then
    reports on the output isotopics.

    Args:
        * `input` (dict or MassStream): If input is present, it set as the component's 
          :attr:`ms_feed <bright.FCComp.ms_feed>`.  If input is a isotopic dictionary (zzaaam keys, float values), this
          dictionary is first converted into a MassStream before being set as :attr:`ms_feed <bright.FCComp.ms_feed>`.

    Returns:
        * `output` (MassStream): :attr:`ms_prod <bright.FCComp.ms_prod>`.


------------------------
Lattice Function Methods
------------------------
The lattice function methods below are all similar in that the calculate and set one of the 
lattice functions, E(F) or F(F), for one of the geometries, Planar or Spherical or Cylindrical.
The lattice functions are required for thermal disadvantage factor calculations.
They all take the equivalent of a fuel radius and a fuel cell pitch as input.  Please see
Lamarsh's `Nuclear Reactor Theory <http://www.amazon.com/Introduction-Nuclear-Reactor-Theory-Lamarsh/dp/0894480405>`_
for more information.


.. method:: Reactor1G.lattice_E_planar(a, b)

    Calculates the lattice function E(F) for planar geometry.  Sets value as :attr:`lattice_E_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].

.. method:: Reactor1G.lattice_F_planar(a, b)

    Calculates the lattice function F(F) for planar geometry.  Sets value as :attr:`lattice_F_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].

.. method:: Reactor1G.lattice_E_spherical(a, b)

    Calculates the lattice function E(F) for spherical geometry.  Sets value as :attr:`lattice_E_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].

.. method:: Reactor1G.lattice_F_spherical(a, b)

    Calculates the lattice function F(F) for spherical geometry.  Sets value as :attr:`lattice_F_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].

.. method:: Reactor1G.lattice_E_cylindrical(a, b)

    Calculates the lattice function E(F) for cylindrical geometry.  Sets value as :attr:`lattice_E_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].

.. method:: Reactor1G.lattice_F_cylindrical(a, b)

    Calculates the lattice function F(F) for cylindrical geometry.  Sets value as :attr:`lattice_F_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].


-----------------------------------
Thermal Disadvantage Factor Methods
-----------------------------------
These functions calculate and set the thermal disadvantage factor :attr:`zeta_F_`.  This is 
ostensibly done via the methodology detailed in 
Lamarsh's `Nuclear Reactor Theory <http://www.amazon.com/Introduction-Nuclear-Reactor-Theory-Lamarsh/dp/0894480405>`_.

Unfortunately, this formulation for the disadvantage factor is **only** valid in the case where ``a << b``!
Often times, modern (thermal) reactors do not satisfy this requirement.
We instead have a 'thin moderator' situation.

To fix this problem properly requires going to a multi-region diffusion/transport calculation.
Doing so is beyond the scope of Bright at this time and certainly beyond the aspirations of a 
one-group methodology.
A strategy that is more in-line with current practice is to use the results of a more sophisticated method,
interpolate over them, and use them here.

Thus in the case where `` 0.1 < VF / VC``, where the fuel is greater than 10% of the coolant, the 
above strategy is what is implemented.
A baseline disadvantage factor is determined from data presented in 
"*Thermal disadvantage factor calculation by the multi-region collision probability method*" by 
B. Ozgener,  and H. A. Ozgener, Institute of Nuclear Energy, Istanbul Technical University 80626 Maslak, Istanbul, Turkey, Received 
28 January 2003;  accepted 20 May 2003.  Available online 6 March 2004.
This baseline value happens to be a function of ``VF/VC``.

The Lamarsh method is then used as a scaling factor on the baseline function to 
account for variations in fuel composition and fluence.

.. method:: Reactor1G.calc_zeta()

    This calculates the thermal disadvantage factor for the geometry specified by :attr:`Lattice`.  The results
    are set to :attr:`zeta_F_`.

.. method:: Reactor1G.calc_zeta_planar()

    This calculates the thermal disadvantage factor for a planar geometry.  The results
    are set to :attr:`zeta_F_`.

.. method:: Reactor1G.calc_zeta_spherical()

    This calculates the thermal disadvantage factor for a spherical geometry.  The results
    are set to :attr:`zeta_F_`.

.. method:: Reactor1G.calc_zeta_cylindrical()

    This calculates the thermal disadvantage factor for a cylindrical geometry.  The results
    are set to :attr:`zeta_F_`.

