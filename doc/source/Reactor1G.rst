***************
Reactor1G Class
***************
Reactors are the most computationally complex of the nuclear fuel cycle components.  Bright handles nuclear reactors
in two distinct object classes: a one neutron energy group methodology (implemented here) and a multi-group algorithm (coming soon).
Moreover, there are several different types of reactors.  Each type has its own characteristic data library associated with it 
and takes on type-specific base case values.  The reactor types that have been fully implemented thus far are a light water 
reactor (LWR) and a fast reactor (FR).  You may read more about these on their own pages.

All one-group (1G) reactors share a common methodological backbone.  This page describes what is fundamentally the same 
about reactor objects via the :class:`bright.Reactor1G` class.  This is a subclass of :class:`BriPy.FCComp`.  Moreover, all 
one-group reactor types have :class:`bright.Reactor1G` as their parent.  The type-specific reactor objects turn out to 
be relatively simple since most of the computational effort is in  :class:`bright.Reactor1G`.

All one energy group reactors are based on a algorithm published by the author in Nuclear Engineering & Design, 
"`A new method for rapid computation of transient fuel cycle material balances <http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6V4D-4W0SSGY-1&_user=10&_coverDate=10/31/2009&_rdoc=1&_fmt=high&_orig=search&_sort=d&_docanchor=&view=c&_acct=C000050221&_version=1&_urlVersion=0&_userid=10&md5=e52ececcd93400f84cc630ba20b01994>`_".

.. currentmodule:: bright
    
.. class:: Reactor1G([reactor_parameters, track_params, name])

    One-Group Reactor Fuel Cycle Component Class.  Daughter of :class:`bright.FCComp` class.

    Args:
        * `reactor_parameters` (:class:`ReactorParameters`): A special data structure that contains information
          on how to setup and run the reactor.
        * `track_params` (string set): A set of strings that represents what parameter data the reactor should 
          store and set.  Different reactor types may have different characteristic parameters that are of interest.
        * `name` (str): The name of the reactor fuel cycle component instance.

    Note that this automatically calls the public :meth:`initialize <bright.Reactor1G.initialize>` C function.

.. note:: 

    Some data members and functions have names that end in ``_F_``.  This indicates that these are a 
    function of fluence, the time integral of the flux.  The ``_Fd_`` suffix implies that the data is 
    evaluated at the discharge fluence.

**Reactor1G Helper & Child Classes**

.. toctree::
    :maxdepth: 1

    Reactor1G_Helpers
    LightWaterReactor1G
    FastReactor1G

.. _Reactor1G_Attributes:

====================
Reactor1G Attributes
====================
As a daughter class of :class:`bright.FCComp`, :class:`Reactor1G` inherits all of the attributes of its parent.  
The following is a listing of the additional attributes specific to this class and those modified from the parent values.

--------------------
Protected Attributes
--------------------
These attributes are used internal to the reactor model and are therefore protected and not accessible from Python.
All other attributes are public and may be retrieved and modified from via the Python bindings.

.. attribute:: Reactor1G.I

    This is an zzaaam-integer isotopic set representing all of the nuclides that are valid inputs to the core.  This
    includes not just the heavy metal, but also coolant and cladding material as well.  :attr:`I` is typically 
    indexed by ``i``.  This is usually a strict subset of :func:`bright.track_isos`.

.. attribute:: Reactor1G.J

    This is an zzaaam-integer isotopic set representing all of the nuclides that are valid outputs from the core.  This
    encompasses all actinides and fission products.  :attr:`J` is typically 
    indexed by ``j``.  This is usually equivalent to :func:`bright.track_isos`.


.. attribute:: Reactor1G.sigma_a_therm

    This dictionary represents static, one-group microscopic thermal absorption cross sections in [barns] as 
    gathered from `KAERI <http://atom.kaeri.re.kr/>`_.  Such data is stored in the ``KaeriData.h5`` database
    that should be present in the ``BRIGHT_DATA`` directory.  This data is only read in if the thermal 
    disadvantage factor is used.

.. attribute:: Reactor1G.sigma_s_therm

    This dictionary represents static, one-group microscopic thermal scattering cross sections in [barns] as 
    gathered from `KAERI <http://atom.kaeri.re.kr/>`_.  Such data is stored in the ``KaeriData.h5`` database
    that should be present in the ``BRIGHT_DATA`` directory.  This data is only read in if the thermal 
    disadvantage factor is used.


----------------------------
Reactor Parameter Attributes
----------------------------
The following parameters are an instance copy of analogous `reactor_parameters` data that was fed to the :class:`Reactor1G` 
constructor.  Storing these functional parameters in each instance allows for several reactors of the same type to act 
in different manners in the same fuel cycle code.

.. attribute:: Reactor1G.B

    This integer is the total number of batches in the fuel management scheme.  :attr:`B` is typically indexed by ``b``.

.. attribute:: Reactor1G.phi

    The nominal flux value (float) that the library for this reactor type was generated with.  Used to correctly
    weight batch-specific fluxes.

.. attribute:: Reactor1G.fuel_chemical_form

    This is the chemical form of fuel as dictionary.  Keys are strings that represent isotopes (mixed form) while 
    values represent the corresponding mass weights.  The heavy metal concentration by the key ``"IHM"``.  This 
    will automatically fill in the nuclides in :attr:`ms_feed <bright.FCComp.ms_feed>` for the ``"IHM"`` weight.  For 
    example, LWRs typically use a UOX fuel form::

        LWR.fuel_chemical_form = {"IHM": 1.0, "O16": 2.0}

.. attribute:: Reactor1G.coolant_chemical_form

    This is the chemical form of coolant as dictionary.  This uses the same notation as :attr:`fuel_chemical_form` except 
    that ``"IHM"`` is no longer a valid key.  The term 'coolant' is used in preference over the term 'moderator' because
    not all reactors moderate neutrons.  For example, LWRs often cool the reactor core with borated water::

        LWR.coolant_chemical_form = {}

        LWR.CoolantchemicalForm["H1"]  = 2.0
        LWR.coolant_chemical_form["O16"] = 1.0
        LWR.coolant_chemical_form["B10"] = 0.199 * 550 * 10.0**-6
        LWR.coolant_chemical_form["B11"] = 0.801 * 550 * 10.0**-6

.. attribute:: Reactor1G.rhoF

    The fuel region density.  A float in units of [g/cm^3].

.. attribute:: Reactor1G.rhoC

    The coolant region density.  A float in units of [g/cm^3].

.. attribute:: Reactor1G.P_NL

    The reactor's non-leakage probability (float).  This is often used as a calibration parameter.

.. attribute:: Reactor1G.target_BU

    The reactor's target discharge burnup (float).  This is given in units of [MWd/kgIHM].  Often the
    actual discharge burnup :attr:`BUd` does not quite hit this value, but comes acceptably close.

.. attribute:: Reactor1G.use_zeta

    A boolean value that determines whether the thermal disadvantage factor is employed or not.  LWRs 
    typically set this as ``True`` while FRs have a ``False`` value.

.. attribute:: Reactor1G.Lattice

    A string flag that represents what lattice type the fuel assemblies are arranged in.  Currently accepted values 
    are ``"Planar"``, ``"Spherical"``, and ``"Cylindrical"``.

.. attribute:: Reactor1G.H_XS_Rescale

    This boolean determines whether the reactor should rescale the Hydrogen-1 destruction rate in the coolant as a
    function of fluence.  The scaling factor is calculated via the following equation:

    .. math:: f(F) = 1.36927 - 0.01119 \cdot BU(F)

    This is typically not done for fast reactors but is a useful correction for LWRs.

.. attribute:: Reactor1G.r

    The radius (float) of the fuel region.  In units of [cm].

.. attribute:: Reactor1G.l

    The pitch or length (float) of the unit fuel pin cell.  In units of [cm].

.. attribute:: Reactor1G.S_O

    The number of slots in a fuel assembly that are open (float).  *I.E.* this is the number of slots that 
    do not contain a fuel pin and are instead filled in by coolant. 

.. attribute:: Reactor1G.S_T

    The total number of fuel pin slots in a fuel assembly.  For a 17x17 bundle, :attr:`S_T` is ``289.0``. 

.. attribute:: Reactor1G.VF

    The relative fuel region volume.  Calculated from above information.

.. attribute:: Reactor1G.VC

    The relative coolant region volume.  Calculated from above information.

-----------------------------
Basic Reactor Data Attributes
-----------------------------
These attributes represent the raw data that is read in from a reactor type library.  These are loaded 
into memory when the :meth:`Reactor1G.loadLib` function is called.

.. attribute:: Reactor1G.libfile

    The path (string) to the reactor data library.  Usually something like ``LWR.h5`` or ``FR.h5``.

.. attribute:: Reactor1G.F

    The fluence points that the reactor library is based on.  This is a vector of floats that have units [n/kb] or 
    [neutrons/kilobarn].  This is read in from :attr:`libfile`.

.. attribute:: Reactor1G.BUi_F_

    The burnup of each initial isotope in the core as a function of fluence.  This is a dictionary whose 
    keys are in :attr:`I` and whose values are vectors of floats.  This data has units of [MWd/kgIHM] and is
    read in from :attr:`libfile`.

.. attribute:: Reactor1G.pi_F_

    The neutron production rate of each initial isotope in the core as a function of fluence.  This is a dictionary whose 
    keys are in :attr:`I` and whose values are vectors of floats.  This data has units of [n/s] or [neutrons/second] and is
    read in from :attr:`libfile`.

.. attribute:: Reactor1G.di_F_

    The neutron destruction rate of each initial isotope in the core as a function of fluence.  This is a dictionary whose 
    keys are in :attr:`I` and whose values are vectors of floats.  This data has units of [n/s] or [neutrons/second] and is
    read in from :attr:`libfile`.

.. attribute:: Reactor1G.Tij_F_

    The transmutation matrix of each initial isotope in the core into daughter nuclides as a function of fluence.  
    This is a dictionary whose 
    keys are in :attr:`I` and whose values also dictionaries.  These nested dictionaries have keys 
    that are members of :attr:`J` and whose values are vectors of floats.  
    This data has units of [kg_i/kgIHM] or kilogram of each initial isotope per kg of total initial heavy metal.
    This matrix is read in from :attr:`libfile`.


----------------------------
Calculated Weight Attributes
----------------------------
This data represents mass weights that are calculated from the initial isotopics :attr:`ms_feed <bright.FCComp.ms_feed>`.
They are assigned appropriate values during the :meth:`Reactor1G.foldMassWeights` execution.

.. attribute:: Reactor1G.A_IHM

    The atomic weight of the initial heavy metal (float).

.. attribute:: Reactor1G.MWF

    The molecular weight of the fuel (float).

.. attribute:: Reactor1G.MWC

    The molecular weight of the coolant (float).

.. attribute:: Reactor1G.niF

    Atomic number weight of the fuel as a function of initial isotope.  Dictionary with zzaaam-integer keys and 
    float values.

.. attribute:: Reactor1G.niC

    Atomic number weight of the coolant as a function of initial isotope.  Dictionary with zzaaam-integer keys and 
    float values.

.. attribute:: Reactor1G.miF

    Mass weight of the fuel as a function of initial isotope.  Dictionary with zzaaam-integer keys and 
    float values.

.. attribute:: Reactor1G.miC

    Mass weight of the coolant as a function of initial isotope.  Dictionary with zzaaam-integer keys and 
    float values.

.. attribute:: Reactor1G.NiF

    Number density of the fuel as a function of initial isotope.  Dictionary with zzaaam-integer keys and 
    float values.

.. attribute:: Reactor1G.NiC

    Number density of the coolant as a function of initial isotope.  Dictionary with zzaaam-integer keys and 
    float values.
    
--------------------------
Calculated Data Attributes
--------------------------
The following represents reactor that are calculated from the initial isotopics :attr:`ms_feed <bright.FCComp.ms_feed>`.
These attributes are assigned appropriate values when :meth:`Reactor1G.foldMassWeights` is called.  Almost all of these
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
    :meth:`calcZeta` is called.


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

.. attribute:: Reactor1G.InU

    The input uranium mass stream, ``Reactor1G.ms_feed.get_u()``.

.. attribute:: Reactor1G.InTRU

    The input transuranic mass stream, ``Reactor1G.ms_feed.get_tru()``.

.. attribute:: Reactor1G.InLAN

    The input lanthanide mass stream, ``Reactor1G.ms_feed.get_lan()``.

.. attribute:: Reactor1G.InACT

    The input actinide mass stream, ``Reactor1G.ms_feed.get_act()``.

.. attribute:: Reactor1G.OutU

    The output uranium mass stream, ``Reactor1G.ms_prod.get_u()``.

.. attribute:: Reactor1G.OutTRU

    The output transuranic mass stream, ``Reactor1G.ms_prod.get_tru()``.

.. attribute:: Reactor1G.OutLAN

    The output lanthanide mass stream, ``Reactor1G.ms_prod.get_lan()``.

.. attribute:: Reactor1G.OutACT

    The output actinide mass stream, ``Reactor1G.ms_prod.get_act()``.


----------------
Other Attributes
----------------

.. attribute:: Reactor1G.deltaR

    The :math:`\delta R` value of the core with ``ms_feed``.  This is equal to the 
    production rate minus the destruction rate at the target burnup::

        deltaR = batchAve(target_BU, "P") - batchAve(target_BU, "D")

    This is computed via the ``Reactor1G.calc_deltaR()`` method.

.. attribute:: Reactor1G.TruCR

    The transuranic conversion ratio of the reactor (float).  This is set via :meth:`calcTruCR`.  


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


.. attribute:: Reactor1G.LatticeE_F_

    The lattice function E(:attr:`F`).

.. attribute:: Reactor1G.LatticeF_F_

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
For instance, every time :attr:`ms_feed <bright.FCComp.ms_feed>` is changed, :meth:`foldMassWeights`
should be called as well.

.. method:: Reactor1G.initialize(reactor_parameters)

    The :meth:`initialize` method for reactors copies all of the reactor specific parameters to this instance.
    Additionally, it calculates and sets the volumes :attr:`VF` and :attr:`VC`.

    Args:
        * `reactor_parameters` (:class:`ReactorParameters`): A special data structure that contains information
          on how to setup and run the reactor.

.. method:: Reactor1G.loadLib([libfile="Reactor.h5"])

    This method finds the HDF5 library for this reactor and extracts the necessary information from it.
    This method is typically called by the constructor of the child reactor type object.  It must be 
    called before attempting to do any real computation.

    Args: 
        * `libfile` (string): Path to the reactor library.


.. method:: Reactor1G.foldMassWeights()

    This method performs the all-important task of doing the isotopically-weighted linear combination of raw data. 
    In a very real sense this is what makes this reactor *this specific reactor*.  The weights are taken 
    as the values of :attr:`ms_feed <bright.FCComp.ms_feed>`.  The raw data must have previously been 
    read in from :meth:`loadLib`.  

    .. warning::

        Anytime any reactor parameter whatsoever (:attr:`ms_feed <bright.FCComp.ms_feed>`, :attr:`P_NL`, *etc*) is 
        altered in any way, the :meth:`foldMassWeights` function must be called to reset all of the resultant data.
        If you are unsure, please call this function anyway to be safe.  There is little harm in calling it twice by accident
        as it is computationally cheap to do so.  


----------------------------
Transmutation Matrix Methods
----------------------------
Computing the transmutation matrix is one of the more expensive :class:`Reactor1G` operations.  This is 
because transmutation is implicitly a function of three parameters: the input isotope ``i``, the output isotope
``j``, and the fluence :attr:`F`.  Therefore, for speedy execution times, the following functions should
only be called as needed (*ie* for generating :attr:`ms_prod <bright.FCComp.ms_prod>`).


.. method:: Reactor1G.mkMj_F_()

    This function calculates and sets the :attr:`Mj_F_` attribute from :attr:`ms_feed <bright.FCComp.ms_feed>` and the 
    raw reactor data :attr:`Tij_F_`.

.. method:: Reactor1G.mkMj_Fd_()

    This function evaluates :attr:`Mj_F_` calculated from :meth:`mkMj_F_` at the discharge fluence :attr:`Fd`.
    The resultant isotopic dictionary is then converted into the :attr:`ms_prod <bright.FCComp.ms_prod>` mass stream
    for this pass through the reactor.  Thus if ever you need to calculate :attr:`ms_prod <bright.FCComp.ms_prod>`
    without going through :meth:`calc`, use this function.


-------------------------
Basic Calculation Methods
-------------------------
The following functions represent basic calculations common to most reactor types.

.. method:: Reactor1G.calcOutIso()

    This is a convenience function that wraps the transmutation matrix methods.  It is equivalent to::

        #Wrapper to calculate discharge isotopics.
        mkMj_F_()
        mkMj_Fd_()

.. method:: Reactor1G.calcSubStreams()

    This sets possibly relevant reactor input and output substreams.  Specifically, it calculates the 
    attributes:

        * :attr:`InU`
        * :attr:`InTRU`
        * :attr:`InLAN`
        * :attr:`InACT`
        * :attr:`OutU`
        * :attr:`OutTRU`
        * :attr:`OutLAN`
        * :attr:`OutACT`

.. method:: Reactor1G.calcTruCR()

    This calculates and sets the transuranic conversion ratio :attr:`TruCR` through the equation:

    .. math:: \mbox{TruCR} = \frac{\mbox{InTRU.mass} - \mbox{OutTRU.mass}}{\frac{\mbox{BUd}}{935.0}}

    Returns:
        * `TruCR` (float): The value of the transuranic conversion ratio just calculated.

.. method:: Reactor1G.calc_deltaR([input])

    Calculates and sets the :math:`\delta R` (:attr:`deltaR`) value of the reactor.  
    This is equal to the production rate minus the destruction rate at the target burnup::

        deltaR = batchAve(target_BU, "P") - batchAve(target_BU, "D")

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

.. method:: Reactor1G.FluenceAtBU(burnup)

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


.. method:: Reactor1G.batchAve(BUd[, PDk_flag = "K"])

    Finds the batch-averaged ``P(F)``, ``D(F)``, or ``k(F)`` when at discharge burnup BUd.
    This function is typically iterated over until a BUd is found such that ``k(F) = 1.0 + err``.

    Args:
        * `BUd` (float): The discharge burnup [MWd/kgIHM] to obtain a batch-averaged value for.
        * `PDk_flag` (string): Flag that determines whether the neutron production rate "P" [n/s], 
          the neutron destruction rate "D" [n/s], or the multiplication factor "K" is reported in the output.

    Returns:
        * `PDk` (float): the batch averaged neutron production rate,
          neutron destruction rate, or the multiplication factor as determined by the input.

.. method:: Reactor1G.batchAveK(BUd)

    Convenience function that calls :meth:`batchAve(BUd, "K") <batchAve>`.

    Args:
        * `BUd` (float): The discharge burnup [MWd/kgIHM] to obtain a batch-averaged value for.

    Returns:
        * `k` (float): the batch averaged multiplication factor.


.. method:: Reactor1G.BUd_BisectionMethod()

    Calculates the maximum discharge burnup via the Bisection Method for a given :attr:`ms_feed <bright.FCComp.ms_feed>`
    in this reactor.  This iterates over values of ``BUd`` to find a batch averaged multiplication factor 
    that is closest to ``1.0``.

    Other root finding methods for determining maximum discharge burnup are certainly possible.  
    However, with Bright's piecewise reactor data, the bisection method was found to return the most reliable results.


.. method:: Reactor1G.Run_PNL(pnl)

    Performs a reactor run for a specific non-leakage probability value.
    This requires that :attr:`ms_feed <bright.FCComp.ms_feed>` be (meaningfully) set and is
    for use with :meth:`Calibrate_PNL_2_BUd`.

    This function amounts to the following code::

        self.P_NL = pnl
        self.foldMassWeights()
        self.BUd_BisectionMethod()

    Args:
        * `pnl` (float): The new non-leakage probability for the reactor.


.. method:: Reactor1G.Calibrate_PNL_2_BUd()

    Often times the non-leakage probability of a reactor is not known, though the input isotopics 
    and the target discharge burnup are.  This function handles that situation by
    calibrating the non-leakage probability of this reactor :attr:`P_NL` to hit its target burnup :attr:`target_BU`.
    Such a calibration proceeds by bisection method as well.  This function is extremely useful for 
    benchmarking calculations.


.. method:: Reactor1G.calc([input])

    Since many other methods provide the computational heavy-lifting of reactor calculations, 
    the :meth:`calc` method is relatively simple::

        self.ms_feed = input
        self.foldMassWeights()
        self.BUd_BisectionMethod()
        self.calcOutIso()
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


.. method:: Reactor1G.LatticeEPlanar(a, b)

    Calculates the lattice function E(F) for planar geometry.  Sets value as :attr:`LatticeE_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].

.. method:: Reactor1G.LatticeFPlanar(a, b)

    Calculates the lattice function F(F) for planar geometry.  Sets value as :attr:`LatticeF_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].

.. method:: Reactor1G.LatticeESpherical(a, b)

    Calculates the lattice function E(F) for spherical geometry.  Sets value as :attr:`LatticeE_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].

.. method:: Reactor1G.LatticeFSpherical(a, b)

    Calculates the lattice function F(F) for spherical geometry.  Sets value as :attr:`LatticeF_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].

.. method:: Reactor1G.LatticeECylindrical(a, b)

    Calculates the lattice function E(F) for cylindrical geometry.  Sets value as :attr:`LatticeE_F_`

    Args:
        * `a` (float): Fuel region radius equivalent [cm].
        * `b` (float): Unit fuel cell pitch length equivalent [cm].

.. method:: Reactor1G.LatticeFCylindrical(a, b)

    Calculates the lattice function F(F) for cylindrical geometry.  Sets value as :attr:`LatticeF_F_`

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

.. method:: Reactor1G.calcZeta()

    This calculates the thermal disadvantage factor for the geometry specified by :attr:`Lattice`.  The results
    are set to :attr:`zeta_F_`.

.. method:: Reactor1G.calcZetaPlanar()

    This calculates the thermal disadvantage factor for a planar geometry.  The results
    are set to :attr:`zeta_F_`.

.. method:: Reactor1G.calcZetaSpherical()

    This calculates the thermal disadvantage factor for a spherical geometry.  The results
    are set to :attr:`zeta_F_`.

.. method:: Reactor1G.calcZetaCylindrical()

    This calculates the thermal disadvantage factor for a cylindrical geometry.  The results
    are set to :attr:`zeta_F_`.

