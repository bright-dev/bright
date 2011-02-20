.. currentmodule:: bright

************************
Reactor1G Helper Classes
************************
The :class:`Reactor1G` class has a couple of helper classes that allow it to function in a more intuitive manner.
The first is a simple :class:`FluencePoint` data structure.  This contains information on where one is in the 
:attr:`Reactor1G.F` vector.  The second object is a :class:`ReactorParamters` class that specifies required data
for the reactor to run.  The use of this class prevents the initialization methods from having 15+ arguments.  
Instead this class is passed to the constructor and related methods.

    
===================
Fluence Point Class
===================

.. class:: FluencePoint()

    Fluence point class with three simple data members.

    .. attribute:: f    

        Index of :attr:`Reactor1G.F` immediately lower than the value of :attr:`F` (int).

    .. attribute:: F    

        Fluence value itself (float).  In units of [n/kb] or [neutrons/kilobarn].

    .. attribute:: m    

        The slope dBU/dF between points f and f+1 (float).  Has the odd units of [MWd kb / kgIHM n]
    

========================
Reactor Parameters Class
========================

.. class:: ReactorParameters()

    This data structure is a set of physical reactor parameters. It may be used to instantiate new reactor objects **OR**
    to define default settings for a reactor type.  The data stored in this class is copied over to 
    a reactor instance in the :meth:`Reactor1G.initialize` method.  However, the attributes of this objects 
    take on more natural names than their :class:`Reactor1G` analogies.  This is because it is this 
    object that Bright users will more often be interacting with. 

    .. attribute:: batches

        This integer is the total number of batches in the fuel management scheme.  This is typically indexed by ``b``.

    .. attribute:: flux

        The nominal flux value (float) that the library for this reactor type was generated with.  Used to correctly
        weight batch-specific fluxes.

    .. attribute:: FuelForm

        This is the chemical form of fuel as dictionary.  Keys are strings that represent isotopes (mixed form) while 
        values represent the corresponding mass weights.  The heavy metal concentration by the key ``"IHM"``.  This 
        will automatically fill in the nuclides in :attr:`ms_feed <bright.FCComp.ms_feed>` for the ``"IHM"`` weight.  For 
        example, LWRs typically use a UOX fuel form::

            ReactorParameters.FuelForm = {"IHM": 1.0, "O16": 2.0}

    .. attribute:: CoolantForm

        This is the chemical form of coolant as dictionary.  This uses the same notation as :attr:`FuelForm` except 
        that ``"IHM"`` is no longer a valid key.  The term 'coolant' is used in preference over the term 'moderator' because
        not all reactors moderate neutrons.  For example, LWRs often cool the reactor core with borated water::

            ReactorParamters.CoolantForm = {}

            ReactorParamters.CoolantForm["H1"]  = 2.0
            ReactorParamters.CoolantForm["O16"] = 1.0
            ReactorParamters.CoolantForm["B10"] = 0.199 * 550 * 10.0**-6
            ReactorParamters.CoolantForm["B11"] = 0.801 * 550 * 10.0**-6

    .. attribute:: FuelDensity

        The fuel region density.  A float in units of [g/cm^3].

    .. attribute:: CoolantDensity

        The coolant region density.  A float in units of [g/cm^3].

    .. attribute:: pnl

        The reactor's non-leakage probability (float).  This is often used as a calibration parameter.

    .. attribute:: BUt

        The reactor's target discharge burnup (float).  This is given in units of [MWd/kgIHM].  Often the
        actual discharge burnup :attr:`BUd` does not quite hit this value, but comes acceptably close.

    .. attribute:: useDisadvantage

        A boolean value that determines whether the thermal disadvantage factor is employed or not.  LWRs 
        typically set this as ``True`` while FRs have a ``False`` value.

    .. attribute:: LatticeType

        A string flag that represents what lattice type the fuel assemblies are arranged in.  Currently accepted values 
        are ``"Planar"``, ``"Spherical"``, and ``"Cylindrical"``.

    .. attribute:: HydrogenRescale

        This boolean determines whether the reactor should rescale the Hydrogen-1 destruction rate in the coolant as a
        function of fluence.  The scaling factor is calculated via the following equation:

        .. math:: f(F) = 1.36927 - 0.01119 \cdot BU(F)

        This is typically not done for fast reactors but is a useful correction for LWRs.

    .. attribute:: Radius

        The radius (float) of the fuel region.  In units of [cm].

    .. attribute:: Length

        The pitch or length (float) of the unit fuel pin cell.  In units of [cm].

    .. attribute:: open_slots

        The number of slots in a fuel assembly that are open (float).  *I.E.* this is the number of slots that 
        do not contain a fuel pin and are instead filled in by coolant. 

    .. attribute:: total_slots

        The total number of fuel pin slots in a fuel assembly.  For a 17x17 bundle, :attr:`S_T` is ``289.0``. 
