*************
Storage Class
*************
The storage class very simply solves the Bateman decay equations on an input mass stream.  For a more complete 
treatment please see Stacey's 
`Nuclear Reactor Physics <http://www.amazon.com/Nuclear-Reactor-Physics-Weston-Stacey/dp/0471391271>`_.
Currently, the storage component uses a recursive method to calculate the decay chains and resultant masses.  
However, a more efficient, more general matrix exponential method is expected soon.  This will probably use a Pade approximation.

To use this class necessary decay information must be available.  This is stored in a ``decay.h5`` database that must be
present within `BRIGHT_DATA` directory.   Instantiation of this class automatically calls the protected 
:meth:`initialize <bright.Storage.initialize>`  C function, which loads the ``decay.h5`` library.

.. currentmodule:: bright
    
.. class:: Storage([name])

    Storage Fuel Cycle Component Class.  Daughter of :class:`bright.FCComp` class.

    Args:
        * `name` (str): The name of the storage fuel cycle component instance.


.. _Storage_Attributes:

==================
Storage Attributes
==================
As a daughter class of :class:`bright.FCComp`, :class:`Storage` inherits all of the attributes of its parent.  
The following is a listing of the additional attributes specific to this class and those modified from the parent.

.. attribute:: Storage.decay 

    This is a protected dictionary that has zzaaam-integer keys and ``FromIsosStruct`` values.  ``FromIsosStruct`` is 
    a C structure that stores information on which nuclides the parent-key nuclide decays into.  This is the 
    algorithmic friendly version of the ``decay.h5`` file.  

.. attribute:: Storage.decay_data

    This is another protected C member that represents the ``decay.h5`` file.  This is a closer approximation
    to how the actual library is structured.  It is an array of ``DecayIso`` types.

.. attribute:: Storage.decay_data_len

    This a protected integer that is equal to the length of :attr:`decay_data`.

.. attribute:: Storage.decay_time

    This the public float (double) attribute that represents how long an input fuel mass should be 
    stored for.  This time in represented in seconds, so be sure to convert to the proper units 
    before using.

.. attribute:: Storage.isochains

    To prevent redundant calculations of decay chains, previously used decays are stored in this 
    protected attribute.  :attr:`isochains` itself is a set of zzaaam-integer vectors.  These vectors 
    represent the progression from mother to daughter to granddaughter to *etc*.

.. attribute:: Storage.params2track

    For :class:`Storage`, the only parameter that is tracked is the aggregate mass.  
    Thus this attribute is automatically set to ``["Mass"]``.

.. _Storage_Methods:

===============
Storage Methods
===============

.. method:: Storage.PrintChain(isochain)

    This is a protected C method that prints the contents of an isotopic chain as a list.  It 
    is primarily used for debugging.

    Args:
        * `isochain` (C-vector): Vector containing zzaaam-nuclides from the parent through the 
          daughters.  Need not end with a stable isotope.

.. method:: Storage.addchains(iso)

    This method adds an isotope's decay chains or a single chain to the :attr:`isochains` set. 
    This functionality is protected and called automatically by :meth:`doCalc`.

    Args:
        * `iso` (int or isochain):  The isotope or chain to add.  If iso is an integer (zzaaam form), 
          then it is treated as the mother nuclide and all chains are added until a stable nuclide has been 
          reached.  If iso is an isotopic chain, then it is added to :attr:`isochains` directly. 

.. method:: Storage.bateman(iso, mass, isochain)

    This is the central part of the storage algorithm as :meth:`bateman` solves the Bateman
    equation for an isotope decaying into the last nuclide in the isochain vector, starting with a given mass.
    Note that this decay is calculated for a :attr:`decay_time` number of seconds.  This method is private 
    and implicitly called by :meth:`doCalc`.

    Args:
        * `iso` (zzaaam int): The parent nuclide.  Must be the same as ``isocahin[0]``.
        * `mass` (double): The initial mass of iso.
        * `isochain` (C-vector): The decay chain for iso into a daughter nuclide, ``isochain[-1]``.

.. method:: Storage.doCalc([input, time])

    As usual, :meth:`doCalc` sets up the Storage component's input stream and calculates the corresponding 
    output :class:`MassStream`.  Here, this amounts to calling :meth:`bateman` for every nuclide in 
    :attr:`IsosIn <bright.FCComp.IsosIn>`, for each chain that ends with a nuclide in :meth:`BriPy.isos2track`.

    This method is public and accessible from Python.

    Args:
        * `input` (dict or MassStream): If input is present, it set as the component's 
          :attr:`IsosIn <bright.FCComp.IsosIn>`.  If input is a isotopic dictionary (zzaaam keys, float values), this
          dictionary is first converted into a MassStream before being set as :attr:`IsosIn <bright.FCComp.IsosIn>`.
        * `time` (float): :attr:`decay_time` is set to the time value here prior to any other calculations.  This
          time has units of seconds.

    Returns:
        * `output` (MassStream): :attr:`IsosOut <bright.FCComp.IsosOut>`.


.. method:: Storage.initialize()

    The :meth:`initialize` function's primary purpose is to load the ``decay.h5`` data library into the 
    appropriate points in memory.  It is protected and called automatically by the :class:`Storage` constructor.


.. method:: Storage.setParams()

    Here the parameters for :class:`Storage` are set.  For storage, this amounts to just
    a "Mass" parameter::

        self.ParamsIn["Mass"]  = self.IsosIn.mass
        self.ParamsOut["Mass"] = self.IsosOut.mass
