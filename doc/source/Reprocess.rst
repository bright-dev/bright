***************
Reprocess Class
***************
The reprocessing object is one of the simplest that Bright offers.  In a technology nonspecific way, it performs
separations on an input mass stream.  Reprocessing facilities have nominal separation efficiencies (SE) for each element 
(and possibly nuclide).  These SE are valued on the range [0,1] and represent the fractional portion of an element that 
is recovered in the reprocessing output (:attr:`ms_prod <bright.FCComp.ms_prod>`).

The SE, if left unspecified, default to full recovery (a value of 1).  Currently, Reprocess does not contain a ``TailsOut`` 
attribute that corresponds to material not recovered in :attr:`ms_prod <bright.FCComp.ms_prod>` (equivalent to 1 - SE).  
However, as such an item is becoming evermore useful, its future implementation is likely.

.. currentmodule:: bright
    
.. class:: Reprocess([sepeff[, name]])

    Reprocess Fuel Cycle Component Class.  Daughter of :class:`bright.FCComp` class.

    Args:
        * `sepeff` (dict): A dictionary containing the separation efficiencies (float) to initialize
          the instance with.  The keys of this dictionary must be strings.  However, the strings may 
          represent either elements or isotopes or both::

                #ssed = string dictionary of separation efficiencies.  
                #Of form {zz: 0.99}, eg 
                ssed = {"92": 0.999, "94": 0.99} 
                #of form {LL: 0.99}, eg 
                ssed = {"U": 0.999, "PU": 0.99} 
                #or of form {mixed: 0.99}, eg 
                ssed = {"U235": 0.9, "922350": 0.999, "94239": 0.99}

        * `name` (str): The name of the reprocessing fuel cycle component instance.

    Note that this automatically calls the public :meth:`initialize` C function.

    .. note::
       The C++ version of the code also allows you to initialize from an int-keyed dictionary (map).
       However, due to a from_python C++ signature ambiguity, you cannot do use this directly in Python.
       Separation efficiencies must therefore be automatically initialized through string dictionaries.
       If you need to initialize via an int dictionary in python, you can always init with an empty
       string dictionary and then manually initialize with an int one.  For example::

            R = Reprocess({}, name)
            R.initialize( {92: 0.99, 942390: 0.9} )

.. _Reprocess_Attributes:

====================
Reprocess Attributes
====================
As a daughter class of :class:`bright.FCComp`, :class:`Reprocess` inherits all of the attributes of its parent.  
The following is a listing of the additional attributes specific to this class.

.. attribute:: Reprocess.sepeff

    This is a dictionary or map representing the separation efficiencies of each isotope in :func:`bright.track_nucs`.
    Therefore it has zzaaam-integer keys and float (double) values.  During initialization, other SE dictionaries are converted 
    to this standard form::

        sepeff = {922350: 0.999, 942390: 0.99}

    This attribute may be accessed and altered directly (public).

.. attribute:: Reprocess.track_params

    For :class:`Reprocess`, the only parameter that is tracked is the aggregate mass.  Thus this attribute
    is automatically set to ``["Mass"]``.

.. _Reprocess_Methods:

=================
Reprocess Methods
=================

.. method:: Reprocess.calc([input])

    This method performs the relatively simply task of multiplying the current input stream by 
    the SE to form a new output stream::

        incomp  = self.ms_feed.mult_by_mass()
        outcomp = {}
        for iso in incomp.keys():
            outcomp[iso] = incomp[iso] * sepeff[iso]
        self.ms_prod = MassStream(outcomp)
        return self.ms_prod

    Args:
        * `input` (dict or MassStream): If input is present, it set as the component's 
          :attr:`ms_feed <bright.FCComp.ms_feed>`.  If input is a isotopic dictionary (zzaaam keys, float values), this
          dictionary is first converted into a MassStream before being set as :attr:`ms_feed <bright.FCComp.ms_feed>`.

    Returns:
        * `output` (MassStream): :attr:`ms_prod <bright.FCComp.ms_prod>`.


.. method:: Reprocess.initialize(sepdict)

    The :meth:`initialize` function calculates the :attr:`sepeff` from an integer-keyed dictionary
    of separation efficiencies.  The difference is that `sepdict` may contain either elemental or
    isotopic keys and need not contain every isotope tracked.  On the other hand, :attr:`sepeff`
    must have only zzaaam keys that match exactly the isotopes in :func:`bright.track_nucs`.

    Args:
        * `sepdict` (dict): Integer valued dictionary of SE to be converted to :attr:`sepeff`.


.. method:: Reprocess.calc_params()

    Here the parameters for :class:`Reprocess` are set.  For reprocessing, this amounts to just
    a "Mass" parameter::

        self.params_prior_calc["Mass"]  = self.ms_feed.mass
        self.params_after_calc["Mass"] = self.ms_prod.mass

