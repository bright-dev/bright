**********************
Fuel Fabrication Class
**********************
Bright deals with fuel fabrication as an explicit fuel cycle component object.
This makes explicit the feedback effect between the Reactor component and what 
fuel is emplaced in the reactor.

A reactor operator orders fuel based on a target FuelFabrication, but the fuel 
fabricator has to know about the sepecifics of the reactor core to be able
to successfully construct an appropriate fuel.

Unlike most fuel cycle components, including the reactor object, fuel 
fabrication facilities may take more than input mass stream.  The fuel 
fabricator adjusts the mass ratios of these streams to calculate the 
required core input mass stream.

Therefore this component does an optimization calculation to determine the
correct reactor input given some set of available fuels. Currently, 
only two streams of all of the available are allowed to vary.  In the future,
a more sophisticated optimizer may be implemented.

.. currentmodule:: bright
    
.. class:: FuelFabrication([mass_streams, mass_weights_in, reactor, track_params, name])

    Fuel Fabrication Fuel Cycle Component Class.  Daughter of :class:`bright.FCComp` class.

    Args:
        * `mass_streams` (dict): A dictionary whose keys are string labels (eg, "U-235", 
          "TRU", "My Fuel") and whose values are mass streams.  For example::

            mass_streams = {
                "U235": MassStream({922350: 1.0}, 1.0, "U-235"),
                "U236": MassStream({922360: 1.0}, 1.0, "U-236"),
                "U238": MassStream({922380: 1.0}, 1.0, "U-238"),
                }

          would be valid for a light water reactor.
        * `mass_weights_in` (dict): A dictionary whose keys are the same as for mass_streams
          and whose values are the associated weight (float) for that stream.  If a stream
          should be allowed to vary (optimized over), specify its weight as a negative number.
          For instance::

            mass_weights_in = {
                "U235": -1.0,
                "U236": 0.005,
                "U238": -1.0,        
                }

          would be valid for a light water reactor with half a percent of U-236 always present.
        * `reactor` (:class:`Reactor1G`): An instance of a Reactor1G class to fabricate fuel for.
        * `track_params` (list of str): Additional parameters to track, if any.        
        * `name` (str): The name of the fuel fabrication fuel cycle component instance.

    Note that this automatically calls the public :meth:`initialize` C function.

.. _FuelFabrication_Attributes:

==========================
FuelFabrication Attributes
==========================
As a daughter class of :class:`bright.FCComp`, :class:`FuelFabrication` inherits all of 
the attributes of its parent.  The following is a listing of the additional 
attributes specific to this class.

.. attribute:: FuelFabrication.mass_streams

    A dictionary of mass streams that are mixed to create a valid fuel form for :attr:`reactor`.

.. attribute:: FuelFabrication.mass_weights_in

    A dictionary representing the initial specification of the mass weights.  Exactly two of these should 
    have negative values.

.. attribute:: FuelFabrication.mass_weights_out

    A dictionary representing the mass weights that are calculated to generate a valid fuel from the :attr:`mass_streams`
    provided. 

.. attribute:: FuelFabrication.deltaRs

    A dictionary representing the :attr:`deltaR <bright.Reactor1G.deltaR>` values of each of the :attr:`mass_streams`.

.. attribute:: FuelFabrication.reactor

    An instance of a :class:`Reactor1G` class that the :attr:`ms_prod <FCComp.ms_prod>` is valid as a fuel for.

.. attribute:: FuelFabrication.track_params

    For :class:`FuelFabrication`, the parameters that are automatically tracked are 
    the keys of :attr:`mass_streams` preceded by ``Weight_`` and ``deltaR_``.  For instance, 
    ``["Weight_U235", "deltaR_U235", "Weight_U238", "deltaR_U238"]``.

.. _FuelFabrication_Methods:

==============================
FuelFabrication Commom Methods
==============================

.. method:: FuelFabrication.calc([mass_streams, mass_weights_in, reactor])

    This method performs an optimization calculation on all input mass streams to determine
    the mass ratios that generate the correct fuel form for the reactor.  It then compiles 
    the fuel and returns the resultant MassStream. 

    Args:
        * `mass_streams` (dict): A dictionary whose keys are string labels and whose values are mass streams.  
        * `mass_weights_in` (dict): A dictionary whose keys are the same as for mass_streams
          and whose values are the associated weight (float) for that stream.
        * `reactor` (:class:`Reactor1G`): An instance of a Reactor1G class to fabricate fuel for.

    Returns:
        * `core_input` (MassStream): :attr:`ms_prod <bright.FCComp.ms_prod>`.


.. method:: FuelFabrication.initialize(mass_streams, mass_weights_in, reactor)

    The :meth:`initialize` function takes the appropriate mass streams, input mass weights,
    a reactor objects and resets the current FuelFabrication instance.

    Args:
        * `mass_streams` (dict): A dictionary whose keys are string labels and whose values are mass streams.  
        * `mass_weights_in` (dict): A dictionary whose keys are the same as for mass_streams
          and whose values are the associated weight (float) for that stream.
        * `reactor` (:class:`Reactor1G`): An instance of a Reactor1G class to fabricate fuel for.


.. method:: FuelFabrication.calc_params()

    Here the parameters for :class:`FuelFabrication` are set.  For example::

        self.params_prior_calc["Weight_U235"]  = self.mass_weights_in["U235"]
        self.Paramsout["Weight_U235"] = self.mass_weights_out["U235"]

        self.params_prior_calc["deltaR_U235"]  = self.deltaRs["U235"]
        self.Paramsout["deltaR_U235"] = self.deltaRs["U235"]

        self.params_prior_calc["Weight_U238"]  = self.mass_weights_in["U238"]
        self.Paramsout["Weight_U238"] = self.mass_weights_out["U238"]

        self.params_prior_calc["deltaR_U238"]  = self.deltaRs["U238"]
        self.Paramsout["deltaR_U238"] = self.deltaRs["U238"]


===================================
FuelFabrication Calculation Methods
===================================

.. method:: FuelFabrication.calc_deltaRs()

    Computes :attr:`deltaRs` for each mass stream.            

.. method:: FuelFabrication.calc_core_input()

    Computes the core input mass stream that becomes ms_prod based on :attr:`mass_streams` and 
    :attr:`mass_weights_out`.

    Returns:
        * `core_input` (MassStream): :attr:`ms_prod <bright.FCComp.ms_prod>`.

.. method:: FuelFabrication.calc_mass_ratios()

    Calculates :attr:`mass_weights_out` by varying the values of the two parameter which had 
    negative values in :attr:`mass_weights_in`.  Therefore, this is the portion of the code that 
    performs the optimization calculation.
