.. _bright_fuel_fabrication:

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

All functionality may be found in the ``fuel_fabrication`` module::

    import bright.fuel_fabrication

.. currentmodule:: bright.fuel_fabrication
    
.. autoclass:: FuelFabrication(mats=None, mws_in=None, r=None, paramtrack=None, n="")
    :members:

    .. attribute:: track_params

        For FuelFabrication, the parameters that are automatically tracked are 
        the keys of :attr:`mass_streams` preceded by ``Weight_`` and ``deltaR_``.  For 
        instance, ``["Weight_U235", "deltaR_U235", "Weight_U238", "deltaR_U238"]``.
