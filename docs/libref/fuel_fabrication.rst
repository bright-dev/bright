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
    
.. autoclass:: FuelFabrication(mass_streams=None, mass_weights_in=None, reactor=None, track_params=None, name="")

    .. autoattribute:: materials
    .. autoattribute:: mass_weights_in
    .. autoattribute:: mass_weights_out
    .. autoattribute:: deltaRs
    .. autoattribute:: reactor
    .. attribute:: track_params

        For FuelFabrication, the parameters that are automatically tracked are 
        the keys of :attr:`mass_streams` preceded by ``Weight_`` and ``deltaR_``.  For 
        instance, ``["Weight_U235", "deltaR_U235", "Weight_U238", "deltaR_U238"]``.

    .. automethod:: initialize(materials, mass_weights_in, reactor)
    .. automethod:: calc(materials=None, mass_weights_in=None, reactor=None)
    .. automethod:: calc_params()
    .. automethod:: calc_deltaRs()
    .. automethod:: calc_core_input()
    .. automethod:: calc_mass_ratios()
