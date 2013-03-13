.. _bright_reactor_parameters:

******************
Reactor Parameters
******************
The reactor classes have a couple of helper classes that allow them to function in a more intuitive manner.
The ReactorParamters class specifies required data for the reactor model to run.  The use of this class 
prevents the initialization methods from having 15+ arguments.  Instead this class is passed to the constructor 
and related methods.

All functionality may be found in the ``reactor_parameters`` module::

    import bright.reactor_parameters

.. currentmodule:: bright.reactor_parameters   

========================
Reactor Parameters Class
========================

.. autoclass:: ReactorParameters()

    .. autoattribute:: batches
    .. autoattribute:: flux
    .. autoattribute:: fuel_form
    .. autoattribute:: cladding_form
    .. autoattribute:: coolant_form
    .. autoattribute:: fuel_density
    .. autoattribute:: cladding_density
    .. autoattribute:: coolant_density
    .. autoattribute:: pnl
    .. autoattribute:: BUt
    .. autoattribute:: specific_power
    .. autoattribute:: burn_regions
    .. autoattribute:: burn_times
    .. autoattribute:: use_disadvantage_factor
    .. autoattribute:: lattice_type
    .. autoattribute:: rescale_hydrogen
    .. autoattribute:: fuel_radius
    .. autoattribute:: void_radius
    .. autoattribute:: clad_radius
    .. autoattribute:: unit_cell_pitch
    .. autoattribute:: open_slots
    .. autoattribute:: total_slots


==========================
Reactor Parameter Defaults
==========================

.. autofunction:: lwr_defaults()

.. autofunction:: fr_defaults()
