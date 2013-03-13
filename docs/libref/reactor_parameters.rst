.. _bright_reactor_parameters:

******************
Reactor Parameters
******************
The reactor classes have a couple of helper classes that allow them to function in 
a more intuitive manner.  The ReactorParamters class specifies required data for the 
reactor model to run.  The use of this class prevents the initialization methods 
from having 15+ arguments.  Instead this class is passed to the constructor and 
related methods.

All functionality may be found in the ``reactor_parameters`` module::

    import bright.reactor_parameters

.. currentmodule:: bright.reactor_parameters   

========================
Reactor Parameters Class
========================

.. autoclass:: ReactorParameters()
    :members:


==========================
Reactor Parameter Defaults
==========================

.. autofunction:: lwr_defaults()

.. autofunction:: fr_defaults()
