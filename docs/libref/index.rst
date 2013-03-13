.. _libref:

=================
Library Reference
=================
This presents information on the Bright API.  The api module is the central place for 
accessing all of Bright's features::

    from bright.api import *

This command implicitly imports all of the underlying extension modules.  While bright is 
primarily a collection of fuel cycle components (subclasses of FCComp), there are a couple 
of package-level functions and attributes that must be set for Bright to run successfully.  
These are explained here.  Meanwhile, because of their potential for complexity, each fuel cycle 
component is given its own page.  


~~~~~~~~~~~
Components
~~~~~~~~~~~

All fuel cycle objects inherit from a common FCComp class.  This takes care of 
all of the bookkeeping, input, and output for all component instances.  Below is 
a diagram of the how all of the comonents inheret from the top-level FCComp object. 
Generally, a bright user will only need to call the bottom-level classes directly.

.. toctree::
    :maxdepth: 1

    fccomp
    enrichment_parameters
    enrichment
    fuel_fabrication
    fluence_point
    reactor_parameters
    reactor1g
    light_water_reactor1g
    fast_reactor1g
    reactormg
    origen_reactormg
    reprocess
    storage


.. inheritance-diagram:: bright.bright_config bright.fccomp bright.enrichment 
    bright.fuel_fabrication bright.fluence_point bright.reactor_parameters 
    bright.reactor1g bright.light_water_reactor1g bright.fast_reactor1g 
    bright.reactormg bright.origen_reactormg bright.reprocess bright.storage
    bright.enrichment_parameters
    :parts: 1

~~~~~~~~~~~
Other APIs
~~~~~~~~~~~

The following parts of Bright aid in the construction of the fuel cycle or 
help provide an API for generating code.

.. toctree::
    :maxdepth: 2

    bright_config
    apigen/index
