************************
Bright Extension Module
************************
The `Bright` extension module contains all of the fuel cycle components that Bright may use and exists 
within the `bright` package.  

The documentation for `bright` is split up based on component object (or class).  This demonstrates the heirarchy 
between classes and allows for greater exploration of each components methodology.

**Components**

.. toctree::
   :maxdepth: 3

   Reprocess
   Storage
   Reactor1G
   FuelFabrication

Currently, all fuel cycle objects inherit from a common `FCComp` class.  This takes care of all of the bookkeeping, 
input, and output for all component instances.  Below is a a diagram of the how all of the comonents inheret 
from the top-level `FCComp` object.  Generally, a Bright user will only need to call the bottom-level classes
directly.


.. inheritance-diagram:: bright
   :parts: 1

Additionally, a multi-component enrichment model also has been developed.  However, it has yet to be brought up to 
Bright standards.   Make a request and it may magically show up!
