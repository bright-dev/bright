.. _bright_storage:

*************
Storage Class
*************
The storage class very simply solves the Bateman decay equations on an input mass stream.  For a more complete 
treatment please see Stacey's 
`Nuclear Reactor Physics <http://www.amazon.com/Nuclear-Reactor-Physics-Weston-Stacey/dp/0471391271>`_.
Currently, the storage component uses a recursive method to calculate the decay chains and resultant masses.  

To use this class the necessary decay information must be available.  This is stored in a ``decay.h5`` database 
that must be present within the `BRIGHT_DATA` directory.   Instantiation of this class automatically calls 
the initialize() method on the C++ level, which loads the ``decay.h5`` library.

All functionality may be found in the ``storage`` module::

    import bright.storage

.. currentmodule:: bright.storage
    
.. autoclass:: Storage(name="")

    .. autoattribute:: decay_time
    .. attribute:: track_params

        For :class:`Storage`, the only parameter that is tracked is the aggregate mass.  
        Thus this attribute is automatically set to ``["Mass"]``.

    .. automethod:: calc(input=None, time=None)
    .. automethod:: calc_params()


