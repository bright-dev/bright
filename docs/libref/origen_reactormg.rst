.. _bright_origen_reactormg:

**********************
Origen ReactorMG Class
**********************
The existing low-level transmutation method for ReactorMG is a matrix exponential method 
which leaves the matrix unadultered and very stiff.  This means that transmutation and burnup
calculations often fail for systems involving a non-trivial number of nuclides.  To avoid this
problem OrigenReactorMG swaps the ReatocMG burnup methods with a thin wrapper around ORIGEN using
the cross-section data from the original reactor.

All functionality may be found in the ``origen_reactormg`` module::

    import bright.reactormg

.. currentmodule:: bright.origen_reactormg
    
.. autoclass:: OrigenReactorMG(tape9=None, reactor_parameters=None, track_params=None, name="")

    .. automethod:: burnup_core()
    .. automethod:: assemble_transmutation_matrices()
    .. automethod:: calc_transmutation()
