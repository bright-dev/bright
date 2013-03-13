.. _bright_fccomp:

************
FCComp Class
************
This is the central class to all bright fuel cycle components.  The premise of FCComp is 
that all components take a material input, perform some operation on the input, and most 
then return a material as output.  This class defines the elements common to all fuel cycle 
components.  

All functionality may be found in the ``fccomp`` module::

    import bright.fccomp

.. currentmodule:: bright.fccomp

.. autoclass:: FCComp(paramlist=None, name="")

    .. autoattribute:: name
    .. autoattribute:: natural_name
    .. autoattribute:: mat_feed
    .. autoattribute:: mat_prod
    .. autoattribute:: params_prior_calc
    .. autoattribute:: params_after_calc
    .. autoattribute:: pass_num
    .. autoattribute:: track_params

    .. automethod:: calc(input=None)
    .. automethod:: calc_params()
    .. automethod:: write()
    .. automethod:: write_text()
    .. automethod:: write_hdf5()
    .. automethod:: write_mat_pass()
    .. automethod:: write_params_pass()
