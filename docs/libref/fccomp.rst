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
    :members:
