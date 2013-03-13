.. _bright_reactormg:

***************
ReactorMG Class
***************
Reactors are the most computationally complex of the nuclear fuel cycle components currently implemented.  
Bright handles nuclear reactors in two distinct object classes: a one neutron energy group methodology 
and a multi-group algorithm (implemented here). 

All multi-group (MG) reactors share a common methodological backbone.  This page describes what is fundamentally 
the same about such reactor objects via the ReactorMG class.  This is a subclass of FCComp.  

The multi-group reactors are based on a algorithm submitted for published by the author to Nuclear 
Engineering & Design and otherwise found in 
`Chapter 5  of the author's dissertation <https://docs.google.com/viewer?a=v&pid=explorer&chrome=true&srcid=0BxUpd34yizZrNTdlZTdlMDYtY2ZhZi00N2RhLWIxM2UtYTFkZDA5NDJhYTEx&hl=en>`_.

**ReactorMG Helper & Child Classes**

.. toctree::
    :maxdepth: 1

    fluence_point
    reactor_parameters
    origen_reactormg

All functionality may be found in the ``reactormg`` module::

    import bright.reactormg

.. currentmodule:: bright.reactormg
    
.. autoclass:: ReactorMG(rp=None, paramtrack=None, n="")
    :members:
