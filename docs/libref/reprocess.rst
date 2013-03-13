.. _bright_reprocess:

***************
Reprocess Class
***************
The reprocessing object is one of the simplest in bright.  In a technology nonspecific way, it performs
separations on an input material.  Reprocessing facilities have nominal separation efficiencies (SE) for 
each element (and possibly nuclide).  These SE are valued on the range [0,1] and represent the fractional 
portion of an element that is recovered in the reprocessing output (mat_prod).

The SE, if left unspecified, default to full recovery (a value of 1).  Currently, Reprocess does not contain 
a mat_tail attribute that corresponds to material not recovered in mat_prod.  

.. currentmodule:: bright.reprocess
    
.. autoclass:: Reprocess(sepeff=None, n="")
    :members:

    .. attribute:: track_params

        For Reprocess, the only parameter that is tracked is the aggregate mass.  
        Thus this attribute is automatically set to ``["Mass"]``.

