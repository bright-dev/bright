.. _bright_reactor1g:

***************
Reactor1G Class
***************
Reactors are the most computationally complex of the nuclear fuel cycle components currently implemented.  
Bright handles nuclear reactors in two distinct object classes: a one neutron energy group methodology 
(implemented here) and a multi-group algorithm. Moreover, there are several different types of reactors.  
Each type has its own characteristic data library associated with it and takes on type-specific base case 
values.  The reactor types that have been fully implemented thus far are a light water reactor (LWR) and a 
fast reactor (FR).  You may read more about these on their own pages.

All one-group (1G) reactors share a common methodological backbone.  This page describes what is fundamentally 
the same about one group reactor objects via the Reactor1G class.  This is a subclass of FCComp.  Moreover, all 
one-group reactor types have :class:`bright.Reactor1G` as their parent.  The type-specific reactor objects turn 
out to be relatively simple since most of the computational effort is in Reactor1G.

All one energy group reactors are based on a algorithm published by the author in Nuclear Engineering & Design, 
"`A new method for rapid computation of transient fuel cycle material balances <http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6V4D-4W0SSGY-1&_user=10&_coverDate=10/31/2009&_rdoc=1&_fmt=high&_orig=search&_sort=d&_docanchor=&view=c&_acct=C000050221&_version=1&_urlVersion=0&_userid=10&md5=e52ececcd93400f84cc630ba20b01994>`_".

**Reactor1G Helper & Child Classes**

.. toctree::
    :maxdepth: 1

    fluence_point
    reactor_parameters
    light_water_reactor1g
    fast_reactor1g

All functionality may be found in the ``reactor1g`` module::

    import bright.reactor1g

.. currentmodule:: bright.reactor1g
    
.. autoclass:: Reactor1G(rp=None, paramtrack=None, n="")
    :members:


