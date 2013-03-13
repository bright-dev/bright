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
    
.. autoclass:: Reactor1G(reactor_parameters=None, track_params=None, name="")

    .. autoattribute:: B
    .. autoattribute:: phi
    .. autoattribute:: fuel_chemical_form
    .. autoattribute:: coolant_chemical_form
    .. autoattribute:: rhoF
    .. autoattribute:: rhoC
    .. autoattribute:: P_NL
    .. autoattribute:: target_BU
    .. autoattribute:: use_zeta
    .. autoattribute:: lattice_flag
    .. autoattribute:: rescale_hydrogen_xs
    .. autoattribute:: r
    .. autoattribute:: l
    .. autoattribute:: S_O
    .. autoattribute:: S_T
    .. autoattribute:: VF
    .. autoattribute:: VC
    .. autoattribute:: libfile
    .. autoattribute:: F
    .. autoattribute:: BUi_F_
    .. autoattribute:: pi_F_
    .. autoattribute:: di_F_
    .. autoattribute:: Tij_F_
    .. autoattribute:: A_IHM
    .. autoattribute:: MWF
    .. autoattribute:: MWC
    .. autoattribute:: niF
    .. autoattribute:: niC
    .. autoattribute:: miF
    .. autoattribute:: miC
    .. autoattribute:: NiF
    .. autoattribute:: NiC
    .. autoattribute:: dF_F_
    .. autoattribute:: dC_F_
    .. autoattribute:: BU_F_
    .. autoattribute:: P_F_
    .. autoattribute:: D_F_
    .. autoattribute:: k_F_
    .. autoattribute:: Mj_F_
    .. autoattribute:: zeta_F_
    .. autoattribute:: fd
    .. autoattribute:: Fd
    .. autoattribute:: BUd
    .. autoattribute:: k
    .. autoattribute:: mat_feed_u
    .. autoattribute:: mat_feed_tru
    .. autoattribute:: mat_feed_lan
    .. autoattribute:: mat_feed_act
    .. autoattribute:: mat_prod_u
    .. autoattribute:: mat_prod_tru
    .. autoattribute:: mat_prod_lan
    .. autoattribute:: mat_prod_act
    .. autoattribute:: deltaR
    .. autoattribute:: tru_cr
    .. autoattribute:: SigmaFa_F_
    .. autoattribute:: SigmaFtr_F_
    .. autoattribute:: kappaF_F_
    .. autoattribute:: SigmaCa_F_
    .. autoattribute:: SigmaCtr_F_
    .. autoattribute:: kappaC_F_
    .. autoattribute:: lattice_E_F_
    .. autoattribute:: lattice_F_F_
    .. automethod:: initialize(reactor_parameters)
    .. automethod:: loadlib(libfile="reactor.h5")
    .. automethod:: fold_mass_weights()
    .. automethod:: calc_Mj_F_()
    .. automethod:: calc_Mj_Fd_()
    .. automethod:: calc_mat_prod()
    .. automethod:: calc_sub_mats()
    .. automethod:: calc_tru_cr()
    .. automethod:: calc_deltaR(input=None)
    .. automethod:: fluence_at_BU(burnup)
    .. automethod:: batch_average(BUd, PDk_flag="K")
    .. automethod:: batch_average_k(BUd)
    .. automethod:: BUd_bisection_method()
    .. automethod:: run_P_NL(pnl)
    .. automethod:: calibrate_P_NL_to_BUd()
    .. automethod:: calc(input=None)
    .. automethod:: lattice_E_planar(a, b)
    .. automethod:: lattice_F_planar(a, b)
    .. automethod:: lattice_E_spherical(a, b)
    .. automethod:: lattice_F_spherical(a, b)
    .. automethod:: lattice_E_cylindrical(a, b)
    .. automethod:: lattice_F_cylindrical(a, b)
    .. automethod:: calc_zeta()
    .. automethod:: calc_zeta_planar()
    .. automethod:: calc_zeta_spherical()
    .. automethod:: calc_zeta_cylindrical()


