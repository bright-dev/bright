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
    
.. autoclass:: ReactorMG(reactor_parameters=None, track_params=None, name="")

    .. autoattribute:: B
    .. autoattribute:: flux
    .. autoattribute:: chemical_form_fuel
    .. autoattribute:: chemical_form_clad
    .. autoattribute:: chemical_form_cool
    .. autoattribute:: rho_fuel
    .. autoattribute:: rho_clad
    .. autoattribute:: rho_cool
    .. autoattribute:: P_NL
    .. autoattribute:: target_BU
    .. autoattribute:: specific_power
    .. autoattribute:: burn_regions
    .. autoattribute:: S
    .. autoattribute:: burn_time
    .. autoattribute:: bt_s
    .. autoattribute:: burn_times
    .. autoattribute:: use_zeta
    .. autoattribute:: lattice_flag
    .. autoattribute:: rescale_hydrogen_xs
    .. autoattribute:: burnup_via_constant
    .. autoattribute:: branch_ratio_cutoff
    .. autoattribute:: r_fuel
    .. autoattribute:: r_void
    .. autoattribute:: r_clad
    .. autoattribute:: pitch
    .. autoattribute:: S_O
    .. autoattribute:: S_T
    .. autoattribute:: V_fuel
    .. autoattribute:: V_clad
    .. autoattribute:: V_cool
    .. autoattribute:: libfile
    .. autoattribute:: I
    .. autoattribute:: J
    .. autoattribute:: K
    .. autoattribute:: K_num
    .. autoattribute:: K_ord
    .. autoattribute:: K_ind
    .. autoattribute:: nperturbations
    .. autoattribute:: perturbed_fields
    .. autoattribute:: G
    .. autoattribute:: E_g
    .. autoattribute:: phi_g
    .. autoattribute:: phi
    .. autoattribute:: Phi
    .. autoattribute:: time0
    .. autoattribute:: BU0
    .. autoattribute:: Ti0
    .. autoattribute:: sigma_t_pg
    .. autoattribute:: sigma_a_pg
    .. autoattribute:: nubar_sigma_f_pg
    .. autoattribute:: chi_pg
    .. autoattribute:: sigma_s_pgh
    .. autoattribute:: sigma_f_pg
    .. autoattribute:: sigma_gamma_pg
    .. autoattribute:: sigma_2n_pg
    .. autoattribute:: sigma_3n_pg
    .. autoattribute:: sigma_alpha_pg
    .. autoattribute:: sigma_proton_pg
    .. autoattribute:: sigma_gamma_x_pg
    .. autoattribute:: sigma_2n_x_pg
    .. autoattribute:: A_HM_t
    .. autoattribute:: MW_fuel_t
    .. autoattribute:: MW_clad_t
    .. autoattribute:: MW_cool_t
    .. autoattribute:: n_fuel_it
    .. autoattribute:: n_clad_it
    .. autoattribute:: n_cool_it
    .. autoattribute:: m_fuel_it
    .. autoattribute:: m_clad_it
    .. autoattribute:: m_cool_it
    .. autoattribute:: N_fuel_it
    .. autoattribute:: N_clad_it
    .. autoattribute:: N_cool_it
    .. autoattribute:: phi_tg
    .. autoattribute:: phi_t
    .. autoattribute:: Phi_t
    .. autoattribute:: BU_t
    .. autoattribute:: zeta_tg
    .. autoattribute:: T_it
    .. autoattribute:: sigma_t_itg
    .. autoattribute:: sigma_a_itg
    .. autoattribute:: nubar_sigma_f_itg
    .. autoattribute:: chi_itg
    .. autoattribute:: sigma_s_itgh
    .. autoattribute:: sigma_f_itg
    .. autoattribute:: sigma_gamma_itg
    .. autoattribute:: sigma_2n_itg
    .. autoattribute:: sigma_3n_itg
    .. autoattribute:: sigma_alpha_itg
    .. autoattribute:: sigma_proton_itg
    .. autoattribute:: sigma_gamma_x_itg
    .. autoattribute:: sigma_2n_x_itg
    .. autoattribute:: Sigma_t_fuel_tg
    .. autoattribute:: Sigma_a_fuel_tg
    .. autoattribute:: nubar_Sigma_f_fuel_tg
    .. autoattribute:: chi_fuel_tg
    .. autoattribute:: Sigma_s_fuel_tgh
    .. autoattribute:: Sigma_f_fuel_tg
    .. autoattribute:: Sigma_gamma_fuel_tg
    .. autoattribute:: Sigma_2n_fuel_tg
    .. autoattribute:: Sigma_3n_fuel_tg
    .. autoattribute:: Sigma_alpha_fuel_tg
    .. autoattribute:: Sigma_proton_fuel_tg
    .. autoattribute:: Sigma_gamma_x_fuel_tg
    .. autoattribute:: Sigma_2n_x_fuel_tg
    .. autoattribute:: kappa_fuel_tg
    .. autoattribute:: Sigma_t_clad_tg
    .. autoattribute:: Sigma_a_clad_tg
    .. autoattribute:: nubar_Sigma_f_clad_tg
    .. autoattribute:: chi_clad_tg
    .. autoattribute:: Sigma_s_clad_tgh
    .. autoattribute:: Sigma_f_clad_tg
    .. autoattribute:: Sigma_gamma_clad_tg
    .. autoattribute:: Sigma_2n_clad_tg
    .. autoattribute:: Sigma_3n_clad_tg
    .. autoattribute:: Sigma_alpha_clad_tg
    .. autoattribute:: Sigma_proton_clad_tg
    .. autoattribute:: Sigma_gamma_x_clad_tg
    .. autoattribute:: Sigma_2n_x_clad_tg
    .. autoattribute:: kappa_clad_tg
    .. autoattribute:: Sigma_t_cool_tg
    .. autoattribute:: Sigma_a_cool_tg
    .. autoattribute:: nubar_Sigma_f_cool_tg
    .. autoattribute:: chi_cool_tg
    .. autoattribute:: Sigma_s_cool_tgh
    .. autoattribute:: Sigma_f_cool_tg
    .. autoattribute:: Sigma_gamma_cool_tg
    .. autoattribute:: Sigma_2n_cool_tg
    .. autoattribute:: Sigma_3n_cool_tg
    .. autoattribute:: Sigma_alpha_cool_tg
    .. autoattribute:: Sigma_proton_cool_tg
    .. autoattribute:: Sigma_gamma_x_cool_tg
    .. autoattribute:: Sigma_2n_x_cool_tg
    .. autoattribute:: kappa_cool_tg
    .. autoattribute:: Sigma_t_tg
    .. autoattribute:: Sigma_a_tg
    .. autoattribute:: nubar_Sigma_f_tg
    .. autoattribute:: chi_tg
    .. autoattribute:: Sigma_s_tgh
    .. autoattribute:: Sigma_f_tg
    .. autoattribute:: Sigma_gamma_tg
    .. autoattribute:: Sigma_2n_tg
    .. autoattribute:: Sigma_3n_tg
    .. autoattribute:: Sigma_alpha_tg
    .. autoattribute:: Sigma_proton_tg
    .. autoattribute:: Sigma_gamma_x_tg
    .. autoattribute:: Sigma_2n_x_tg
    .. autoattribute:: A_tgh
    .. autoattribute:: F_tgh
    .. autoattribute:: A_inv_tgh
    .. autoattribute:: A_inv_F_tgh
    .. autoattribute:: nearest_neighbors
    .. autoattribute:: k_t
    .. autoattribute:: td_n
    .. autoattribute:: td
    .. autoattribute:: BUd
    .. autoattribute:: Phid
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
    .. automethod:: initialize(reactor_parameters)
    .. automethod:: loadlib(libfile="reactor.h5")
    .. automethod:: interpolate_cross_sections()
    .. automethod:: calc_mass_weights()
    .. automethod:: fold_mass_weights()
    .. automethod:: assemble_multigroup_matrices()
    .. automethod:: assemble_transmutation_matrices()
    .. automethod:: calc_criticality()
    .. automethod:: calc_transmutation()
    .. automethod:: init_core()
    .. automethod:: burnup_core()
    .. automethod:: calc_nearest_neighbors()
    .. automethod:: calc_T_itd()
    .. automethod:: calc_mat_prod()
    .. automethod:: calc_sub_mats()
    .. automethod:: calc_tru_cr()
    .. automethod:: fluence_at_BU(burnup)
    .. automethod:: batch_average_k(BUd)
    .. automethod:: BUd_bisection_method()
    .. automethod:: run_P_NL(pnl)
    .. automethod:: calibrate_P_NL_to_BUd()
    .. automethod:: calc(input=None)
