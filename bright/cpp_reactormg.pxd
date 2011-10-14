"""Cython header for  library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

from pyne cimport std
from pyne cimport cpp_material
from pyne cimport material

from bright.cpp_fccomp cimport FCComp
from bright.cpp_reactor_parameters cimport ReactorParameters
from bright.cpp_fluence_point cimport FluencePoint

cdef extern from "reactormg.h" namespace "bright":

    cdef cppclass ReactorMG(FCComp):
        # Constructors        
        ReactorMG() except +
        ReactorMG(std.string) except +
        ReactorMG(set[std.string], std.string) except +
        ReactorMG(ReactorParameters, std.string) except +
        ReactorMG(ReactorParameters, set[std.string], std.string) except +

        # Attributes
        int B
        double flux

        map[std.string, double] chemical_form_fuel
        map[std.string, double] chemical_form_clad
        map[std.string, double] chemical_form_cool

        double rho_fuel
        double rho_clad
        double rho_cool

        double P_NL
        double target_BU
        double specific_power
        int burn_regions
        int S
        double burn_time
        int bt_s
        vector[double] burn_times

        bint use_zeta
        std.string lattice_flag
        bint rescale_hydrogen_xs
        std.string burnup_via_constant
        double branch_ratio_cutoff

        double r_fuel
        double r_void
        double r_clad
        double pitch

        double S_O
        double S_T
        double V_fuel
        double V_clad
        double V_cool

        std.string libfile

        set[int] I
        set[int] J
        set[int] K

        int K_num
        vector[int] K_ord
        map[int, int] K_ind

#        vector[vector[double]] decay_matrix
        vector[vector[double]] thermal_yield_matrix
        vector[vector[double]] fast_yield_matrix
        vector[vector[vector[double]]] fission_product_yield_matrix

        # Perturbation table goes here
        int nperturbations
        map[std.string, vector[double]] perturbed_fields

        int G
        vector[double] E_g
        vector[vector[double]] phi_g
        vector[double] phi
        vector[double] Phi
        vector[double] time0
        vector[double] BU0

        map[int, vector[double]] Ti0
        map[int, vector[vector[double]]] sigma_t_pg
        map[int, vector[vector[double]]] sigma_a_pg
        map[int, vector[vector[double]]] nubar_sigma_f_pg
        map[int, vector[vector[double]]] chi_pg
        map[int, vector[vector[vector[double]]]] sigma_s_pgh
        map[int, vector[vector[double]]] sigma_f_pg
        map[int, vector[vector[double]]] sigma_gamma_pg
        map[int, vector[vector[double]]] sigma_2n_pg
        map[int, vector[vector[double]]] sigma_3n_pg
        map[int, vector[vector[double]]] sigma_alpha_pg
        map[int, vector[vector[double]]] sigma_proton_pg
        map[int, vector[vector[double]]] sigma_gamma_x_pg
        map[int, vector[vector[double]]] sigma_2n_x_pg

        vector[double] A_HM_t
        vector[double] MW_fuel_t
        vector[double] MW_clad_t
        vector[double] MW_cool_t
        map[int, vector[double]] n_fuel_it
        map[int, vector[double]] n_clad_it
        map[int, vector[double]] n_cool_it
        map[int, vector[double]] m_fuel_it
        map[int, vector[double]] m_clad_it
        map[int, vector[double]] m_cool_it
        map[int, vector[double]] N_fuel_it
        map[int, vector[double]] N_clad_it
        map[int, vector[double]] N_cool_it

        vector[vector[double]] phi_tg
        vector[double] phi_t
        vector[double] Phi_t
        vector[double] BU_t

        vector[vector[double]] zeta_tg
        vector[vector[double]] lattice_E_tg
        vector[vector[double]] lattice_F_tg

        map[int, vector[double]] T_it
        map[int, vector[vector[double]]] sigma_t_itg
        map[int, vector[vector[double]]] sigma_a_itg
        map[int, vector[vector[double]]] nubar_sigma_f_itg
        map[int, vector[vector[double]]] chi_itg
        map[int, vector[vector[vector[double]]]] sigma_s_itgh
        map[int, vector[vector[double]]] sigma_f_itg
        map[int, vector[vector[double]]] sigma_gamma_itg
        map[int, vector[vector[double]]] sigma_2n_itg
        map[int, vector[vector[double]]] sigma_3n_itg
        map[int, vector[vector[double]]] sigma_alpha_itg
        map[int, vector[vector[double]]] sigma_proton_itg
        map[int, vector[vector[double]]] sigma_gamma_x_itg
        map[int, vector[vector[double]]] sigma_2n_x_itg

        vector[vector[double]] Sigma_t_fuel_tg
        vector[vector[double]] Sigma_a_fuel_tg
        vector[vector[double]] nubar_Sigma_f_fuel_tg
        vector[vector[double]] chi_fuel_tg
        vector[vector[vector[double]]] Sigma_s_fuel_tgh
        vector[vector[double]] Sigma_f_fuel_tg
        vector[vector[double]] Sigma_gamma_fuel_tg
        vector[vector[double]] Sigma_2n_fuel_tg
        vector[vector[double]] Sigma_3n_fuel_tg
        vector[vector[double]] Sigma_alpha_fuel_tg
        vector[vector[double]] Sigma_proton_fuel_tg
        vector[vector[double]] Sigma_gamma_x_fuel_tg
        vector[vector[double]] Sigma_2n_x_fuel_tg
        vector[vector[double]] kappa_fuel_tg

        vector[vector[double]] Sigma_t_clad_tg
        vector[vector[double]] Sigma_a_clad_tg
        vector[vector[double]] nubar_Sigma_f_clad_tg
        vector[vector[double]] chi_clad_tg
        vector[vector[vector[double]]] Sigma_s_clad_tgh
        vector[vector[double]] Sigma_f_clad_tg
        vector[vector[double]] Sigma_gamma_clad_tg
        vector[vector[double]] Sigma_2n_clad_tg
        vector[vector[double]] Sigma_3n_clad_tg
        vector[vector[double]] Sigma_alpha_clad_tg
        vector[vector[double]] Sigma_proton_clad_tg
        vector[vector[double]] Sigma_gamma_x_clad_tg
        vector[vector[double]] Sigma_2n_x_clad_tg
        vector[vector[double]] kappa_clad_tg

        vector[vector[double]] Sigma_t_cool_tg
        vector[vector[double]] Sigma_a_cool_tg
        vector[vector[double]] nubar_Sigma_f_cool_tg
        vector[vector[double]] chi_cool_tg
        vector[vector[vector[double]]] Sigma_s_cool_tgh
        vector[vector[double]] Sigma_f_cool_tg
        vector[vector[double]] Sigma_gamma_cool_tg
        vector[vector[double]] Sigma_2n_cool_tg
        vector[vector[double]] Sigma_3n_cool_tg
        vector[vector[double]] Sigma_alpha_cool_tg
        vector[vector[double]] Sigma_proton_cool_tg
        vector[vector[double]] Sigma_gamma_x_cool_tg
        vector[vector[double]] Sigma_2n_x_cool_tg
        vector[vector[double]] kappa_cool_tg

        vector[vector[double]] Sigma_t_tg
        vector[vector[double]] Sigma_a_tg
        vector[vector[double]] nubar_Sigma_f_tg
        vector[vector[double]] chi_tg
        vector[vector[vector[double]]] Sigma_s_tgh
        vector[vector[double]] Sigma_f_tg
        vector[vector[double]] Sigma_gamma_tg
        vector[vector[double]] Sigma_2n_tg
        vector[vector[double]] Sigma_3n_tg
        vector[vector[double]] Sigma_alpha_tg
        vector[vector[double]] Sigma_proton_tg
        vector[vector[double]] Sigma_gamma_x_tg
        vector[vector[double]] Sigma_2n_x_tg

        vector[vector[vector[double]]] A_tgh
        vector[vector[vector[double]]] F_tgh
        vector[vector[vector[double]]] A_inv_tgh
        vector[vector[vector[double]]] A_inv_F_tgh
#        vector[vector[vector[double]]] T_int_tij
#        vector[vector[vector[double]]] M_tij

        vector[int] nearest_neighbors

        vector[double] k_t

        int td_n
        double td
        double BUd
        double Phid
        double k

        cpp_material.Material mat_feed_u
        cpp_material.Material mat_feed_tru
        cpp_material.Material mat_feed_lan
        cpp_material.Material mat_feed_act
        cpp_material.Material mat_prod_u
        cpp_material.Material mat_prod_tru
        cpp_material.Material mat_prod_lan
        cpp_material.Material mat_prod_act

        double deltaR
        double tru_cr

        # Methods
        void initialize(ReactorParameters) except +
        void loadlib(std.string) except +
        void interpolate_cross_sections() except +
        void calc_mass_weights() except +
        void fold_mass_weights() except +
        void assemble_multigroup_matrices() except +
        void assemble_transmutation_matrices() except +
        void calc_criticality() except +
        void calc_transmutation() except +

        void init_core() except +
        void burnup_core() except +

        void calc_nearest_neighbors() except +

        void calc_T_itd() except +

        void calc_mat_prod() except +
        void calcSubStreams() except +
        double calc_tru_cr() except +

        FluencePoint fluence_at_BU(double) except +
        double batch_average_k(double) except +
        void BUd_bisection_method() except +
        void run_P_NL(double) except +
        void calibrate_P_NL_to_BUd() except +

        cpp_material.Material calc() except +
        cpp_material.Material calc(map[int, double]) except +
        cpp_material.Material calc(cpp_material.Material) except +
