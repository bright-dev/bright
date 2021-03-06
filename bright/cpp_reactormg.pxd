################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################


from bright cimport cpp_fccomp
from bright cimport cpp_fluence_point
from bright cimport cpp_reactor_parameters
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as cpp_vector
from pyne cimport cpp_material

cdef extern from "reactormg.h" namespace "bright":

    cdef cppclass ReactorMG(cpp_fccomp.FCComp):
        # constructors
        ReactorMG() except +
        ReactorMG(std_string) except +
        ReactorMG(cpp_set[std_string]) except +
        ReactorMG(cpp_set[std_string], std_string) except +
        ReactorMG(cpp_reactor_parameters.ReactorParameters) except +
        ReactorMG(cpp_reactor_parameters.ReactorParameters, std_string) except +
        ReactorMG(cpp_reactor_parameters.ReactorParameters, cpp_set[std_string]) except +
        ReactorMG(cpp_reactor_parameters.ReactorParameters, cpp_set[std_string], std_string) except +

        # attributes
        cpp_vector[double] A_HM_t
        int B
        cpp_vector[double] BU0
        cpp_vector[double] BU_t
        double BUd
        cpp_vector[double] E_g
        int G
        cpp_set[int] I
        cpp_set[int] J
        cpp_set[int] K
        cpp_map[int, int] K_ind
        int K_num
        cpp_vector[int] K_ord
        cpp_vector[double] MW_clad_t
        cpp_vector[double] MW_cool_t
        cpp_vector[double] MW_fuel_t
        cpp_map[int, cpp_vector[double]] N_clad_it
        cpp_map[int, cpp_vector[double]] N_cool_it
        cpp_map[int, cpp_vector[double]] N_fuel_it
        double P_NL
        cpp_vector[double] Phi
        cpp_vector[double] Phi_t
        double Phid
        int S
        double S_O
        double S_T
        cpp_map[int, cpp_vector[double]] T_it
        cpp_map[int, cpp_vector[double]] Ti0
        double V_clad
        double V_cool
        double V_fuel
        double branch_ratio_cutoff
        int bt_s
        int burn_regions
        double burn_time
        cpp_vector[double] burn_times
        std_string burnup_via_constant
        cpp_map[std_string, double] chemical_form_clad
        cpp_map[std_string, double] chemical_form_cool
        cpp_map[std_string, double] chemical_form_fuel
        double deltaR
        double flux
        double k
        cpp_vector[double] k_t
        std_string lattice_flag
        std_string libfile
        cpp_map[int, cpp_vector[double]] m_clad_it
        cpp_map[int, cpp_vector[double]] m_cool_it
        cpp_map[int, cpp_vector[double]] m_fuel_it
        cpp_material.Material mat_feed_act
        cpp_material.Material mat_feed_lan
        cpp_material.Material mat_feed_tru
        cpp_material.Material mat_feed_u
        cpp_material.Material mat_prod_act
        cpp_material.Material mat_prod_lan
        cpp_material.Material mat_prod_tru
        cpp_material.Material mat_prod_u
        cpp_map[int, cpp_vector[double]] n_clad_it
        cpp_map[int, cpp_vector[double]] n_cool_it
        cpp_map[int, cpp_vector[double]] n_fuel_it
        cpp_vector[int] nearest_neighbors
        int nperturbations
        cpp_map[std_string, cpp_vector[double]] perturbed_fields
        cpp_vector[double] phi
        cpp_vector[double] phi_t
        double pitch
        double r_clad
        double r_fuel
        double r_void
        bint rescale_hydrogen_xs
        double rho_clad
        double rho_cool
        double rho_fuel
        double specific_power
        double target_BU
        double td
        int td_n
        cpp_vector[double] time0
        cpp_vector[double] trans_consts
        double tru_cr
        bint use_zeta

        # methods
        void BUd_bisection_method() except +
        void add_transmutation_chains(cpp_vector[int]) except +
        void assemble_multigroup_matrices() except +
        void assemble_transmutation_matrices() except +
        double batch_average_k(double) except +
        double bateman(int, int, double) except +
        double bateman_chain(int, int, int, double) except +
        void burnup_core() except +
        cpp_material.Material calc() except +
        cpp_material.Material calc(cpp_map[int, double]) except +
        cpp_material.Material calc(cpp_material.Material) except +
        void calc_T_itd() except +
        void calc_criticality() except +
        void calc_mass_weights() except +
        void calc_mat_prod() except +
        void calc_nearest_neighbors() except +
        void calc_sub_mats() except +
        void calc_transmutation() except +
        double calc_tru_cr() except +
        void calc_zeta() except +
        void calibrate_P_NL_to_BUd() except +
        cpp_fluence_point.FluencePoint fluence_at_BU(double) except +
        void fold_mass_weights() except +
        void init_core() except +
        void initialize(cpp_reactor_parameters.ReactorParameters) except +
        void interpolate_cross_sections() except +
        void lattice_E_cylindrical(double, double) except +
        void lattice_E_planar(double, double) except +
        void lattice_E_spherical(double, double) except +
        void lattice_F_cylindrical(double, double) except +
        void lattice_F_planar(double, double) except +
        void lattice_F_spherical(double, double) except +
        void loadlib() except +
        void loadlib(std_string) except +
        void run_P_NL(double) except +
        pass




