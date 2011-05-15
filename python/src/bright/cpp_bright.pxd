"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

cimport std
cimport cpp_mass_stream
cimport mass_stream

cdef extern from "../../../cpp/bright.h" namespace "bright":
    std.string BRIGHT_DATA

    void bright_start() except +


cdef extern from "../../../cpp/FCComps.h" namespace "FCComps":
    set[int] track_isos
    vector[int] track_isos_order

    void load_track_isos_hdf5(std.string, std.string, bint) except +
    void load_track_isos_text(std.string, bint) except +

    void sort_track_isos()

    int verbosity
    bint write_hdf5
    bint write_text

    std.string output_filename


cdef extern from "../../../cpp/FCComp.h":
    cdef cppclass FCComp:
        # Constructors
        FCComp() except +
        FCComp(std.string) except +
        FCComp(set[std.string], std.string) except +

        # Attributes
        std.string name 
        std.string natural_name

        cpp_mass_stream.MassStream ms_feed
        cpp_mass_stream.MassStream ms_prod

        map[std.string, double] params_prior_calc
        map[std.string, double] params_after_calc

        int pass_num
        set[std.string] track_params

        # Methods
        void calc_params() except +
        void write_ms_pass() except +
        void write_params_pass() except +
        void write_text() except +
        void write_hdf5() except +
        void write() except +
        cpp_mass_stream.MassStream calc() except +


cdef extern from "../../../cpp/Enrichment.h":

    cdef cppclass EnrichmentParameters:
        # Constructors
        EnrichmentParameters() except +

        # Attributes
        double alpha_0
        double Mstar_0

        int j
        int k

        double N0
        double M0

        double xP_j
        double xW_j

    EnrichmentParameters fillUraniumEnrichmentDefaults() except +

    cdef cppclass Enrichment(FCComp): 
        # Constructors
        Enrichment() except +
        Enrichment(std.string) except +
        Enrichment(EnrichmentParameters, std.string) except +

        # Attributes
        double alpha_0
        double Mstar_0
        double Mstar
        cpp_mass_stream.MassStream ms_tail

        int j
        int k
        double xP_j
        double xW_j

        double N
        double M
        double N0
        double M0

        double TotalPerFeed
        double SWUperFeed
        double SWUperProduct

        # Methods
        void initialize(EnrichmentParameters) except +
        void calc_params () except +
        cpp_mass_stream.MassStream calc () except +
        cpp_mass_stream.MassStream calc (map[int, double]) except +
        cpp_mass_stream.MassStream calc (cpp_mass_stream.MassStream) except +

        double PoverF (double, double, double) except +
        double WoverF (double, double, double) except +

        # The following are methods I am too lazy to expose to Python
        # FIXME
        #double get_alphastar_i (double)

        #double get_Ei (double)
        #double get_Si (double)
        #void FindNM()

        #double xP_i(int)
        #double xW_i(int)
        #void SolveNM()
        #void Comp2UnitySecant()
        #void Comp2UnityOther()
        #double deltaU_i_OverG(int)
        #void LoverF()
        #void MstarOptimize()



cdef extern from "../../../cpp/Reprocess.h":

    cdef cppclass Reprocess(FCComp):
        # Constructors
        Reprocess() except +
        Reprocess(map[int, double], std.string) except +

        # Attributes
        map[int, double] sepeff

        # Methods
        void initialize(map[int, double]) except +
        void calc_params() except +
        cpp_mass_stream.MassStream calc() except +
        cpp_mass_stream.MassStream calc(map[int, double]) except +
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream) except +


cdef extern from "../../../cpp/Storage.h":

    cdef cppclass Storage(FCComp):
        # Constructors
        Storage() except +
        Storage(std.string) except +

        # Attributes
        double decay_time

        # Methods
        void calc_params() except +
        cpp_mass_stream.MassStream calc() except +
        cpp_mass_stream.MassStream calc(map[int, double]) except +
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream) except +
        cpp_mass_stream.MassStream calc(double) except +
        cpp_mass_stream.MassStream calc(map[int, double], double) except +
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream, double) except +



cdef extern from "../FluencePoint.h":

    cdef cppclass FluencePoint:
        # Constructors        
        FluencePoint() except +

        # Attributes
        int f
        double F
        double m


cdef extern from "../../../cpp/ReactorParameters.h":

    cdef cppclass ReactorParameters:
        # Constructors        
        ReactorParameters() except +

        # Attributes
        int batches
        double flux

        map[std.string, double] fuel_form
        map[std.string, double] cladding_form
        map[std.string, double] coolant_form

        double fuel_density
        double cladding_density
        double coolant_density

        double pnl
        double BUt
        double specific_power
        int burn_regions
        vector[double] burn_times

        bint use_disadvantage_factor
        std.string lattice_type
        bint rescale_hydrogen

        double fuel_radius
        double void_radius
        double clad_radius
        double unit_cell_pitch

        double open_slots
        double total_slots

    ReactorParameters fill_lwr_defaults() except +

    ReactorParameters fill_fr_defaults() except +


cdef extern from "../../../cpp/Reactor1G.h":

    cdef cppclass Reactor1G(FCComp):
        # Constructors        
        Reactor1G() except +
        Reactor1G(std.string) except +
        Reactor1G(set[std.string], std.string) except +
        Reactor1G(ReactorParameters, std.string) except +
        Reactor1G(ReactorParameters, set[std.string], std.string) except +

        # Attributes
        int B
        double phi
        map[std.string, double] fuel_chemical_form
        map[std.string, double] coolant_chemical_form
        double rhoF
        double rhoC
        double P_NL
        double target_BU
        bint use_zeta
        std.string lattice_flag
        bint rescale_hydrogen_xs

        double r
        double l
        double S_O
        double S_T
        double VF
        double VC

        std.string libfile
        vector[double] F
        map[int, vector[double]] BUi_F_
        map[int, vector[double]] pi_F_
        map[int, vector[double]] di_F_
        map[int, map[int, vector[double]]] Tij_F_

        double A_IHM
        double MWF
        double MWC
        map[int, double] niF
        map[int, double] niC
        map[int, double] miF
        map[int, double] miC
        map[int, double] NiF
        map[int, double] NiC

        vector[double] dF_F_
        vector[double] dC_F_
        vector[double] BU_F_
        vector[double] P_F_
        vector[double] D_F_
        vector[double] k_F_
        map[int, vector[double]] Mj_F_
        vector[double] zeta_F_

        int fd
        double Fd
        double BUd
        double k

        cpp_mass_stream.MassStream ms_feed_u
        cpp_mass_stream.MassStream ms_feed_tru
        cpp_mass_stream.MassStream ms_feed_lan
        cpp_mass_stream.MassStream ms_feed_act
        cpp_mass_stream.MassStream ms_prod_u
        cpp_mass_stream.MassStream ms_prod_tru
        cpp_mass_stream.MassStream ms_prod_lan
        cpp_mass_stream.MassStream ms_prod_act

        double deltaR
        double tru_cr

        vector[double] SigmaFa_F_
        vector[double] SigmaFtr_F_
        vector[double] kappaF_F_

        vector[double] SigmaCa_F_
        vector[double] SigmaCtr_F_
        vector[double] kappaC_F_

        vector[double] lattice_E_F_
        vector[double] lattice_F_F_

        # Methods
        void initialize(ReactorParameters) except +
        void loadlib(std.string) except +
        void fold_mass_weights() except +

        void calc_Mj_F_() except +
        void calc_Mj_Fd_() except +

        void calc_ms_prod() except +
        void calcSubStreams() except +
        double calc_tru_cr() except +

        double calc_deltaR() except +
        double calc_deltaR(map[int, double]) except +
        double calc_deltaR(cpp_mass_stream.MassStream) except +

        FluencePoint fluence_at_BU(double) except +
        double batch_average(double, std.string) except +
        double batch_average_k(double) except +
        void BUd_bisection_method() except +
        void run_P_NL(double) except +
        void calibrate_P_NL_to_BUd() except +

        cpp_mass_stream.MassStream calc() except +
        cpp_mass_stream.MassStream calc(map[int, double]) except +
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream) except +

        void lattice_E_planar(double, double) except +
        void lattice_F_planar(double, double) except +
        void lattice_E_spherical(double, double) except +
        void lattice_F_spherical(double, double) except +
        void lattice_E_cylindrical(double, double) except +
        void lattice_F_cylindrical(double, double) except +

        void calc_zeta() except +
        void calc_zeta_planar() except +
        void calc_zeta_spherical() except +
        void calc_zeta_cylindrical() except +




cdef extern from "../../../cpp/LightWaterReactor1G.h":

    cdef cppclass LightWaterReactor1G(Reactor1G):
        # Constructors        
        LightWaterReactor1G() except +
        LightWaterReactor1G(std.string, std.string) except +
        LightWaterReactor1G(ReactorParameters, std.string) except +
        LightWaterReactor1G(std.string, ReactorParameters, std.string) except +

        # Methods
        void calc_params() except +




cdef extern from "../../../cpp/FastReactor1G.h":

    cdef cppclass FastReactor1G(Reactor1G):
        # Constructors        
        FastReactor1G() except +
        FastReactor1G(std.string, std.string) except +
        FastReactor1G(ReactorParameters, std.string) except +
        FastReactor1G(std.string, ReactorParameters, std.string) except +

        # Methods
        void calc_params() except +



cdef extern from "../../../cpp/FuelFabrication.h":

    cdef cppclass FuelFabrication(FCComp):
        # Constructors        
        FuelFabrication() except +
        FuelFabrication(std.string) except +
        FuelFabrication(set[std.string], std.string) except +
        FuelFabrication(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G, std.string) except +
        FuelFabrication(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G, set[std.string], std.string) except +

        # Attributes
        map[std.string, mass_stream.msp] mass_streams
        map[std.string, double] mass_weights_in
        map[std.string, double] mass_weights_out
        map[std.string, double] deltaRs

        Reactor1G reactor

        # Methods
        void initialize(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G) except +
        void calc_params() except +

        void calc_deltaRs() except +
        cpp_mass_stream.MassStream calc_core_input() except +
        void calc_mass_ratios() except +

        cpp_mass_stream.MassStream calc() except +
        cpp_mass_stream.MassStream calc(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G) except +



cdef extern from "../../../cpp/ReactorMG.h":

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
        vector[int] J_order
        map[int, int] J_index

        vector[vector[double]] decay_matrix
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

        map[int, vector[double]] T_it
        map[int, vector[vector[double]]] sigma_t_itg
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

        vector[vector[double]] Sigma_t_tg
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
        vector[vector[vector[double]]] T_int_tij
        vector[vector[vector[double]]] M_tij

        vector[int] nearest_neighbors

        vector[double] k_t

        int td_n
        double td
        double BUd
        double Phid
        double k

        cpp_mass_stream.MassStream ms_feed_u
        cpp_mass_stream.MassStream ms_feed_tru
        cpp_mass_stream.MassStream ms_feed_lan
        cpp_mass_stream.MassStream ms_feed_act
        cpp_mass_stream.MassStream ms_prod_u
        cpp_mass_stream.MassStream ms_prod_tru
        cpp_mass_stream.MassStream ms_prod_lan
        cpp_mass_stream.MassStream ms_prod_act

        double deltaR
        double tru_cr

        # Methods
        void initialize(ReactorParameters) except +
        void loadlib(std.string) except +
        void interpolate_cross_sections() except +
        void calc_mass_weights() except +
        void fold_mass_weights() except +
        void assemble_multigroup_matrices() except +
        void calc_criticality() except +

        void burnup_core() except +

        void calc_nearest_neighbors() except +

        void calc_T_itd() except +

        void calc_ms_prod() except +
        void calcSubStreams() except +
        double calc_tru_cr() except +

        FluencePoint fluence_at_BU(double) except +
        double batch_average_k(double) except +
        void BUd_bisection_method() except +
        void run_P_NL(double) except +
        void calibrate_P_NL_to_BUd() except +

        cpp_mass_stream.MassStream calc() except +
        cpp_mass_stream.MassStream calc(map[int, double]) except +
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream) except +
