// ReactorMG.h
// Header for general Multi-Group Reactor Model

#if !defined(_BRIGHT_REACTORMG_)
#define _BRIGHT_REACTORMG_


// Boost
#include <boost/math/special_functions/bessel.hpp>

// Bright Libs
#include "fccomp.h"
#include "fluence_point.h"
#include "reactor_parameters.h"

namespace bright {

  /***********************************************/
  /*** ReactorMG Component Class and Functions ***/
  /***********************************************/

  typedef std::vector<double> time_data;

  typedef std::vector<time_data> time_g;

  typedef std::vector<time_g> time_gh;

  typedef std::map<int, std::vector<double> > iso_time_map;
  typedef iso_time_map::iterator iso_time_iter;

  typedef std::map<int, std::vector<time_data> > iso_time_g;
  typedef iso_time_g::iterator iso_time_g_iter;

  typedef std::map<int, std::vector< std::vector<time_data> > > iso_time_gh;
  typedef iso_time_gh::iterator iso_time_gh_iter;

  typedef std::vector<double> pert_data;
  typedef std::vector< std::vector<double> > pert_data_g; 
  typedef std::vector< std::vector< std::vector<double> > > pert_data_gh; 

  typedef std::vector<int> iso_vec;
  typedef std::map<int, int> iso_map;


  class ReactorMG : public FCComp
  {
  /** Reactor class
   *  Basic One-Group Reactor Model.  Computes one Burnup Calculationn with the option of computing output isotopics.
   *  Specific reactor types inherit this class and change base parameters.
   */
  protected:
    // FIXME these should be public eventually, once python bindings are written
    bright::SparseMatrix<double> decay_matrix;
    bright::SparseMatrix<double> thermal_yield_matrix;
    bright::SparseMatrix<double> fast_yield_matrix;
    std::vector< bright::SparseMatrix<double> > fission_product_yield_matrix;

    pert_data_g phi_g;        // Group fluxes

    std::map<int, std::map<int, std::vector< std::vector<int> > > > transmutation_chains;

    std::map<int, pert_data_g> sigma_t_pg;       // Total cross section from data library
    std::map<int, pert_data_g> sigma_a_pg;       // Absorption cross section from data library
    std::map<int, pert_data_g> nubar_sigma_f_pg;   // Neutrons per fission times Fission cross section from data library
    std::map<int, pert_data_g> chi_pg;           // Fission energy spectrum from data library
    std::map<int, pert_data_gh> sigma_s_pgh;     // Group to group scattering cross section from data library
    std::map<int, pert_data_g> sigma_f_pg;       // Fission cross section from data library
    std::map<int, pert_data_g> sigma_gamma_pg;     // Capture cross section from data library
    std::map<int, pert_data_g> sigma_2n_pg;      // (n, 2n) cross section from data library
    std::map<int, pert_data_g> sigma_3n_pg;      // (n, 3n) cross section from data library
    std::map<int, pert_data_g> sigma_alpha_pg;     // (n, alpha) cross section from data library
    std::map<int, pert_data_g> sigma_proton_pg;    // (n, proton) cross section from data library
    std::map<int, pert_data_g> sigma_gamma_x_pg;   // Capture cross section (excited) from data library
    std::map<int, pert_data_g> sigma_2n_x_pg;    // (n, 2n *) cross section from data library

    std::vector< std::vector<double> > branch_ratios;

    h5wrap::HomogenousTypeTable<double> perturbations;  // Load perturbation table

    time_g phi_tg;    // Group fluxes as a function of time
    time_g lattice_E_tg;  // Lattice function E
    time_g lattice_F_tg;  // Lattuce function F

    time_g zeta_tg;     // Group disadvantage factors
    iso_time_g sigma_t_itg;        // Total cross section as a function of nuclide and burn_time
    iso_time_g sigma_a_itg;        // Absorption cross section as a function of nuclide and burn_time
    iso_time_g nubar_sigma_f_itg;    // Neutrons per fission times Fission cross section as a function of nuclide and burn_time
    iso_time_g chi_itg;            // Fission neutron energy spectrum as a function of nuclide and burn_time
    iso_time_gh sigma_s_itgh;      // Group to group scattering cross section as a function of nuclide and burn_time
    iso_time_g sigma_f_itg;        // Fission cross section as a function of nuclide and burn_time
    iso_time_g sigma_gamma_itg;    // Capture cross section as a function of nuclide and burn_time
    iso_time_g sigma_2n_itg;       // (n, 2n) cross section as a function of nuclide and burn_time
    iso_time_g sigma_3n_itg;       // (n, 3n) cross section as a function of nuclide and burn_time
    iso_time_g sigma_alpha_itg;    // (n, alpha) cross section as a function of nuclide and burn_time
    iso_time_g sigma_proton_itg;     // (n, proton) cross section as a function of nuclide and burn_time
    iso_time_g sigma_gamma_x_itg;    // Capture cross section (excited) as a function of nuclide and burn_time
    iso_time_g sigma_2n_x_itg;     // (n, 2n *) cross section as a function of nuclide and burn_time

    time_g Sigma_t_fuel_tg;        // Core-average Macroscopic total cross-section as a function of time and energy group
    time_g Sigma_a_fuel_tg;        // Core-average Macroscopic absorption cross section as a function of time and energy group
    time_g nubar_Sigma_f_fuel_tg;    // Core-average nubar times the Macroscopic fission cross-section as a function of time and energy group
    time_g chi_fuel_tg;            // Core-average Fission neutron energy spectrum as a function of time and energy group
    time_gh Sigma_s_fuel_tgh;      // Core-average Macroscopic scattering kernel cross-section as a function of time
    time_g Sigma_f_fuel_tg;        // Core-average Macroscopic fission cross-section as a function of time and energy group
    time_g Sigma_gamma_fuel_tg;    // Core-average Macroscopic capture cross-section as a function of time and energy group
    time_g Sigma_2n_fuel_tg;       // Core-average Macroscopic (n, 2n) cross-section as a function of time and energy group
    time_g Sigma_3n_fuel_tg;       // Core-average Macroscopic (n, 3n) cross-section as a function of time and energy group
    time_g Sigma_alpha_fuel_tg;    // Core-average Macroscopic (n, alpha) cross-section as a function of time and energy group
    time_g Sigma_proton_fuel_tg;     // Core-average Macroscopic (n, proton) cross-section as a function of time and energy group
    time_g Sigma_gamma_x_fuel_tg;    // Core-average Macroscopic capture cross-section (excited) as a function of time and energy group
    time_g Sigma_2n_x_fuel_tg;     // Core-average Macroscopic (n, 2n *) cross-section as a function of time and energy group
    time_g kappa_fuel_tg;          // Inverse of the diffusion coefficent

    time_g Sigma_t_clad_tg;        // Core-average Macroscopic total cross-section as a function of time and energy group
    time_g Sigma_a_clad_tg;        // Core-average Macroscopic absorption cross section as a function of time and energy group
    time_g nubar_Sigma_f_clad_tg;    // Core-average nubar times the Macroscopic fission cross-section as a function of time and energy group
    time_g chi_clad_tg;            // Core-average Fission neutron energy spectrum as a function of time and energy group
    time_gh Sigma_s_clad_tgh;      // Core-average Macroscopic scattering kernel cross-section as a function of time
    time_g Sigma_f_clad_tg;        // Core-average Macroscopic fission cross-section as a function of time and energy group
    time_g Sigma_gamma_clad_tg;    // Core-average Macroscopic capture cross-section as a function of time and energy group
    time_g Sigma_2n_clad_tg;       // Core-average Macroscopic (n, 2n) cross-section as a function of time and energy group
    time_g Sigma_3n_clad_tg;       // Core-average Macroscopic (n, 3n) cross-section as a function of time and energy group
    time_g Sigma_alpha_clad_tg;    // Core-average Macroscopic (n, alpha) cross-section as a function of time and energy group
    time_g Sigma_proton_clad_tg;     // Core-average Macroscopic (n, proton) cross-section as a function of time and energy group
    time_g Sigma_gamma_x_clad_tg;    // Core-average Macroscopic capture cross-section (excited) as a function of time and energy group
    time_g Sigma_2n_x_clad_tg;     // Core-average Macroscopic (n, 2n *) cross-section as a function of time and energy group
    time_g kappa_clad_tg;          // Inverse of the diffusion coefficent

    time_g Sigma_t_cool_tg;        // Core-average Macroscopic total cross-section as a function of time and energy group
    time_g Sigma_a_cool_tg;        // Core-average Macroscopic absorption cross section as a function of time and energy group
    time_g nubar_Sigma_f_cool_tg;    // Core-average nubar times the Macroscopic fission cross-section as a function of time and energy group
    time_g chi_cool_tg;            // Core-average Fission neutron energy spectrum as a function of time and energy group
    time_gh Sigma_s_cool_tgh;      // Core-average Macroscopic scattering kernel cross-section as a function of time
    time_g Sigma_f_cool_tg;        // Core-average Macroscopic fission cross-section as a function of time and energy group
    time_g Sigma_gamma_cool_tg;    // Core-average Macroscopic capture cross-section as a function of time and energy group
    time_g Sigma_2n_cool_tg;       // Core-average Macroscopic (n, 2n) cross-section as a function of time and energy group
    time_g Sigma_3n_cool_tg;       // Core-average Macroscopic (n, 3n) cross-section as a function of time and energy group
    time_g Sigma_alpha_cool_tg;    // Core-average Macroscopic (n, alpha) cross-section as a function of time and energy group
    time_g Sigma_proton_cool_tg;     // Core-average Macroscopic (n, proton) cross-section as a function of time and energy group
    time_g Sigma_gamma_x_cool_tg;    // Core-average Macroscopic capture cross-section (excited) as a function of time and energy group
    time_g Sigma_2n_x_cool_tg;     // Core-average Macroscopic (n, 2n *) cross-section as a function of time and energy group
    time_g kappa_cool_tg;          // Inverse of the diffusion coefficent

    time_g Sigma_t_tg;        // Core-average Macroscopic total cross-section as a function of time and energy group
    time_g Sigma_a_tg;        // Core-average Macroscopic absorption cross section as a function of time and energy group
    time_g nubar_Sigma_f_tg;    // Core-average nubar times the Macroscopic fission cross-section as a function of time and energy group
    time_g chi_tg;            // Core-average Fission neutron energy spectrum as a function of time and energy group
    time_gh Sigma_s_tgh;      // Core-average Macroscopic scattering kernel cross-section as a function of time
    time_g Sigma_f_tg;        // Core-average Macroscopic fission cross-section as a function of time and energy group
    time_g Sigma_gamma_tg;    // Core-average Macroscopic capture cross-section as a function of time and energy group
    time_g Sigma_2n_tg;       // Core-average Macroscopic (n, 2n) cross-section as a function of time and energy group
    time_g Sigma_3n_tg;       // Core-average Macroscopic (n, 3n) cross-section as a function of time and energy group
    time_g Sigma_alpha_tg;    // Core-average Macroscopic (n, alpha) cross-section as a function of time and energy group
    time_g Sigma_proton_tg;     // Core-average Macroscopic (n, proton) cross-section as a function of time and energy group
    time_g Sigma_gamma_x_tg;    // Core-average Macroscopic capture cross-section (excited) as a function of time and energy group
    time_g Sigma_2n_x_tg;     // Core-average Macroscopic (n, 2n *) cross-section as a function of time and energy group

    time_gh A_fuel_tgh;     // Absorprion Matrix, as a function of time
    time_gh F_fuel_tgh;     // Fission Matrix, as a function of time
    time_gh A_inv_fuel_tgh;   // Inverse of Absorprion Matrix, as a function of time
    time_gh A_inv_F_fuel_tgh; // Inverse of Absorprion Matrix mult by the Fission Matrix, as a function of time

    time_gh A_clad_tgh;     // Absorprion Matrix, as a function of time
    time_gh A_inv_clad_tgh;   // Inverse of Absorprion Matrix, as a function of time

    time_gh A_cool_tgh;     // Absorprion Matrix, as a function of time
    time_gh A_inv_cool_tgh;   // Inverse of Absorprion Matrix, as a function of time

    time_gh A_tgh;     // Absorprion Matrix, as a function of time
    time_gh F_tgh;     // Fission Matrix, as a function of time
    time_gh A_inv_tgh;   // Inverse of Absorprion Matrix, as a function of time
    time_gh A_inv_F_tgh; // Inverse of Absorprion Matrix mult by the Fission Matrix, as a function of time

    std::vector< bright::SparseMatrix<double> > T_int_tij;   // Energy Integral of the Transmutation Matrix, as a function of time
    std::vector< bright::SparseMatrix<double> > M_tij;     // Burnup Matrix, T_int Matrix plus the Decay Matrix, as a function of time

  public:
    // ReactorMG Constructors
    ReactorMG(std::string n="");
    ReactorMG(std::set<std::string> paramtrack, std::string n="");
    ReactorMG(ReactorParameters rp, std::string n="");
    ReactorMG(ReactorParameters rp, std::set<std::string> paramtrack, std::string n="");
    ~ReactorMG();
  
    // Public data

    // Attributes from ReactorParameters (initialization)
    int B;        // Total number of fuel loading batches
    double flux;    // Flux used for Fluence

    std::map<std::string, double> chemical_form_fuel; // Chemical form of Fuel as Dictionary.  Keys are elements or nuclides while values represent mass weights.  Denote heavy metal by key "IHM".
    std::map<std::string, double> chemical_form_clad; // Same a fuel chemical form but for cladding.  Should not have "IHM"
    std::map<std::string, double> chemical_form_cool; // Same a fuel chemical form but for coolant.  Should not have "IHM"

    double rho_fuel;    // Fuel Density
    double rho_clad;    // Cladding Density
    double rho_cool;    // Coolant Density

    double P_NL;          // Non-Leakage Probability
    double target_BU;     // Target Discharge Burnup, only used for graphing inside of this component
    double specific_power;  // The specific power of the fuel
    int burn_regions;     // Number of burn regions for this reactor 
    int S;                // Number of burnup time steps.
    double burn_time;     // Curent burnup time.
    int bt_s;             // Curent burnup time index.  burn_time == burn_times[bt_s]
    time_data burn_times;   // A non-negative, monotonically increasing vector of burnup steps

    bool use_zeta;            // Boolean value on whether or not the disadvantage factor should be used
    std::string lattice_flag;   // lattice_flag Type (Planar || Spherical || Cylindrical)
    bool rescale_hydrogen_xs;   // Rescale the Hydrogen-1 XS?
    std::string burnup_via_constant;
    double branch_ratio_cutoff; // Cut-off for bateman chains

    double r_fuel;  // Fuel region radius
    double r_void;  // Void region radius
    double r_clad;  // Cladding region radius
    double pitch;   // Unit cell side length

    // Attributes calculated from initialization
    double S_O; // Number of open slots in fuel assembly
    double S_T; // Total number of Fuel assembly slots.
    double V_fuel; 	// Fuel Volume
    double V_clad;  // Cladding Volume
    double V_cool;  // Coolant Volume

    // Path where the reactor's HDF5 library
    std::string libfile;

    // Attributes read in from data library
    nuc_set I;   // Set of nuclides that may be in mat_feed.
    nuc_set J;   // Set of nuclides that may be in mat_prod.
    nuc_set K;   // Set of nuclides that is the union of all isos in mat_feed and all isos in nuc_data

    int K_num;
    iso_vec K_ord; // Lowest-to-highest order of J.
    iso_map K_ind; // Lowest-to-highest map of J into matrix position.

    std::vector<double> trans_consts;

    int nperturbations; // number of rows in the pertubtaion table
    std::map<std::string, std::vector<double> > perturbed_fields;  // {field_name: [min, max, delta

    int G;                    // number of energu bins
    std::vector<double> E_g;    // Energy bin boundaries
    pert_data phi;            // Total fluxes
    pert_data Phi;            // Fluence
    pert_data time0;          // initial time steps
    pert_data BU0;            // initial burnups

    std::map<int, pert_data> Ti0;                // Data library's transmutation vector


    // Attributes calculated from fold_mass_weights()
    time_data A_HM_t;     // Atomic weight of IHM
    time_data MW_fuel_t;    // Fuel Molecular Weight
    time_data MW_clad_t;    // Cladding Molecular Weight
    time_data MW_cool_t;    // Coolant Molecular Weight
    iso_time_map n_fuel_it; // Fuel Atom Number Weight
    iso_time_map n_clad_it; // Cladding Atom Number Weight
    iso_time_map n_cool_it; // Coolant Atom Number Weight
    iso_time_map m_fuel_it; // Fuel Mass Weight
    iso_time_map m_clad_it; // Cladding Mass Weight
    iso_time_map m_cool_it; // Coolant Mass Weight
    iso_time_map N_fuel_it; // Fuel Number Density
    iso_time_map N_clad_it; // Cladding Number Density
    iso_time_map N_cool_it; // Coolant Number Density


    // Attributes caluclated from burnup_core()
    time_data phi_t;    // Group fluxes as a function of time
    time_data Phi_t;    // Fluence in [n/kb]
    time_data BU_t;		// Burnup [MWd/kgIHM]


    iso_time_map T_it;             // Transformation Matrix [kg_i/kgIHM]


    // Attribute that denotes the indices of the perturbation table the 
    // are closest to the current state of the reactor (ie densities, burn_time, etc.)
    std::vector<int> nearest_neighbors;


    // applying mass weights to burnup core    
    time_data k_t;    // k(t) --- P(t) / D(t)


    // Attributes that are results of BUd_bisection_method() 
    int td_n;    // Lower index of discharge time
    double td;   // Discharge time
    double BUd;  // Discharge Burnup
    double Phid; // Discharge Fluence
    double k;    // Multiplication factor 


    // Results of calc_sub_streams()
    pyne::Material mat_feed_u;   // Input Uranium pyne::Material
    pyne::Material mat_feed_tru; // Input Transuranic pyne::Material
    pyne::Material mat_feed_lan; // Input Lanthinide pyne::Material
    pyne::Material mat_feed_act; // Input Actinide pyne::Material
    pyne::Material mat_prod_u;   // Output Uranium pyne::Material
    pyne::Material mat_prod_tru; // Output Transuranic pyne::Material
    pyne::Material mat_prod_lan; // Output Lanthinide pyne::Material
    pyne::Material mat_prod_act; // Output Actinide pyne::Material


    // The production rate subtracted by the destruction rate at target_BU
    double deltaR;

    // Transuranic Conversion Ratio
    double tru_cr;


    // Public access functions
    void initialize(ReactorParameters rp);
    void loadlib(std::string lib="Reactor.h5");
    void interpolate_cross_sections();
    void calc_mass_weights();
    void fold_mass_weights();
    void assemble_multigroup_matrices();
    void assemble_transmutation_matrices();
    void calc_criticality();
    void calc_transmutation();

    void init_core();
    void burnup_core();

    void calc_nearest_neighbors();

    void       calc_T_itd();
    void       calc_mat_prod();
    void       calc_sub_mats();
    double     calc_tru_cr();

    FluencePoint fluence_at_BU(double BU);
    double     batch_average_k(double BUd);
    void       BUd_bisection_method();
    void       run_P_NL(double temp_pnl);
    void       calibrate_P_NL_to_BUd();

    pyne::Material calc();
    pyne::Material calc(pyne::comp_map incomp);
    pyne::Material calc(pyne::Material mat);	

    // Lattice functions
    void lattice_E_planar(double a, double b);
    void lattice_F_planar(double a, double b);
    void lattice_E_spherical(double a, double b);
    void lattice_F_spherical(double a, double b);
    void lattice_E_cylindrical(double a, double b);
    void lattice_F_cylindrical(double a, double b);
    void calc_zeta();

    // Transmutation chain functions
    void add_transmutation_chains(std::vector<int> tc);
    double bateman_chain(int i, int j, int c, double t);
    double bateman(int i, int j, double t);
  };

// end bright
};

#endif
