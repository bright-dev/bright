// ReactorMG.h
// Header for general Multi-Group Reactor Model

#if !defined(_Bright_ReactorMG_)
#define _Bright_ReactorMG_

//C++ stdlib
#include <map>
#include <set>
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <algorithm>

//HDF5
//#include "hdf5.h"
//#include "hdf5_hl.h"
#include "H5Cpp.h"
#include "h5wrap.h"

//Boost
#include <boost/math/special_functions/bessel.hpp>

//Bright Libs
#include "FCComp.h"
#include "MassStream.h"
#include "isoname.h"

#include "FluencePoint.h"
#include "ReactorParameters.h"

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

typedef std::set<int> iso_set;
typedef iso_set::iterator iso_iter;

typedef std::vector<int> iso_vec;
typedef std::map<int, int> iso_map;


class ReactorMG : public FCComp
{
/** Reactor class
 *  Basic One-Group Reactor Model.  Computes one Burnup Calculationn with the option of computing output isotopics.
 *  Specific reactor types inherit this class and change base parameters.
 */

public:
    // ReactorMG Constructors
    ReactorMG ();
    ReactorMG (std::string);
    ReactorMG (std::set<std::string>, std::string = "");
    ReactorMG (ReactorParameters, std::string = "");
    ReactorMG (ReactorParameters, std::set<std::string>, std::string = "");
    ~ReactorMG ();
    
    // Public data

    // Attributes from ReactorParameters (initialization)
    int B;          // Total number of fuel loading batches
    double flux;    // Flux used for Fluence

    std::map<std::string, double> chemical_form_fuel; // Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
    std::map<std::string, double> chemical_form_clad; // Same a fuel chemical form but for cladding.  Should not have "IHM"
    std::map<std::string, double> chemical_form_cool; // Same a fuel chemical form but for coolant.  Should not have "IHM"

    double rho_fuel;    // Fuel Density
    double rho_clad;    // Cladding Density
    double rho_cool;    // Coolant Density

    double P_NL;            // Non-Leakage Probability
    double target_BU;       // Target Discharge Burnup, only used for graphing inside of this component
    double specific_power;  // The specific power of the fuel
    int burn_regions;       // Number of burn regions for this reactor 
    int S;                  // Number of burnup time steps.
    double burn_time;       // Curent burnup time.
    int bt_s;               // Curent burnup time index.  burn_time == burn_times[bt_s]
    time_data burn_times;   // A non-negative, monotonically increasing vector of burnup steps

    bool use_zeta;              // Boolean value on whether or not the disadvantage factor should be used
    std::string lattice_flag;   // lattice_flag Type (Planar || Spherical || Cylindrical)
    bool rescale_hydrogen_xs;   // Rescale the Hydrogen-1 XS?

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
    iso_set I;   // Set of isotopes that may be in ms_feed.
    iso_set J;   // Set of isotopes that may be in ms_prod.

    int J_size;
    iso_vec J_order; // Lowest-to-highest order of J.
    iso_map J_index; // Lowest-to-highest map of J into matrix position.


    std::vector< std::vector<double> > decay_matrix;
    std::vector< std::vector<double> > thermal_yield_matrix;
    std::vector< std::vector<double> > fast_yield_matrix;
    std::vector< std::vector< std::vector<double> > > fission_product_yield_matrix;


    h5wrap::HomogenousTypeTable<double> perturbations;  // Load perturbation table
    int nperturbations; // number of rows in the pertubtaion table
    std::map<std::string, std::vector<double> > perturbed_fields;  // {field_name: [min, max, delta

    int G;                      // number of energu bins
    std::vector<double> E_g;    // Energy bin boundaries
    pert_data_g phi_g;          // Group fluxes
    pert_data phi;              // Total fluxes
    pert_data Phi;              // Fluence
    pert_data time0;            // initial time steps
    pert_data BU0;              // initial burnups

    std::map<int, pert_data> Ti0;                  // Data library's transmutation vector
    std::map<int, pert_data_g> sigma_t_pg;         // Total cross section from data library
    std::map<int, pert_data_g> nubar_sigma_f_pg;   // Neutrons per fission times Fission cross section from data library
    std::map<int, pert_data_g> chi_pg;             // Fission energy spectrum from data library
    std::map<int, pert_data_gh> sigma_s_pgh;       // Group to group scattering cross section from data library
    std::map<int, pert_data_g> sigma_f_pg;         // Fission cross section from data library
    std::map<int, pert_data_g> sigma_gamma_pg;     // Capture cross section from data library
    std::map<int, pert_data_g> sigma_2n_pg;        // (n, 2n) cross section from data library
    std::map<int, pert_data_g> sigma_3n_pg;        // (n, 3n) cross section from data library
    std::map<int, pert_data_g> sigma_alpha_pg;     // (n, alpha) cross section from data library
    std::map<int, pert_data_g> sigma_proton_pg;    // (n, proton) cross section from data library
    std::map<int, pert_data_g> sigma_gamma_x_pg;   // Capture cross section (excited) from data library
    std::map<int, pert_data_g> sigma_2n_x_pg;      // (n, 2n *) cross section from data library


    // Attributes calculated from fold_mass_weights()
    time_data A_HM_t;       // Atomic weight of IHM
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
    time_g phi_tg;      // Group fluxes as a function of time
    time_data phi_t;    // Group fluxes as a function of time
    time_data Phi_t;    // Fluence in [n/kb]
    time_data BU_t;		// Burnup [MWd/kgIHM]

    iso_time_map T_it;               // Transformation Matrix [kg_i/kgIHM]
    iso_time_g sigma_t_itg;          // Total cross section as a function of isotope and burn_time
    iso_time_g nubar_sigma_f_itg;    // Neutrons per fission times Fission cross section as a function of isotope and burn_time
    iso_time_g chi_itg;              // Fission neutron energy spectrum as a function of isotope and burn_time
    iso_time_gh sigma_s_itgh;        // Group to group scattering cross section as a function of isotope and burn_time
    iso_time_g sigma_f_itg;          // Fission cross section as a function of isotope and burn_time
    iso_time_g sigma_gamma_itg;      // Capture cross section as a function of isotope and burn_time
    iso_time_g sigma_2n_itg;         // (n, 2n) cross section as a function of isotope and burn_time
    iso_time_g sigma_3n_itg;         // (n, 3n) cross section as a function of isotope and burn_time
    iso_time_g sigma_alpha_itg;      // (n, alpha) cross section as a function of isotope and burn_time
    iso_time_g sigma_proton_itg;     // (n, proton) cross section as a function of isotope and burn_time
    iso_time_g sigma_gamma_x_itg;    // Capture cross section (excited) as a function of isotope and burn_time
    iso_time_g sigma_2n_x_itg;       // (n, 2n *) cross section as a function of isotope and burn_time


    time_g Sigma_t_tg;          // Core-average Macroscopic total cross-section as a function of time and energy group
    time_g nubar_Sigma_f_tg;    // Core-average nubar times the Macroscopic fission cross-section as a function of time and energy group
    time_g chi_tg;              // Core-average Fission neutron energy spectrum as a function of time and energy group
    time_gh Sigma_s_tgh;        // Core-average Macroscopic scattering kernel cross-section as a function of time
    time_g Sigma_f_tg;          // Core-average Macroscopic fission cross-section as a function of time and energy group
    time_g Sigma_gamma_tg;      // Core-average Macroscopic capture cross-section as a function of time and energy group
    time_g Sigma_2n_tg;         // Core-average Macroscopic (n, 2n) cross-section as a function of time and energy group
    time_g Sigma_3n_tg;         // Core-average Macroscopic (n, 3n) cross-section as a function of time and energy group
    time_g Sigma_alpha_tg;      // Core-average Macroscopic (n, alpha) cross-section as a function of time and energy group
    time_g Sigma_proton_tg;     // Core-average Macroscopic (n, proton) cross-section as a function of time and energy group
    time_g Sigma_gamma_x_tg;    // Core-average Macroscopic capture cross-section (excited) as a function of time and energy group
    time_g Sigma_2n_x_tg;       // Core-average Macroscopic (n, 2n *) cross-section as a function of time and energy group

    time_gh A_tgh;       // Absorprion Matrix, as a function of time
    time_gh F_tgh;       // Fission Matrix, as a function of time
    time_gh A_inv_tgh;   // Inverse of Absorprion Matrix, as a function of time
    time_gh A_inv_F_tgh; // Inverse of Absorprion Matrix mult by the Fission Matrix, as a function of time
    time_gh T_int_tij;   // Energy Integral of the Transmutation Matrix, as a function of time
    time_gh M_tij;       // Burnup Matrix, T_int Matrix plus the Decay Matrix, as a function of time




    // Attribute that denotes the indices of the perturbation table the 
    // are closest to the current state of the reactor (ie densities, burn_time, etc.)
    std::vector<int> nearest_neighbors;


    // applying mass weights to burnup core    
    time_data P_t;      // Full-core Production Rate P(t)
    time_data D_t;      // Full-Core Destruction Rate D(t)
    time_data k_t;      // k(t) --- P(t) / D(t)
    time_data zeta_t;   // Disadvantage Factor zeta(t)


    // Attributes that are results of BUd_bisection_method() 
    int td_n;   // Lower index of discharge time
    double td;  // Discharge time
    double BUd; // Discharge Burnup
    double k;   // Multiplication factor 


    // Results of calc_sub_streams()
    MassStream ms_feed_u;   // Input Uranium MassStream
    MassStream ms_feed_tru; // Input Transuranic MassStream
    MassStream ms_feed_lan; // Input Lanthinide MassStream
    MassStream ms_feed_act; // Input Actinide MassStream
    MassStream ms_prod_u;   // Output Uranium MassStream
    MassStream ms_prod_tru; // Output Transuranic MassStream
    MassStream ms_prod_lan; // Output Lanthinide MassStream
    MassStream ms_prod_act; // Output Actinide MassStream


    // The production rate subtracted by the destruction rate at target_BU
    double deltaR;


    // Transuranic Conversion Ratio
    double tru_cr;


    // Attributes related to disadvantage factor calculation
    time_data SigmaF_at;    // Fuel Macro Absorption XS, Sigma^F_a(F)
    time_data SigmaF_trt;   // Fuel Macro Transport XS, Sigma^F_tr(F)
    time_data kappaF_t;     // Fuel kappa, kappa^F(F)

    time_data SigmaC_at;    // Coolant Macro Absorption XS, Sigma^C_a(F)
    time_data SigmaC_trt;   // Coolant Macro Transport XS, Sigma^C_tr(F)
    time_data kappaC_t;     // Coolant kappa, kappa^C(F)

    time_data lattice_E_t;  // Values for lattice function E(F)
    time_data lattice_F_t;  // Values for lattice function F(F)


    // Public access functions
    void initialize(ReactorParameters);
    void loadlib(std::string libfile = "Reactor.h5");
    void interpolate_cross_sections();
    void calc_mass_weights();
    void fold_mass_weights();
    void assemble_multigroup_matrices();
    void calc_criticality();
    void calc_transmutation();

    void burnup_core();

    void calc_nearest_neighbors();

    void         calc_T_itd();
    void         calc_ms_prod();
    void         calcSubStreams();
    double       calc_tru_cr();

    double       calc_deltaR();
    double       calc_deltaR(CompDict);
    double       calc_deltaR(MassStream);

    FluencePoint fluence_at_BU(double);
    double       batch_average(double, std::string = "K");
    double       batch_average_k(double);
    void         BUd_bisection_method();
    void         run_P_NL(double);
    void         calibrate_P_NL_to_BUd();

    MassStream   calc ();
    MassStream   calc (CompDict);
    MassStream   calc (MassStream);	

    void lattice_E_planar(double, double);
    void lattice_F_planar(double, double);
    void lattice_E_spherical(double, double);
    void lattice_F_spherical(double, double);
    void lattice_E_cylindrical(double, double);
    void lattice_F_cylindrical(double, double);

    void         calc_zeta();
    void         calc_zeta_planar();
    void         calc_zeta_spherical();
    void         calc_zeta_cylindrical();

};


#endif
