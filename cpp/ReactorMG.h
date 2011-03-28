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


class ReactorMG : public FCComp
{
/** Reactor class
 *  Basic One-Group Reactor Model.  Computes one Burnup Calculationn with the option of computing output isotopics.
 *  Specific reactor types inherit this class and change base parameters.
 */
protected:
    //Protected data

    //Thermal XS data is read in from static KAERI Data
    //Only read in if the disadvantage factor will be used.
    std::map<int, double> sigma_a_therm;			//Microscopic Thermal Absorption XS 
    std::map<int, double> sigma_s_therm;			//Microscopic Thermal Scattering XS 

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

    std::map<std::string, double> fuel_chemical_form;       // Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
    std::map<std::string, double> coolant_chemical_form;    // Same a fuel chemical form but for coolant.  Should not have "IHM"

    double rho_fuel;    // Fuel Density
    double rho_clad;    // Cladding Density
    double rho_cool;    // Coolant Density

    double P_NL;            // Non-Leakage Probability
    double target_BU;       // Target Discharge Burnup, only used for graphing inside of this component
    double specific_power;  // The specific power of the fuel
    int burn_regions;       // Number of burn regions for this reactor 
    int S;                  // Number of burnup time steps.
    double burn_time;       // Cuurent burnup time.
    int bt_s;               // Cuurent burnup time index.  burn_time == burn_times[bt_s]
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
    double VF; 	// Fuel Volume
    double VC;  // Coolant Volume


    // Path where the reactor's HDF5 library
    std::string libfile;


    // Attributes read in from data library
    iso_set I;   // Set of isotopes that may be in ms_feed.
    iso_set J;   // Set of isotopes that may be in ms_prod.

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

    std::map<int, pert_data> Ti0;               // Data library's transmutation vector
    std::map<int, pert_data_g> sigma_a;         // Absorption cross section from data library
    std::map<int, pert_data_g> sigma_s;         // Scattering cross section from data library
    std::map<int, pert_data_g> sigma_f;         // Fission cross section from data library
    std::map<int, pert_data_g> nubar_sigma_f;   // Neutrons per fission times Fission cross section from data library
    std::map<int, pert_data_g> nubar;           // Neutrons per fission from data library
    std::map<int, pert_data_gh> sigma_s_gh;     // Group to group scattering cross section from data library


    // Attributes calculated from fold_mass_weights()
    double A_IHM;   // Atomic weight of IHM
    double MWF;     // Fuel Molecular Weight
    double MWC;     // Coolant Molecular Weight
    CompDict niF;   // Fuel Atom Number Weight
    CompDict niC;   // Coolant Atom Number Weight
    CompDict miF;   // Fuel Mass Weight
    CompDict miC;   // Coolant Mass Weight
    CompDict NiF;   // Fuel Number Density
    CompDict NiC;   // Coolant Number Density


    // Attributes caluclated from burnup_core()
    time_data Phi_t;    // Fluence in [n/kb]
    time_data BU_t;		// Burnup [MWd/kgIHM]
    time_data pF_t;     // Production rate of the fuel [n/s]
    time_data dF_t;     // Destruction rate of the fuel [n/s]
    time_data dC_t;     // Destruction rate of the coolant [n/s]

    iso_time_map T_it;              // Transformation Matrix [kg_i/kgIHM]
    iso_time_g sigma_a_it;          // Absorption cross section as a function of isotope and burn_time
    iso_time_g sigma_s_it;          // Scattering cross section as a function of isotope and burn_time
    iso_time_g sigma_f_it;          // Fission cross section  as a function of isotope and burn_time
    iso_time_g nubar_sigma_f_it;    // Neutrons per fission times Fission cross section as a function of isotope and burn_time
    iso_time_g nubar_it;            // Neutrons per fission from data library
    iso_time_gh sigma_s_it_gh;      // Group to group scattering cross section as a function of isotope and burn_time


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
    void fold_mass_weights();
    void interpolate_cross_sections();
    void burnup_core();
    void old_burnup_core();

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
