// Reactor1G.h
// Header for general One Group Reactor Model

#if !defined(_BRIGHT_REACTOR1G_)
#define _BRIGHT_REACTOR1G_

// Boost
#include "boost/math/special_functions/bessel.hpp"

// Bright Libs
#include "fccomp.h"
#include "fluence_point.h"
#include "reactor_parameters.h"

namespace bright {

  /***********************************************/
  /*** Reactor1G Component Class and Functions ***/
  /***********************************************/

  typedef std::map<int, std::vector<double> > nuc_fluence_dict;
  typedef nuc_fluence_dict::iterator nuc_fluence_iter;

  typedef std::vector<double> data_F_;

  class Reactor1G : public FCComp
  {
  /** Reactor class
   *  Basic One-Group Reactor Model.  Computes one Burnup Calculationn with the option of computing output isotopics.
   *  Specific reactor types inherit this class and change base parameters.
   */
  protected:
    // Protected data
    nuc_set I;  // Set of isotopes that may be in mat_feed.
    nuc_set J;  // Set of isotopes that may be in mat_prod.

    // Thermal XS data is read in from static KAERI Data
    // Only read in if the disadvantage factor will be used.
    std::map<int, double> sigma_a_therm;  // Microscopic Thermal Absorption XS 
    std::map<int, double> sigma_s_therm;  // Microscopic Thermal Scattering XS 

  public:
    // Reactor1G Constructors
    Reactor1G ();
    Reactor1G (std::string);
    Reactor1G (std::set<std::string>, std::string = "");
    Reactor1G (ReactorParameters, std::string = "");
    Reactor1G (ReactorParameters, std::set<std::string>, std::string = "");
    ~Reactor1G ();
    
    // Public data
    int B;      // Total number of fuel loading batches
    double phi; // Flux used for Fluence
    std::map<std::string, double> fuel_chemical_form;     // Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
    std::map<std::string, double> coolant_chemical_form;  // Same a fuel chemical form but for coolant.  Should not have "IHM"
    double rhoF;              // Fuel Density
    double rhoC;              // Coolant Density
    double P_NL; 		          // Non-Leakage Probability
    double target_BU;         // Target Discharge Burnup, only used for graphing inside of this component
    bool use_zeta;            // Boolean value on whether or not the disadvantage factor should be used
    std::string lattice_flag; // lattice_flag Type (Planar || Spherical || Cylindrical)
    bool rescale_hydrogen_xs; // Rescale the Hydrogen-1 XS?

    double r;   // Fuel region radius
    double l;   // Unit cell side length
    double S_O; // Number of open slots in fuel assembly
    double S_T; // Total number of Fuel assembly slots.
    double VF;  // Fuel Volume
    double VC;  // Coolant Volume

    std::string libfile;    // Path where the reactor's HDF5 library
    std::vector<double> F;  // Fluence in [n/kb]
    nuc_fluence_dict BUi_F_;  // Burnup [MWd/kgIHM]
    nuc_fluence_dict pi_F_;   // Production rate [n/s]
    nuc_fluence_dict di_F_;   // Destruction rate [n/s]
    std::map<int, nuc_fluence_dict> Tij_F_; //T ransformation Matrix [kg_i/kgIHM]

    double A_IHM; // Atomic weight of IHM
    double MWF;   // Fuel Molecular Weight
    double MWC;   // Coolant Molecular Weight
    pyne::comp_map niF; // Fuel Atom Number Weight
    pyne::comp_map niC; // Coolant Atom Number Weight
    pyne::comp_map miF; // Fuel Mass Weight
    pyne::comp_map miC; // Coolant Mass Weight
    pyne::comp_map NiF; // Fuel Number Density
    pyne::comp_map NiC; // Coolant Number Density

    data_F_ dF_F_;  // Fuel Destuction Rate d^F(F)
    data_F_ dC_F_;  // Coolant Destuction Rate d^C(F)
    data_F_ BU_F_;  // Burnup BU(F)
    data_F_ P_F_;   // Production Rate P(F)
    data_F_ D_F_;   // Destruction Rate P(F)
    data_F_ k_F_;   // k(F) -- almost meaningless
    std::map<int, data_F_ > Mj_F_;  // Transmuted to Matrix Mj(F)
    data_F_ zeta_F_;  // Disadvantage Factor zeta(F)

    int fd;			// Lower index of discharge fluence
    double Fd;	// Discharge Fluence
    double BUd; // Discharge Burnup
    double k;   // Multiplication factor 

    pyne::Material mat_feed_u;    // Input Uranium pyne::Material
    pyne::Material mat_feed_tru;  // Input Transuranic pyne::Material
    pyne::Material mat_feed_lan; 	// Input Lanthinide pyne::Material
    pyne::Material mat_feed_act;  // Input Actinide pyne::Material
    pyne::Material mat_prod_u;    // Output Uranium pyne::Material
    pyne::Material mat_prod_tru;  // Output Transuranic pyne::Material
    pyne::Material mat_prod_lan;  // Output Lanthinide pyne::Material
    pyne::Material mat_prod_act; 	// Output Actinide pyne::Material

    double deltaR;  // The production rate subtracted by the destruction rate at target_BU
    double tru_cr;  // Transuranic Conversion Ratio

    data_F_ SigmaFa_F_;   // Fuel Macro Absorption XS, Sigma^F_a(F)
    data_F_ SigmaFtr_F_;  // Fuel Macro Transport XS, Sigma^F_tr(F)
    data_F_ kappaF_F_;    // Fuel kappa, kappa^F(F)

    data_F_ SigmaCa_F_;   // Coolant Macro Absorption XS, Sigma^C_a(F)
    data_F_ SigmaCtr_F_;  // Coolant Macro Transport XS, Sigma^C_tr(F)
    data_F_ kappaC_F_;    // Coolant kappa, kappa^C(F)

    data_F_ lattice_E_F_;	// Values for lattice function E(F)
    data_F_ lattice_F_F_; // Values for lattice function F(F)


    // Public access functions
    void initialize(ReactorParameters);
    void loadlib(std::string libfile = "Reactor.h5");
    void fold_mass_weights();

    void calc_Mj_F_();
    void calc_Mj_Fd_();

    void   calc_mat_prod();
    void   calcSubStreams();
    double calc_tru_cr();

    double calc_deltaR();
    double calc_deltaR(pyne::comp_map);
    double calc_deltaR(pyne::Material);

    FluencePoint fluence_at_BU(double);
    double batch_average(double, std::string = "K");
    double batch_average_k(double);
    void BUd_bisection_method();
    void run_P_NL(double);
    void calibrate_P_NL_to_BUd();

    pyne::Material calc();
    pyne::Material calc(pyne::comp_map);
    pyne::Material calc(pyne::Material);	

    void lattice_E_planar(double, double);
    void lattice_F_planar(double, double);
    void lattice_E_spherical(double, double);
    void lattice_F_spherical(double, double);
    void lattice_E_cylindrical(double, double);
    void lattice_F_cylindrical(double, double);

    void  calc_zeta();
    void  calc_zeta_planar();
    void  calc_zeta_spherical();
    void  calc_zeta_cylindrical();
  };

// end bright
};

#endif
