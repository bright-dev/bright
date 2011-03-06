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

typedef std::map<int, std::vector<double> > IsoFluenceDict;
typedef IsoFluenceDict::iterator IsoFluenceIter;

typedef std::vector<double> Data_F_;
typedef std::vector<double> pert_data;
typedef std::vector< std::vector<double> > pert_data_g; 
typedef std::vector< std::vector< std::vector<double> > > pert_data_gh; 

typedef std::set<int> IsoSet;
typedef IsoSet::iterator IsoIter;


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
    //ReactorMG Constructors
    ReactorMG ();
    ReactorMG (std::string);
    ReactorMG (std::set<std::string>, std::string = "");
    ReactorMG (ReactorParameters, std::string = "");
    ReactorMG (ReactorParameters, std::set<std::string>, std::string = "");
    ~ReactorMG ();
    
    //Public data
    int B; 								//Total number of fuel loading batches
    double flux;							//Flux used for Fluence
    std::map<std::string, double> fuel_chemical_form;			//Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
    std::map<std::string, double> coolant_chemical_form;		//Same a fuel chemical form but for coolant.  Should not have "IHM"
    double rho_fuel;							//Fuel Density
    double rho_clad; 							//Coolant Density
    double rho_cool; 							//Coolant Density
    double P_NL; 							//Non-Leakage Probability
    double target_BU; 						//Target Discharge Burnup, only used for graphing inside of this component
    bool use_zeta; 							//Boolean value on whether or not the disadvantage factor should be used
    std::string lattice_flag;						//lattice_flag Type (Planar || Spherical || Cylindrical)
    bool rescale_hydrogen_xs;						//Rescale the Hydrogen-1 XS?

    double r_fuel; 							//Fuel region radius
    double r_void; 							//Void region radius
    double r_clad; 							//Cladding region radius
    double pitch; 							//Unit cell side length
    double S_O; 							//Number of open slots in fuel assembly
    double S_T; 							//Total number of Fuel assembly slots.
    double VF; 							//Fuel Volume
    double VC;							//Coolant Volume

    std::string libfile;						//Path where the reactor's HDF5 library
    std::vector<double> F;   					//Fluence in [n/kb]
    IsoFluenceDict BUi_F_;						//Burnup [MWd/kgIHM]
    IsoFluenceDict pi_F_; 						//Production rate [n/s]
    IsoFluenceDict di_F_;		      	   			//Destruction rate [n/s]
    std::map<int, IsoFluenceDict> Tij_F_;				//Transformation Matrix [kg_i/kgIHM]

    IsoSet I;						//Set of isotopes that may be in ms_feed.
    IsoSet J;						//Set of isotopes that may be in ms_prod.

    h5wrap::HomogenousTypeTable<double> perturbations;  // Load perturbation table
    int nperturbations;

    int G;                      // number of energu bins
    std::vector<double> E_g;    // Energy bin boundaries
    pert_data_g phi_g;          // Group fluxes
    pert_data phi;              // Total fluxes
    pert_data Phi;              // Fluence
    pert_data time0;            // initial time steps
    pert_data BU0;              // initial burnups

    std::map<int, pert_data> Ti0;
    std::map<int, pert_data_g> sigma_a;
    std::map<int, pert_data_g> sigma_s;
    std::map<int, pert_data_g> sigma_f;
    std::map<int, pert_data_g> nubar_sigma_f;
    std::map<int, pert_data_g> nubar;
    std::map<int, pert_data_gh> sigma_s_gh;

    double A_IHM;							//Atomic weight of IHM
    double MWF;							//Fuel Molecular Weight
    double MWC;							//Coolant Molecular Weight
    CompDict niF;							//Fuel Atom Number Weight
    CompDict niC;							//Coolant Atom Number Weight
    CompDict miF;							//Fuel Mass Weight
    CompDict miC;							//Coolant Mass Weight
    CompDict NiF;							//Fuel Number Density
    CompDict NiC;							//Coolant Number Density

    Data_F_ dF_F_;							//Fuel Destuction Rate d^F(F)
    Data_F_ dC_F_;							//Coolant Destuction Rate d^C(F)
    Data_F_ BU_F_;							//Burnup BU(F)
    Data_F_ P_F_;							//Production Rate P(F)
    Data_F_ D_F_;							//Destruction Rate P(F)
    Data_F_ k_F_;							//k(F) -- almost meaningless
    std::map<int, Data_F_ > Mj_F_;					//Transmuted to Matrix Mj(F)
    Data_F_ zeta_F_;						//Disadvantage Factor zeta(F)

    int fd;								//Lower index of discharge fluence
    double Fd;							//Discharge Fluence
    double BUd;							//Discharge Burnup
    double k;							//Multiplication factor 

    MassStream ms_feed_u; 						//Input Uranium MassStream
    MassStream ms_feed_tru; 						//Input Transuranic MassStream
    MassStream ms_feed_lan; 						//Input Lanthinide MassStream
    MassStream ms_feed_act; 						//Input Actinide MassStream
    MassStream ms_prod_u; 						//Output Uranium MassStream
    MassStream ms_prod_tru; 						//Output Transuranic MassStream
    MassStream ms_prod_lan; 						//Output Lanthinide MassStream
    MassStream ms_prod_act; 						//Output Actinide MassStream

    double deltaR;                          // The production rate subtracted by the destruction rate at target_BU
    double tru_cr;							//Transuranic Conversion Ratio

    Data_F_ SigmaFa_F_;						//Fuel Macro Absorption XS, Sigma^F_a(F)
    Data_F_ SigmaFtr_F_;						//Fuel Macro Transport XS, Sigma^F_tr(F)
    Data_F_ kappaF_F_;						//Fuel kappa, kappa^F(F)

    Data_F_ SigmaCa_F_;						//Coolant Macro Absorption XS, Sigma^C_a(F)
    Data_F_ SigmaCtr_F_;						//Coolant Macro Transport XS, Sigma^C_tr(F)
    Data_F_ kappaC_F_;						//Coolant kappa, kappa^C(F)

    Data_F_ lattice_E_F_;						//Values for lattice function E(F)
    Data_F_ lattice_F_F_;						//Values for lattice function F(F)


    //Public access functions
    void              initialize(ReactorParameters);
    void              loadlib(std::string libfile = "Reactor.h5");
    void              fold_mass_weights();
    std::vector<int>  nearest_neighbors();

    void         calc_Mj_F_();
    void         calc_Mj_Fd_();

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
