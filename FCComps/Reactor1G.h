// Reactor1G.h
// Header for general One Group Reactor Model

#if !defined(_Bright_Reactor1G_)
#define _Bright_Reactor1G_

//C++ stdlib
#include <map>
#include <set>
#include <vector>
#include <string>
#include <exception>
#include <iostream>

//HDF5
#include "hdf5.h"
#include "hdf5_hl.h"

//Boost
#include <boost/math/special_functions/bessel.hpp>

//Bright Libs
#include "FCComp.h"
#include "MassStream.h"
#include "isoname.h"

/***********************************************/
/*** Reactor1G Component Class and Functions ***/
/***********************************************/

typedef std::map<int, std::vector<float> > IsoFluenceDict;
typedef IsoFluenceDict::iterator IsoFluenceIter;

typedef std::vector<double> Data_F_;

typedef std::set<int> IsoSet;
typedef IsoSet::iterator IsoIter;

struct FluencePoint
{
    int f;
    double F;
    double m;

    FluencePoint()
    {
        f = 0;
        F = 0.0;
        m = 0.0;
    };
};

struct ReactorParameters
{
        /** Set of physical reactor parameters.
         *  May be used to instantiate new reactor objects, -or-
         *  to define default settings for a reactor type.
         */

    //General Info
        int batches;					//Total number of fuel loading batches
        double flux;					//Flux used for Fluence
        std::map<std::string, double> FuelForm;		//Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent$
        std::map<std::string, double> CoolantForm;	//Same a fuel chemical form but for coolant.  Should not have "IHM"
        double FuelDensity;				//Fuel Density
        double CoolantDensity;				//Coolant Density
        double pnl;					//Non-Leakage Probability
        double BUt;					//Target Discharge Burnup, only used for graphing inside of this component
        bool useDisadvantage;				//Boolean value on whether or not the disadvantage factor should be used
        std::string LatticeType;			//Lattice Type (Planar || Spherical || Cylindrical)
        bool HydrogenRescale;				//Rescale the Hydrogen-1 XS?

        //Volumetric Info
        double Radius;					//Fuel region radius
        double Length;					//Unit cell side length
        double open_slots;				//Number of open slots in fuel assembly
        double total_slots;				//Total number of Fuel assembly slots.


    //Get Functions
    int                           get_batches()         const {return batches;};
    double                        get_flux()            const {return flux;};
    std::map<std::string, double> get_FuelForm()        const {return FuelForm;};
    std::map<std::string, double> get_CoolantForm()     const {return CoolantForm;};
    double                        get_FuelDensity()     const {return FuelDensity;};
    double                        get_CoolantDensity()  const {return CoolantDensity;};
    double                        get_pnl()             const {return pnl;};
    double                        get_BUt()             const {return BUt;};
    bool                          get_useDisadvantage() const {return useDisadvantage;};
    std::string                   get_LatticeType()     const {return LatticeType;};
    bool                          get_HydrogenRescale() const {return HydrogenRescale;};

    double                        get_Radius()          const {return Radius;};
    double                        get_Length()          const {return Length;};
    double                        get_open_slots()      const {return open_slots;};
    double                        get_total_slots()     const {return total_slots;};

    //Set Functions
    void set_batches(int b)                                {batches = b;};
    void set_flux(double f)                                {flux = f;};
    void set_FuelForm(std::map<std::string, double> ff)    {FuelForm = ff;};
    void set_CoolantForm(std::map<std::string, double> cf) {CoolantForm = cf;};
    void set_FuelDensity(double fd)                        {FuelDensity = fd;};
    void set_CoolantDensity(double cd)                     {CoolantDensity = cd;};
    void set_pnl(double prob_nl)                           {pnl = prob_nl;};
    void set_BUt(double buT)                               {BUt = buT;};
    void set_useDisadvantage(bool uD)                      {useDisadvantage = uD;};
    void set_LatticeType(std::string lt)                   {LatticeType = lt;};
    void set_HydrogenRescale(bool hr)                      {HydrogenRescale = hr;};

    void set_Radius(double r)                              {Radius = r;};
    void set_Length(double l)                              {Length = l;};
    void set_open_slots(double os)                         {open_slots = os;};
    void set_total_slots(double ts)                        {total_slots = ts;};
};

class Reactor1G : public FCComp
{
/** Reactor class
 *  Basic One-Group Reactor Model.  Computes one Burnup Calculationn with the option of computing output isotopics.
 *  Specific reactor types inherit this class and change base parameters.
 */
protected:
    //Protected data
    IsoSet I;						//Set of isotopes that may be in IsosIn.
    IsoSet J;						//Set of isotopes that may be in IsosOut.

    //Thermal XS data is read in from static KAERI Data
    //Only read in if the disadvantage factor will be used.
    std::map<int, double> sigma_a_therm;			//Microscopic Thermal Absorption XS 
    std::map<int, double> sigma_s_therm;			//Microscopic Thermal Scattering XS 

public:
    //Public data
    int B; 								//Total number of fuel loading batches
    double phi;							//Flux used for Fluence
    std::map<std::string, double> FuelChemicalForm;			//Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
    std::map<std::string, double> CoolantChemicalForm;		//Same a fuel chemical form but for coolant.  Should not have "IHM"
    double rhoF;							//Fuel Density
    double rhoC; 							//Coolant Density
    double P_NL; 							//Non-Leakage Probability
    double TargetBU; 						//Target Discharge Burnup, only used for graphing inside of this component
    bool useZeta; 							//Boolean value on whether or not the disadvantage factor should be used
    std::string Lattice;						//Lattice Type (Planar || Spherical || Cylindrical)
    bool H_XS_Rescale;						//Rescale the Hydrogen-1 XS?

    double r; 							//Fuel region radius
    double l; 							//Unit cell side length
    double S_O; 							//Number of open slots in fuel assembly
    double S_T; 							//Total number of Fuel assembly slots.
    double VF; 							//Fuel Volume
    double VC;							//Coolant Volume

    std::string libfile;						//Path where the reactor's HDF5 library
    std::vector<float> F;   					//Fluence in [n/kb]
    IsoFluenceDict BUi_F_;						//Burnup [MWd/kgIHM]
    IsoFluenceDict pi_F_; 						//Production rate [n/s]
        IsoFluenceDict di_F_;		      	   			//Destruction rate [n/s]
    std::map<int, IsoFluenceDict> Tij_F_;				//Transformation Matrix [kg_i/kgIHM]

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

    MassStream InU; 						//Input Uranium MassStream
    MassStream InTRU; 						//Input Transuranic MassStream
    MassStream InLAN; 						//Input Lanthinide MassStream
    MassStream InACT; 						//Input Actinide MassStream
    MassStream OutU; 						//Output Uranium MassStream
    MassStream OutTRU; 						//Output Transuranic MassStream
    MassStream OutLAN; 						//Output Lanthinide MassStream
    MassStream OutACT; 						//Output Actinide MassStream

    double TruCR;							//Transuranic Conversion Ratio

    Data_F_ SigmaFa_F_;						//Fuel Macro Absorption XS, Sigma^F_a(F)
    Data_F_ SigmaFtr_F_;						//Fuel Macro Transport XS, Sigma^F_tr(F)
    Data_F_ kappaF_F_;						//Fuel kappa, kappa^F(F)

    Data_F_ SigmaCa_F_;						//Coolant Macro Absorption XS, Sigma^C_a(F)
    Data_F_ SigmaCtr_F_;						//Coolant Macro Transport XS, Sigma^C_tr(F)
    Data_F_ kappaC_F_;						//Coolant kappa, kappa^C(F)

    Data_F_ LatticeE_F_;						//Values for lattice function E(F)
    Data_F_ LatticeF_F_;						//Values for lattice function F(F)

    //Reactor1G Constructors
    Reactor1G ();
    Reactor1G (std::string);
    Reactor1G (std::set<std::string>, std::string = "");
    Reactor1G (ReactorParameters, std::string = "");
    Reactor1G (ReactorParameters, std::set<std::string>, std::string = "");
    ~Reactor1G ();
    
    //Get Functions
    int                           get_B()                   const {return B;};
    double                        get_phi()                 const {return phi;};
    std::map<std::string, double> get_FuelChemicalForm()    const {return FuelChemicalForm;};
    std::map<std::string, double> get_CoolantChemicalForm() const {return CoolantChemicalForm;};
    double                        get_rhoF()                const {return rhoF;};
    double                        get_rhoC()                const {return rhoC;};
    double                        get_P_NL()                const {return P_NL;};
    double                        get_TargetBU()            const {return TargetBU;};
    bool                          get_useZeta()             const {return useZeta;};
    std::string                   get_Lattice()             const {return Lattice;};
    bool                          get_H_XS_Rescale()        const {return H_XS_Rescale;};

    double                        get_r()                   const {return r;};
    double                        get_l()                   const {return l;};
    double                        get_S_O()                 const {return S_O;}; 
    double                        get_S_T()                 const {return S_T;};
    double                        get_VF()                  const {return VF;};
    double                        get_VC()                  const {return VC;};

    std::string                   get_libfile()             const {return libfile;};
    std::vector<float>            get_F()                   const {return F;};
    IsoFluenceDict                get_BUi_F_()              const {return BUi_F_;};
    IsoFluenceDict                get_pi_F_()               const {return pi_F_;};
    IsoFluenceDict                get_di_F_()               const {return di_F_;};
    std::map<int, IsoFluenceDict> get_Tij_F_()              const {return Tij_F_;};

    IsoSet get_I() const {return I;};
    IsoSet get_J() const {return J;};
    std::map<int, double> get_sigma_a_therm() const {return sigma_a_therm;};
    std::map<int, double> get_sigma_s_therm() const {return sigma_s_therm;};

    double                        get_A_IHM()               const {return A_IHM;};
    double                        get_MWF()                 const {return MWF;};
    double                        get_MWC()                 const {return MWC;};
    CompDict                      get_niF()                 const {return niF;};
    CompDict                      get_niC()                 const {return niC;};
    CompDict                      get_miF()                 const {return miF;};
    CompDict                      get_miC()                 const {return miC;};
    CompDict                      get_NiF()                 const {return NiF;};
    CompDict                      get_NiC()                 const {return NiC;};

    Data_F_                       get_dF_F_()               const {return dF_F_;};
    Data_F_                       get_dC_F_()               const {return dC_F_;};
    Data_F_                       get_BU_F_()               const {return BU_F_;};
    Data_F_                       get_P_F_()                const {return P_F_;};
    Data_F_                       get_D_F_()                const {return D_F_;};
    Data_F_                       get_k_F_()                const {return k_F_;};
    std::map<int, Data_F_ >       get_Mj_F_()               const {return Mj_F_;};
    Data_F_                       get_zeta_F_()             const {return zeta_F_;};

    int                           get_fd()                  const {return fd;};
    double                        get_Fd()                  const {return Fd;};
    double                        get_BUd()                 const {return BUd;};
    double                        get_k()                   const {return k;};

    MassStream                    get_InU()                 const {return InU;};
    MassStream                    get_InTRU()               const {return InTRU;};
    MassStream                    get_InLAN()               const {return InLAN;};
    MassStream                    get_InACT()               const {return InACT;};
    MassStream                    get_OutU()                const {return OutU;};
    MassStream                    get_OutTRU()              const {return OutTRU;};
    MassStream                    get_OutLAN()              const {return OutLAN;};
    MassStream                    get_OutACT()              const {return OutACT;};

    double                        get_TruCR()               const {return TruCR;};

    Data_F_                       get_SigmaFa_F_()          const {return SigmaFa_F_;};
    Data_F_                       get_SigmaFtr_F_()         const {return SigmaFtr_F_;};
    Data_F_                       get_kappaF_F_()           const {return kappaF_F_;};

    Data_F_                       get_SigmaCa_F_()          const {return SigmaCa_F_;};
    Data_F_                       get_SigmaCtr_F_()         const {return SigmaCtr_F_;};
    Data_F_                       get_kappaC_F_()           const {return kappaC_F_;};

    Data_F_                       get_LatticeE_F_()         const {return LatticeE_F_;};
    Data_F_                       get_LatticeF_F_()         const {return LatticeF_F_;};

    //Set Functions
    void set_B(int b)                                               {B = b;};
    void set_phi(double p)                                          {phi = p;};
    void set_FuelChemicalForm(std::map<std::string, double> fcf)    {FuelChemicalForm = fcf;};
    void set_CoolantChemicalForm(std::map<std::string, double> ccf) {CoolantChemicalForm = ccf;};
    void set_rhoF(double rF)                                        {rhoF = rF;};
    void set_rhoC(double rC) 					{rhoC = rC;};
    void set_P_NL(double p_nl)                                      {P_NL= p_nl;};
    void set_TargetBU(double tbu)                                   {TargetBU = tbu;};
    void set_useZeta(bool uZ)                                       {useZeta = uZ;}; 	
    void set_Lattice(std::string lt)                                {Lattice = lt;};
    void set_H_XS_Rescale(bool hxs)                                 {H_XS_Rescale = hxs;}; 	

    void set_r(double R)                                            {r = R;};
    void set_l(double L)                                            {l = L;};
    void set_S_O(double s_o)                                        {S_O = s_o;};
    void set_S_T(double s_t)                                        {S_T = s_t;};
    void set_VF(double vf)                                          {VF = vf;};
    void set_VC(double vc)                                          {VC = vc;};

    void set_libfile(std::string lf)                                {libfile = lf;};

    void set_A_IHM(double a_ihm)                                    {A_IHM = a_ihm;};
    void set_MWF(double mwf)                                        {MWF = mwf;};
    void set_MWC(double mwc)                                        {MWF = mwc;};
    void set_niF(CompDict NIf)                                      {niF = NIf;};
    void set_niC(CompDict NIc)                                      {niC = NIc;};
    void set_miF(CompDict MIf)                                      {miF = MIf;};
    void set_miC(CompDict MIc)                                      {miC = MIc;};
    void set_NiF(CompDict nIf)                                      {NiF = nIf;};
    void set_NiC(CompDict nIc)                                      {NiC = nIc;};

    void set_fd(int FD)                                             {fd = FD;};
    void set_Fd(double fD)                                          {Fd = fD;};
    void set_BUd(double buD)                                        {BUd = buD;};
    void set_k(double K)                                            {k = K;};

    void set_InU(MassStream ms)                                     {InU = ms;};
    void set_InTRU(MassStream ms)                                   {InTRU = ms;};
    void set_InLAN(MassStream ms)                                   {InLAN = ms;};
    void set_InACT(MassStream ms)                                   {InACT = ms;};
    void set_OutU(MassStream ms)                                    {OutU = ms;};
    void set_OutTRU(MassStream ms)                                  {OutTRU = ms;};
    void set_OutLAN(MassStream ms)                                  {OutLAN = ms;};
    void set_OutACT(MassStream ms)                                  {OutACT = ms;};

    void set_TruCR(double tcr)                                      {TruCR = tcr;};

    //Public access functions
    void         initialize(ReactorParameters);
    void         loadLib(std::string libfile = "Reactor.h5");
    void         foldMassWeights();

    void         mkMj_F_();
    void         mkMj_Fd_();

    void         calcOutIso();
    void         calcSubStreams();
    double       calcTruCR();

    FluencePoint FluenceAtBU(double);
    double       batchAve(double, std::string = "K");
    double       batchAveK(double);
    void         BUd_BisectionMethod();
    void         Run_PNL(double);
    void         Calibrate_PNL_2_BUd();

    MassStream   doCalc ();
    MassStream   doCalc (CompDict);
    MassStream   doCalc (MassStream);	

    void LatticeEPlanar(double, double);
    void LatticeFPlanar(double, double);
    void LatticeESpherical(double, double);
    void LatticeFSpherical(double, double);
    void LatticeECylindrical(double, double);
    void LatticeFCylindrical(double, double);

    void         calcZeta();
    void         calcZetaPlanar();
    void         calcZetaSpherical();
    void         calcZetaCylindrical();
};

/******************/
/*** Exceptions ***/
/******************/

class BadFuelForm : public std::exception
{
//Exception for valid fuel form.
public:
    BadFuelForm () {};
    ~BadFuelForm () throw () {};

    static char * name ()
    {
        return (char *) "BadFuelForm";
    };

    virtual const char* what() const throw()
    {
        if (1 < FCComps::verbosity)
        {
            std::cout << "\n";
            std::cout << "***********************************************\n";
            std::cout << "* WARNING: FUEL COMPOSITION NOT COMPUTABLE!!! *\n";
            std::cout << "***********************************************\n";
            std::cout << "\n";
        };
        std::string BFFstr ("FUEL COMPOSITION NOT COMPUTABLE!");
        return (const char *) BFFstr.c_str();
    };
};

class BisectionMethodNotPerformed : public std::exception
{
//Exception for when the bisection method is not calculated.
public:
    BisectionMethodNotPerformed ()
    {
        errstr = "Bisection method was not performed.";
    };
    BisectionMethodNotPerformed (std::string calctype)
    {
        errstr = "Bisection method durring " + calctype + " calculation was not performed.";
    };
    ~BisectionMethodNotPerformed () throw () {};

    static char * name ()
    {
        return (char *) "BisectionMethodNotPerformed";
    };

    virtual const char* what() const throw()
    {
        if (1 < FCComps::verbosity)
        {
            std::cout << "\n";
            std::cout << "**************\n";
            std::cout << "* WARNING!!! *\n";
            std::cout << "**************\n";
            std::cout << "\n";
            std::cout << errstr << "\n";
            std::cout << "\n";
        };
        return (const char *) errstr.c_str();
    };
private:
    std::string errstr;
};

#endif
