// Enrichment.h
// Header for general Fuel Cycle Component Objects

#if !defined(_BRIGHT_ENRICHMENT_)
#define _BRIGHT_ENRICHMENT_
#include <map>
#include <set>
#include <string>
#include <vector>
#include <exception>
#include <math.h>
#include "FCComp.h"
#include "MassStream.h"
#include "isoname.h"

/************************************************/
/*** Enrichment Component Class and Functions ***/
/************************************************/

static std::string enr_p2t [] = {"MassFeed", "MassProduct", "MassTails", \
    "N", "M", "Mstar", "TotalPerFeed", "SWUperFeed", "SWUperProduct"};
static std::set<std::string> enr_p2track (enr_p2t, enr_p2t+9);


struct EnrichmentParameters 
{
    /** Set of physical parameters used to specify an enrichment cascade.
     */

    double alpha_0; //Initial stage separation factor
    double Mstar_0; //Initial guess for mass separation factor

    int j; //Component to enrich (U-235), zzaaam form
    int k; //Component to de-enrich, or strip (U-238), zzaaam form

    double N0; //Initial guess for the number of enriching stages
    double M0; //Initial guess for the number of stripping stages

    double xP_j; //Enrichment of the jth isotope in the product stream
    double xW_j; //Enrichment of the jth isotope in the waste (tails) stream

    //Get Functions
    double get_alpha_0() const {return alpha_0;};
    double get_Mstar_0() const {return Mstar_0;};
    int    get_j()       const {return j;};
    int    get_k()       const {return k;};
    double get_N0()      const {return N0;};
    double get_M0()      const {return M0;};
    double get_xP_j()    const {return xP_j;};
    double get_xW_j()    const {return xW_j;};

    //Set Functions
    void set_alpha_0(double a0) {alpha_0 = a0;};
    void set_Mstar_0(double m0) {Mstar_0 = m0;};
    void set_j(int jso)         {j = jso;};
    void set_k(int kso)         {k = kso;};
    void set_N0(double n0)      {N0 = n0;};
    void set_M0(double m0)      {M0 = m0;};
    void set_xP_j(double xpj)   {xP_j = xpj;};
    void set_xW_j(double xwj)   {xW_j = xwj;};
};

EnrichmentParameters fillUraniumEnrichmentDefaults();
extern EnrichmentParameters UraniumEnrichmentDefaults;

class Enrichment : public FCComp
{
//Reprocessing class
protected:
    //Protected functions

public:
    //Public data
    double alpha_0;         //specify on init.
    double Mstar_0;         //specify on init.
    double Mstar;           //Current Mstar
    MassStream IsosTail;    //Waste Stream

    //key isotopic info
    int j;          //The jth isotope is the key, in zzaaam form, must be in IsosIn.
    int k;          //The kth isotope is the other key to separate j away from.
    double xP_j;    //Product enrichment of jth isotope
    double xW_j;    //Waste/Tails enrichment of the jth isotope

    //Stage info
    double N;       //N Enriching Stages
    double M;       //M Stripping Stages
    double N0;      //initial guess of N-stages
    double M0;      //initial guess of M-stages

    //Flow Rates
    double TotalPerFeed;    //Total flow rate per feed rate.
    double SWUperFeed;      //This is the SWU for 1 kg of Feed material.
    double SWUperProduct;   //This is the SWU for 1 kg of Product material.

    //Reprocessing Constructors
    Enrichment ();
    Enrichment (std::string);
    Enrichment (EnrichmentParameters, std::string = "");
    
    //Get Functions
    double get_alpha_0() const {return alpha_0;};
    double get_Mstar_0() const {return Mstar_0;};
    int    get_j()       const {return j;};
    int    get_k()       const {return k;};
    double get_N0()      const {return N0;};
    double get_M0()      const {return M0;};
    double get_xP_j()    const {return xP_j;};
    double get_xW_j()    const {return xW_j;};

    double     get_Mstar()    const {return Mstar;};
    MassStream get_IsosTail() const {return IsosTail;};
    double     get_N()        const {return N;};
    double     get_M()        const {return M;};

    double get_TotalPerFeed()  const {return TotalPerFeed;};
    double get_SWUperFeed()    const {return SWUperFeed;};
    double get_SWUperProduct() const {return SWUperProduct;};

    //Set Functions
    void set_alpha_0(double a0) {alpha_0 = a0;};
    void set_Mstar_0(double m0) {Mstar_0 = m0;};
    void set_j(int jso)         {j = jso;};
    void set_k(int kso)         {k = kso;};
    void set_N0(double n0)      {N0 = n0;};
    void set_M0(double m0)      {M0 = m0;};
    void set_xP_j(double xpj)   {xP_j = xpj;};
    void set_xW_j(double xwj)   {xW_j = xwj;};

    void set_Mstar(double m)         {Mstar = m;};
    void set_IsosTail(MassStream it) {IsosTail = it;};
    void set_N(double n)             {N = n;};
    void set_M(double m)             {M = m;};

    void set_TotalPerFeed(double tpf)    {TotalPerFeed = tpf;};
    void set_SWUperFeed(double swupf)    {SWUperFeed = swupf;};
    void set_SWUperProduct(double swupp) {SWUperProduct = swupp;};

    //Public access functions
    void initialize(EnrichmentParameters);		//Initializes the constructors.
    void setParams ();
    MassStream doCalc ();
    MassStream doCalc (CompDict);
    MassStream doCalc (MassStream);	

    double PoverF (double, double, double);
    double WoverF (double, double, double);

    double get_alphastar_i (double);

    double get_Ei (double);
    double get_Si (double);
    void FindNM();

    double xP_i(int);
    double xW_i(int);
    void SolveNM();
    void Comp2UnitySecant();
    void Comp2UnityOther();
    double deltaU_i_OverG(int);
    void LoverF();
    void MstarOptimize();
};


//Exceptions
class EnrichmentInfiniteLoopError: public std::exception
{
    virtual const char* what() const throw()
    {
        return "Inifinite loop found while calculating enrichment cascade!  Breaking...";
    };
};

class EnrichmentIterationLimit: public std::exception
{
    virtual const char* what() const throw()
    {
        return "Iteration limit hit durring enrichment calculation!  Breaking...";
    };
};

class EnrichmentIterationNaN: public std::exception
{
    virtual const char* what() const throw()
    {
        return "Iteration has hit a point where some values are not-a-number!  Breaking...";
    };
};


#endif


