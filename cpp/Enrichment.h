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
#include "pyne::Material.h"
#include "isoname.h"

/************************************************/
/*** Enrichment Component Class and Functions ***/
/************************************************/

//static std::string enr_p2t [] = {"MassFeed", "MassProduct", "MassTails", \
//    "N", "M", "Mstar", "TotalPerFeed", "SWUperFeed", "SWUperProduct"};
//static std::set<std::string> enr_p2track (enr_p2t, enr_p2t+9);
extern std::string enr_p2t [];
extern std::set<std::string> enr_p2track;


class EnrichmentParameters 
{
    /** Set of physical parameters used to specify an enrichment cascade.
     */

public:
    // Constructors
    EnrichmentParameters();
    ~EnrichmentParameters();

    // Attributes
    double alpha_0; //Initial stage separation factor
    double Mstar_0; //Initial guess for mass separation factor

    int j; //Component to enrich (U-235), zzaaam form
    int k; //Component to de-enrich, or strip (U-238), zzaaam form

    double N0; //Initial guess for the number of enriching stages
    double M0; //Initial guess for the number of stripping stages

    double xP_j; //Enrichment of the jth isotope in the product stream
    double xW_j; //Enrichment of the jth isotope in the waste (tails) stream
};

EnrichmentParameters fillUraniumEnrichmentDefaults();
extern EnrichmentParameters UraniumEnrichmentDefaults;

class Enrichment : public FCComp
{
//Reprocessing class
protected:
    //Protected functions

public:
    //Reprocessing Constructors
    Enrichment ();
    Enrichment (std::string);
    Enrichment (EnrichmentParameters, std::string = "");
    ~Enrichment ();

    //Public data
    double alpha_0;         //specify on init.
    double Mstar_0;         //specify on init.
    double Mstar;           //Current Mstar
    pyne::Material mat_tail;    //Waste Stream

    //key isotopic info
    int j;          //The jth isotope is the key, in zzaaam form, must be in mat_feed.
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


    //Public access functions
    void initialize(EnrichmentParameters);		//Initializes the constructors.
    void calc_params ();
    pyne::Material calc ();
    pyne::Material calc (pyne::comp_map);
    pyne::Material calc (pyne::Material);	

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


