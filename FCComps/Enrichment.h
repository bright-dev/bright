// Enrichment.h
// Header for general Fuel Cycle Component Objects

#if !defined(_BRIGHT_ENRICHMENT_)
#define _BRIGHT_ENRICHMENT_
#include <map>
#include <set>
#include <string>
#include <list>
#include <exception>
#include <math.h>
#include "FCComp.h"
#include "MassStream.h"
#include "isoname.h"

/************************************************/
/*** Enrichment Component Class and Functions ***/
/************************************************/

static std::string enr_p2t [] = {"Mass"};
static std::set<std::string> enr_p2track (enr_p2t, enr_p2t+1);


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
    MassStream IsosIn;      //Feed stream
    MassStream IsosOut;     //Product Stream
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
    Enrichment (std::string = "");
    
    //Get Functions
//    SepEffDict get_sepeff() const {return sepeff;};

    //Set Functions
//    void set_sepeff(SepEffDict sed) {sepeff = sed;};

    //Public access functions
    void initialize();		//Initializes the constructors.
    void setParams ();
    MassStream doCalc ();
    MassStream doCalc (CompDict);
    MassStream doCalc (MassStream);	

    double PoverF (double, double, double);
    double WoverF (double, double, double);

    double get_alphastari (double);

    double get_Ei (double);
    double get_Si (double);
    void FindNM():

    double xP_i();
    double xW_i();
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


#endif


