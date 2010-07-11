// Enrichment.h
// Header for general Fuel Cycle Component Objects

#if !defined(_BRIGHT_ENRICHMENT_)
#define _BRIGHT_ENRICHMENT_
#include <map>
#include <set>
#include <string>
#include <exception>
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
    double alpha_0; //specify on init.
    MassStream IsosIn;
    MassStream IsosOut;
    MassStream IsosTail;

    //key isotopic info
    int j;          //The jth isotope is the key, in zzaaam form, must be in IsosIn.
    double xP_j;    //Product enrichment of jth isotope
    double xW_j;    //Waste/Tails enrichment of the jth isotope

    //Stage info
    double N;       //N-stages
    double M;       //M-stages
    double N0;      //initial guess of N-stages
    double M0;      //initial guess of M-stages


    //Reprocessing Constructors
    Enrichment ();
    Enrichment (std::string = "");
//    Reprocess (std::map<std::string, double>, std::string = "");
    
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

    double get_alphastari (double, double);

    double get_Ei (double, double);
    double get_Si (double, double);
    void FindNM(double):

    def xiP(comp, N, M, alpha0, Mstar, compF, j,  xjP, xjW):
    def xiW(comp, N, M, alpha0, Mstar, compF, j,  xjP, xjW):
    def SolveNM(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0):
    def Comp2UnitySecant(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0):
def Comp2UnityOther(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0):
def eq28denom(comp, alpha0, Mstar, j):
def LoverF(alpha0, Mstar, compF, j,  xjP, xjW, N0, M0, k):
def MstarOptimize(alpha0, Mstar0, compF, j,  xjP, xjW, N0, M0, k):
};

#endif


