// Reactor1G.h
// Header for general One Group Reactor Model

#if !defined(_Bright_FuelFabrication_)
#define _Bright_FuelFabrication_

//C++ stdlib
#include <map>
#include <set>
#include <vector>
#include <string>
#include <exception>
#include <iostream>

//Bright Libs
#include "Reactor1G.h"
#include "MassStream.h"
#include "isoname.h"

/***********************************************/
/*** Reactor1G Component Class and Functions ***/
/***********************************************/

typedef std::map<std::string, MassStream *> MassStreams;
typedef std::map<std::string, double> MassWeights;
typedef std::map<std::string, double> MassDeltaRs;

namespace bright 
{
    static std::set<std::string> empty_string_set;
    std::set<std::string> make_fuel_fab_params_set(MassStreams *, std::set<std::string> = empty_string_set);
};

class FuelFabrication : public FCComp
{
/** Fule Fabrication class
 *  Computes the value of different MassStreams inside of a reactor.  
 *  From here an optimum fuel for this reactor may be found.
 */
protected:

public:
    //Reactor1G Constructors
    FuelFabrication ();
    FuelFabrication (std::string);
    FuelFabrication (std::set<std::string>, std::string = "");
    FuelFabrication (MassStreams, MassWeights, Reactor1G, std::string = "");
    FuelFabrication (MassStreams, MassWeights, Reactor1G, std::set<std::string>, std::string = "");
    ~FuelFabrication ();
    
    //Public data
    MassStreams mass_streams;
    MassWeights mass_weights_in;
    MassWeights mass_weights_out;
    MassDeltaRs deltaRs;

    Reactor1G reactor;

    //Public access functions
    void initialize(MassStreams, MassWeights, Reactor1G);
    void setParams ();

    void calc_deltaRs();
    MassStream calc_core_input();
    void calc_mass_ratios();

    MassStream   doCalc ();
    MassStream   doCalc (MassStreams, MassWeights, Reactor1G);

};


class BadFuelWeights : public std::exception
{
//Exception for valid fuel weights.
public:
    double k_a;
    double k_b;
    
    BadFuelWeights () {};
    ~BadFuelWeights () throw () {};

    BadFuelWeights (double ka, double kb)
    {
        k_a = ka;
        k_b = kb;
    };    
    
    static char * name ()
    {
        return (char *) "BadFuelWeights";
    };
    
    virtual const char* what() const throw()
    {
        std::string BFWstr ("BadFuelWeights: Multiplication Factor Opperates on Range [");
        BFWstr += k_a; 
        BFWstr += ", "; 
        BFWstr += k_b;
        BFWstr += "]";
        return (const char * ) BFWstr.c_str();
    };
};


#endif
