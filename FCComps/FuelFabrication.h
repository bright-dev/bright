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
    std::set<std::string> empty_string_set;
    std::set<std::string> make_fuel_fab_params_set(MassStreams *, std::set<std::string> = empty_string_set);
};

class FuelFabrication : public FCComp
{
/** Fule Fabrication class
 *  Computes the value of different MMassStreams inside of a reactor.  
 *  From here an optimum fuel for this reactor may be found.
 */
protected:

public:
    //Public data
    MassStreams mass_streams;
    MassWeights mass_weights;
    MassDeltaRs mass_deltaRs;

    Reactor1G * reactor;

    //Reactor1G Constructors
    FuelFabrication ();
    FuelFabrication (std::string);
    FuelFabrication (std::set<std::string>, std::string = "");
    FuelFabrication (MassStreams, MassWeights, Reactor1G *, std::string = "");
    FuelFabrication (MassStreams, MassWeights, Reactor1G *, std::set<std::string>, std::string = "");
    ~FuelFabrication ();
    
    //Get Functions
    MassStreams get_mass_streams() const {return mass_streams;};
    MassWeights get_mass_weights() const {return mass_weights;};
    MassDeltaRs get_mass_deltaRs() const {return mass_deltaRs;};

    Reactor1G * get_reactor() const {return reactor;};

    //Set Functions
    void set_mass_streams(MassStreams mss) {mass_streams = mss;};
    void set_mass_weights(MassWeights mws) {mass_weights = mws;};
    void set_mass_deltaRs(MassDeltaRs drs) {mass_deltaRs = drs;};

    void set_reactor(Reactor1G * r) {reactor = r;};


    //Public access functions
    void initialize(MassStreams, MassWeights, Reactor1G *);
};

#endif
