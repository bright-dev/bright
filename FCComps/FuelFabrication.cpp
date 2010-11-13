// Fuel Fabrication Class

#include "FuelFabrication.h"

/*** Common Function ***/

std::set<std::string> bright::make_fuel_fab_params_set(MassStreams * mss, std::set<std::string> orig_set)
{
    std::set<std::string> new_set (orig_set);

    for ( MassStreams::iterator ms=(*mss).begin() ; ms != (*mss).end(); ms++ )
    {
        new_set.insert( "Weight_" + (*ms).first );
        new_set.insert( "deltaR_" + (*ms).first );
    };

    return new_set;
};

/*****************************************************/
/*** FuelFabrication Component Class and Functions ***/
/*****************************************************/

FuelFabrication::FuelFabrication()
{
};

FuelFabrication::FuelFabrication(std::string n) : FCComp(n)
{
};

FuelFabrication::FuelFabrication(std::set<std::string> paramtrack, std::string n) : FCComp(paramtrack, n)
{
};

FuelFabrication::FuelFabrication(MassStreams mss, MassWeights mws, Reactor1G r, std::string n) : \
    FCComp(bright::make_fuel_fab_params_set(&mss), n)
{
    initialize(mss, mws, &r);
};

FuelFabrication::FuelFabrication(MassStreams mss, MassWeights mws, Reactor1G * r, std::string n) : \
    FCComp(bright::make_fuel_fab_params_set(&mss), n)
{
    initialize(mss, mws, r);
};

FuelFabrication::FuelFabrication(MassStreams mss, MassWeights mws, Reactor1G r, std::set<std::string> paramtrack, std::string n) : \
    FCComp(bright::make_fuel_fab_params_set(&mss, paramtrack), n)
{
    initialize(mss, mws, &r);
};

FuelFabrication::FuelFabrication(MassStreams mss, MassWeights mws, Reactor1G * r, std::set<std::string> paramtrack, std::string n) : \
    FCComp(bright::make_fuel_fab_params_set(&mss, paramtrack), n)
{
    initialize(mss, mws, r);
};

FuelFabrication::~FuelFabrication()
{
};

void FuelFabrication::initialize(MassStreams mss, MassWeights mws, Reactor1G * r)
{
    /** Sets the fuel fabrication specific parameters.
     *  Must be done once at the beginning of fuel fabrication object life.
     */

    mass_streams = mss;
    mass_weights = mws;

    reactor = r;

};

