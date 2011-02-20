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

FuelFabrication::FuelFabrication(MassStreams mss, MassWeights mws_in, Reactor1G r, std::string n) : \
    FCComp(bright::make_fuel_fab_params_set(&mss), n)
{
    initialize(mss, mws_in, r);
};


FuelFabrication::FuelFabrication(MassStreams mss, MassWeights mws_in, Reactor1G r, std::set<std::string> paramtrack, std::string n) : \
    FCComp(bright::make_fuel_fab_params_set(&mss, paramtrack), n)
{
    initialize(mss, mws_in, r);
};


FuelFabrication::~FuelFabrication()
{
};



void FuelFabrication::initialize(MassStreams mss, MassWeights mws_in, Reactor1G r)
{
    /** Sets the fuel fabrication specific parameters.
     *  Must be done once at the beginning of fuel fabrication object life.
     */

    mass_streams = mss;
    mass_weights_in = mws_in;

    reactor = r;

    calc_deltaRs();
};



void FuelFabrication::calc_deltaRs()
{
    // Calculates the deltaR for each mass stream
    for (MassStreams::iterator mss = mass_streams.begin(); mss != mass_streams.end(); mss++)
    {
        MassStream ms = (*(*mss).second);
        ms.normalize();

        double dR = reactor.calc_deltaR(ms);

        deltaRs[(*mss).first] = dR;
        ms.mass = mass_weights_in[(*mss).first];
    };
};


MassStream FuelFabrication::calc_core_input()
{
    CompDict cd;
    MassStream core_input (cd, 0.0);

    for (MassWeights::iterator mws = mass_weights_out.begin(); mws != mass_weights_out.end(); mws++)
    {
        core_input = core_input + ( (*mass_streams[(*mws).first]) * (*mws).second );
    };    

    core_input.normalize();
    core_input.name = "CoreInput";
    return core_input;
};


void FuelFabrication::calc_mass_ratios()
{
    // Find the key mass stremas based on which two have weights that are less than zero
    std::string key_a;
    std::string key_b;

    bool found_a = false;

    for (MassWeights::iterator mws = mass_weights_in.begin(); mws != mass_weights_in.end(); mws++)
    {
        if ((*mws).second < 0.0)
        {
            if (!found_a)
            {
                key_a = (*mws).first;
                found_a = true;
            }
            else
            {
                key_b = (*mws).first;
                break;
            };
        };
    };

    // deltaR for key a
    MassStream ms_a = *mass_streams[key_a];
    ms_a.normalize();
    double dR_a = reactor.calc_deltaR( ms_a );

    // deltaR for key b
    MassStream ms_b = *mass_streams[key_b];
    ms_b.normalize();
    double dR_b = reactor.calc_deltaR( ms_b );

    //First Guess for key mass stream masses; each get half of the remaining mass space.
    double top_up_mass_space = 1.0; 
    for (MassWeights::iterator mws = mass_weights_in.begin(); mws != mass_weights_in.end(); mws++)
    {
        if (0.0 <= (*mws).second)
            top_up_mass_space = top_up_mass_space - (*mws).second;
    };


    // Initialize mass_weights_out as a copy of mass_weights_in
    for (MassWeights::iterator mws = mass_weights_in.begin(); mws != mass_weights_in.end(); mws++)
    {
        mass_weights_out[(*mws).first] = (*mws).second;
    };

    double dR_guess;
    MassStream core_input;

    double k_a, k_b;
    double sign_a, sign_b;


    // Find bound for All Mass Stream A
    mass_weights_out[key_a] = top_up_mass_space;
    mass_weights_out[key_b] = 0.0;
    core_input = calc_core_input();
    dR_guess = reactor.calc_deltaR( core_input );
    k_a = reactor.batchAveK( reactor.TargetBU );
    sign_a = (1.0 - k_a) / fabs(1.0 - k_a);

    // Find bound for All Mass Stream B
    mass_weights_out[key_a] = 0.0;
    mass_weights_out[key_b] = top_up_mass_space;
    core_input = calc_core_input();
    dR_guess = reactor.calc_deltaR( core_input );
    k_b = reactor.batchAveK( reactor.TargetBU );
    sign_b = (1.0 - k_b) / fabs(1.0 - k_b);

    // Ensure calculation is possible
    if (sign_a == sign_b)
        throw BadFuelWeights(k_a, k_b);

    //Continue nomrally
    mass_weights_out[key_a] = top_up_mass_space * 0.5;
    mass_weights_out[key_b] = top_up_mass_space - mass_weights_out[key_a];

    //Calculate delta R for the Guess
    core_input = calc_core_input();
    dR_guess = reactor.calc_deltaR( core_input );

    int n;
    double k;
    double dMass;
    
    k = reactor.batchAveK( reactor.TargetBU );
    n = 0;
    if (0 < FCComps::verbosity)
        std::cout << n << ") " << k << " "; 

    while (0.001 < fabs(1.0 - k) && n < 10)
    {
        //Adjust Masses based on pertubation guess.
        dMass = - dR_guess / (dR_a - dR_b);
        mass_weights_out[key_a] = mass_weights_out[key_a] + dMass;
        mass_weights_out[key_b] = mass_weights_out[key_b] - dMass;

        //Recalculate core parameters for new masses guess
        core_input = calc_core_input();
        dR_guess = reactor.calc_deltaR( core_input );
        k = reactor.batchAveK( reactor.TargetBU );
        n = n+1;
        if (0 < FCComps::verbosity)
            std::cout << k << " ";
    };

    if (0 < FCComps::verbosity)
        std::cout << "\n\n";
    
};


MassStream FuelFabrication::doCalc()
{
    calc_mass_ratios();
    ms_prod = calc_core_input();
    return ms_prod;
};


MassStream FuelFabrication::doCalc(MassStreams mss, MassWeights mws_in, Reactor1G r)
{
    initialize(mss, mws_in, r);
    return doCalc();
};


void FuelFabrication::setParams ()
{
    for ( MassWeights::iterator mws = mass_weights_in.begin(); mws != mass_weights_in.end(); mws++ )
    {
        std::string mw_key ("Weight_" + (*mws).first);
        std::string dR_key ("deltaR_" + (*mws).first);

        params_prior_calc[mw_key]  = mass_weights_in[(*mws).first];
        params_after_calc[mw_key] = mass_weights_out[(*mws).first];

        params_prior_calc[dR_key]  = deltaRs[(*mws).first];
        params_after_calc[dR_key] = deltaRs[(*mws).first];
    };
};

