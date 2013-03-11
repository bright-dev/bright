// Fuel Fabrication Class

#include "fuel_fabrication.h"

/*** Common Function ***/

std::set<std::string> bright::make_fuel_fab_params_set(material_dict * mats, std::set<std::string> orig_set)
{
  std::set<std::string> new_set (orig_set);

  for ( material_dict::iterator mat=(*mats).begin() ; mat != (*mats).end(); mat++ )
  {
    new_set.insert( "Weight_" + (*mat).first );
    new_set.insert( "deltaR_" + (*mat).first );
  };

  return new_set;
};


/*****************************************************/
/*** FuelFabrication Component Class and Functions ***/
/*****************************************************/

bright::FuelFabrication::FuelFabrication(std::string n) : bright::FCComp(n)
{
};


bright::FuelFabrication::FuelFabrication(std::set<std::string> paramtrack, std::string n) : bright::FCComp(paramtrack, n)
{
};

bright::FuelFabrication::FuelFabrication(material_dict mats, mass_weight_dict mws_in, Reactor1G r, std::string n) : bright::FCComp(bright::make_fuel_fab_params_set(&mats), n)
{
  initialize(mats, mws_in, r);
};


bright::FuelFabrication::FuelFabrication(material_dict mats, mass_weight_dict mws_in, Reactor1G r, std::set<std::string> paramtrack, std::string n) : bright::FCComp(bright::make_fuel_fab_params_set(&mats, paramtrack), n)
{
  initialize(mats, mws_in, r);
};


bright::FuelFabrication::~FuelFabrication()
{
};



void bright::FuelFabrication::initialize(material_dict mats, mass_weight_dict mws_in, Reactor1G r)
{
  /** Sets the fuel fabrication specific parameters.
   *  Must be done once at the beginning of fuel fabrication object life.
   */

  materials = mats;
  mass_weights_in = mws_in;

  reactor = r;

  calc_deltaRs();
};



void bright::FuelFabrication::calc_deltaRs()
{
  // Calculates the deltaR for each mass stream
  for (material_dict::iterator mat = materials.begin(); mat != materials.end(); mat++)
  {
    pyne::Material m = (*(*mat).second);
    m.normalize();

    double dR = reactor.calc_deltaR(m);

    deltaRs[(*mat).first] = dR;
    m.mass = mass_weights_in[(*mat).first];
  };
};


pyne::Material bright::FuelFabrication::calc_core_input()
{
  pyne::comp_map cm;
  pyne::Material core_input (cm, 0.0);

  for (mass_weight_dict::iterator mws = mass_weights_out.begin(); mws != mass_weights_out.end(); mws++)
    core_input = core_input + ( (*materials[(*mws).first]) * (*mws).second );

  core_input.normalize();
  return core_input;
};


void bright::FuelFabrication::calc_mass_ratios()
{
  // Find the key mass stremas based on which two have weights that are less than zero
  std::string key_a;
  std::string key_b;

  bool found_a = false;

  for (mass_weight_dict::iterator mws = mass_weights_in.begin(); mws != mass_weights_in.end(); mws++)
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
  pyne::Material mat_a = *materials[key_a];
  mat_a.normalize();
  double dR_a = reactor.calc_deltaR( mat_a );

  // deltaR for key b
  pyne::Material mat_b = *materials[key_b];
  mat_b.normalize();
  double dR_b = reactor.calc_deltaR( mat_b );

  //First Guess for key mass stream masses; each get half of the remaining mass space.
  double top_up_mass_space = 1.0; 
  for (mass_weight_dict::iterator mws = mass_weights_in.begin(); mws != mass_weights_in.end(); mws++)
  {
    if (0.0 <= (*mws).second)
      top_up_mass_space = top_up_mass_space - (*mws).second;
  };


  // Initialize mass_weights_out as a copy of mass_weights_in
  for (mass_weight_dict::iterator mws = mass_weights_in.begin(); mws != mass_weights_in.end(); mws++)
    mass_weights_out[(*mws).first] = (*mws).second;

  double dR_guess;
  pyne::Material core_input;

  double k_a, k_b;
  double sign_a, sign_b;


  // Find bound for All Mass Stream A
  mass_weights_out[key_a] = top_up_mass_space;
  mass_weights_out[key_b] = 0.0;
  core_input = calc_core_input();
  dR_guess = reactor.calc_deltaR( core_input );
  k_a = reactor.batch_average_k( reactor.target_BU );
  sign_a = (1.0 - k_a) / fabs(1.0 - k_a);

  // Find bound for All Mass Stream B
  mass_weights_out[key_a] = 0.0;
  mass_weights_out[key_b] = top_up_mass_space;
  core_input = calc_core_input();
  dR_guess = reactor.calc_deltaR( core_input );
  k_b = reactor.batch_average_k( reactor.target_BU );
  sign_b = (1.0 - k_b) / fabs(1.0 - k_b);

  // Ensure calculation is possible
  if (sign_a == sign_b)
    throw BadFuelWeights(k_a, k_b);

  // Continue nomrally
  mass_weights_out[key_a] = top_up_mass_space * 0.5;
  mass_weights_out[key_b] = top_up_mass_space - mass_weights_out[key_a];

  // Calculate delta R for the Guess
  core_input = calc_core_input();
  dR_guess = reactor.calc_deltaR( core_input );

  int n;
  double k;
  double dMass;
  
  k = reactor.batch_average_k( reactor.target_BU );
  n = 0;
  if (0 < bright::verbosity)
    std::cout << n << ") " << k << " "; 

  while (0.001 < fabs(1.0 - k) && n < 10)
  {
    // Adjust Masses based on pertubation guess.
    dMass = - dR_guess / (dR_a - dR_b);
    mass_weights_out[key_a] = mass_weights_out[key_a] + dMass;
    mass_weights_out[key_b] = mass_weights_out[key_b] - dMass;

    // Recalculate core parameters for new masses guess
    core_input = calc_core_input();
    dR_guess = reactor.calc_deltaR( core_input );
    k = reactor.batch_average_k( reactor.target_BU );
    n = n+1;
    if (0 < bright::verbosity)
      std::cout << k << " ";
  };

  if (0 < bright::verbosity)
    std::cout << "\n\n";
};


pyne::Material bright::FuelFabrication::calc()
{
  calc_mass_ratios();
  mat_prod = calc_core_input();
  return mat_prod;
};


pyne::Material bright::FuelFabrication::calc(material_dict mats, mass_weight_dict mws_in, Reactor1G r)
{
  initialize(mats, mws_in, r);
  return calc();
};


void bright::FuelFabrication::calc_params ()
{
  for (mass_weight_dict::iterator mws = mass_weights_in.begin(); mws != mass_weights_in.end(); mws++)
  {
    std::string mw_key ("Weight_" + (*mws).first);
    std::string dR_key ("deltaR_" + (*mws).first);

    params_prior_calc[mw_key]  = mass_weights_in[(*mws).first];
    params_after_calc[mw_key] = mass_weights_out[(*mws).first];

    params_prior_calc[dR_key]  = deltaRs[(*mws).first];
    params_after_calc[dR_key] = deltaRs[(*mws).first];
  };
};

