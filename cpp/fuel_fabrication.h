// Reactor1G.h
// Header for general One Group Reactor Model

#if !defined(_BRIGHT_FUEL_FABRICATION_)
#define _BRIGHT_FUEL_FABRICATION_

#include "reactor1g.h"


namespace bright {

  /***********************************************/
  /*** Reactor1G Component Class and Functions ***/
  /***********************************************/

  typedef std::map<std::string, pyne::Material *> material_dict;
  typedef std::map<std::string, double> mass_weight_dict;
  typedef std::map<std::string, double> deltaR_dict;

  static std::set<std::string> empty_string_set;
  std::set<std::string> make_fuel_fab_params_set(material_dict *, std::set<std::string> = empty_string_set);

  class FuelFabrication : public FCComp
  {
  /** Fule Fabrication class
   *  Computes the value of different pyne::Materials inside of a reactor.  
   *  From here an optimum fuel for this reactor may be found.
   */
  public:
    // Reactor1G Constructors
    FuelFabrication ();
    FuelFabrication (std::string);
    FuelFabrication (std::set<std::string>, std::string = "");
    FuelFabrication (material_dict, mass_weight_dict, Reactor1G, std::string = "");
    FuelFabrication (material_dict, mass_weight_dict, Reactor1G, std::set<std::string>, std::string = "");
    ~FuelFabrication ();
    
    //Public data
    material_dict materials;
    mass_weight_dict mass_weights_in;
    mass_weight_dict mass_weights_out;
    deltaR_dict deltaRs;

    Reactor1G reactor;

    //Public access functions
    void initialize(material_dict, mass_weight_dict, Reactor1G);
    void calc_params ();

    void calc_deltaRs();
    pyne::Material calc_core_input();
    void calc_mass_ratios();

    pyne::Material calc();
    pyne::Material calc(material_dict, mass_weight_dict, Reactor1G);
  };


  class BadFuelWeights : public std::exception
  {
  // Exception for valid fuel weights.
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

// end bright
};

#endif
