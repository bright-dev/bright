// Storage.h
// Header for general Fuel Cycle Component Objects

#if !defined(_BRIGHT_STORAGE_)
#define _BRIGHT_STORAGE_

#include "fccomp.h"

namespace bright {

  /*********************************************/
  /*** Storage Component Class and Functions ***/
  /*********************************************/
  typedef struct decay_nuc
  {
    int    fromiso;
    double halflife;
    double decayconst;
    int    toiso;
    double branchratio;
  } decay_nuc;

  // Sexy decay library form
  typedef std::map<int, double *> to_nuc_dict;
  typedef to_nuc_dict::iterator to_nuc_iter;
  typedef struct from_nuc_struct
  {
    double * halflife;
    double * decayconst;
    to_nuc_dict toiso;
  } from_nuc_struct;
  typedef std::map<int, from_nuc_struct> decay_dict;
  typedef decay_dict::iterator decay_iter;

  typedef std::vector<int> nuc_chain;
  typedef nuc_chain::iterator nuc_chain_iter;
  typedef std::set<nuc_chain> nuc_chain_set;
  typedef nuc_chain_set::iterator nuc_chain_set_iter;

  static std::string stor_p2t [] = {"Mass"};
  static std::set<std::string> stor_p2track (stor_p2t, stor_p2t+1);


  class Storage : public FCComp
  {
  // Storage/Cooling/Decay Fuel Cycle Component.
  protected:
    // Protected Data
    nuc_chain_set nucchains;
    decay_nuc * decay_data;
    int decay_data_len;
    decay_dict decay;

    // Protected functions
    void initialize ();                       // Initializes the constructors.
    double get_decay ();
    double bateman(int nuc, double mass, nuc_chain nucchain);  // Solves the Bateman Decay equation.
    void addchains(nuc_chain nc);               // Get the decay chain for a mother isotope
    void addchains(int nuc);                     // Get the decay chain for a mother isotope

    void print_chain (nuc_chain nc);

  public:
    // Storage Constructors	
    Storage(std::string n="");
    ~Storage();

    //Public data
    double decay_time;			//time to decay for

    //Public Functions
    void calc_params();
    pyne::Material calc();
    pyne::Material calc(pyne::comp_map incomp);
    pyne::Material calc(pyne::Material mat);
    pyne::Material calc(double t);
    pyne::Material calc(pyne::comp_map, double t);
    pyne::Material calc(pyne::Material mat, double t);
  };

// end bright
}; 

#endif
