// One-Group Fast Reactor Component Class

#include "fast_reactor1g.h"

bright::FastReactor1G::FastReactor1G() : bright::Reactor1G(fr_defaults, fr_p2track)
{
};


bright::FastReactor1G::FastReactor1G(std::string lib, std::string n) : bright::Reactor1G(fr_defaults, fr_p2track, n)
{
  libfile = lib;
  loadlib(lib);
};

bright::FastReactor1G::FastReactor1G(ReactorParameters rp, std::string n) : bright::Reactor1G(rp, fr_p2track, n)
{
};

bright::FastReactor1G::FastReactor1G(std::string lib, ReactorParameters rp, std::string n) : bright::Reactor1G(rp, fr_p2track, n)
{
  libfile = lib;
  loadlib(lib);
};

bright::FastReactor1G::~FastReactor1G() 
{
};


void bright::FastReactor1G::calc_params()
{
  /** Sets relevent FR parameters.
   *  Overwrites standard, do-nothing calc_params() function.
   */

  calc_sub_mats();

  params_prior_calc["BUd"]  = 0.0;
  params_after_calc["BUd"] = BUd;

  params_prior_calc["TRUCR"]  = 0.0;
  params_after_calc["TRUCR"] = calc_tru_cr();

  params_prior_calc["P_NL"]  = 0.0;
  params_after_calc["P_NL"] = P_NL;

  params_prior_calc["U"]  = mat_feed_u.mass;
  params_after_calc["U"] = mat_prod_u.mass;

  params_prior_calc["TRU"]  = mat_feed_tru.mass;
  params_after_calc["TRU"] = mat_prod_tru.mass;

  params_prior_calc["ACT"]  = mat_feed_act.mass;
  params_after_calc["ACT"] = mat_prod_act.mass;

  params_prior_calc["LAN"]  = mat_feed_lan.mass;
  params_after_calc["LAN"] = mat_prod_lan.mass;

  params_prior_calc["FP"]  = 1.0 - mat_feed_act.mass  - mat_feed_lan.mass;
  params_after_calc["FP"] = 1.0 - mat_prod_act.mass - mat_prod_lan.mass;

  return;
};
