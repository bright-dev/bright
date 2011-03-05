// One-Group Light Water Reactor Component Class

#include "LightWaterReactor1G.h"
LightWaterReactor1G::LightWaterReactor1G() : Reactor1G(lwr_defaults, lwr_p2track)
{
};

LightWaterReactor1G::LightWaterReactor1G(std::string h5lib, std::string n) : Reactor1G(lwr_defaults, lwr_p2track, n)
{
    libfile = h5lib;
    loadlib(h5lib);
};

LightWaterReactor1G::LightWaterReactor1G(ReactorParameters rp, std::string n) : Reactor1G(rp, lwr_p2track, n)
{
};

LightWaterReactor1G::LightWaterReactor1G(std::string h5lib, ReactorParameters rp, std::string n) : Reactor1G(rp, lwr_p2track, n)
{
    libfile = h5lib;
    loadlib(h5lib);
};

LightWaterReactor1G::~LightWaterReactor1G() 
{
};

void LightWaterReactor1G::calc_params()
{
    /** Sets relevent LWR parameters.
     *  Overwrites standard, do-nothing calc_params() function.
     */

    calcSubStreams();

    params_prior_calc["BUd"]  = 0.0;
    params_after_calc["BUd"] = BUd;

    params_prior_calc["U"]  = ms_feed_u.mass;
    params_after_calc["U"] = ms_prod_u.mass;

    params_prior_calc["TRU"]  = ms_feed_tru.mass;
    params_after_calc["TRU"] = ms_prod_tru.mass;

    params_prior_calc["ACT"]  = ms_feed_act.mass;
    params_after_calc["ACT"] = ms_prod_act.mass;

    params_prior_calc["LAN"]  = ms_feed_lan.mass;
    params_after_calc["LAN"] = ms_prod_lan.mass;

    params_prior_calc["FP"]  = 1.0 - ms_feed_act.mass  - ms_feed_lan.mass;
    params_after_calc["FP"] = 1.0 - ms_prod_act.mass - ms_prod_lan.mass;

    return;
};
