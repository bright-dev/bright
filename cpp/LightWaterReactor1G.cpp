// One-Group Light Water Reactor Component Class

#include "LightWaterReactor1G.h"

ReactorParameters fillLWRDefaults ()
{
    //Default LWR physical parameters
    ReactorParameters lwrd;

    lwrd.batches = 3;
    lwrd.flux = 4.0 * pow(10.0, 14);

    lwrd.FuelForm["IHM"] = 1.0;
    lwrd.FuelForm["O16"] = 2.0;
    lwrd.CoolantForm["O16"] = 1.0;
    lwrd.CoolantForm["H1"]  = 2.0;
    lwrd.CoolantForm["B10"] = 0.199 * 550 * pow(10.0, -6);
    lwrd.CoolantForm["B11"] = 0.801 * 550 * pow(10.0, -6);

    lwrd.FuelDensity = 10.7;
    lwrd.CoolantDensity = 0.73;

    lwrd.pnl = 0.98;
    lwrd.BUt = 0.0;

    lwrd.useDisadvantage = true;
    lwrd.LatticeType = "Cylindrical";
    lwrd.HydrogenRescale = true;

    lwrd.Radius = 0.412;
    lwrd.Length = 1.33;
    lwrd.open_slots = 25.0;
    lwrd.total_slots = 289.0;

    return lwrd;
};
ReactorParameters LWRDefaults (fillLWRDefaults());

LightWaterReactor1G::LightWaterReactor1G() : Reactor1G(LWRDefaults, lwr_p2track)
{
};

LightWaterReactor1G::LightWaterReactor1G(std::string h5lib, std::string n) : Reactor1G(LWRDefaults, lwr_p2track, n)
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
