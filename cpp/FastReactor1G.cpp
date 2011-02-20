// One-Group Fast Reactor Component Class

#include "FastReactor1G.h"

ReactorParameters fillFRDefaults ()
{
    //Default FR physical parameters
    ReactorParameters frd;

    frd.batches = 3;
    frd.flux = 2.0 * pow(10.0, 15);

    frd.FuelForm["IHM"] = 1.0;
    frd.CoolantForm["NA23"] = 1.0;

    frd.FuelDensity = 18.0;
    frd.CoolantDensity = 0.927;

    frd.pnl = 0.65;
    frd.BUt = 0.0;

    frd.useDisadvantage = false;
    frd.LatticeType = "Cylindrical";
    frd.HydrogenRescale = false;
    
    frd.Radius = 0.3115;
    frd.Length = 0.956;
    frd.open_slots = 19.0;
    frd.total_slots = 163.0;

    return frd;
};
ReactorParameters FRDefaults (fillFRDefaults());

FastReactor1G::FastReactor1G() : Reactor1G(FRDefaults, fr_p2track)
{
};

FastReactor1G::FastReactor1G(std::string h5lib, std::string n) : Reactor1G(FRDefaults, fr_p2track, n)
{
    libfile = h5lib;
    loadLib(h5lib);
};

FastReactor1G::FastReactor1G(ReactorParameters rp, std::string n) : Reactor1G(rp, fr_p2track, n)
{
};

FastReactor1G::FastReactor1G(std::string h5lib, ReactorParameters rp, std::string n) : Reactor1G(rp, fr_p2track, n)
{
    libfile = h5lib;
    loadLib(h5lib);
};

FastReactor1G::~FastReactor1G() 
{
};

void FastReactor1G::calc_params()
{
    /** Sets relevent FR parameters.
     *  Overwrites standard, do-nothing calc_params() function.
     */

    calcSubStreams();

    params_prior_calc["BUd"]  = 0.0;
    params_after_calc["BUd"] = BUd;

    params_prior_calc["TRUCR"]  = 0.0;
    params_after_calc["TRUCR"] = calc_tru_cr();

    params_prior_calc["P_NL"]  = 0.0;
    params_after_calc["P_NL"] = P_NL;

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
