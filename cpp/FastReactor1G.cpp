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
    params_after_calc["TRUCR"] = calcTruCR();

    params_prior_calc["P_NL"]  = 0.0;
    params_after_calc["P_NL"] = P_NL;

    params_prior_calc["U"]  = InU.mass;
    params_after_calc["U"] = OutU.mass;

    params_prior_calc["TRU"]  = InTRU.mass;
    params_after_calc["TRU"] = OutTRU.mass;

    params_prior_calc["ACT"]  = InACT.mass;
    params_after_calc["ACT"] = OutACT.mass;

    params_prior_calc["LAN"]  = InLAN.mass;
    params_after_calc["LAN"] = OutLAN.mass;

    params_prior_calc["FP"]  = 1.0 - InACT.mass  - InLAN.mass;
    params_after_calc["FP"] = 1.0 - OutACT.mass - OutLAN.mass;

    return;
};
