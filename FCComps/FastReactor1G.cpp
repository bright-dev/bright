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

void FastReactor1G::setParams()
{
    /** Sets relevent FR parameters.
     *  Overwrites standard, do-nothing setParams() function.
     */

    calcSubStreams();

    ParamsIn["BUd"]  = 0.0;
    ParamsOut["BUd"] = BUd;

    ParamsIn["TRUCR"]  = 0.0;
    ParamsOut["TRUCR"] = calcTruCR();

    ParamsIn["P_NL"]  = 0.0;
    ParamsOut["P_NL"] = P_NL;

    ParamsIn["U"]  = InU.mass;
    ParamsOut["U"] = OutU.mass;

    ParamsIn["TRU"]  = InTRU.mass;
    ParamsOut["TRU"] = OutTRU.mass;

    ParamsIn["ACT"]  = InACT.mass;
    ParamsOut["ACT"] = OutACT.mass;

    ParamsIn["LAN"]  = InLAN.mass;
    ParamsOut["LAN"] = OutLAN.mass;

    ParamsIn["FP"]  = 1.0 - InACT.mass  - InLAN.mass;
    ParamsOut["FP"] = 1.0 - OutACT.mass - OutLAN.mass;

    return;
};
