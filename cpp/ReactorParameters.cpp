// ReactorParameters Reactor Helper Class

#include "ReactorParameters.h"

/*******************************/
/*** ReactorParameters Class ***/
/*******************************/

ReactorParameters::ReactorParameters()
{
    batches = 0;
    flux = 0.0;

    fuel_form = std::map<std::string, double>();
    cladding_form = std::map<std::string, double>();
    coolant_form = std::map<std::string, double>();

    fuel_density = 0.0;
    cladding_density = 0.0;
    coolant_density = 0.0;

    pnl = 0.0;
    BUt = 0.0;
    specific_power = 0.0;
    burn_regions = 0;
    burn_times = std::vector<double>();

    use_disadvantage_factor = false;
    lattice_type = std::string();
    rescale_hydrogen = false;
    burnup_via_constant = "";

    fuel_radius = 0.0;
    void_radius = 0.0;
    clad_radius = 0.0;
    unit_cell_pitch = 0.0;

    open_slots = 0.0;
    total_slots = 0.0;

    branch_ratio_cutoff = 0.0;
};


ReactorParameters::~ReactorParameters()
{
};



/************/
/* Defaults */
/************/

// Light Water Reactor Defaults

ReactorParameters fill_lwr_defaults ()
{
    //Default LWR physical parameters
    ReactorParameters lwrd;

    lwrd.batches = 3;
    lwrd.flux = 4.0 * pow(10.0, 14);

    lwrd.fuel_form["IHM"] = 1.0;
    lwrd.fuel_form["O16"] = 2.0;

    // Default zircaloy
    // Natural Zirconium
    lwrd.cladding_form["ZR90"] = 0.98135 * 0.5145;
    lwrd.cladding_form["ZR91"] = 0.98135 * 0.1122;
    lwrd.cladding_form["ZR92"] = 0.98135 * 0.1715;
    lwrd.cladding_form["ZR94"] = 0.98135 * 0.1738;
    lwrd.cladding_form["ZR96"] = 0.98135 * 0.0280;
    // The plastic is all melted and the natural Chromium too..
    lwrd.cladding_form["CR50"] = 0.00100 * 0.04345;
    lwrd.cladding_form["CR52"] = 0.00100 * 0.83789;
    lwrd.cladding_form["CR53"] = 0.00100 * 0.09501;
    lwrd.cladding_form["CR54"] = 0.00100 * 0.02365;
    // Natural Iron
    lwrd.cladding_form["FE54"] = 0.00135 * 0.05845;
    lwrd.cladding_form["FE56"] = 0.00135 * 0.91754;
    lwrd.cladding_form["FE57"] = 0.00135 * 0.02119;
    lwrd.cladding_form["FE58"] = 0.00135 * 0.00282;
    // Natural Nickel
    lwrd.cladding_form["NI58"] = 0.00055 * 0.68077;
    lwrd.cladding_form["NI60"] = 0.00055 * 0.26223;
    lwrd.cladding_form["NI61"] = 0.00055 * 0.01140;
    lwrd.cladding_form["NI62"] = 0.00055 * 0.03634;
    lwrd.cladding_form["NI64"] = 0.00055 * 0.00926;
    // Natural Tin
    lwrd.cladding_form["SN112"] = 0.01450 * 0.0097;
    lwrd.cladding_form["SN114"] = 0.01450 * 0.0065;
    lwrd.cladding_form["SN115"] = 0.01450 * 0.0034;
    lwrd.cladding_form["SN116"] = 0.01450 * 0.1454;
    lwrd.cladding_form["SN117"] = 0.01450 * 0.0768;
    lwrd.cladding_form["SN118"] = 0.01450 * 0.2422;
    lwrd.cladding_form["SN119"] = 0.01450 * 0.0858;
    lwrd.cladding_form["SN120"] = 0.01450 * 0.3259;
    lwrd.cladding_form["SN122"] = 0.01450 * 0.0463;
    lwrd.cladding_form["SN124"] = 0.01450 * 0.0579;
    // We Need Oxygen!
    lwrd.cladding_form["O16"] = 0.00125;

    lwrd.coolant_form["O16"] = 1.0;
    lwrd.coolant_form["H1"]  = 2.0;
    lwrd.coolant_form["B10"] = 0.199 * 550 * pow(10.0, -6);
    lwrd.coolant_form["B11"] = 0.801 * 550 * pow(10.0, -6);

    lwrd.fuel_density = 10.7;
    lwrd.cladding_density = 5.87;
    lwrd.coolant_density = 0.73;

    lwrd.pnl = 0.98;
    lwrd.BUt = 0.0;
    lwrd.specific_power = 0.04;
    lwrd.burn_regions = 1;
    lwrd.burn_times = std::vector<double>(144, -1.0);
    for (int t = 0; t < 144; t++)
    {
        lwrd.burn_times[t] = 30.0 * t;
    };


    lwrd.use_disadvantage_factor = true;
    lwrd.lattice_type = "Cylindrical";
    lwrd.rescale_hydrogen = true;
    lwrd.burnup_via_constant = "power";

    lwrd.fuel_radius = 0.412;
    lwrd.void_radius = 0.4205;
    lwrd.clad_radius = 0.475;
    lwrd.unit_cell_pitch = 1.33;

    lwrd.open_slots = 25.0;
    lwrd.total_slots = 289.0;

    lwrd.branch_ratio_cutoff = 1E-20;

    return lwrd;
};

ReactorParameters lwr_defaults (fill_lwr_defaults());



// Fast Reactor Defaults

ReactorParameters fill_fr_defaults ()
{
    //Default FR physical parameters
    ReactorParameters frd;

    frd.batches = 3;
    frd.flux = 2.0 * pow(10.0, 15);

    frd.fuel_form["IHM"] = 1.0;
    frd.coolant_form["NA23"] = 1.0;

    frd.fuel_density = 18.0;
    frd.coolant_density = 0.927;

    frd.pnl = 0.65;
    frd.BUt = 0.0;
    frd.specific_power = 0.04;
    frd.burn_regions = 1;
    frd.burn_times = std::vector<double>(144, -1.0);
    for (int t = 0; t < 144; t++)
    {
        frd.burn_times[t] = 31.0 * t;
    };

    frd.use_disadvantage_factor = false;
    frd.lattice_type = "Cylindrical";
    frd.rescale_hydrogen = false;
    frd.burnup_via_constant = "power";

    frd.fuel_radius = 0.3115;
    frd.void_radius = 0.3115;
    frd.clad_radius = 0.3115;
    frd.unit_cell_pitch = 0.956;

    frd.open_slots = 19.0;
    frd.total_slots = 163.0;

    frd.branch_ratio_cutoff = 1E-20;

    return frd;
};

ReactorParameters fr_defaults (fill_fr_defaults());

