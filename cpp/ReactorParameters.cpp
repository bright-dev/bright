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
    coolant_form = std::map<std::string, double>();

    fuel_density = 0.0;
    cladding_density = 0.0;
    coolant_density = 0.0;

    pnl = 0.0;
    BUt = 0.0;
    use_disadvantage_factor = false;
    lattice_type = std::string();
    rescale_hydrogen = false;

    fuel_radius = 0.0;
    void_radius = 0.0;
    clad_radius = 0.0;
    unit_cell_pitch = 0.0;

    open_slots = 0.0;
    total_slots = 0.0;
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

    lwrd.coolant_form["O16"] = 1.0;
    lwrd.coolant_form["H1"]  = 2.0;
    lwrd.coolant_form["B10"] = 0.199 * 550 * pow(10.0, -6);
    lwrd.coolant_form["B11"] = 0.801 * 550 * pow(10.0, -6);

    lwrd.fuel_density = 10.7;
    lwrd.cladding_density = 5.87;
    lwrd.coolant_density = 0.73;

    lwrd.pnl = 0.98;
    lwrd.BUt = 0.0;

    lwrd.use_disadvantage_factor = true;
    lwrd.lattice_type = "Cylindrical";
    lwrd.rescale_hydrogen = true;

    lwrd.fuel_radius = 0.412;
    lwrd.void_radius = 0.4205;
    lwrd.clad_radius = 0.475;
    lwrd.unit_cell_pitch = 1.33;

    lwrd.open_slots = 25.0;
    lwrd.total_slots = 289.0;

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

    frd.use_disadvantage_factor = false;
    frd.lattice_type = "Cylindrical";
    frd.rescale_hydrogen = false;

    frd.fuel_radius = 0.3115;
    frd.void_radius = 0.3115;
    frd.clad_radius = 0.3115;
    frd.unit_cell_pitch = 0.956;

    frd.open_slots = 19.0;
    frd.total_slots = 163.0;

    return frd;
};

ReactorParameters fr_defaults (fill_fr_defaults());

