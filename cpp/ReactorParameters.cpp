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
    coolant_density = 0.0;
    pnl = 0.0;
    BUt = 0.0;
    use_disadvantage_factor = false;
    lattice_type = std::string();
    rescale_hydrogen = false;
    radius = 0.0;
    pitch = 0.0;
    open_slots = 0.0;
    total_slots = 0.0;
};


ReactorParameters::~ReactorParameters()
{
};
