// ReactorParameters.h
// Header for ReactorParameters reactor helper class

#if !defined(_BRIGHT_REACTOR_PARAMETERS_)
#define _BRIGHT_REACTOR_PARAMETERS_

#include <math.h>
#include <map>
#include <vector>
#include <string>

namespace bright {

  class ReactorParameters
  {
  public:
    // Constructors
    ReactorParameters();
    ~ReactorParameters();

    // Attributes
    int batches;
    double flux;

    std::map<std::string, double> fuel_form;
    std::map<std::string, double> cladding_form;
    std::map<std::string, double> coolant_form;

    double fuel_density;
    double cladding_density;
    double coolant_density;

    double pnl;
    double BUt;
    double specific_power;
    int burn_regions;
    std::vector<double> burn_times;

    bool use_disadvantage_factor;
    std::string lattice_type;
    bool rescale_hydrogen;
    std::string burnup_via_constant;

    double fuel_radius;
    double void_radius;
    double clad_radius;
    double unit_cell_pitch;

    double open_slots;
    double total_slots;

    double branch_ratio_cutoff;
  };


  // Preset Defaults
  ReactorParameters fill_lwr_defaults();
  extern ReactorParameters lwr_defaults;

  ReactorParameters fill_fr_defaults();
  extern ReactorParameters fr_defaults;

// end bright
};

#endif
