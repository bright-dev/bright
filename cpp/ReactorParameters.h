// ReactorParameters.h
// Header for ReactorParameters reactor helper class

#if !defined(_Bright_ReactorParameters_)
#define _Bright_ReactorParameters_

#include <map>
#include <string>

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
    std::map<std::string, double> coolant_form;
    double fuel_density;
    double coolant_density;
    double pnl;
    double BUt;
    bool use_disadvantage_factor;
    std::string lattice_type;
    bool rescale_hydrogen;
    double radius;
    double pitch;
    double open_slots;
    double total_slots;

};

#endif
