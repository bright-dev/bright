// FR_1G.h
// Header for One-Group Light Water Reactor Model

#if !defined(_Bright_LWR_1G_)
#define _Bright_LWR_1G_

//C++ stdlib

//HDF5

//Bright Libs
#include "Reactor1G.h"

static std::string lwr_p2t [] = {"BUd", "U", "TRU", "ACT", "LAN", "FP"};
static std::set<std::string> lwr_p2track (lwr_p2t, lwr_p2t+6);

ReactorParameters fillLWRDefaults();
extern ReactorParameters LWRDefaults;

/************************************************************/
/*** Light Water Reactor 1G Component Class and Functions ***/
/************************************************************/

class LightWaterReactor1G : public Reactor1G
{
/** One-Group Light Water Reactor Model Class.
 *  Has default, but overridable values, for light water reactors that are loaded on initiation.
 */
protected:
public:
        //Public data

    //LightWaterReactor1G Constructors
    LightWaterReactor1G();
    LightWaterReactor1G(std::string, std::string = "");
    LightWaterReactor1G(ReactorParameters, std::string = "");
    LightWaterReactor1G(std::string, ReactorParameters, std::string = "");
    ~LightWaterReactor1G();

    //Get Functions

    //Set Functions

    //Public access functions
    void calc_params();
};

#endif
