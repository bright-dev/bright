// FR_1G.h
// Header for One-Group Fast Reactor Model

#if !defined(_Bright_FR_1G_)
#define _Bright_FR_1G_

//C++ stdlib

//HDF5

//Bright Libs
#include "Reactor1G.h"

static std::string fr_p2t [] = {"BUd", "TRUCR", "P_NL", "U", "TRU", "ACT", "LAN", "FP"};
static std::set<std::string> fr_p2track (fr_p2t, fr_p2t+8);

ReactorParameters fillFRDefaults();
extern ReactorParameters FRDefaults;

/*****************************************************/
/*** Fast Reactor 1G Component Class and Functions ***/
/*****************************************************/

class FastReactor1G : public Reactor1G
{
/** One-Group Fast Reactor Model Class.
 *  Has default, but overridable values, for a fast reactor parameters that are loaded on initiation.
 */
protected:
public:
        //Public data

    //FastReactor1G Constructors
    FastReactor1G();
    FastReactor1G(std::string, std::string = "");
        FastReactor1G(ReactorParameters, std::string = "");
        FastReactor1G(std::string, ReactorParameters, std::string = "");
    ~FastReactor1G();

    //Get Functions

    //Set Functions

    //Public access functions
    void setParams();
};

#endif
