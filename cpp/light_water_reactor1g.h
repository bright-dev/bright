// Header for One-Group Light Water Reactor Model

#if !defined(_BRIGHT_LWR_1G_)
#define _BRIGHT_LWR_1G_


//Bright Libs
#include "reactor1g.h"

namespace bright{

  static std::string lwr_p2t [] = {"BUd", "U", "TRU", "ACT", "LAN", "FP"};
  static std::set<std::string> lwr_p2track (lwr_p2t, lwr_p2t+6);

  /************************************************************/
  /*** Light Water Reactor 1G Component Class and Functions ***/
  /************************************************************/

  class LightWaterReactor1G : public Reactor1G
  {
  /** One-Group Light Water Reactor Model Class.
  *  Has default, but overridable values, for light water reactors that are loaded on initiation.
  */
  public:
    // LightWaterReactor1G Constructors
    LightWaterReactor1G();
    LightWaterReactor1G(std::string lib, std::string n="");
    LightWaterReactor1G(ReactorParameters rp, std::string n="");
    LightWaterReactor1G(std::string lib, ReactorParameters rp, std::string n="");
    ~LightWaterReactor1G();

    // Public access functions
    void calc_params();
  };

// end bright
};
#endif
