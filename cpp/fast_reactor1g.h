// Header for One-Group Fast Reactor Model

#if !defined(_BRIGHT_FR_1G_)
#define _BRIGHT_FR_1G_

//Bright Libs
#include "reactor1g.h"

namespace bright {

  static std::string fr_p2t [] = {"BUd", "TRUCR", "P_NL", "U", "TRU", "ACT", "LAN", "FP"};
  static std::set<std::string> fr_p2track (fr_p2t, fr_p2t+8);

  /*****************************************************/
  /*** Fast Reactor 1G Component Class and Functions ***/
  /*****************************************************/

  class FastReactor1G : public Reactor1G
  {
  /** One-Group Fast Reactor Model Class.
  *  Has default, but overridable values, for a fast reactor parameters that are loaded on initiation.
  */
  public:
    // FastReactor1G Constructors
    FastReactor1G();
    FastReactor1G(std::string lib, std::string n="");
    FastReactor1G(ReactorParameters rp, std::string n="");
    FastReactor1G(std::string lib, ReactorParameters rp, std::string n="");
    ~FastReactor1G();

    //Public access functions
    void calc_params();
  };

// end bright
};
#endif
