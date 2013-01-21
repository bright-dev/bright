// Teast header!
// Toaster.h

#if !defined(_BRIGHT_TOASTER_)
#define _BRIGHT_TOASTER_

#include "../../../cpp/fccomp.h"

/*********************************************/
/*** Toaster Component Class and Functions ***/
/*********************************************/

namespace bright {

  class Toaster : public FCComp
  {
  // Toaster class
  public:
    // Toaster Constructors
    Toaster();
    ~Toaster();
    
    // Public data
    std::string toastiness;
    unsigned int nslices;
    float rate;

    // Public access functions
    int make_toast(std::string when, unsigned int nslices=1);
  };

// end namespace
};

#endif
