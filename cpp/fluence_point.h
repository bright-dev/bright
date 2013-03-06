// FluencePoint.h
// Header for Flunece point class that helps out reactors

#if !defined(_BRIGHT_FLUENCE_POINT_)
#define _BRIGHT_FLUENCE_POINT_

namespace bright {

  class FluencePoint
  {
  public:
    // Constructors
    FluencePoint(int f_=0, double F_=0.0, double m_=0.0);
    ~FluencePoint();

    // Attributes
    int f;
    double F;
    double m;
  };

// end bright
};

#endif 
