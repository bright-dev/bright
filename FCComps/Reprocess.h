// Reprocess.h
// Header for general Fuel Cycle Component Objects

#if !defined(_Bright_Reprocess_)
#define _Bright_Reprocess_
#include <map>
#include <set>
#include <string>
#include <exception>
#include "FCComp.h"
#include "MassStream.h"
#include "isoname.h"

/**************************************************/
/*** Reprocessing Component Class and Functions ***/
/**************************************************/

typedef std::map<int, double> SepEffDict;
typedef SepEffDict::iterator SepEffIter;

static std::string rep_p2t [] = {"Mass"};
static std::set<std::string> rep_p2track (rep_p2t, rep_p2t+1);


class Reprocess : public FCComp
{
//Reprocessing class
protected:
    //Protected functions

public:
    //Public data
    SepEffDict sepeff;			//separation efficiency dictionary

    //Reprocessing Constructors
    Reprocess ();
    Reprocess (SepEffDict, std::string = "");
    Reprocess (std::map<std::string, double>, std::string = "");
    
    //Get Functions
    SepEffDict get_sepeff() const {return sepeff;};

    //Set Functions
    void set_sepeff(SepEffDict sed) {sepeff = sed;};

    //Public access functions
    void initialize(SepEffDict);		//Initializes the constructors.
    void setParams ();
    MassStream doCalc ();
    MassStream doCalc (CompDict);
    MassStream doCalc (MassStream);	
};

#endif
