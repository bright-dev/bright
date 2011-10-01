// Reprocess.h
// Header for general Fuel Cycle Component Objects

#if !defined(_Bright_Reprocess_)
#define _Bright_Reprocess_
#include <map>
#include <set>
#include <string>
#include <exception>
#include "FCComp.h"
#include "pyne::Material.h"
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
    //Reprocessing Constructors
    Reprocess ();
    Reprocess (SepEffDict, std::string = "");
    Reprocess (std::map<std::string, double>, std::string = "");
    ~Reprocess ();
    
    //Public data
    SepEffDict sepeff;			//separation efficiency dictionary

    //Public access functions
    void initialize(SepEffDict);		//Initializes the constructors.
    void calc_params ();
    pyne::Material calc ();
    pyne::Material calc (pyne::comp_map);
    pyne::Material calc (pyne::Material);	
};

#endif
