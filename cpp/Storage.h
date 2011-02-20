// Storage.h
// Header for general Fuel Cycle Component Objects

#if !defined(_Bright_Storage_)
#define _Bright_Storage_
#include <map>
#include <set>
#include <vector>
#include <math.h>
#include <string.h>
#include "isoname.h"
#include "MassStream.h"
#include "FCComp.h"
#include "hdf5.h"

/*********************************************/
/*** Storage Component Class and Functions ***/
/*********************************************/

typedef struct DecayIso
{
    int 	fromiso;
    double	halflife;
    double	decayconst;
    int 	toiso;
    double	branchratio;
} DecayIso;

//Sexy decay library form
typedef std::map<int, double *> ToIsoDict;
typedef ToIsoDict::iterator ToIsoIter;
typedef struct FromIsoStuct
{
    double * halflife;
    double * decayconst;
    ToIsoDict toiso;
} FromIsoStruct;
typedef std::map<int, FromIsoStruct> DecayDict;
typedef DecayDict::iterator DecayIter;

typedef std::vector<int> IsoChain;
typedef IsoChain::iterator IsoChainIter;
typedef std::set<IsoChain> IsoChainSet;
typedef IsoChainSet::iterator IsoChainSetIter;

static std::string stor_p2t [] = {"Mass"};
static std::set<std::string> stor_p2track (stor_p2t, stor_p2t+1);

class Storage : public FCComp
{
//Storage/Cooling/Decay Fuel Cycle Component.
protected:
    //Protected Data
    IsoChainSet isochains;
    DecayIso * decay_data;
    int decay_data_len;
    DecayDict decay;

    //Protected functions
    void initialize ();						//Initializes the constructors.
    double getDecay ();
    double bateman (int, double, IsoChain);				//Solves the Bateman Decay equation.
    void addchains (IsoChain);				//Get the decay chain for a mother isotope
    void addchains (int);						//Get the decay chain for a mother isotope

    void PrintChain (IsoChain);

public:
    //Storage Constructors	
    Storage ();
    Storage(std::string);
    ~Storage ();

    //Public data
    double decay_time;			//time to decay for

    //Public Functions
    void setParams();
    MassStream doCalc ();
    MassStream doCalc (CompDict);
    MassStream doCalc (MassStream);
    MassStream doCalc (double);
    MassStream doCalc (CompDict, double);
    MassStream doCalc (MassStream, double);
};

#endif
