// FCComp.h 
// Header for general Fuel Cycle Component Objects

#if !defined(_Bright_FCComp_)
#define _Bright_FCComp_
#include <math.h>
#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <set>
#include <map>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Cpp.h"

#include "bright.h"
#include "isoname.h"
#include "MassStream.h"

//Declare Global Fuel Cycle Components in 
//Thier own namespace!
namespace FCComps {
    extern std::set<int> isos2track;	            //Set of isotopes to track for all components.
    extern void load_isos2track_hdf5(std::string, std::string = "", bool = false);   //Load isotopic tracking list from HDF5 file.
    extern void load_isos2track_text(std::string, int = 0,          bool = false);   //Load isotopic tracking list from text file.

    extern int verbosity;			//How much should the components talk to us? 0 = None, 1 = a little, 2 = a lot!, etc.
    extern int write_text;
    extern int write_hdf5;

    extern std::set<int> get_isos2track();
    extern int           get_verbosity();
    extern int           get_write_text();
    extern int           get_write_hdf5();

    extern void set_isos2track(std::set<int> i2t);
    extern void set_verbosity(int v);
    extern void set_write_text(int wt);
    extern void set_write_hdf5(int wh);

};

/**************************************************/
/*** Fuel Cycle Component Classes And Functions ***/
/**************************************************/

typedef std::map<std::string, double> ParamDict; 
typedef ParamDict::iterator ParamDictIter;
    
class FCComp
{
//Parent class for all fuel cycle components.
protected:
    //Protected access data
    std::set<std::string> params2track;	//Set of Parameters to track for this component

    //Protected function data
    void initialize (std::set<std::string>, std::string = "");	//initializes empty variables
    void initialize_Text ();	                                //initializes Text output files
    void initialize_HDF5 ();	                                //initializes HDF5 output files

    void appendHDF5array(H5::H5File *, std::string, double *, const int *, hsize_t [], hsize_t [], hsize_t []);

public:
    //Public access data
    std::string name;			//Component name
    MassStream IsosIn;			//Nuclides flowing into the component.
    MassStream IsosOut;			//Nuclides flowing out of the component.
    ParamDict ParamsIn;			//Input paramater values.
    ParamDict ParamsOut;			//Output parameter values.
    int PassNum;				//Cycle Number currently on [int].

    //FCComp Constructors
    FCComp ();
    FCComp (std::string);
    FCComp (std::set<std::string>, std::string = "");

    //Get Functions
    std::set<std::string> get_params2track() const {return params2track;};

    std::string get_name()      const {return name;};
    MassStream  get_IsosIn()    const {return IsosIn;};
    MassStream  get_IsosOut()   const {return IsosOut;};
    ParamDict   get_ParamsIn()  const {return ParamsIn;};
    ParamDict   get_ParamsOut() const {return ParamsOut;};
    int         get_PassNum()   const {return PassNum;};

    //Set Functions
    void set_params2track(std::set<std::string> sset) {params2track = sset;};

    void set_name(std::string s)     {name      = s;};
    void set_IsosIn(MassStream ms)   {IsosIn    = ms;};
    void set_IsosOut(MassStream ms)  {IsosOut   = ms;};
    void set_ParamsIn(ParamDict pd)  {ParamsIn  = pd;};
    void set_ParamsOut(ParamDict pd) {ParamsOut = pd;};
    void set_PassNum(int i)          {PassNum   = i;};

    //Public access functions
    virtual void setParams ();
    void writeIsoPass ();
    void writeParamPass ();
    void writeText ();
    void writeHDF5 ();
    void writeout ();
    virtual MassStream doCalc ();
    virtual MassStream doCalc (CompDict);
    virtual MassStream doCalc (MassStream);
};

#endif
