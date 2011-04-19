// FCComp.h 
// Header for general Fuel Cycle Component Objects

#if !defined(_Bright_FCComp_)
#define _Bright_FCComp_
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <set>
#include <map>
#include <vector>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Cpp.h"

#include "bright.h"
#include "isoname.h"
#include "MassStream.h"
#include "FCComps.h"

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

    //Protected function data
    void initialize (std::set<std::string>, std::string = "");	//initializes empty variables
    void initialize_Text ();	                                //initializes Text output files
    void initialize_HDF5 ();	                                //initializes HDF5 output files

    void appendHDF5array(H5::H5File *, std::string, double *, const int *, hsize_t [], hsize_t [], hsize_t []);

public:
    //FCComp Constructors
    FCComp ();
    FCComp (std::string);
    FCComp (std::set<std::string>, std::string = "");
    ~FCComp ();

    //Public access data
    std::string name;			        //Component name
    std::string natural_name;           //Component natural name
    MassStream ms_feed;			        //Nuclides flowing into the component.
    MassStream ms_prod;		        	//Nuclides flowing out of the component.
    ParamDict params_prior_calc;			        //Input paramater values.
    ParamDict params_after_calc;		        //Output parameter values.
    int pass_num;			        	//Cycle Number currently on [int].
    std::set<std::string> track_params;	//Set of Parameters to track for this component

    //Public access functions
    virtual void calc_params ();
    void write_ms_pass ();
    void write_params_pass ();
    void write_text ();
    void write_hdf5 ();
    void write ();
    virtual MassStream calc ();
    virtual MassStream calc (CompDict);
    virtual MassStream calc (MassStream);
};

#endif
