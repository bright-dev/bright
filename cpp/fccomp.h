// FCComp.h 
// Header for general Fuel Cycle Component Objects

#if !defined(_BRIGHT_FCCOMP_)
#define _BRIGHT_FCCOMP_

#include "utils.h"

/**************************************************/
/*** Fuel Cycle Component Classes And Functions ***/
/**************************************************/

namespace bright {

  typedef std::map<std::string, double> param_dict; 
  typedef param_dict::iterator param_dict_iter;
    
  class FCComp
  {
  // Parent class for all fuel cycle components.
  protected:
    // Protected access data

    // Protected function data
    void initialize(std::set<std::string> paramtrack, std::string n=""); // initializes empty variables
    void initialize_text();	                                  // initializes Text output files
    void initialize_hdf5();	                                  // initializes HDF5 output files

    void appendHDF5array(H5::H5File * dbFile, std::string set_name, \
      double * append_value, const int * rank, hsize_t dims [], hsize_t offset[], \
      hsize_t extend_size[]);

  public:
    // FCComp Constructors
    FCComp (std::string n="");
    FCComp (std::set<std::string> paramtrack, std::string n="");
    ~FCComp ();

    // Public access data
    std::string name;			              // Component name
    std::string natural_name;           // Component natural name
    pyne::Material mat_feed;			      // Nuclides flowing into the component.
    pyne::Material mat_prod;            // Nuclides flowing out of the component.
    param_dict params_prior_calc;			  // Input paramater values.
    param_dict params_after_calc;		    // Output parameter values.
    int pass_num;			        	        // Cycle Number currently on [int].
    std::set<std::string> track_params; // Set of Parameters to track for this component

    // Public access functions
    virtual void calc_params();
    void write_mat_pass();
    void write_params_pass();
    void write_text();
    void write_hdf5();
    void write();
    virtual pyne::Material calc();
    virtual pyne::Material calc(pyne::comp_map incomp);
    virtual pyne::Material calc(pyne::Material mat);
  };

// end bright
};

#endif
