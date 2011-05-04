// FCComp.h 
// Header for general Fuel Cycle Component Namespace

#if !defined(_Bright_FCComps_)
#define _Bright_FCComps_

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <stdint.h>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Cpp.h"

#include "bright.h"
#include "isoname.h"
#include "MassStream.h"

// Declare Global Fuel Cycle Components in 
// Their own namespace!

namespace FCComps {
    extern std::set<int> track_isos;	      // Set of isotopes to track for all components.
    extern std::vector<int> track_isos_order; // Vector of isotopes to track for all components.

    extern void load_track_isos_hdf5(std::string, std::string = "", bool = false);  //Load isotopic tracking list from HDF5 file.
    extern void load_track_isos_text(std::string, bool = false);                    //Load isotopic tracking list from text file.

    extern void sort_track_isos(); // Sets the isotopic tracking by zzaaam from lowest to highest and stores it in track_isos_order

    extern int verbosity;			//How much should the components talk to us? 0 = None, 1 = a little, 2 = a lot!, etc.
    extern int write_text;
    extern int write_hdf5;

    extern std::string output_filename;


    // Some HDF5 helpers
    extern H5::StrType iso_LL_type;

    extern hsize_t cinder_g_dims [1];
    extern H5::ArrayType cinder_g_type;


    // Isotopic decay HDF5 interface
    typedef struct decay_iso_struct {
        char from_iso_LL[6];
        int from_iso_zz;

        double half_life;
        double decay_const;

        char to_iso_LL[6];
        int to_iso_zz;

        double branch_ratio;
    } decay_iso_stuct;

    H5::CompType make_decay_iso_desc();
    extern H5::CompType decay_iso_desc;


    // Fission Product HDF5 interface
    typedef struct fission_struct {
        char iso_LL[6];
        int iso_zz;

        int8_t thermal_yield;
        int8_t fast_yield;
        int8_t high_energy_yield;

        double xs [63];
    } fission_struct;

    H5::CompType make_fission_desc();
    extern H5::CompType fission_desc;


    typedef struct fission_product_yields_struct {
        int16_t index;

        char from_iso_LL[6];
        int from_iso_zz;

        char to_iso_LL[6];
        int to_iso_zz;

        double mass_frac;
    } fission_product_yields_struct;

    H5::CompType make_fission_product_yields_desc();
    extern H5::CompType fission_product_yields_desc;




};

#endif
