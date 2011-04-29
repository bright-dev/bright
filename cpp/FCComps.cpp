// Fuel Cycle Component Namespace

#include "FCComps.h"

/*********************************************/
/*** Fuel Cycle Components Namespace Stuff ***/
/*********************************************/

#ifdef _WIN32
    int null_set [1] = {922350};
    std::set<int> FCComps::track_isos (null_set, null_set+1);
    std::vector<int> FCComps::track_isos_order (null_set, null_set+1);
#else
    int null_set [0] = {};
    std::set<int> FCComps::track_isos (null_set, null_set+0);
    std::vector<int> FCComps::track_isos_order (null_set, null_set+0);
#endif

int FCComps::verbosity  = 0;
int FCComps::write_text = 1;
int FCComps::write_hdf5 = 0;

std::string FCComps::output_filename = "fuel_cycle.h5";


void FCComps::sort_track_isos()
{
    track_isos_order = std::vector<int> (track_isos.begin(), track_isos.end());
    std::sort(track_isos_order.begin(), track_isos_order.end());
};



void FCComps::load_track_isos_hdf5(std::string filename, std::string datasetname, bool clear_prev)
{
    // Check that the file is there
    if (!bright::FileExists(filename))
        throw bright::FileNotFound(filename);

    //Load values into track_isos from an hdf5 file.
    //If the dataspace name is not given, try some defaults.
    int dslen = 14;
    std::string defaultsets [14] = {
        "/track_isos",
        "/Isos2Track",
        "/isostrack",   
        "/IsosTrack",
        "/isotrack",   
        "/IsoTrack",    
        "/ToIso",
        "/ToIsos",
        "/ToIso_zz",
        "/ToIso_MCNP",  
        "/FromIso",  
        "/FromIsos",  
        "/FromIso_zz",
        "/FromIso_MCNP"
        };

    //Open file
    H5::H5File isofile(filename, H5F_ACC_RDONLY);

    //Open dataset
    H5::DataSet isoset;
    if (datasetname.length() != 0)
        isoset = isofile.openDataSet(datasetname);
    else
    {
        //Try to grab one of the default data sets
        //...while suppressing HDF5 errors.
        H5::FileIException not_found_error;
        not_found_error.dontPrint();

        int n = 0;
        while (n < dslen)
        {
            try {  
                isoset = isofile.openDataSet(defaultsets[n]);
                break;
            }
            catch( H5::FileIException not_found_error )
            {
            };
            n++;
        };
        if (n == dslen)
            throw H5::FileIException("load_track_isos", "Dataset not found!");
    };

    //Read in isos from dataset.
    H5::DataSpace isospace = isoset.getSpace();
    hsize_t isolen[1];
    int isodim = isospace.getSimpleExtentDims(isolen, NULL);

    //Try native int data type first
    #ifdef _WIN32
        // Windows VC++ doesn't accept variable length arrays!
        // So let's make the read in array greater than the number of known isotopes...
        int         iso_out_int [5000];
    #else
        int         iso_out_int [isolen[0]];
    #endif
    isoset.read(iso_out_int, H5::PredType::NATIVE_INT);
    //Maybe add other data types in the future... 

    //Clear previous entries
    if (clear_prev)
        track_isos.clear();

    //load into track_isos
    for(int n = 0; n < isolen[0]; n++)
    {
        track_isos.insert(isoname::mixed_2_zzaaam(iso_out_int[n]));
    };

    // Sort the results
    sort_track_isos();
};

void FCComps::load_track_isos_text(std::string filename, bool clear_prev)
{
    // Check that the file is there
    if (!bright::FileExists(filename))
        throw bright::FileNotFound(filename);

    //Clear previous entries
    if (clear_prev)
        track_isos.clear();

    //open file
    std::fstream isofile;
    isofile.open(filename.c_str(), std::fstream::in);

    char isoraw [20];
    std::string isostr;
    while (!isofile.eof())
    {
        isofile.width(20);
        isofile >> isoraw;
        isostr.assign(isoraw);
        isostr = bright::MultiStrip(isostr, "()[],.;{}!#|");
        track_isos.insert(isoname::mixed_2_zzaaam(isostr));
    };

    //close file
    isofile.close();

    // Sort the results
    sort_track_isos();
};




//
// Some HDF5 helpera
//

hsize_t FCComps::iso_LL_dims [1] = {6};
H5::ArrayType FCComps::iso_LL_type = H5::ArrayType(H5::PredType::NATIVE_CHAR, 1, FCComps::iso_LL_dims);

hsize_t FCComps::cinder_g_dims [1] = {63};
H5::ArrayType FCComps::cinder_g_type = H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 1, FCComps::cinder_g_dims);



//
// Isotopic decay HDF5 interface
//

H5::CompType FCComps::make_decay_iso_desc()
{
    //  Makes a decay isotope compound datatype
    H5::CompType didesc( sizeof(decay_iso_struct) );

    didesc.insertMember( "from_iso_LL", HOFFSET(decay_iso_struct, from_iso_LL), FCComps::iso_LL_type);
    didesc.insertMember( "from_iso_zz", HOFFSET(decay_iso_struct, from_iso_zz), H5::PredType::NATIVE_INT);

    didesc.insertMember( "half_life", HOFFSET(decay_iso_struct, half_life), H5::PredType::NATIVE_DOUBLE);
    didesc.insertMember( "decay_const", HOFFSET(decay_iso_struct, decay_const), H5::PredType::NATIVE_DOUBLE);

    didesc.insertMember( "to_iso_LL", HOFFSET(decay_iso_struct, to_iso_LL), FCComps::iso_LL_type);
    didesc.insertMember( "to_iso_zz", HOFFSET(decay_iso_struct, to_iso_zz), H5::PredType::NATIVE_INT);

    didesc.insertMember( "branch_ratio", HOFFSET(decay_iso_struct, branch_ratio), H5::PredType::NATIVE_DOUBLE);
 
    return didesc;
};

H5::CompType FCComps::decay_iso_desc = FCComps::make_decay_iso_desc();



//
// Fission Product HDF5 interface
//

H5::CompType FCComps::make_fission_desc()
{
    //  Makes a fission compound datatype
    H5::CompType fdesc( sizeof(fission_struct) );

    fdesc.insertMember( "iso_LL", HOFFSET(fission_struct, iso_LL), FCComps::iso_LL_type);
    fdesc.insertMember( "iso_zz", HOFFSET(fission_struct, iso_zz), H5::PredType::NATIVE_INT);

    fdesc.insertMember( "thermal_yield", HOFFSET(fission_struct, thermal_yield), H5::PredType::NATIVE_INT8);
    fdesc.insertMember( "fast_yield", HOFFSET(fission_struct, fast_yield), H5::PredType::NATIVE_INT8);
    fdesc.insertMember( "high_energy_yield", HOFFSET(fission_struct, high_energy_yield), H5::PredType::NATIVE_INT8);

    fdesc.insertMember( "xs", HOFFSET(fission_struct, xs), FCComps::cinder_g_type);
 
    return fdesc;
};

H5::CompType FCComps::fission_desc = FCComps::make_fission_desc();



H5::CompType FCComps::make_fission_product_yields_desc()
{
    //  Makes a decay isotope compound datatype
    H5::CompType fpydesc( sizeof(fission_product_yields_struct) );

    fpydesc.insertMember( "index", HOFFSET(fission_product_yields_struct, index), H5::PredType::NATIVE_INT16);

    fpydesc.insertMember( "from_iso_LL", HOFFSET(fission_product_yields_struct, from_iso_LL), FCComps::iso_LL_type);
    fpydesc.insertMember( "from_iso_zz", HOFFSET(fission_product_yields_struct, from_iso_zz), H5::PredType::NATIVE_INT);

    fpydesc.insertMember( "to_iso_LL", HOFFSET(fission_product_yields_struct, to_iso_LL), FCComps::iso_LL_type);
    fpydesc.insertMember( "to_iso_zz", HOFFSET(fission_product_yields_struct, to_iso_zz), H5::PredType::NATIVE_INT);

    fpydesc.insertMember( "mass_frac", HOFFSET(fission_product_yields_struct, mass_frac), H5::PredType::NATIVE_DOUBLE);
 
    return fpydesc;
};

H5::CompType FCComps::fission_product_yields_desc = FCComps::make_fission_product_yields_desc();
