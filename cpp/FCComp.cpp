// Fuel Cycle Component Parent Class

#include "FCComp.h"

/**************************************************/
/*** Fuel Cycle Component Classes And Functions ***/
/**************************************************/

/***************************/
/*** Protected Functions ***/
/***************************/

void FCComp::initialize (std::set<std::string> ptrack, std::string n)
{
    //Protected Variables
    track_params = ptrack;

    //Public Variables
    name = n;
    natural_name = bright::natural_naming(n);
    if (natural_name.length() == 0)
        natural_name = "this_is_not_a_name";
    
    pass_num = 0;

    if (bright::write_text)
        initialize_Text();

    if (bright::write_hdf5)
        initialize_HDF5();
}
    
void FCComp::initialize_Text ()
{
    //Initialize the Isotopic tracking file
    if (!bright::track_isos.empty())
    {
        std::ofstream isofile ( (name + "Isos.txt").c_str() );
        isofile << "Isotope\n";
        for (std::set<int>::iterator iso = bright::track_isos.begin(); iso != bright::track_isos.end(); iso++)
        {
            isofile << pyne::nucname::zzaaam_2_LLAAAM(*iso) << "\n"; 
        }
        isofile.close();
    }

    //Initialize the Parameter tracking file.
    if (!track_params.empty())
    {	
        std::ofstream paramfile ( (name + "Params.txt").c_str() );
        paramfile << "Param\n";
        for ( std::set<std::string>::iterator p = track_params.begin(); p != track_params.end(); p++)
        {
                paramfile << *p + "\n";
        }
        paramfile.close();
    }

}

void FCComp::initialize_HDF5 ()
{
    // Turn off annoying HDF5 errors
    H5::Exception::dontPrint();

    //Create 1D arrays of doubles.
    const int     RANK = 1;
    hsize_t dims[1]    = {0};               // dataset dimensions at creation
    hsize_t maxdims[1] = {H5S_UNLIMITED};
    H5::DataSpace ext_1D_space(RANK, dims, maxdims);

    //Create new/open datafile.
    H5::H5File dbFile;
    if (bright::FileExists(bright::output_filename))
        dbFile = H5::H5File(bright::output_filename, H5F_ACC_RDWR);
    else
        dbFile = H5::H5File(bright::output_filename, H5F_ACC_TRUNC);

    //Modify dataset creation properties.
    H5::DSetCreatPropList double_params;

    hsize_t chunk_dims[1] = {10};
    double_params.setChunk(RANK, chunk_dims);

    double fill_val = -1.0;
    double_params.setFillValue(H5::PredType::NATIVE_DOUBLE, &fill_val);

    /***
     NOTE!  All of the stupid try/catch, open/create stuff is simply to test 
            if a group or dataset already exist.  If it does, open() passes.
            If it doesn't exist, open() fails and a default group/dataset
            is created in its place.
     ***/

    // Open/Create group for this FCComp
    std::string comp_path ("/" + natural_name);
    H5::Group gFCComp;
    try 
        { gFCComp = dbFile.openGroup(comp_path); }
    catch (H5::Exception fgerror) 
        { gFCComp = dbFile.createGroup(comp_path); }

    //Initialize the IsoStreams 
    if (!bright::track_isos.empty())
    {
        // Open/Create ms_feed group
        H5::Group gms_feed;
        try
          { gms_feed = dbFile.openGroup(comp_path + "/ms_feed"); }
        catch (H5::Exception fgerror)
          { gms_feed = dbFile.createGroup(comp_path + "/ms_feed"); }

        // Open/Create ms_prod group
        H5::Group gms_prod;
        try
            { gms_prod = dbFile.openGroup(comp_path + "/ms_prod"); }
        catch (H5::Exception fgerror)
            { gms_prod = dbFile.createGroup(comp_path + "/ms_prod"); }

        // Open/Create /ms_feed/Mass Dataset
        H5::DataSet dsms_feedMass;
        try 
            { dsms_feedMass = dbFile.openDataSet(comp_path + "/ms_feed/Mass"); }
        catch (H5::Exception fgerror)
            { dsms_feedMass = dbFile.createDataSet(comp_path + "/ms_feed/Mass",  H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }

        // Open/Create /ms_prod/Mass Dataset
        H5::DataSet dsms_prodMass;
        try
            { dsms_prodMass = dbFile.openDataSet(comp_path + "/ms_prod/Mass"); }
        catch (H5::Exception fgerror)
            { dsms_prodMass = dbFile.createDataSet(comp_path + "/ms_prod/Mass", H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }

        // Open/Create /Isos[In|Out]/iso Datasets
        H5::DataSet dsms_feedIso;
        H5::DataSet dsms_prodIso;
        for (std::set<int>::iterator iso = bright::track_isos.begin(); iso != bright::track_isos.end(); iso++)
        {
            std::string isoLL = pyne::nucname::zzaaam_2_LLAAAM(*iso);

            try
                { dsms_feedIso = dbFile.openDataSet(comp_path + "/ms_feed/" + isoLL); }
            catch (H5::Exception fgerror)
                { dsms_feedIso = dbFile.createDataSet(comp_path + "/ms_feed/" + isoLL, H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }

            try
                { dsms_prodIso = dbFile.openDataSet(comp_path + "/ms_prod/" + isoLL); }
            catch (H5::Exception fgerror)
                { dsms_prodIso = dbFile.createDataSet(comp_path + "/ms_prod/" + isoLL, H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }
        }
    }

    //Initiallize the Parameters
    if (!track_params.empty())
    {	
        // Open/Create params_prior_calc group
        H5::Group gparams_prior_calc;
        try
          { gparams_prior_calc = dbFile.openGroup(comp_path + "/params_prior_calc"); }
        catch (H5::Exception fgerror)
          { gparams_prior_calc = dbFile.createGroup(comp_path + "/params_prior_calc"); }

        // Open/Create params_after_calc group
        H5::Group gparams_after_calc;
        try
            { gparams_after_calc = dbFile.openGroup(comp_path + "/params_after_calc"); }
        catch (H5::Exception fgerror)
            { gparams_after_calc = dbFile.createGroup(comp_path + "/params_after_calc"); }

        // Open/Create /Params[In|Out]/param Datasets
        H5::DataSet dsparams_prior_calcParam;
        H5::DataSet dsparams_after_calcParam;
        for ( std::set<std::string>::iterator p = track_params.begin(); p != track_params.end(); p++)
        {
            try
                { dsparams_prior_calcParam = dbFile.openDataSet(comp_path + "/params_prior_calc/"  + *p); }
            catch (H5::Exception fgerror)
                { dsparams_prior_calcParam = dbFile.createDataSet(comp_path + "/params_prior_calc/"  + *p, H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }

            try
                { dsparams_after_calcParam = dbFile.openDataSet(comp_path + "/params_after_calc/" + *p); }
            catch (H5::Exception fgerror)
                { dsparams_after_calcParam = dbFile.createDataSet(comp_path + "/params_after_calc/" + *p, H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }
        }
    }

    //Close out thr HDF5 database file
    dbFile.close();
}

/***************************/
/*** FCComp Constructors ***/
/***************************/
    
FCComp::FCComp ()
{
    //Parent class for all fuel cycle components.
    std::set<std::string> emptystrset;
    initialize(emptystrset);
}

FCComp::FCComp (std::string n)
{
    //Parent class for all fuel cycle components.
    std::set<std::string> emptystrset;
    initialize(emptystrset, n);
}

FCComp::FCComp (std::set<std::string> ptrack, std::string n)
{
    //Parent class for all fuel cycle components.
    initialize(ptrack, n);	
}

FCComp::~FCComp ()
{
}


/*******************************/
/*** Public Access Functions ***/
/*******************************/

void FCComp::calc_params ()
{
    //Placeholder function that sets the states of params_prior_calc and params_after_calc.
    for ( std::set<std::string>::iterator p2t = track_params.begin(); p2t != track_params.end(); p2t++)
    {
        params_prior_calc[*p2t]  = 0.0;
        params_after_calc[*p2t] = 0.0;
    }
}

void FCComp::write_ms_pass ()
{
    //Writes a single pass to the isotopic tracking file.
    std::ifstream isofilein  ( (name + "Isos.txt").c_str() );
    std::stringstream isobuf;
    isobuf.precision(6);
    isobuf << std::scientific << std::uppercase;

    //normalize ms_prod per kgms_feed
    double outmassfrac = ms_prod.mass / ms_feed.mass;

    while (!isofilein.eof() )
    {
        char line [3000];
        isofilein.getline(line, 3000);
        isobuf << line;

        std::string isoflag = bright::getFlag(line, 10);
        if (isoflag == "Isotope")
                isobuf << "\t" << bright::to_str(pass_num) << "in\t\t" << bright::to_str(pass_num) << "out\t";
        else
        {
            try
            {
                int isoInLine = pyne::nucname::LLAAAM_2_zzaaam(isoflag);
                if (0 < ms_feed.comp.count(isoInLine) )
                    isobuf << "\t" << ms_feed.comp[isoInLine];
            
                else 
                    isobuf << "\t" << 0.0;
    
                if (0 < ms_prod.comp.count(isoInLine) )
                    isobuf << "\t" << ms_prod.comp[isoInLine] * outmassfrac;
                else 
                    isobuf << "\t" << 0.0;
            }
            catch (std::exception& e)
            {
                continue;
            }
        }
        isobuf << "\n";
    }
    isofilein.close();
    std::ofstream isofileout ( (name + "Isos.txt").c_str() );
    isofileout <<  isobuf.rdbuf();
    isofileout.close();
}

void FCComp::write_params_pass ()
{
    //Writes a single pass to the parameter tracking file.
    std::ifstream paramfilein  ( (name + "Params.txt").c_str() );
    std::stringstream parambuf;
    parambuf.precision(6);
    parambuf << std::scientific << std::uppercase;

    while (!paramfilein.eof() )
    {
        char line [3000];
        paramfilein.getline(line, 3000);
        parambuf << line;

        std::string paramflag = bright::getFlag(line, 10);
        if (paramflag == "Param")
            parambuf << "\t" << bright::to_str(pass_num) << "in\t\t" << bright::to_str(pass_num) << "out\t";
        else if (0 < params_prior_calc.count(paramflag))
            parambuf << "\t" << params_prior_calc[paramflag] << "\t" << params_after_calc[paramflag];
        parambuf << "\n";
    }
    paramfilein.close();
    std::ofstream paramfileout ( (name + "Params.txt").c_str() );
    paramfileout <<  parambuf.rdbuf();
    paramfileout.close();
}

void FCComp::write_text()
{
    //Write the isotopic streams
    write_ms_pass();

    //Write the parameters if they are there to write!
    if (!track_params.empty()) 
        write_params_pass();
}

void FCComp::appendHDF5array(H5::H5File *dbFile, std::string set_name, double *append_value, const int *rank, \
hsize_t dims[], hsize_t offset[], hsize_t extend_size[])
{
    H5::DataSet array_set = (*dbFile).openDataSet(set_name);
    array_set.extend(extend_size);

    H5::DataSpace append_space = array_set.getSpace();
    append_space.selectHyperslab(H5S_SELECT_SET, dims, offset);

    H5::DataSpace value_space(*rank, dims);

    double data_value[1] = { *append_value };
    array_set.write(data_value, H5::PredType::NATIVE_DOUBLE, value_space, append_space);

    array_set.close();
}

void FCComp::write_hdf5 ()
{
    //Writes the fuel cycle component to an HDF5 file
    const int    RANK   = 1;
    hsize_t dims[1]     = {1};
    hsize_t offset[1]   = {pass_num - 1};
    hsize_t ext_size[1] = {pass_num};
        
    //Open the HDF5 file
    H5::H5File dbFile (bright::output_filename, H5F_ACC_RDWR);
    std::string comp_path ("/" + natural_name);

    //Write the isotopic component input and output streams
    if (!bright::track_isos.empty())
    {
        appendHDF5array(&dbFile, comp_path + "/ms_feed/Mass",  &(ms_feed.mass),  &RANK, dims, offset, ext_size);
        appendHDF5array(&dbFile, comp_path + "/ms_prod/Mass", &(ms_prod.mass), &RANK, dims, offset, ext_size);

        for (std::set<int>::iterator iso = bright::track_isos.begin(); iso != bright::track_isos.end(); iso++)
        {
            std::string isoLL = pyne::nucname::zzaaam_2_LLAAAM(*iso);
            appendHDF5array(&dbFile, comp_path + "/ms_feed/"  + isoLL, &(ms_feed.comp[*iso]),  &RANK, dims, offset, ext_size);
            appendHDF5array(&dbFile, comp_path + "/ms_prod/" + isoLL, &(ms_prod.comp[*iso]), &RANK, dims, offset, ext_size);
        }
    }    

    //Write the parameter tracking
    if (!track_params.empty())
    {
        for ( std::set<std::string>::iterator p = track_params.begin(); p != track_params.end(); p++)
        {
            appendHDF5array(&dbFile, comp_path + "/params_prior_calc/"  + (*p), &(params_prior_calc[*p]),  &RANK, dims, offset, ext_size);
            appendHDF5array(&dbFile, comp_path + "/params_after_calc/" + (*p), &(params_after_calc[*p]), &RANK, dims, offset, ext_size);
        }
    }

    //close the HDF5 File
    dbFile.close();   
}

void FCComp::write ()
{
    //Now that we are ready to start writing out data, let's update the pass number that we are on.
    pass_num++;

    //Set the parameters for this pass.
    if (!track_params.empty())
        calc_params();

    //Writes the output table files.
    if (bright::write_text)
        write_text();
    
    if (bright::write_hdf5)
        write_hdf5();
}

MassStream FCComp::calc ()
{
    //Placehodler function for the calculation of all relevant isotopes and parameters.
    //Returns an empty MassStream object.
    return MassStream ();
}

MassStream FCComp::calc (CompDict cd)
{
    //Placehodler function for the calculation of all relevant isotopes and parameters.
    //Returns an empty MassStream object.
    return MassStream ();
}

MassStream FCComp::calc (MassStream ms)
{
    //Placehodler function for the calculation of all relevant isotopes and parameters.
    //Returns an empty MassStream object.
    return MassStream ();
}
