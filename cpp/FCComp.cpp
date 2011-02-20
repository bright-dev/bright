// Fuel Cycle Component Parent Class

#include "FCComp.h"

/*********************************************/
/*** Fuel Cycle Components Namespace Stuff ***/
/*********************************************/

#ifdef _WIN32
    int null_set [1] = {922350};
    std::set<int> FCComps::isos2track (null_set, null_set+1);
#else
    int null_set [0] = {};
    std::set<int> FCComps::isos2track (null_set, null_set+0);
#endif

int FCComps::verbosity  = 0;
int FCComps::write_text = 1;
int FCComps::write_hdf5 = 0;

std::string FCComps::output_filename = "fuel_cycle.h5";

void FCComps::load_isos2track_hdf5(std::string filename, std::string datasetname, bool clear_prev)
{
    // Check that the file is there
    if (!bright::FileExists(filename))
        throw bright::FileNotFound(filename);

    //Load values into isos2track from an hdf5 file.
    //If the dataspace name is not given, try some defaults.
    int dslen = 14;
    std::string defaultsets [14] = {
        "/isos2track",
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
            throw H5::FileIException("load_isos2track", "Dataset not found!");
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
        isos2track.clear();

    //load into isos2track
    for(int n = 0; n < isolen[0]; n++)
    {
        isos2track.insert(isoname::mixed_2_zzaaam(iso_out_int[n]));
    };
    
};

void FCComps::load_isos2track_text(std::string filename, bool clear_prev)
{
    // Check that the file is there
    if (!bright::FileExists(filename))
        throw bright::FileNotFound(filename);

    //Clear previous entries
    if (clear_prev)
        isos2track.clear();

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
        isos2track.insert(isoname::mixed_2_zzaaam(isostr));
    };

    //close file
    isofile.close();
};

/**************************************************/
/*** Fuel Cycle Component Classes And Functions ***/
/**************************************************/

/***************************/
/*** Protected Functions ***/
/***************************/

void FCComp::initialize (std::set<std::string> ptrack, std::string n)
{
    //Protected Variables
    params2track = ptrack;

    //Public Variables
    name = n;
    natural_name = bright::natural_naming(n);
    if (natural_name.length() == 0)
        natural_name = "this_is_not_a_name";
    
    PassNum = 0;

    if (FCComps::write_text)
        initialize_Text();

    if (FCComps::write_hdf5)
        initialize_HDF5();
}
    
void FCComp::initialize_Text ()
{
    //Initialize the Isotopic tracking file
    if (!FCComps::isos2track.empty())
    {
        std::ofstream isofile ( (name + "Isos.txt").c_str() );
        isofile << "Isotope\n";
        for (std::set<int>::iterator iso = FCComps::isos2track.begin(); iso != FCComps::isos2track.end(); iso++)
        {
            isofile << isoname::zzaaam_2_LLAAAM(*iso) << "\n"; 
        }
        isofile.close();
    }

    //Initialize the Parameter tracking file.
    if (!params2track.empty())
    {	
        std::ofstream paramfile ( (name + "Params.txt").c_str() );
        paramfile << "Param\n";
        for ( std::set<std::string>::iterator p = params2track.begin(); p != params2track.end(); p++)
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
    if (bright::FileExists(FCComps::output_filename))
        dbFile = H5::H5File(FCComps::output_filename, H5F_ACC_RDWR);
    else
        dbFile = H5::H5File(FCComps::output_filename, H5F_ACC_TRUNC);

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
    if (!FCComps::isos2track.empty())
    {
        // Open/Create IsosIn group
        H5::Group gIsosIn;
        try
          { gIsosIn = dbFile.openGroup(comp_path + "/IsosIn"); }
        catch (H5::Exception fgerror)
          { gIsosIn = dbFile.createGroup(comp_path + "/IsosIn"); }

        // Open/Create IsosOut group
        H5::Group gIsosOut;
        try
            { gIsosOut = dbFile.openGroup(comp_path + "/IsosOut"); }
        catch (H5::Exception fgerror)
            { gIsosOut = dbFile.createGroup(comp_path + "/IsosOut"); }

        // Open/Create /IsosIn/Mass Dataset
        H5::DataSet dsIsosInMass;
        try 
            { dsIsosInMass = dbFile.openDataSet(comp_path + "/IsosIn/Mass"); }
        catch (H5::Exception fgerror)
            { dsIsosInMass = dbFile.createDataSet(comp_path + "/IsosIn/Mass",  H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }

        // Open/Create /IsosOut/Mass Dataset
        H5::DataSet dsIsosOutMass;
        try
            { dsIsosOutMass = dbFile.openDataSet(comp_path + "/IsosOut/Mass"); }
        catch (H5::Exception fgerror)
            { dsIsosOutMass = dbFile.createDataSet(comp_path + "/IsosOut/Mass", H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }

        // Open/Create /Isos[In|Out]/iso Datasets
        H5::DataSet dsIsosInIso;
        H5::DataSet dsIsosOutIso;
        for (std::set<int>::iterator iso = FCComps::isos2track.begin(); iso != FCComps::isos2track.end(); iso++)
        {
            std::string isoLL = isoname::zzaaam_2_LLAAAM(*iso);

            try
                { dsIsosInIso = dbFile.openDataSet(comp_path + "/IsosIn/" + isoLL); }
            catch (H5::Exception fgerror)
                { dsIsosInIso = dbFile.createDataSet(comp_path + "/IsosIn/" + isoLL, H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }

            try
                { dsIsosOutIso = dbFile.openDataSet(comp_path + "/IsosOut/" + isoLL); }
            catch (H5::Exception fgerror)
                { dsIsosOutIso = dbFile.createDataSet(comp_path + "/IsosOut/" + isoLL, H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }
        }
    }

    //Initiallize the Parameters
    if (!params2track.empty())
    {	
        // Open/Create ParamsIn group
        H5::Group gParamsIn;
        try
          { gParamsIn = dbFile.openGroup(comp_path + "/ParamsIn"); }
        catch (H5::Exception fgerror)
          { gParamsIn = dbFile.createGroup(comp_path + "/ParamsIn"); }

        // Open/Create ParamsOut group
        H5::Group gParamsOut;
        try
            { gParamsOut = dbFile.openGroup(comp_path + "/ParamsOut"); }
        catch (H5::Exception fgerror)
            { gParamsOut = dbFile.createGroup(comp_path + "/ParamsOut"); }

        // Open/Create /Params[In|Out]/param Datasets
        H5::DataSet dsParamsInParam;
        H5::DataSet dsParamsOutParam;
        for ( std::set<std::string>::iterator p = params2track.begin(); p != params2track.end(); p++)
        {
            try
                { dsParamsInParam = dbFile.openDataSet(comp_path + "/ParamsIn/"  + *p); }
            catch (H5::Exception fgerror)
                { dsParamsInParam = dbFile.createDataSet(comp_path + "/ParamsIn/"  + *p, H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }

            try
                { dsParamsOutParam = dbFile.openDataSet(comp_path + "/ParamsOut/" + *p); }
            catch (H5::Exception fgerror)
                { dsParamsOutParam = dbFile.createDataSet(comp_path + "/ParamsOut/" + *p, H5::PredType::NATIVE_DOUBLE, ext_1D_space, double_params); }
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

void FCComp::setParams ()
{
    //Placeholder function that sets the states of ParamsIn and ParamsOut.
    for ( std::set<std::string>::iterator p2t = params2track.begin(); p2t != params2track.end(); p2t++)
    {
        ParamsIn[*p2t]  = 0.0;
        ParamsOut[*p2t] = 0.0;
    }
}

void FCComp::writeIsoPass ()
{
    //Writes a single pass to the isotopic tracking file.
    std::ifstream isofilein  ( (name + "Isos.txt").c_str() );
    std::stringstream isobuf;
    isobuf.precision(6);
    isobuf << std::scientific << std::uppercase;

    //normalize IsosOut per kgIsosIn
    double outmassfrac = IsosOut.mass / IsosIn.mass;

    while (!isofilein.eof() )
    {
        char line [3000];
        isofilein.getline(line, 3000);
        isobuf << line;

        std::string isoflag = bright::getFlag(line, 10);
        if (isoflag == "Isotope")
                isobuf << "\t" << bright::to_str(PassNum) << "in\t\t" << bright::to_str(PassNum) << "out\t";
        else
        {
            try
            {
                int isoInLine = isoname::LLAAAM_2_zzaaam(isoflag);
                if (0 < IsosIn.comp.count(isoInLine) )
                    isobuf << "\t" << IsosIn.comp[isoInLine];
            
                else 
                    isobuf << "\t" << 0.0;
    
                if (0 < IsosOut.comp.count(isoInLine) )
                    isobuf << "\t" << IsosOut.comp[isoInLine] * outmassfrac;
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

void FCComp::writeParamPass ()
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
            parambuf << "\t" << bright::to_str(PassNum) << "in\t\t" << bright::to_str(PassNum) << "out\t";
        else if (0 < ParamsIn.count(paramflag))
            parambuf << "\t" << ParamsIn[paramflag] << "\t" << ParamsOut[paramflag];
        parambuf << "\n";
    }
    paramfilein.close();
    std::ofstream paramfileout ( (name + "Params.txt").c_str() );
    paramfileout <<  parambuf.rdbuf();
    paramfileout.close();
}

void FCComp::writeText()
{
    //Write the isotopic streams
    writeIsoPass();

    //Write the parameters if they are there to write!
    if (!params2track.empty()) 
        writeParamPass();
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

void FCComp::writeHDF5 ()
{
    //Writes the fuel cycle component to an HDF5 file
    const int    RANK   = 1;
    hsize_t dims[1]     = {1};
    hsize_t offset[1]   = {PassNum - 1};
    hsize_t ext_size[1] = {PassNum};
        
    //Open the HDF5 file
    H5::H5File dbFile (FCComps::output_filename, H5F_ACC_RDWR);
    std::string comp_path ("/" + natural_name);

    //Write the isotopic component input and output streams
    if (!FCComps::isos2track.empty())
    {
        appendHDF5array(&dbFile, comp_path + "/IsosIn/Mass",  &(IsosIn.mass),  &RANK, dims, offset, ext_size);
        appendHDF5array(&dbFile, comp_path + "/IsosOut/Mass", &(IsosOut.mass), &RANK, dims, offset, ext_size);

        for (std::set<int>::iterator iso = FCComps::isos2track.begin(); iso != FCComps::isos2track.end(); iso++)
        {
            std::string isoLL = isoname::zzaaam_2_LLAAAM(*iso);
            appendHDF5array(&dbFile, comp_path + "/IsosIn/"  + isoLL, &(IsosIn.comp[*iso]),  &RANK, dims, offset, ext_size);
            appendHDF5array(&dbFile, comp_path + "/IsosOut/" + isoLL, &(IsosOut.comp[*iso]), &RANK, dims, offset, ext_size);
        }
    }    

    //Write the parameter tracking
    if (!params2track.empty())
    {
        for ( std::set<std::string>::iterator p = params2track.begin(); p != params2track.end(); p++)
        {
            appendHDF5array(&dbFile, comp_path + "/ParamsIn/"  + (*p), &(ParamsIn[*p]),  &RANK, dims, offset, ext_size);
            appendHDF5array(&dbFile, comp_path + "/ParamsOut/" + (*p), &(ParamsOut[*p]), &RANK, dims, offset, ext_size);
        }
    }

    //close the HDF5 File
    dbFile.close();   
}

void FCComp::writeout ()
{
    //Now that we are ready to start writing out data, let's update the pass number that we are on.
    PassNum++;

    //Set the parameters for this pass.
    if (!params2track.empty())
        setParams();

    //Writes the output table files.
    if (FCComps::write_text)
        writeText();
    
    if (FCComps::write_hdf5)
        writeHDF5();
}

MassStream FCComp::doCalc ()
{
    //Placehodler function for the calculation of all relevant isotopes and parameters.
    //Returns an empty MassStream object.
    return MassStream ();
}

MassStream FCComp::doCalc (CompDict cd)
{
    //Placehodler function for the calculation of all relevant isotopes and parameters.
    //Returns an empty MassStream object.
    return MassStream ();
}

MassStream FCComp::doCalc (MassStream ms)
{
    //Placehodler function for the calculation of all relevant isotopes and parameters.
    //Returns an empty MassStream object.
    return MassStream ();
}
