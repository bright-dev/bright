//The Bright Python wrapper for "MassStream"

//standard library headers

//Needed Boost Modules
#include <boost/python.hpp>
#include <boost/python/overloads.hpp>

//My Files
#include "BriPyHelper.h"

#ifdef _WIN32
    #include "../FCComps/FCComp.h"
    #include "../FCComps/Reprocess.h"
    #include "../FCComps/Storage.h"
    #include "../FCComps/Enrichment.h"
    #include "../FCComps/Reactor1G.h"
    #include "../FCComps/FastReactor1G.h"
    #include "../FCComps/LightWaterReactor1G.h"
    #include "../FCComps/FuelFabrication.h"
#else
    #include "src/FCComp.h"
    #include "src/Reprocess.h"
    #include "src/Storage.h"
    #include "src/Enrichment.h"
    #include "src/Reactor1G.h"
    #include "src/FastReactor1G.h"
    #include "src/LightWaterReactor1G.h"
    #include "src/FuelFabrication.h"
#endif

//Need to keep Boost::Python in its own namespace.
//	NO: using namespace boost::python;
//possible because both Boost::Python and HDF5 define
//a 'ssize_t' variable that colides in an abiguous 
//definition at global namespce.  Unfortunately, HDF5
//does not encapsulate itself in a namespce.
//Thus to use these two codes together we have to use
//Boost::Python's namespace explicitly.
//Renamed here for ease of use.
namespace bp = boost::python;

//Have to define FCComps Globals to allow Python to import them.
//int null_set [0] = {};
//std::set<int> FCComps::isos2track (null_set, null_set+0);
//int FCComps::verbosity = 0;


//Exception Hack
//static PyObject * Py_BadFuelForm = PyErr_NewException((char *) "FCComps.BadFuelForm", 0, 0); 
//void BFF_Trans(BadFuelForm const& e)
//{
//    PyErr_SetString(Py_BadFuelForm, e.what());
//};

//static PyObject * Py_BisectionMethodNotPerformed = PyErr_NewException((char *) "FCComps.BisectionMethodNotPerformed", 0, 0); 
//void BMNP_Trans(BisectionMethodNotPerformed const& e)
//{
//    PyErr_SetString(Py_BisectionMethodNotPerformed, e.what());
//};


//Default Argument Function Overload Macros
BOOST_PYTHON_FUNCTION_OVERLOADS(load_isos2track_hdf5_overloads, FCComps::load_isos2track_hdf5, 1, 3)
BOOST_PYTHON_FUNCTION_OVERLOADS(load_isos2track_text_overloads, FCComps::load_isos2track_text, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(batchAve_overloads, batchAve, 1, 2)

BOOST_PYTHON_MODULE(FCComps)
{
    //sets the current scope...
    bp::scope current;
    current.attr("__doc__") = "Fuel Cycle Component Classes & Module Wrapper";

    bp::def("BrightStart", &bright::BrightStart, "Initializes Bright");

    char isos2track__doc__ []                   = "Gets/Sets the iso2track set.";
    std::set<int> (*isos2track1)()              = &FCComps::get_isos2track;
    void          (*isos2track2)(std::set<int>) = &FCComps::set_isos2track;
    bp::def("isos2track", isos2track1, isos2track__doc__);
    bp::def("isos2track", isos2track2, isos2track__doc__);

    bp::def("load_isos2track_hdf5", &FCComps::load_isos2track_hdf5, load_isos2track_hdf5_overloads("Loads isos2track from an HDF5 file."));
    bp::def("load_isos2track_text", &FCComps::load_isos2track_text, load_isos2track_text_overloads("Loads isos2track from a text file."));

    char verbosity__doc__ [] = "Gets/Sets the verbosity level.";
    int  (*verbosity1)()     = &FCComps::get_verbosity;
    void (*verbosity2)(int)  = &FCComps::set_verbosity;
    bp::def("verbosity", verbosity1, verbosity__doc__);
    bp::def("verbosity", verbosity2, verbosity__doc__);

    char write_text__doc__ [] = "Gets/Sets the write-to-text output files flag.";
    int  (*write_text1)()     = &FCComps::get_write_text;
    void (*write_text2)(int)  = &FCComps::set_write_text;
    bp::def("write_text", write_text1, write_text__doc__);
    bp::def("write_text", write_text2, write_text__doc__);

    char write_hdf5__doc__ [] = "Gets/Sets the write-to-HDF5 output files flag.";
    int  (*write_hdf51)()     = &FCComps::get_write_hdf5;
    void (*write_hdf52)(int)  = &FCComps::set_write_hdf5;
    bp::def("write_hdf5", write_hdf51, write_hdf5__doc__);
    bp::def("write_hdf5", write_hdf52, write_hdf5__doc__);

    char output_filename__doc__ []               = "Gets/Sets the fuel cycle output file name for HDF5.";
    std::string (*output_filename1)()            = &FCComps::get_output_filename;
    void        (*output_filename2)(std::string) = &FCComps::set_output_filename;
    bp::def("output_filename", output_filename1, output_filename__doc__);
    bp::def("output_filename", output_filename2, output_filename__doc__);

    //Basic to- and from-converters
    dict2map<Iso, Weight>();
    bp::to_python_converter< CompDict, map2dict<Iso, Weight> >();

    //FCComp to- and from-converters
    dict2map<std::string, double>();
    bp::to_python_converter< ParamDict, map2dict<std::string, double> >();

    py_list2c_set<int>();
    bp::to_python_converter< std::set<int> , c_set2py_list<int> >();

    py_list2c_set<std::string>();
    bp::to_python_converter< std::set<std::string> , c_set2py_list<std::string> >();

    //FCComp Method Function Overloads
    MassStream (FCComp::*FCComp_doCalc_NA)()           = &FCComp::doCalc;
    MassStream (FCComp::*FCComp_doCalc_CD)(CompDict)   = &FCComp::doCalc;
    MassStream (FCComp::*FCComp_doCalc_MS)(MassStream) = &FCComp::doCalc;
        
    bp::class_< FCComp >("FCComp", "Fuel Cycle Component Base Class",  bp::init<>("Empty Fuel Cycle Component.") )
        //Basic FCComp Properties
        .add_property("params2track", &FCComp::get_params2track, &FCComp::set_params2track)

        .add_property("name",         &FCComp::get_name,         &FCComp::set_name)
        .add_property("natural_name", &FCComp::get_natural_name, &FCComp::set_natural_name)
        .add_property("IsosIn",       &FCComp::get_IsosIn,       &FCComp::set_IsosIn)
        .add_property("IsosOut",      &FCComp::get_IsosOut,      &FCComp::set_IsosOut)
        .add_property("ParamsIn",     &FCComp::get_ParamsIn,     &FCComp::set_ParamsIn)
        .add_property("ParamsOut",    &FCComp::get_ParamsOut,    &FCComp::set_ParamsOut)
        .add_property("PassNum",      &FCComp::get_PassNum,      &FCComp::set_PassNum)

        //FCComp Constructor Overloads
        .def(bp::init< bp::optional<std::string> >("Fuel Cycle Component tracking isotopes."))
        .def(bp::init< std::set<std::string>, bp::optional<std::string> >("Fuel Cycle Component tracking isotopes and parameters."))

        //Useful Functions
        .def("setParams",      &FCComp::setParams)
        .def("writeIsoPass",   &FCComp::writeIsoPass)
        .def("writeParamPass", &FCComp::writeParamPass)
        .def("writeHDF5",      &FCComp::writeHDF5)
        .def("writeText",      &FCComp::writeText)
        .def("writeout",       &FCComp::writeout)
        .def("doCalc", FCComp_doCalc_NA)
        .def("doCalc", FCComp_doCalc_CD)
        .def("doCalc", FCComp_doCalc_MS)
    ;

    //Reprocess to- and from-converters
    dict2map<int, double>();
    bp::to_python_converter< SepEffDict, map2dict<int, double> >();

    dict2map<std::string, double>();
    bp::to_python_converter< std::map<std::string, double>, map2dict<std::string, double> >();

    //Reprocess Method Function Overloads
    MassStream (Reprocess::*Reprocess_doCalc_NA)()           = &Reprocess::doCalc;
    MassStream (Reprocess::*Reprocess_doCalc_CD)(CompDict)   = &Reprocess::doCalc;
    MassStream (Reprocess::*Reprocess_doCalc_MS)(MassStream) = &Reprocess::doCalc;
        
    bp::class_< Reprocess, bp::bases<FCComp> >("Reprocess", "Reprocessing Facility Fuel Cycle Component",  bp::init<>() )
        //Basic Reprocess Properties
        .add_property("sepeff", &Reprocess::get_sepeff, &Reprocess::set_sepeff)

        //Reprocess Constructor Overloads
        .def(bp::init< std::map<std::string, double>, bp::optional<std::string> >())
        //NOTE! Due to from_pyton C++ signature ambiguity, you cannot overload both non-trival constructors.
        //Separation efficiencies must therefore be automatically initilized through string ditionaries.
        //Since SepEffDict is an int dictionary DO NOT add a line like the following...
        //	.def(init< SepEffDict, std::set<int>, optional<std::string> >())
        //If you need to initizilze via an int dictionary in python, you can always init with an empty
        //string dictionary and then manually initilaize with an int one.  For example,
        //	R = FCComps.Reprocess({}, isotrack, name)
        //	R.initialize( {92: 0.99, 942390: 0.9, ...} )

        //Useful Functions
        .def("initialize", &Reprocess::initialize)
        .def("setParams",  &Reprocess::setParams)
        .def("doCalc", Reprocess_doCalc_NA)
        .def("doCalc", Reprocess_doCalc_CD)
        .def("doCalc", Reprocess_doCalc_MS)
    ;

    //Storage to- and from-converters

    //Storage Method Function Overloads
    MassStream (Storage::*Storage_doCalc_NA)()                      = &Storage::doCalc;
    MassStream (Storage::*Storage_doCalc_CD)(CompDict)              = &Storage::doCalc;
    MassStream (Storage::*Storage_doCalc_MS)(MassStream)            = &Storage::doCalc;
    MassStream (Storage::*Storage_doCalc_DB)(double)                = &Storage::doCalc;
    MassStream (Storage::*Storage_doCalc_CD_DB)(CompDict, double)   = &Storage::doCalc;
    MassStream (Storage::*Storage_doCalc_MS_DB)(MassStream, double) = &Storage::doCalc;
        
    bp::class_< Storage, bp::bases<FCComp> >("Storage", "Storage/Decay Facility Fuel Cycle Component",  bp::init<>() )
        //Basic Storage Properties
        .add_property("decay_time", &Storage::get_decay_time, &Storage::set_decay_time)

        //Storage Constructor Overloads
        .def(bp::init< bp::optional<std::string> >())

        //Useful Functions
        .def("setParams", &Storage::setParams)
        .def("doCalc", Storage_doCalc_NA)
        .def("doCalc", Storage_doCalc_CD)
        .def("doCalc", Storage_doCalc_MS)
        .def("doCalc", Storage_doCalc_DB)
        .def("doCalc", Storage_doCalc_CD_DB)
        .def("doCalc", Storage_doCalc_MS_DB)
    ;



    /******************/
    /*** Enrichment ***/
    /******************/

    //EnrichmentParameters Class...Used in Enrichment
    bp::class_< EnrichmentParameters >("EnrichmentParameters", "Set of physical parameters that define an enrichment cascade. This is used to instantiate new enrichment objects.", bp::init<>() )
        .add_property("alpha_0", &EnrichmentParameters::get_alpha_0, &EnrichmentParameters::set_alpha_0)
        .add_property("Mstar_0", &EnrichmentParameters::get_Mstar_0, &EnrichmentParameters::set_Mstar_0)
        .add_property("j",       &EnrichmentParameters::get_j,       &EnrichmentParameters::set_j)
        .add_property("k",       &EnrichmentParameters::get_k,       &EnrichmentParameters::set_k)
        .add_property("N0",      &EnrichmentParameters::get_N0,      &EnrichmentParameters::set_N0)
        .add_property("M0",      &EnrichmentParameters::get_M0,      &EnrichmentParameters::set_M0)
        .add_property("xP_j",    &EnrichmentParameters::get_xP_j,    &EnrichmentParameters::set_xP_j)
        .add_property("xW_j",    &EnrichmentParameters::get_xW_j,    &EnrichmentParameters::set_xW_j)
    ;

    //Grabs a copy of the urnaium enrichment default Settings
    bp::def("UraniumEnrichmentDefaults", &fillUraniumEnrichmentDefaults);

    //Enrichment Method Function Overloads
    MassStream (Enrichment::*Enrichment_doCalc_NA)()           = &Enrichment::doCalc;
    MassStream (Enrichment::*Enrichment_doCalc_CD)(CompDict)   = &Enrichment::doCalc;
    MassStream (Enrichment::*Enrichment_doCalc_MS)(MassStream) = &Enrichment::doCalc;
        
    bp::class_< Enrichment, bp::bases<FCComp> >("Enrichment", "Enrichmenting Facility Fuel Cycle Component",  bp::init< >() )
        //Basic Enrichment Properties
        .add_property("alpha_0", &Enrichment::get_alpha_0, &Enrichment::set_alpha_0)
        .add_property("Mstar_0", &Enrichment::get_Mstar_0, &Enrichment::set_Mstar_0)
        .add_property("j",       &Enrichment::get_j,       &Enrichment::set_j)
        .add_property("k",       &Enrichment::get_k,       &Enrichment::set_k)
        .add_property("N0",      &Enrichment::get_N0,      &Enrichment::set_N0)
        .add_property("M0",      &Enrichment::get_M0,      &Enrichment::set_M0)
        .add_property("xP_j",    &Enrichment::get_xP_j,    &Enrichment::set_xP_j)
        .add_property("xW_j",    &Enrichment::get_xW_j,    &Enrichment::set_xW_j)

        .add_property("Mstar",    &Enrichment::get_Mstar,    &Enrichment::set_Mstar)
        .add_property("IsosTail", &Enrichment::get_IsosTail, &Enrichment::set_IsosTail)
        .add_property("N",        &Enrichment::get_N,        &Enrichment::set_N)
        .add_property("M",        &Enrichment::get_M,        &Enrichment::set_M)

        .add_property("TotalPerFeed",  &Enrichment::get_TotalPerFeed,  &Enrichment::set_TotalPerFeed)
        .add_property("SWUperFeed",    &Enrichment::get_SWUperFeed,    &Enrichment::set_SWUperFeed)
        .add_property("SWUperProduct", &Enrichment::get_SWUperProduct, &Enrichment::set_SWUperProduct)

        //Enrichment Constructor Overloads
        .def(bp::init<std::string>())
        .def(bp::init<EnrichmentParameters,  bp::optional<std::string> >())

        //Useful Functions
        .def("initialize", &Enrichment::initialize)
        .def("setParams",  &Enrichment::setParams)
        .def("doCalc", Enrichment_doCalc_NA)
        .def("doCalc", Enrichment_doCalc_CD)
        .def("doCalc", Enrichment_doCalc_MS)
    ;




    /*****************/
    /*** Reactor1G ***/
    /*****************/

    //FleneceIndex Class...Used in Reactor1G
    bp::class_< FluencePoint >("FluencePoint", "Structure for Fluence Point",  bp::init<>() )
        .def_readwrite("f", &FluencePoint::f)
        .def_readwrite("F", &FluencePoint::F)
        .def_readwrite("m", &FluencePoint::m)
    ;

    //ReactorParameters Class...Used in Reactor1G
    bp::class_< ReactorParameters >("ReactorParameters", "Set of physical reactor parameters. May be used to instantiate new reactor objects, -or- to define default settings for a reactor type.", bp::init<>() )
        .add_property("batches",         &ReactorParameters::get_batches,         &ReactorParameters::set_batches)
        .add_property("flux",            &ReactorParameters::get_flux,            &ReactorParameters::set_flux)
        .add_property("FuelForm",        &ReactorParameters::get_FuelForm,        &ReactorParameters::set_FuelForm)
        .add_property("CoolantForm",     &ReactorParameters::get_CoolantForm,     &ReactorParameters::set_CoolantForm)
        .add_property("FuelDensity",     &ReactorParameters::get_FuelDensity,     &ReactorParameters::set_FuelDensity)
        .add_property("CoolantDensity",  &ReactorParameters::get_CoolantDensity,  &ReactorParameters::set_CoolantDensity)
        .add_property("pnl",             &ReactorParameters::get_pnl,             &ReactorParameters::set_pnl)
        .add_property("BUt",             &ReactorParameters::get_BUt,             &ReactorParameters::set_BUt)
        .add_property("useDisadvantage", &ReactorParameters::get_useDisadvantage, &ReactorParameters::set_useDisadvantage)
        .add_property("LatticeType",     &ReactorParameters::get_LatticeType,     &ReactorParameters::set_LatticeType)
        .add_property("HydrogenRescale", &ReactorParameters::get_HydrogenRescale, &ReactorParameters::set_HydrogenRescale)

        .add_property("Radius",          &ReactorParameters::get_Radius,          &ReactorParameters::set_Radius)
        .add_property("Length",          &ReactorParameters::get_Length,          &ReactorParameters::set_Length)
        .add_property("open_slots",      &ReactorParameters::get_open_slots,      &ReactorParameters::set_open_slots)
        .add_property("total_slots",     &ReactorParameters::get_total_slots,     &ReactorParameters::set_total_slots)
    ;

    //Reactor1G Errors
    //bp::register_exception_translator<BadFuelForm>( BFF_Trans );
    //bp::register_exception_translator<BisectionMethodNotPerformed>( BMNP_Trans );
    //Add to-python converter for exceptions to get the next line to work and expose them!
    //current.attr("BadFuelForm") = &Py_BadFuelForm;
    //current.attr("BisectionMethodNotPerformed") = &Py_BisectionMethodNotPerformed;
    bp::register_exception_translator<BadFuelForm>( BriPy_Exception_Translator<BadFuelForm> );
    bp::register_exception_translator<BisectionMethodNotPerformed>( BriPy_Exception_Translator<BisectionMethodNotPerformed> );

    //Reactor1G to- and from-converters
    bp::to_python_converter< std::vector<float> , c_vector2py_list<float> >();
    bp::to_python_converter< std::map<int, std::vector<float> >, map2dict<int, std::vector<float> > >();
    bp::to_python_converter< std::map<int, std::map<int, std::vector<float> > >, map2dict<int, std::map<int, std::vector<float> > > >();

    bp::to_python_converter< std::vector<double> , c_vector2py_list<double> >();
    bp::to_python_converter< std::map<int, std::vector<double> >, map2dict<int, std::vector<double> > >();
    bp::to_python_converter< std::map<int, std::map<int, std::vector<double> > >, map2dict<int, std::map<int, std::vector<double> > > >();

    //Reactor1G Method Function Overloads
    MassStream (Reactor1G::*Reactor1G_doCalc_NA)()           = &Reactor1G::doCalc;
    MassStream (Reactor1G::*Reactor1G_doCalc_CD)(CompDict)   = &Reactor1G::doCalc;
    MassStream (Reactor1G::*Reactor1G_doCalc_MS)(MassStream) = &Reactor1G::doCalc;
        
    bp::class_< Reactor1G, bp::bases<FCComp> >("Reactor1G", "One-Group Nuclear Reactor Fuel Cycle Component",  bp::init<>() )
        //Basic MassStream Properties
        .add_property("B",                   &Reactor1G::get_B,                   &Reactor1G::set_B)
        .add_property("phi",                 &Reactor1G::get_phi,                 &Reactor1G::set_phi)
        .add_property("FuelChemicalForm",    &Reactor1G::get_FuelChemicalForm,    &Reactor1G::set_FuelChemicalForm)
        .add_property("CoolantChemicalForm", &Reactor1G::get_CoolantChemicalForm, &Reactor1G::set_CoolantChemicalForm)
        .add_property("rhoF",                &Reactor1G::get_rhoF,                &Reactor1G::set_rhoF)
        .add_property("rhoC",                &Reactor1G::get_rhoC,                &Reactor1G::set_rhoC)
        .add_property("P_NL",                &Reactor1G::get_P_NL,                &Reactor1G::set_P_NL)
        .add_property("TargetBU",            &Reactor1G::get_TargetBU,            &Reactor1G::set_TargetBU)
        .add_property("useZeta",             &Reactor1G::get_useZeta,             &Reactor1G::set_useZeta)
        .add_property("Lattice",             &Reactor1G::get_Lattice,             &Reactor1G::set_Lattice)
        .add_property("H_XS_Rescale",        &Reactor1G::get_H_XS_Rescale,        &Reactor1G::set_H_XS_Rescale)

        .add_property("r",                   &Reactor1G::get_r,                   &Reactor1G::set_r)
        .add_property("l",                   &Reactor1G::get_l,                   &Reactor1G::set_l)
        .add_property("S_O",                 &Reactor1G::get_S_O,                 &Reactor1G::set_S_O)
        .add_property("S_T",                 &Reactor1G::get_S_T,                 &Reactor1G::set_S_T)
        .add_property("VF",                  &Reactor1G::get_VF,                  &Reactor1G::set_VF)
        .add_property("VC",                  &Reactor1G::get_VC,                  &Reactor1G::set_VC)

        .add_property("libfile",             &Reactor1G::get_libfile,             &Reactor1G::set_libfile)
        //The following attributes are read only!
        .add_property("F",                   &Reactor1G::get_F)
        .add_property("BUi_F_",              &Reactor1G::get_BUi_F_)
        .add_property("pi_F_",               &Reactor1G::get_pi_F_)
        .add_property("di_F_",               &Reactor1G::get_di_F_)
        .add_property("Tij_F_",              &Reactor1G::get_Tij_F_)

        //The following attributes are read only!
        .add_property("I", &Reactor1G::get_I)
        .add_property("J", &Reactor1G::get_J)
        .add_property("sigma_a_therm", &Reactor1G::get_sigma_a_therm)
        .add_property("sigma_s_therm", &Reactor1G::get_sigma_s_therm)


        .add_property("A_IHM",               &Reactor1G::get_A_IHM,               &Reactor1G::set_A_IHM)
        .add_property("MWF",                 &Reactor1G::get_MWF,                 &Reactor1G::set_MWF)
        .add_property("MWC",                 &Reactor1G::get_MWC,                 &Reactor1G::set_MWC)
        .add_property("niF",                 &Reactor1G::get_niF,                 &Reactor1G::set_niF)
        .add_property("niC",                 &Reactor1G::get_niC,                 &Reactor1G::set_niC)
        .add_property("miF",                 &Reactor1G::get_miF,                 &Reactor1G::set_miF)
        .add_property("miC",                 &Reactor1G::get_miC,                 &Reactor1G::set_miC)
        .add_property("NiF",                 &Reactor1G::get_NiF,                 &Reactor1G::set_NiF)
        .add_property("NiC",                 &Reactor1G::get_NiC,                 &Reactor1G::set_NiC)

        //The following attributes are read only!
        .add_property("dF_F_",               &Reactor1G::get_dF_F_)
        .add_property("dC_F_",               &Reactor1G::get_dC_F_)
        .add_property("BU_F_",               &Reactor1G::get_BU_F_)
        .add_property("P_F_",                &Reactor1G::get_P_F_)
        .add_property("D_F_",                &Reactor1G::get_D_F_)
        .add_property("k_F_",                &Reactor1G::get_k_F_)
        .add_property("Mj_F_",               &Reactor1G::get_Mj_F_)
        .add_property("zeta_F_",             &Reactor1G::get_zeta_F_)

        .add_property("fd",                  &Reactor1G::get_fd,                 &Reactor1G::set_fd)
        .add_property("Fd",                  &Reactor1G::get_Fd,                 &Reactor1G::set_Fd)
        .add_property("BUd",                 &Reactor1G::get_BUd,                &Reactor1G::set_BUd)
        .add_property("k",                   &Reactor1G::get_k,                  &Reactor1G::set_k)

        .add_property("InU",                 &Reactor1G::get_InU,                &Reactor1G::set_InU)
        .add_property("InTRU",               &Reactor1G::get_InTRU,              &Reactor1G::set_InTRU)
        .add_property("InLAN",               &Reactor1G::get_InLAN,              &Reactor1G::set_InLAN)
        .add_property("InACT",               &Reactor1G::get_InACT,              &Reactor1G::set_InACT)
        .add_property("OutU",                &Reactor1G::get_OutU,               &Reactor1G::set_OutU)
        .add_property("OutTRU",              &Reactor1G::get_OutTRU,             &Reactor1G::set_OutTRU)
        .add_property("OutLAN",              &Reactor1G::get_OutLAN,             &Reactor1G::set_OutLAN)
        .add_property("OutACT",              &Reactor1G::get_OutACT,             &Reactor1G::set_OutACT)

        .add_property("TruCR",               &Reactor1G::get_TruCR,              &Reactor1G::set_TruCR)

        //The following attributes are read only!
        .add_property("SigmaFa_F_",          &Reactor1G::get_SigmaFa_F_)
        .add_property("SigmaFtr_F_",         &Reactor1G::get_SigmaFtr_F_)
        .add_property("kappaF_F_",           &Reactor1G::get_kappaF_F_)

        .add_property("SigmaCa_F_",          &Reactor1G::get_SigmaCa_F_)
        .add_property("SigmaCtr_F_",         &Reactor1G::get_SigmaCtr_F_)
        .add_property("kappaC_F_",           &Reactor1G::get_kappaC_F_)

        .add_property("LatticeE_F_",         &Reactor1G::get_LatticeE_F_)
        .add_property("LatticeF_F_",         &Reactor1G::get_LatticeF_F_)

        //Reactor1G Constructor Overloads
        .def(bp::init< bp::optional<std::string> >())
        .def(bp::init< std::set<std::string>, bp::optional<std::string> >())
        .def(bp::init< ReactorParameters, bp::optional<std::string> >() )
        .def(bp::init< ReactorParameters, std::set<std::string>, bp::optional<std::string> >() )

        //Useful Functions
        .def("loadLib",             &Reactor1G::loadLib)
        .def("initialize",          &Reactor1G::initialize)
        .def("foldMassWeights",     &Reactor1G::foldMassWeights)
        .def("mkMj_F_",             &Reactor1G::mkMj_F_)
        .def("mkMj_Fd_",            &Reactor1G::mkMj_Fd_)
        .def("calcOutIso",          &Reactor1G::calcOutIso)
        .def("calcSubStreams",      &Reactor1G::calcSubStreams)
        .def("calcTruCR",           &Reactor1G::calcTruCR)
        .def("FluenceAtBU",         &Reactor1G::FluenceAtBU)
        .def("batchAve",            &Reactor1G::batchAve, batchAve_overloads("Obtains batch averaged value of P, D, or k."))
        .def("batchAveK",           &Reactor1G::batchAveK)
        .def("BUd_BisectionMethod", &Reactor1G::BUd_BisectionMethod)
        .def("Run_PNL",             &Reactor1G::Run_PNL)
        .def("Calibrate_PNL_2_BUd", &Reactor1G::Calibrate_PNL_2_BUd)

        .def("doCalc", Reactor1G_doCalc_NA)
        .def("doCalc", Reactor1G_doCalc_CD)
        .def("doCalc", Reactor1G_doCalc_MS)
    ;

    //Grabs a copy of the Fast Reactor Default Settings
    bp::def("FRDefaults", &fillFRDefaults);
    
    bp::class_< FastReactor1G, bp::bases<Reactor1G> >("FastReactor1G", "One-Group Fast Reactor Model",  bp::init<>() )
        //FastReactor1G Constructor Overloads
        .def(bp::init< std::string, bp::optional<std::string> >())
        .def(bp::init< ReactorParameters, bp::optional<std::string> >())
        .def(bp::init< std::string, ReactorParameters, bp::optional<std::string> >())

        //Useful Functions
        .def("setParams",  &FastReactor1G::setParams)
    ;

    //Grabs a copy of the Light Water Reactor Default Settings
    bp::def("LWRDefaults", &fillLWRDefaults);
    
    bp::class_< LightWaterReactor1G, bp::bases<Reactor1G> >("LightWaterReactor1G", "One-Group Light Water Reactor Model",  bp::init<>() )
        //LightWaterReactor1G Constructor Overloads
        .def(bp::init< std::string, bp::optional<std::string> >())
        .def(bp::init< ReactorParameters, bp::optional<std::string> >())
        .def(bp::init< std::string, ReactorParameters, bp::optional<std::string> >())

        //Useful Functions
        .def("setParams",  &LightWaterReactor1G::setParams)
    ;



    // Fuel Fabrication to- and from-converters
    dict2map<std::string, MassStream *>();
    bp::to_python_converter< MassStreams, map2dict<std::string, MassStream *> >();

    // Fuel Fabrication Facility
    bp::class_< FuelFabrication, bp::bases<FCComp> >("FuelFabrication", "Fuel Fabrication Facility", bp::init<>() )
        // Class Attributes
        .add_property("mass_streams", &FuelFabrication::get_mass_streams, &FuelFabrication::set_mass_streams)
        .add_property("mass_weights", &FuelFabrication::get_mass_weights, &FuelFabrication::set_mass_weights)
        .add_property("mass_deltaRs", &FuelFabrication::get_mass_deltaRs, &FuelFabrication::set_mass_deltaRs)

        .add_property("reactor", &FuelFabrication::get_reactor, &FuelFabrication::set_reactor)

        // Fuel Fabrication Component Constructor
        .def(bp::init< std::string >())
        .def(bp::init< std::set<std::string>, bp::optional<std::string> >())
        .def(bp::init< MassStreams, MassWeights, Reactor1G, bp::optional< std::set<std::string>, std::string > >())

        // Useful Functions
        .def("initialize", &FuelFabrication::initialize)
    ;

};
