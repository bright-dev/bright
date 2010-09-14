//The Bright Python wrapper for "isoname"

//standard library headers

//Needed Boost Modules 
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

//My Files
#include "BriPyHelper.h"

#ifdef _WIN32
    #include "../FCComps/isoname.h"
#else
    #include "src/isoname.h"
#endif

//Need to keep Boost::Python in its own namespace.
//      NO: using namespace boost::python;
//possible because both Boost::Python and HDF5 define
//a 'ssize_t' variable that colides in an abiguous 
//definition at global namespce.  Unfortunately, HDF5
//does not encapsulate itself in a namespce.
//Thus to use these two codes together we have to use     
//Boost::Python's namespace explicitly.
//Renamed here for ease of use.
namespace bp = boost::python;

BOOST_PYTHON_MODULE(isoname)
{
    //sets the current scope...
    bp::scope current;
    current.attr("__doc__") = "Nuclide naming convention definitions and functions.";

    //to- and from-converters
    dict2map<isoname::LL, isoname::zz>();
    bp::to_python_converter< isoname::LLzzType, map2dict<isoname::LL, isoname::zz> >();

    dict2map<isoname::zz, isoname::LL>();
    bp::to_python_converter< isoname::zzLLType, map2dict<isoname::zz, isoname::LL> >();

    py_list2c_set<isoname::LL_Name>();
    bp::to_python_converter< isoname::LL_Group, c_set2py_list<isoname::LL_Name> >();

    py_list2c_set<isoname::zz_Name>();
    bp::to_python_converter< isoname::zz_Group, c_set2py_list<isoname::zz_Name> >();

    //Element dictionaries
    current.attr("LLzz") = &isoname::LLzz;
    current.attr("zzLL") = &isoname::zzLL;

    //Element groups
    current.attr("LAN") = &isoname::LAN;
    current.attr("ACT") = &isoname::ACT;
    current.attr("TRU") = &isoname::TRU;
    current.attr("MA")  = &isoname::MA;
    current.attr("FP")  = &isoname::FP;

    current.attr("lan") = &isoname::lan;
    current.attr("act") = &isoname::act;
    current.attr("tru") = &isoname::tru;
    current.attr("ma")  = &isoname::ma;
    current.attr("fp")  = &isoname::fp;

    //CurrentForm Wrapper
    char CurrentForm__doc__ [] = "Determines the current form of a nuclide from [LLAAAM, zzaaam, MCNP].";
    std::string (*CurrentForm1)(std::string) = &isoname::CurrentForm;
    std::string (*CurrentForm2)(int)         = &isoname::CurrentForm;
    bp::def("CurrentForm", CurrentForm1, (bp::arg("nuc")), CurrentForm__doc__);
    bp::def("CurrentForm", CurrentForm2, (bp::arg("nuc")), CurrentForm__doc__);

    //LLAAAM_2_* Functions
    bp::def("LLAAAM_2_zzaaam", &isoname::LLAAAM_2_zzaaam, "Converts nuclide from LLAAAM -> zzaaam.");
    bp::def("LLAAAM_2_MCNP",   &isoname::LLAAAM_2_MCNP,   "Converts nuclide from LLAAAM -> MCNP."  );

    //zzaaam_2_* Functions
    bp::def("zzaaam_2_LLAAAM", &isoname::zzaaam_2_LLAAAM, "Converts nuclide from zzaaam -> LLAAAM.");
    bp::def("zzaaam_2_MCNP",   &isoname::zzaaam_2_MCNP,   "Converts nuclide from zzaaam -> MCNP."  );

    //MCNP_2_* Functions
    bp::def("MCNP_2_zzaaam", &isoname::MCNP_2_zzaaam, "Converts nuclide from MCNP -> zzaaam.");
    bp::def("MCNP_2_LLAAAM", &isoname::MCNP_2_LLAAAM, "Converts nuclide from MCNP -> LLAAAM.");

    //mixed_2_* Functions
    char mixed_2_zzaaam__doc__ [] = "Converts a nuclide to zzaaam form.";
    int (*mixed_2_zzaaam1)(std::string) = &isoname::mixed_2_zzaaam;
    int (*mixed_2_zzaaam2)(int)         = &isoname::mixed_2_zzaaam;
    bp::def("mixed_2_zzaaam", mixed_2_zzaaam1, mixed_2_zzaaam__doc__);
    bp::def("mixed_2_zzaaam", mixed_2_zzaaam2, mixed_2_zzaaam__doc__);

    char mixed_2_LLAAAM__doc__ [] = "Converts a nuclide to LLAAAM form.";
    std::string (*mixed_2_LLAAAM1)(std::string) = &isoname::mixed_2_LLAAAM;
    std::string (*mixed_2_LLAAAM2)(int)         = &isoname::mixed_2_LLAAAM;
    bp::def("mixed_2_LLAAAM", mixed_2_LLAAAM1, mixed_2_LLAAAM__doc__);
    bp::def("mixed_2_LLAAAM", mixed_2_LLAAAM2, mixed_2_LLAAAM__doc__);

    char mixed_2_MCNP__doc__ [] = "Converts a nuclide to MCNP form.";
    int (*mixed_2_MCNP1)(std::string) = &isoname::mixed_2_MCNP;
    int (*mixed_2_MCNP2)(int)         = &isoname::mixed_2_MCNP;
    bp::def("mixed_2_MCNP", mixed_2_MCNP1, mixed_2_MCNP__doc__);
    bp::def("mixed_2_MCNP", mixed_2_MCNP2, mixed_2_MCNP__doc__);
};
