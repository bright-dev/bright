//The Bright Python wrapper for "MassStream"

//standard library headers

//Needed Boost Modules
#include <boost/python.hpp>
#include <boost/python/overloads.hpp>

//My Files
#include "BriPyHelper.h"
#include "src/MassStream.h"

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

//Default Argument Function Overload Macros
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(gSS_overloads, getSubStream, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getU_overloads,   getU,   0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getPU_overloads,  getPU,  0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getLAN_overloads, getLAN, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getACT_overloads, getACT, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getTRU_overloads, getTRU, 0, 1) 
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getMA_overloads,  getMA,  0, 1)  
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getFP_overloads,  getFP,  0, 1)  

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(lfhdf5_overloads, load_from_hdf5, 2, 3)

BOOST_PYTHON_MODULE(MassStream)
{
    //sets the current scope...
    bp::scope current;
    current.attr("__doc__") = "Mass Stream Class & Module Wrapper";

    //to- and from-converters
    dict2map<Iso, Weight>();
    bp::to_python_converter< CompDict, map2dict<Iso, Weight> >();

    py_list2c_set<int>();
    bp::to_python_converter< std::set<int> , c_set2py_list<int> >();

    py_list2c_set<std::string>();
    bp::to_python_converter< std::set<std::string> , c_set2py_list<std::string> >();

    //Function Pointers
    MassStream (MassStream::*getSubStreamInt)(std::set<int>, std::string) = &MassStream::getSubStream;
    MassStream (MassStream::*getSubStreamString)(std::set<std::string>, std::string) = &MassStream::getSubStream;
        
    bp::class_< MassStream >("MassStream", "Bright/Python Mass Stream Class",  bp::init<>("Empty MassStream Constructor.") )
        //Basic MassStream Properties
        .add_property("comp", &MassStream::get_comp, &MassStream::set_comp)
        .add_property("mass", &MassStream::get_mass, &MassStream::set_mass)
        .add_property("name", &MassStream::get_name, &MassStream::set_name)

        //MassStream Constructors
        .def(bp::init<CompDict,    bp::optional<double, std::string> >("MassStream Isotopic Component Dictionary Constructor."))
        .def(bp::init<char *,      bp::optional<double, std::string> >("MassStream Filename (Character Array) Constructor."))
        .def(bp::init<std::string, bp::optional<double, std::string> >("MassStream Filename (String) Constructor."))

        .def("load_from_hdf5", &MassStream::load_from_hdf5, lfhdf5_overloads("Loads data into the MassStream from an HDF5 file."))
        .def("load_from_text", &MassStream::load_from_text, "Loads data into the MassStream from a text file.")

        //Useful Functions
        .def("Print",      &MassStream::Print, "Print MassStream to std::out.")
        .def("Normalize",  &MassStream::Normalize, "Normalize MassStream by setting mass = 1.")
        .def("multByMass", &MassStream::multByMass, "Return Isotopic Component Dictionary whose values = comp[iso] * mass.")

        //Sub-Stream Generators
        .def("getSubStreamInt", getSubStreamInt, gSS_overloads( "Obtains a sub-stream from this MassStream.\nSub-stream contains only the isotopes (int) in the list passed.\nUsage:\n\tthis_ms.getSubStream([92, 942390, ...])" ))
        .def("getSubStreamStr", getSubStreamString, gSS_overloads("Obtains a sub-stream from this MassStream.\nSub-stream contains only the isotopes (string) in the list passed.\nUsage: this_ms.getSubStream([92, 942390, ...])") )
        .def("getU",   &MassStream::getU,   getU_overloads("Obtains a sub-stream of Uranium from this MassStream.") )
        .def("getPU",  &MassStream::getPU,  getPU_overloads("Obtains a sub-stream of Plutonium from this MassStream.") )
        .def("getLAN", &MassStream::getLAN, getLAN_overloads("Obtains a sub-stream of Lanthanides from this MassStream.") )
        .def("getACT", &MassStream::getACT, getACT_overloads("Obtains a sub-stream of Actinides from this MassStream.") )
        .def("getTRU", &MassStream::getTRU, getTRU_overloads("Obtains a sub-stream of Transuranics from this MassStream.") )
        .def("getMA",  &MassStream::getMA,  getMA_overloads("Obtains a sub-stream of Minor Actinides from this MassStream.") )
        .def("getFP",  &MassStream::getFP,  getFP_overloads("Obtains a sub-stream of Fission Products from this MassStream.") )

        //Opperator Overloads
        .def(bp::self_ns::str(bp::self))

        .def(bp::self + double())
        .def(double() + bp::self)
        .def(bp::self + bp::self)

        .def(bp::self * double())
        .def(double() * bp::self)
        .def(bp::self / double())
    ;
};
