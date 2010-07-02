// Header file for Helper functions for the Bright Python Wrapper

#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <boost/python/str.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python/extract.hpp>

#include <boost/python/converter/registered.hpp>
#include <boost/python/exception_translator.hpp>

#include <map>
#include <set>
#include <vector>

#include <iostream>

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

template <class T, class U>
struct map2dict
{
    static PyObject* convert(const std::map<T, U>& m)
    {
        bp::dict d;
        for (typename std::map<T, U>::const_iterator i = m.begin(); i != m.end(); i++)
        {
            d[(i->first)] = (i->second);
        };
        return bp::incref(d.ptr());
    }
};

template <class T, class U>
struct dict2map
{
    //dict2map constuctor
    dict2map()
    {
        bp::converter::registry::push_back( &convertible, &construct, bp::type_id< std::map<T, U> >() );
    }

    // Determine if obj_ptr can be converted into a map
    static void* convertible(PyObject* obj_ptr)
    {
        if (!PyDict_Check(obj_ptr)) 
            return 0;
        return obj_ptr;
    };

    // Convert obj_ptr into a map
    static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data)
    {
        // Extract the character data from the python string
        std::map<T, U> m;
        PyObject * pyd_items = PyDict_Items( obj_ptr );
        for (int i = 0; i < PyDict_Size(obj_ptr); i++)
        {
            PyObject * k_i = PyList_GetItem( pyd_items, i);
            T newkey = bp::extract<T>( PyTuple_GetItem(k_i, 0) ); 
            U newval = bp::extract<U>( PyTuple_GetItem(k_i, 1) ); 
            m[ newkey ] = newval;
        }

        // Verify that obj_ptr is a dictionary (should be ensured by convertible())
//        #ifndef _WIN32
//            assert(m);
//        #endif

        // Grab pointer to memory into which to construct the new map
        void* storage = ( (bp::converter::rvalue_from_python_storage< std::map<T, U> >*) data)->storage.bytes;

        // in-place construct the new map using the data
        // extraced from the python object
        new (storage) std::map<T, U> (m);

        // Stash the memory chunk pointer for later use by boost.python
         data->convertible = storage;
    };
};

template <class T>
struct c_set2py_list 
{
    static PyObject* convert(const std::set<T>& s )
    {
            bp::list l;
        	for (typename std::set<T>::const_iterator i = s.begin(); i != s.end(); i++)
            {
        	        l.append( *i );
            };
        	return bp::incref(l.ptr());
    }
};

template <class T>
struct py_list2c_set
{
    //py_list2c_set constuctor
    py_list2c_set()
    {
        bp::converter::registry::push_back( &convertible, &construct, bp::type_id< std::set<T> >() );
    }

    // Determine if obj_ptr can be converted into a map
    static void* convertible(PyObject* obj_ptr)
    {
        if (!PyList_Check(obj_ptr)) 
            return 0;
        return obj_ptr;
    };

    // Convert obj_ptr into a map
    static void construct(PyObject* obj_ptr, bp::converter::rvalue_from_python_stage1_data* data)
    {
        // Extract the character data from the python string
        std::set<T> s;
        for (int i = 0; i < PyList_Size(obj_ptr); i++)
        {
            T newval = bp::extract<T>( PyList_GetItem(obj_ptr, i) ); 
            s.insert( newval );
        }

        // Verify that obj_ptr is a py_list (should be ensured by convertible())
//        #ifndef _WIN32
//            assert(s);
//        #endif

        // Grab pointer to memory into which to construct the new map
        void* storage = ( (bp::converter::rvalue_from_python_storage< std::set<T> >*) data)->storage.bytes;

        // in-place construct the new map using the data
        // extraced from the python object
        new (storage) std::set<T> (s);

        // Stash the memory chunk pointer for later use by boost.python
         data->convertible = storage;
    };
};

template <class T>
struct c_vector2py_list 
{
    static PyObject* convert(const std::vector<T>& v )
    {
            bp::list l;
        	for (typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); i++)
            {
        	        l.append( *i );
            };
        	return bp::incref(l.ptr());
    }
};


//Failed Exception Template
//template <class T>
//void BriPy_Exception_Translator(T const& e)
//{
    //class T is the C exception to be translated
    //class U is a charater string that is the name 
    //	of the exception, in both Python and C.
//	static PyObject * NewErr; 
//	NewErr = PyErr_NewException(e.name(), PyExc_StandardError, 0);
//	PyErr_SetString(NewErr, e.what());
    //PyErr_SetString(&PyErr_NewException(e.name(), NULL, NULL), e.what());
//};

template <class E>
void BriPy_Exception_Translator(E const& e)
{
    PyErr_SetString(PyExc_RuntimeError, ((std::string) e.name() + ": " + (std::string) e.what()).c_str() );
};
