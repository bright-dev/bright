"""Python wrapper for fccomp."""
# Cython imports
from libcpp.utility cimport pair as cpp_pair
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free

cimport numpy as np
import numpy as np

from pyne cimport std
from pyne cimport nucname

from pyne cimport stlconverters as conv
from pyne import stlconverters as conv

cimport pyne.cpp_material
cimport pyne.material
import pyne.material

cimport cpp_fccomp

cimport fccomp


####################
### FCComp Class ###
####################

cdef class FCComp:
    """Base Fuel Cycle Component Class.

    Args:
        * paramlist (sequence of str): A set of parameter names (str) that the component will track.
        * name (str): The name of the fuel cycle component instance.

    Note that this automatically calls the protected :meth:`initialize` C function.
    """

    def __cinit__(self, params=None, char * name="", *args, **kwargs):
        cdef cpp_set[std.string] param_set

        if params is None:
            self._inst = new cpp_fccomp.FCComp(std.string(name))
        else:
            param_set = conv.py_to_cpp_set_str(params)
            self._inst = new cpp_fccomp.FCComp(param_set, std.string(name))

        self._free_inst = True

        # property defaults
        self._params_prior_calc = None
        self._params_after_calc = None


    def __dealloc__(self):
        if self._free_inst:
            free(self._inst)


    #
    # Class Attributes
    #

    property name:
        def __get__(self):
            cdef std.string n = (<cpp_fccomp.FCComp *> self._inst).name
            return n.c_str()

        def __set__(self, char * n):
            (<cpp_fccomp.FCComp *> self._inst).name = std.string(n)


    property natural_name:
        def __get__(self):
            cdef std.string n = (<cpp_fccomp.FCComp *> self._inst).natural_name
            return n.c_str()

        def __set__(self, char * n):
            (<cpp_fccomp.FCComp *> self._inst).natural_name = std.string(n)


    property mat_feed:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_fccomp.FCComp *> self._inst).mat_feed
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_fccomp.FCComp *> self._inst).mat_feed = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_fccomp.FCComp *> self._inst).mat_prod
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_fccomp.FCComp *> self._inst).mat_prod = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property params_prior_calc:
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._params_prior_calc is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_fccomp.FCComp *> self._inst).params_prior_calc
                self._params_prior_calc = proxy

            return self._params_prior_calc

            #return conv.map_to_dict_str_dbl((<cpp_fccomp.FCComp *> self._inst).params_prior_calc)

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_fccomp.FCComp *> self._inst).params_prior_calc = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_fccomp.FCComp *> self._inst).params_prior_calc = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_fccomp.FCComp *> self._inst).params_prior_calc = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._params_prior_calc = None
            #(<cpp_fccomp.FCComp *> self._inst).params_prior_calc = conv.dict_to_map_str_dbl(pi)


    property params_after_calc:
        def __get__(self):
            cdef conv._MapStrDouble proxy

            if self._params_after_calc is None:
                proxy = conv.MapStrDouble(False, False)
                proxy.map_ptr = &(<cpp_fccomp.FCComp *> self._inst).params_after_calc
                self._params_after_calc = proxy

            return self._params_after_calc

            #return conv.map_to_dict_str_dbl((<cpp_fccomp.FCComp *> self._inst).params_after_calc)

        def __set__(self, value):
            cdef std.string s
            cdef cpp_pair[std.string, double] item
            cdef cpp_map[std.string, double]  m

            if isinstance(value, conv._MapStrDouble):
                (<cpp_fccomp.FCComp *> self._inst).params_after_calc = deref((<conv._MapStrDouble> value).map_ptr)
            elif hasattr(value, 'items'):
                m = cpp_map[std.string, double]()
                for k, v in value.items():
                    s = std.string(k)
                    item = cpp_pair[std.string, double](s, v)
                    m.insert(item)
                (<cpp_fccomp.FCComp *> self._inst).params_after_calc = m
            elif hasattr(value, '__len__'):
                m = cpp_map[std.string, double]()
                for i in value:
                    s = std.string(i[0])
                    item = cpp_pair[std.string, double](s, i[1])
                    m.insert(item)
                (<cpp_fccomp.FCComp *> self._inst).params_after_calc = m
            else:
                raise TypeError('{0} cannot be converted to a map.'.format(type(value)))

            self._params_after_calc = None

            #(<cpp_fccomp.FCComp *> self._inst).params_after_calc = conv.dict_to_map_str_dbl(po)


    property pass_num:
        def __get__(self):
            return (<cpp_fccomp.FCComp *> self._inst).pass_num

        def __set__(self, int pn):
            (<cpp_fccomp.FCComp *> self._inst).pass_num = pn


    property track_params:
        def __get__(self):
            return conv.cpp_to_py_set_str((<cpp_fccomp.FCComp *> self._inst).track_params)

        def __set__(self, set p2t):
            (<cpp_fccomp.FCComp *> self._inst).track_params = conv.py_to_cpp_set_str(p2t)

    #
    # Class Methods
    #

    def write_mat_pass(self):
        """This method is responsible for adding a new pass to the output text file 
        "{FCComp.name}Isos.txt" for this component.  Further calculations should
        not be performed after :meth:`write_mat_pass` has been called.

        This function has one very important subtlety: it does not write out mass streams data.
        Rather, input columns are given as normalized isotopic vectors.
        As weight fractions, input columns are in units of [kgInIso/kgmat_feed.mass].
        Moreover, the output columns are given in terms relative to the mass of the input mass, 
        [kgOutIso/kgmat_feed.mass].  These are calculated via the following expressions.

        .. math::

            \mbox{inpcol[iso]} = \mbox{mat_feed.comp[iso]}

            \mbox{outcol[iso]} = \mbox{mat_prod.comp[iso]} \times \frac{\mbox{mat_prod.mass}}{\mbox{mat_feed.mass}}

        Because of the units of these two columns, total mass flow data may often only be recovered via the 
        a "Mass" parameter in the "{FCComp.name}Paramat.txt" file.  Here is a sample LWRIsos.txt file for a
        light water reactor for the first pass::

            Isotope 1in             1out    
            H1      0.000000E+00    0.000000E+00
            H3      0.000000E+00    8.568522E-08
            HE4     0.000000E+00    4.421615E-07
            B10     0.000000E+00    0.000000E+00
            B11     0.000000E+00    0.000000E+00
            C14     0.000000E+00    4.015091E-11
            O16     0.000000E+00    0.000000E+00
            SR90    0.000000E+00    8.221283E-04
            TC99    0.000000E+00    1.112580E-03
            CS137   0.000000E+00    1.821226E-03
            U234    0.000000E+00    2.807466E-06
            U235    4.773292E-02    8.951725E-03
            U236    0.000000E+00    6.155297E-03
            U237    0.000000E+00    1.719458E-05
            U238    9.522671E-01    9.211956E-01
            U239    0.000000E+00    6.953862E-07
            NP237   0.000000E+00    8.057270E-04
            PU238   0.000000E+00    2.842232E-04
            PU239   0.000000E+00    5.353362E-03
            PU240   0.000000E+00    2.114728E-03

        """
        (<cpp_fccomp.FCComp *> self._inst).write_mat_pass()


    def write_params_pass(self):
        """What write_ms_pass() does for a component's input and output isotopics, 
        this function does for the components parameters.  To ensure that meaningful 
        data is available, write_params_pass() first must have calc_params()
        called elsewhere in the program.  Note that to get the pass numbering correct, 
        pass_num should always be incremented prior to this method.  The 
        following is an example of "{FCComp.name}Paramat.txt" for a light water 
        reactor spent fuel reprocessing facility::

            Param   1in             1out    
            Mass    9.985828E-01    9.975915E-01

        """
        (<cpp_fccomp.FCComp *> self._inst).write_params_pass()
        

    def write_text(self):
        """This method calls write_ms_pass() and then, if available, calls 
        write_params_pass().  This is convience function for producing 
        text-based output.  However, using write() is recommended.
        """
        (<cpp_fccomp.FCComp *> self._inst).write_text()


    def write_hdf5(self):
        """This method writes out the isotopic pass data to an HDF5 file. 
        Then, if available, it also writes parameter data as well.  
        Using write() instead is recommended.
        """
        (<cpp_fccomp.FCComp *> self._inst).write_hdf5()


    def write(self):
        """This is a convenience function that first increments up pass_num.
        Then, it checks to see if there are any parameters for this component.
        If there are, it sets the current values using :meth:`calc_params`.

        If bright.write_hdf5 is set, then write_hdf5() is called.

        If bright.write_text is set, then write_text() is called.

        This is what is most often used to write Bright output.  Therefore it is
        seen as the last step for every component in each pass.
        """
        (<cpp_fccomp.FCComp *> self._inst).write()


    # Virtual methods

    def calc_params(self):
        """By calling this method, all parameter values are calculated and set for the fuel cycle component.
        This should be done following a calc() calculation but before data is written out.
        If a component has important parameters associated with it, this function must be overridden and called.

        Note that this is called first thing when write_params_pass() is called.  For example, reprocessing only 
        has a "Mass" parameter.  Translated into Python, calc_params() here looks like the following::

            def calc_params(self):
                self.params_prior_calc["Mass"]  = self.mat_feed.mass
                self.params_after_calc["Mass"] = self.mat_prod.mass
                return
        """
        (<cpp_fccomp.FCComp *> self._inst).calc_params()


    def calc(self):
        """This method is used to determine a component's output isotopics from its input isotopics.
        Therefore, this is typically where the bulk of a fuel cycle component's algorithm lies.
        As each component type has a distinct methodology, the calc() method  needs 
        to be overridden child classes.

        This method should return mat_prod so that component calculations may be easily 
        daisy-chained together.

        Args:
            * input (dict or Material): If input is present, it set as the component's 
              mat_feed.  If input is a isotopic dictionary (zzaaam keys, float values), this
              dictionary is first converted into a Material before being set as mat_feed.

        Returns:
            * output (Material): mat_prod.

        """
        cdef pyne.material._Material pymat = pyne.material.Material()
        pymat.mat_pointer[0] = (<cpp_fccomp.FCComp *> self._inst).calc()
        return pymat

