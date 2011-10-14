"""Python wrapper for RMG."""
# Cython imports
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

cimport pyne.cpp_material
cimport pyne.material
import pyne.material

cimport cpp_fccomp
cimport cpp_reactor_parameters
cimport cpp_fluence_point
cimport cpp_reactormg

from bright.reactor_parameters cimport ReactorParameters
from bright.reactor_parameters import ReactorParameters

from bright.fluence_point cimport FluencePoint
from bright.fluence_point import FluencePoint

cimport fccomp
import fccomp


#######################
### ReactorMG Class ###
#######################

cdef class ReactorMG(fccomp.FCComp):
    """Multi-Group Reactor Fuel Cycle Component Class.  Daughter of bright.FCComp class.

    Args:
        * reactor_parameters (ReactorParameters): A special data structure that contains information
          on how to setup and run the reactor.
        * track_params (string set): A set of strings that represents what parameter data the reactor should 
          store and set.  Different reactor types may have different characteristic parameters that are of interest.
        * name (str): The name of the reactor fuel cycle component instance.

    Note that this automatically calls the public initialize() C function.

    .. note:: 

        Some data members and functions have names that end in '_F_'.  This indicates that these are a 
        function of fluence, the time integral of the flux.  The '_Fd_' suffix implies that the data is 
        evaluated at the discharge fluence.
    """

    def __cinit__(self, reactor_parameters=None, track_params=None, char * name="", *args, **kwargs):
        cdef ReactorParameters rp
        cdef std.string cpp_name = std.string(name)

        if (reactor_parameters is None) and (track_params is None):
            self._inst = new cpp_reactormg.ReactorMG(cpp_name)

        elif (reactor_parameters is None) and isinstance(track_params, set):
            self._inst = new cpp_reactormg.ReactorMG(conv.py_to_cpp_set_str(track_params), cpp_name)

        elif isinstance(reactor_parameters, ReactorParameters) and (track_params is None):
            rp = reactor_parameters
            self._inst = new cpp_reactormg.ReactorMG(<cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0], cpp_name)

        elif isinstance(reactor_parameters, ReactorParameters) and isinstance(track_params, set):
            rp = reactor_parameters
            self._inst = new cpp_reactormg.ReactorMG(<cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0], conv.py_to_cpp_set_str(track_params), cpp_name)

        else:
            if reactor_parameters is not None:
                raise TypeError("The reactor_parameters keyword must be an instance of the ReactorParameters class or None.  Got " + str(type(reactor_parameters)))

            if track_params is not None:
                raise TypeError("The track_params keyword must be a set of strings or None.  Got " + str(type(track_params)))

    #
    # Class Attributes
    #

    # ReactorMG attributes

    property B:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).B

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).B = value


    property flux:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).flux

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).flux = value





    property chemical_form_fuel:
        def __get__(self):
            return conv.map_to_dict_str_dbl((<cpp_reactormg.ReactorMG *> self._inst).chemical_form_fuel)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_fuel = conv.dict_to_map_str_dbl(value)


    property chemical_form_clad:
        def __get__(self):
            return conv.map_to_dict_str_dbl((<cpp_reactormg.ReactorMG *> self._inst).chemical_form_clad)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_clad = conv.dict_to_map_str_dbl(value)


    property chemical_form_cool:
        def __get__(self):
            return conv.map_to_dict_str_dbl((<cpp_reactormg.ReactorMG *> self._inst).chemical_form_cool)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chemical_form_cool = conv.dict_to_map_str_dbl(value)





    property rho_fuel:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).rho_fuel

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).rho_fuel = value


    property rho_clad:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).rho_clad

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).rho_clad = value


    property rho_cool:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).rho_cool

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).rho_cool = value





    property P_NL:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).P_NL

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).P_NL = value


    property target_BU:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).target_BU

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).target_BU = value


    property specific_power:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).specific_power

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).specific_power = value


    property burn_regions:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).burn_regions

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).burn_regions = value


    property S:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).S

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).S = value


    property burn_time:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).burn_time

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).burn_time = value


    property bt_s:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).bt_s

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).bt_s = value


    property burn_times:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).burn_times)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).burn_times = conv.array_to_vector_1d_dbl(value)





    property use_zeta:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).use_zeta

        def __set__(self, bint value):
            (<cpp_reactormg.ReactorMG *> self._inst).use_zeta = value


    property lattice_flag:
        def __get__(self):
            cdef std.string value = (<cpp_reactormg.ReactorMG *> self._inst).lattice_flag
            return value.c_str()

        def __set__(self, char * value):
            (<cpp_reactormg.ReactorMG *> self._inst).lattice_flag = std.string(value)


    property rescale_hydrogen_xs:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).rescale_hydrogen_xs

        def __set__(self, bint value):
            (<cpp_reactormg.ReactorMG *> self._inst).rescale_hydrogen_xs = value


    property burnup_via_constant:
        def __get__(self):
            cdef std.string value = (<cpp_reactormg.ReactorMG *> self._inst).burnup_via_constant
            return value.c_str()

        def __set__(self, char * value):
            (<cpp_reactormg.ReactorMG *> self._inst).burnup_via_constant = std.string(value)


    property branch_ratio_cutoff:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).branch_ratio_cutoff

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).branch_ratio_cutoff = value





    property r_fuel:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).r_fuel

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).r_fuel = value


    property r_void:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).r_void

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).r_void = value


    property r_clad:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).r_clad

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).r_clad = value


    property pitch:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).pitch

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).pitch = value





    property S_O:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).S_O

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).S_O = value


    property S_T:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).S_T

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).S_T = value


    property V_fuel:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).V_fuel

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).V_fuel = value


    property V_clad:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).V_clad

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).V_clad = value


    property V_cool:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).V_cool

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).V_cool = value







    property libfile:
        def __get__(self):
            cdef std.string value = (<cpp_reactormg.ReactorMG *> self._inst).libfile
            return value.c_str()

        def __set__(self, char * value):
            (<cpp_reactormg.ReactorMG *> self._inst).libfile = std.string(value)






    property I:
        def __get__(self):
            return conv.cpp_to_py_set_int((<cpp_reactormg.ReactorMG *> self._inst).I)

        def __set__(self, set value):
            (<cpp_reactormg.ReactorMG *> self._inst).I = conv.py_to_cpp_set_int(value)


    property J:
        def __get__(self):
            return conv.cpp_to_py_set_int((<cpp_reactormg.ReactorMG *> self._inst).J)

        def __set__(self, set value):
            (<cpp_reactormg.ReactorMG *> self._inst).J = conv.py_to_cpp_set_int(value)


    property K:
        def __get__(self):
            return conv.cpp_to_py_set_int((<cpp_reactormg.ReactorMG *> self._inst).K)

        def __set__(self, set value):
            (<cpp_reactormg.ReactorMG *> self._inst).K = conv.py_to_cpp_set_int(value)







    property K_num:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).K_num

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).K_num = value


    property K_ord:
        def __get__(self):
            return conv.vector_to_array_1d_int((<cpp_reactormg.ReactorMG *> self._inst).K_ord)

        def __set__(self, set value):
            (<cpp_reactormg.ReactorMG *> self._inst).K_ord = conv.array_to_vector_1d_int(value)


    property K_ind:
        def __get__(self):
            return conv.map_to_dict_int_int((<cpp_reactormg.ReactorMG *> self._inst).K_ind)

        def __set__(self, set value):
            (<cpp_reactormg.ReactorMG *> self._inst).K_ind = conv.dict_to_map_int_int(value)





    """\
    property decay_matrix:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).decay_matrix)

        def __set__(self, value):
            (<cpp_reactormg.ReactorMG *> self._inst).decay_matrix = conv.array_to_vector_2d_dbl(value)


    property thermal_yield_matrix:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).thermal_yield_matrix)

        def __set__(self, value):
            (<cpp_reactormg.ReactorMG *> self._inst).thermal_yield_matrix = conv.array_to_vector_2d_dbl(value)


    property fast_yield_matrix:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).fast_yield_matrix)

        def __set__(self, value):
            (<cpp_reactormg.ReactorMG *> self._inst).fast_yield_matrix = conv.array_to_vector_2d_dbl(value)


    property fission_product_yield_matrix:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).fission_product_yield_matrix)

        def __set__(self, value):
            (<cpp_reactormg.ReactorMG *> self._inst).fission_product_yield_matrix = conv.array_to_vector_3d_dbl(value)

    """






    # Perturbation table goes here


    property nperturbations:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).nperturbations

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).nperturbations = value


    property perturbed_fields:
        def __get__(self):
            return conv.map_to_dict_str_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).perturbed_fields)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).perturbed_fields = conv.dict_to_map_str_array_to_vector_1d_dbl(value)








    property G:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).G

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).G = value


    property E_g:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).E_g)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).E_g = conv.array_to_vector_1d_dbl(value)


    property phi_g:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).phi_g)

        def __set__(self, np.ndarray[np.float64_t, ndim=2] value):
            (<cpp_reactormg.ReactorMG *> self._inst).phi_g = conv.array_to_vector_2d_dbl(value)


    property phi:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).phi)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).phi = conv.array_to_vector_1d_dbl(value)


    property Phi:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Phi)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).Phi = conv.array_to_vector_1d_dbl(value)


    property time0:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).time0)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).time0 = conv.array_to_vector_1d_dbl(value)


    property BU0:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).BU0)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).BU0 = conv.array_to_vector_1d_dbl(value)






    property Ti0:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Ti0)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Ti0 = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property sigma_t_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_t_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_t_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)
    
    
    property sigma_a_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_a_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_a_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)
    
    
    property nubar_sigma_f_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_sigma_f_pg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_sigma_f_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property chi_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_s_pgh:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_s_pgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_s_pgh = conv.dict_to_map_int_array_to_vector_3d_dbl(value)


    property sigma_f_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_f_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_f_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_gamma_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)
    

    property sigma_2n_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_3n_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_3n_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_3n_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_alpha_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_alpha_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_alpha_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_proton_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_proton_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_proton_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_gamma_x_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_x_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_x_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)
    

    property sigma_2n_x_pg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_x_pg)
    
        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_x_pg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)






    property A_HM_t:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).A_HM_t)

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).A_HM_t = conv.array_to_vector_1d_dbl(value)


    property MW_fuel_t:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).MW_fuel_t)

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).MW_fuel_t = conv.array_to_vector_1d_dbl(value)


    property MW_clad_t:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).MW_clad_t)

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).MW_clad_t = conv.array_to_vector_1d_dbl(value)


    property MW_cool_t:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).MW_cool_t)

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).MW_cool_t = conv.array_to_vector_1d_dbl(value)


    property n_fuel_it:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).n_fuel_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).n_fuel_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property n_clad_it:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).n_clad_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).n_clad_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)

    
    property n_cool_it:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).n_cool_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).n_cool_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property m_fuel_it:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).m_fuel_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).m_fuel_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property m_clad_it:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).m_clad_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).m_clad_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)

    
    property m_cool_it:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).m_cool_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).m_cool_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property N_fuel_it:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).N_fuel_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).N_fuel_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property N_clad_it:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).N_clad_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).N_clad_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)

    
    property N_cool_it:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).N_cool_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).N_cool_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)









    property zeta_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).zeta_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).zeta_tg = conv.array_to_vector_2d_dbl(value)


    property phi_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).phi_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).phi_tg = conv.array_to_vector_2d_dbl(value)


    property phi_t:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).phi_t)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).phi_t = conv.array_to_vector_1d_dbl(value)


    property Phi_t:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Phi_t)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).Phi_t = conv.array_to_vector_1d_dbl(value)


    property BU_t:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).BU_t)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).BU_t = conv.array_to_vector_1d_dbl(value)







    property T_it:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).T_it)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).T_it = conv.dict_to_map_int_array_to_vector_1d_dbl(value)


    property sigma_t_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_t_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_t_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_a_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_a_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_a_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property nubar_sigma_f_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_sigma_f_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_sigma_f_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property chi_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_s_itgh:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_s_itgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_s_itgh = conv.dict_to_map_int_array_to_vector_3d_dbl(value)


    property sigma_f_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_f_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_f_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_gamma_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_2n_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_3n_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_3n_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_3n_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_alpha_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_alpha_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_alpha_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_proton_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_proton_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_proton_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_gamma_x_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_x_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_gamma_x_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)


    property sigma_2n_x_itg:
        def __get__(self):
            return conv.map_to_dict_int_vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_x_itg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).sigma_2n_x_itg = conv.dict_to_map_int_array_to_vector_2d_dbl(value)









    property Sigma_t_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_a_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property nubar_Sigma_f_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property chi_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_s_fuel_tgh:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_fuel_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_fuel_tgh = conv.array_to_vector_3d_dbl(value)


    property Sigma_f_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_3n_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_alpha_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_proton_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_x_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_x_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_fuel_tg = conv.array_to_vector_2d_dbl(value)


    property kappa_fuel_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).kappa_fuel_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).kappa_fuel_tg = conv.array_to_vector_2d_dbl(value)









    property Sigma_t_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_a_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_clad_tg = conv.array_to_vector_2d_dbl(value)


    property nubar_Sigma_f_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_clad_tg = conv.array_to_vector_2d_dbl(value)


    property chi_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_s_clad_tgh:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_clad_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_clad_tgh = conv.array_to_vector_3d_dbl(value)


    property Sigma_f_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_3n_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_alpha_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_proton_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_x_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_clad_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_x_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_clad_tg = conv.array_to_vector_2d_dbl(value)


    property kappa_clad_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).kappa_clad_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).kappa_clad_tg = conv.array_to_vector_2d_dbl(value)










    property Sigma_t_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_a_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_cool_tg = conv.array_to_vector_2d_dbl(value)


    property nubar_Sigma_f_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_cool_tg = conv.array_to_vector_2d_dbl(value)


    property chi_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_s_cool_tgh:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_cool_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_cool_tgh = conv.array_to_vector_3d_dbl(value)


    property Sigma_f_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_3n_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_alpha_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_proton_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_x_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_cool_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_x_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_cool_tg = conv.array_to_vector_2d_dbl(value)


    property kappa_cool_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).kappa_cool_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).kappa_cool_tg = conv.array_to_vector_2d_dbl(value)












    property Sigma_t_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_t_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_a_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_a_tg = conv.array_to_vector_2d_dbl(value)


    property nubar_Sigma_f_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).nubar_Sigma_f_tg = conv.array_to_vector_2d_dbl(value)


    property chi_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).chi_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).chi_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_s_tgh:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_s_tgh = conv.array_to_vector_3d_dbl(value)


    property Sigma_f_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_f_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_3n_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_3n_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_alpha_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_alpha_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_proton_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_proton_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_gamma_x_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_gamma_x_tg = conv.array_to_vector_2d_dbl(value)


    property Sigma_2n_x_tg:
        def __get__(self):
            return conv.vector_to_array_2d_dbl((<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_tg)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).Sigma_2n_x_tg = conv.array_to_vector_2d_dbl(value)









    property A_tgh:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).A_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).A_tgh = conv.array_to_vector_3d_dbl(value)


    property F_tgh:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).F_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).F_tgh = conv.array_to_vector_3d_dbl(value)


    property A_inv_tgh:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).A_inv_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).A_inv_tgh = conv.array_to_vector_3d_dbl(value)


    property A_inv_F_tgh:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).A_inv_F_tgh)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).A_inv_F_tgh = conv.array_to_vector_3d_dbl(value)

    """\
    property T_int_tij:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).T_int_tij)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).T_int_tij = conv.array_to_vector_3d_dbl(value)


    property M_tij:
        def __get__(self):
            return conv.vector_to_array_3d_dbl((<cpp_reactormg.ReactorMG *> self._inst).M_tij)

        def __set__(self, dict value):
            (<cpp_reactormg.ReactorMG *> self._inst).M_tij = conv.array_to_vector_3d_dbl(value)

    """










    property nearest_neighbors:
        def __get__(self):
            return conv.vector_to_array_1d_int((<cpp_reactormg.ReactorMG *> self._inst).nearest_neighbors)

        def __set__(self, np.ndarray[np.int32_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).nearest_neighbors = conv.array_to_vector_1d_int(value)






    property k_t:
        def __get__(self):
            return conv.vector_to_array_1d_dbl((<cpp_reactormg.ReactorMG *> self._inst).k_t)

        def __set__(self, np.ndarray[np.float64_t, ndim=1] value):
            (<cpp_reactormg.ReactorMG *> self._inst).k_t = conv.array_to_vector_1d_dbl(value)





    property td_n:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).td_n

        def __set__(self, int value):
            (<cpp_reactormg.ReactorMG *> self._inst).td_n = value


    property td:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).td

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).td = value


    property BUd:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).BUd

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).BUd = value


    property Phid:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).Phid

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).Phid = value


    property k:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).k

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).k = value






    property mat_feed_u:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_u
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_u = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_feed_tru:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_tru
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_tru = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_feed_lan:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_lan
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_lan = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_feed_act:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_act
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_feed_act = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_u:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_u
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_u = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_tru:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_tru
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_tru = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_lan:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_lan
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_lan = <pyne.cpp_material.Material> mat.mat_pointer[0]


    property mat_prod_act:
        def __get__(self):
            cdef pyne.material._Material pymat = pyne.material.Material()
            pymat.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_act
            return pymat

        def __set__(self, pyne.material._Material mat):
            (<cpp_reactormg.ReactorMG *> self._inst).mat_prod_act = <pyne.cpp_material.Material> mat.mat_pointer[0]





    property deltaR:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).deltaR

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).deltaR = value


    property tru_cr:
        def __get__(self):
            return (<cpp_reactormg.ReactorMG *> self._inst).tru_cr

        def __set__(self, double value):
            (<cpp_reactormg.ReactorMG *> self._inst).tru_cr = value





    #
    # Class Methods
    # 

    def initialize(self, ReactorParameters reactor_parameters):
        """The initialize() method for reactors copies all of the reactor specific parameters to this instance.
        Additionally, it calculates and sets the volumes VF and VC.

        Args:
            * reactor_parameters (ReactorParameters): A special data structure that contains information
              on how to setup and run the reactor.
        """
        cdef ReactorParameters rp = reactor_parameters
        (<cpp_reactormg.ReactorMG *> self._inst).initialize(<cpp_reactor_parameters.ReactorParameters> rp.rp_pointer[0])


    def loadlib(self, char * libfile="reactor.h5"):
        """This method finds the HDF5 library for this reactor and extracts the necessary information from it.
        This method is typically called by the constructor of the child reactor type object.  It must be 
        called before attempting to do any real computation.

        Args: 
            * libfile (str): Path to the reactor library.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).loadlib(std.string(libfile))


    def interpolate_cross_sections(self):
        """This method iterpolates the isotopic, time-dependent cross-sections based on the current 
        state of the burn_time, bt_s, and nearest_neighbors attributes.  It is prudent to call 
        the calc_nearest_neighbors() method before this one.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).interpolate_cross_sections()


    def calc_mass_weights(self):
        """Calculates the mass weights for this time step.  Needed for fold_mass_weights() method.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_mass_weights()


    def fold_mass_weights(self):
        """This method performs the all-important task of doing the isotopically-weighted linear combination of raw data. 
        In a very real sense this is what makes this reactor *this specific reactor*.  The weights are taken 
        as the values of mat_feed.  The raw data must have previously been read in from loadlib().  

        .. warning::

            Anytime any reactor parameter whatsoever (mat_feed, P_NL, *etc*) is altered in any way, 
            the fold_mass_weights() function must be called to reset all of the resultant data.
            If you are unsure, please call this function anyway to be safe.  There is little 
            harm in calling it twice by accident.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).fold_mass_weights()


    def assemble_multigroup_matrices(self):
        """Folds mass weight in with cross-sections for current time step.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).assemble_multigroup_matrices()


    def assemble_transmutation_matrices(self):
        """Calculates the transmutation matrices for the current time step
        """
        (<cpp_reactormg.ReactorMG *> self._inst).assemble_transmutation_matrices()


    def calc_criticality(self):
        """Performs the criticality calculation to find k for this time step.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_criticality()


    def calc_transmutation(self):
        """Burns up the fuel using a matrix exponential method for this time step.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_transmutation()








    def init_core(self):
        """This method generates a time-dependent parameters from an reactor's initial conditions.
        This includes all burnup and criticality calculations.  These time-dependent data
        are then used to determine discharge compositions and other parameters.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).init_core()


    def burnup_core(self):
        """This method generates a time-dependent parameters from an reactor's initial conditions.
        This includes all burnup and criticality calculations.  These time-dependent data
        are then used to determine discharge compositions and other parameters.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).burnup_core()






    def calc_nearest_neighbors(self):
        """Calculates a sorted array that indexes the nearest neighbors of the 
        perturbations based off of the current state of the reactor.  The results may
        be found in the neareest_neighbors attribute.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_nearest_neighbors()







    def calc_T_itd(self):
        """This function evaluates transmutation matrix at the discharge time td.
        The resultant isotopic dictionary is then converted into the mat_prod mass stream
        for this pass through the reactor.  Thus if ever you need to calculate mat_prod
        without going through calc(), use this function.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_T_itd()






    def calc_mat_prod(self):
        """This is a convenience function that wraps the transmutation matrix methods.  
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calc_mat_prod()

    def calcSubStreams(self):
        """This sets possibly relevant reactor input and output substreams.  Specifically, it calculates the 
        attributes:

            * mat_feed_u
            * mat_feed_tru
            * mat_feed_lan
            * mat_feed_act
            * mat_prod_u
            * mat_prod_tru
            * mat_prod_lan
            * mat_prod_act

        """
        (<cpp_reactormg.ReactorMG *> self._inst).calcSubStreams()


    def calc_tru_cr(self):
        """This calculates and sets the transuranic conversion ratio tru_cr through the equation:

        .. math:: \mbox{tru_cr} = \frac{\mbox{mat_feed_tru.mass} - \mbox{mat_prod_tru.mass}}{\frac{\mbox{BUd}}{935.0}}

        Returns:
            * tru_cr (float): The value of the transuranic conversion ratio just calculated.
        """
        return (<cpp_reactormg.ReactorMG *> self._inst).calc_tru_cr()







    def fluence_at_BU(self, double burnup):
        """This function takes a burnup value  and returns a special fluence point object.  
        The fluence point is an amalgamation of data where the at which the burnup occurs.
        This object instance FP contains three pieces of information::
    
            FP.f    #Index immediately lower than where BU achieved (int)
            FP.F    #Fluence value itself (float)
            FP.m    #Slope dBU/dF between points f and f+1 (double)

        Args:
            * burnup (float): Burnup [MWd/kgIHM] at which to calculate the corresponding fluence.

        Returns:
            * fp (FluencePoint): A class containing fluence information.
        """
        cdef FluencePoint fp = FluencePoint()
        fp.fp_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).fluence_at_BU(burnup)
        return fp


    def batch_average_k(self, double BUd):
        """Finds the batch average k(F) when at discharge burnup BUd.
        This function is typically iterated over until a BUd is found such that k(F) = 1.0 + err.

        Args:
            * BUd (float): The discharge burnup [MWd/kgIHM] to obtain a batch-averaged value for.

        Returns:
            * k (float): the batch averaged multiplication factor.
        """
        cdef double k = (<cpp_reactormg.ReactorMG *> self._inst).batch_average_k(BUd)
        return k


    def BUd_bisection_method(self):
        """Calculates the maximum discharge burnup via the Bisection Method for a given mat_feed
        in this reactor.  This iterates over values of BUd to find a batch averaged multiplication factor 
        that is closest to 1.0.

        Other root finding methods for determining maximum discharge burnup are certainly possible.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).BUd_bisection_method()


    def run_P_NL(self, double pnl):
        """Performs a reactor run for a specific non-leakage probability value.
        This requires that mat_feed be (meaningfully) set and is for use with calibrate_P_NL_to_BUd().

        This function amounts to the following code::

            self.P_NL = pnl
            self.fold_mass_weights()
            self.BUd_bisection_method()

        Args:
            * pnl (float): The new non-leakage probability for the reactor.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).run_P_NL(pnl)
    

    def calibrate_P_NL_to_BUd(self):
        """Often times the non-leakage probability of a reactor is not known, though the input isotopics 
        and the target discharge burnup are.  This function handles that situation by
        calibrating the non-leakage probability of this reactor P_NL to hit its target burnup target_BU.
        Such a calibration proceeds by bisection method as well.  This function is extremely useful for 
        benchmarking calculations.
        """
        (<cpp_reactormg.ReactorMG *> self._inst).calibrate_P_NL_to_BUd()




    def calc(self, input=None):
        """Since many other methods provide the computational heavy-lifting of reactor calculations, 
        the calc() method is relatively simple::

            self.mat_feed = input
            self.fold_mass_weights()
            self.BUd_bisection_method()
            self.calc_mat_prod()
            return self.mat_prod

        As you can see, all this function does is set burn an input stream to its maximum 
        discharge burnup and then reports on the output isotopics.

        Args:
            * input (dict or Material): If input is present, it set as the component's 
              mat_feed.  If input is a isotopic dictionary (zzaaam keys, float values), this
              dictionary is first converted into a Material before being set as mat_feed.

        Returns:
            * output (Material): mat_prod.
        """
        cdef pyne.material._Material in_mat 
        cdef pyne.material._Material output = pyne.material.Material()

        if input is None:
            output.mat_pointer[0] = (<cpp_reactormg.FCComp *> self._inst).calc()
        elif isinstance(input, dict):
            output.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).calc(conv.dict_to_map_int_dbl(input))
        elif isinstance(input, pyne.material._Material):
            in_mat = input
            output.mat_pointer[0] = (<cpp_reactormg.ReactorMG *> self._inst).calc(<pyne.cpp_material.Material> in_mat.mat_pointer[0])

        return output
