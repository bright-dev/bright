################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################


cimport pyne.stlcontainers
from bright cimport cpp_fccomp
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from libcpp.string cimport string as std_string
from pyne cimport cpp_material
from pyne cimport material



cdef class FCComp:
    cdef void * _inst
    cdef public bint _free_inst
    cdef public material._Material _mat_feed
    cdef public material._Material _mat_prod
    cdef public pyne.stlcontainers._MapStrDouble _params_after_calc
    cdef public pyne.stlcontainers._MapStrDouble _params_prior_calc
    cdef public pyne.stlcontainers._SetStr _track_params
    pass    




