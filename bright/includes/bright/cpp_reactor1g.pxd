"""Cython header for reactor1g library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

from pyne cimport std
from pyne cimport cpp_material
from pyne cimport material

from bright.cpp_fccomp cimport FCComp
from bright.cpp_reactor_parameters cimport ReactorParameters
from bright.cpp_fluence_point cimport FluencePoint


cdef extern from "reactor1g.h" namespace "bright":

    cdef cppclass Reactor1G(FCComp):
        # Constructors        
        Reactor1G() except +
        Reactor1G(std.string) except +
        Reactor1G(set[std.string], std.string) except +
        Reactor1G(ReactorParameters, std.string) except +
        Reactor1G(ReactorParameters, set[std.string], std.string) except +

        # Attributes
        int B
        double phi
        map[std.string, double] fuel_chemical_form
        map[std.string, double] coolant_chemical_form
        double rhoF
        double rhoC
        double P_NL
        double target_BU
        bint use_zeta
        std.string lattice_flag
        bint rescale_hydrogen_xs

        double r
        double l
        double S_O
        double S_T
        double VF
        double VC

        std.string libfile
        vector[double] F
        map[int, vector[double]] BUi_F_
        map[int, vector[double]] pi_F_
        map[int, vector[double]] di_F_
        map[int, map[int, vector[double]]] Tij_F_

        double A_IHM
        double MWF
        double MWC
        map[int, double] niF
        map[int, double] niC
        map[int, double] miF
        map[int, double] miC
        map[int, double] NiF
        map[int, double] NiC

        vector[double] dF_F_
        vector[double] dC_F_
        vector[double] BU_F_
        vector[double] P_F_
        vector[double] D_F_
        vector[double] k_F_
        map[int, vector[double]] Mj_F_
        vector[double] zeta_F_

        int fd
        double Fd
        double BUd
        double k

        cpp_material.Material mat_feed_u
        cpp_material.Material mat_feed_tru
        cpp_material.Material mat_feed_lan
        cpp_material.Material mat_feed_act
        cpp_material.Material mat_prod_u
        cpp_material.Material mat_prod_tru
        cpp_material.Material mat_prod_lan
        cpp_material.Material mat_prod_act

        double deltaR
        double tru_cr

        vector[double] SigmaFa_F_
        vector[double] SigmaFtr_F_
        vector[double] kappaF_F_

        vector[double] SigmaCa_F_
        vector[double] SigmaCtr_F_
        vector[double] kappaC_F_

        vector[double] lattice_E_F_
        vector[double] lattice_F_F_

        # Methods
        void initialize(ReactorParameters) except +
        void loadlib(std.string) except +
        void fold_mass_weights() except +

        void calc_Mj_F_() except +
        void calc_Mj_Fd_() except +

        void calc_mat_prod() except +
        void calc_sub_mats() except +
        double calc_tru_cr() except +

        double calc_deltaR() except +
        double calc_deltaR(map[int, double]) except +
        double calc_deltaR(cpp_material.Material) except +

        FluencePoint fluence_at_BU(double) except +
        double batch_average(double, std.string) except +
        double batch_average_k(double) except +
        void BUd_bisection_method() except +
        void run_P_NL(double) except +
        void calibrate_P_NL_to_BUd() except +

        cpp_material.Material calc() except +
        cpp_material.Material calc(map[int, double]) except +
        cpp_material.Material calc(cpp_material.Material) except +

        void lattice_E_planar(double, double) except +
        void lattice_F_planar(double, double) except +
        void lattice_E_spherical(double, double) except +
        void lattice_F_spherical(double, double) except +
        void lattice_E_cylindrical(double, double) except +
        void lattice_F_cylindrical(double, double) except +

        void calc_zeta() except +
        void calc_zeta_planar() except +
        void calc_zeta_spherical() except +
        void calc_zeta_cylindrical() except +


