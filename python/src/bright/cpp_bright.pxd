"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

cimport std
cimport cpp_mass_stream
cimport mass_stream

cdef extern from "../bright.h" namespace "bright":
    std.string BRIGHT_DATA

    void bright_start()


cdef extern from "../FCComp.h" namespace "FCComps":
    set[int] track_isos
    void load_track_isos_hdf5(std.string, std.string, bint)
    void load_track_isos_text(std.string, bint)

    int verbosity
    bint write_hdf5
    bint write_text

    std.string output_filename


cdef extern from "../FCComp.h":
    cdef cppclass FCComp:
        # Constructors
        FCComp()
        FCComp(std.string)
        FCComp(set[std.string], std.string)

        # Attributes
        std.string name 
        std.string natural_name

        cpp_mass_stream.MassStream ms_feed
        cpp_mass_stream.MassStream ms_prod

        map[std.string, double] params_prior_calc
        map[std.string, double] params_after_calc

        int pass_num
        set[std.string] track_params

        # Methods
        void calc_params()
        void write_ms_pass()
        void write_params_pass()
        void write_text()
        void write_hdf5()
        void write()
        cpp_mass_stream.MassStream calc()


cdef extern from "../Enrichment.h":

    cdef cppclass EnrichmentParameters:
        # Constructors
        EnrichmentParameters()

        # Attributes
        double alpha_0
        double Mstar_0

        int j
        int k

        double N0
        double M0

        double xP_j
        double xW_j

    EnrichmentParameters fillUraniumEnrichmentDefaults()

    cdef cppclass Enrichment(FCComp): 
        # Constructors
        Enrichment()
        Enrichment(std.string)
        Enrichment(EnrichmentParameters, std.string)

        # Attributes
        double alpha_0
        double Mstar_0
        double Mstar
        cpp_mass_stream.MassStream ms_tail

        int j
        int k
        double xP_j
        double xW_j

        double N
        double M
        double N0
        double M0

        double TotalPerFeed
        double SWUperFeed
        double SWUperProduct

        # Methods
        void initialize(EnrichmentParameters)
        void calc_params ()
        cpp_mass_stream.MassStream calc ()
        cpp_mass_stream.MassStream calc (map[int, double])
        cpp_mass_stream.MassStream calc (cpp_mass_stream.MassStream)

        double PoverF (double, double, double)
        double WoverF (double, double, double)

        # The following are methods I am too lazy to expose to Python
        # FIXME
        #double get_alphastar_i (double)

        #double get_Ei (double)
        #double get_Si (double)
        #void FindNM()

        #double xP_i(int)
        #double xW_i(int)
        #void SolveNM()
        #void Comp2UnitySecant()
        #void Comp2UnityOther()
        #double deltaU_i_OverG(int)
        #void LoverF()
        #void MstarOptimize()



cdef extern from "../Reprocess.h":

    cdef cppclass Reprocess(FCComp):
        # Constructors
        Reprocess()
        Reprocess(map[int, double], std.string)

        # Attributes
        map[int, double] sepeff

        # Methods
        void initialize(map[int, double])
        void calc_params()
        cpp_mass_stream.MassStream calc()
        cpp_mass_stream.MassStream calc(map[int, double])
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream)


cdef extern from "../Storage.h":

    cdef cppclass Storage(FCComp):
        # Constructors
        Storage()
        Storage(std.string)

        # Attributes
        double decay_time

        # Methods
        void calc_params()
        cpp_mass_stream.MassStream calc()
        cpp_mass_stream.MassStream calc(map[int, double])
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream)
        cpp_mass_stream.MassStream calc(double)
        cpp_mass_stream.MassStream calc(map[int, double], double)
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream, double)



cdef extern from "../Reactor1G.h":

    cdef cppclass FluencePoint:
        # Constructors        
        FluencePoint()

        # Attributes
        int f
        double F
        double m


    cdef cppclass ReactorParameters:
        # Constructors        
        ReactorParameters()

        # Attributes
        int batches
        double flux
        map[std.string, double] fuel_form
        map[std.string, double] coolant_form
        double fuel_density
        double coolant_density
        double pnl
        double BUt
        bint use_disadvantage_factor
        std.string lattice_type
        bint rescale_hydrogen
        double radius
        double pitch
        double open_slots
        double total_slots


    cdef cppclass Reactor1G(FCComp):
        # Constructors        
        Reactor1G()
        Reactor1G(std.string)
        Reactor1G(set[std.string], std.string)
        Reactor1G(ReactorParameters, std.string)
        Reactor1G(ReactorParameters, set[std.string], std.string)

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

        cpp_mass_stream.MassStream ms_feed_u
        cpp_mass_stream.MassStream ms_feed_tru
        cpp_mass_stream.MassStream ms_feed_lan
        cpp_mass_stream.MassStream ms_feed_act
        cpp_mass_stream.MassStream ms_prod_u
        cpp_mass_stream.MassStream ms_prod_tru
        cpp_mass_stream.MassStream ms_prod_lan
        cpp_mass_stream.MassStream ms_prod_act

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
        void initialize(ReactorParameters)
        void loadlib(std.string)
        void fold_mass_weights()

        void calc_Mj_F_()
        void calc_Mj_Fd_()

        void calc_ms_prod()
        void calcSubStreams()
        double calc_tru_cr()

        double calc_deltaR()
        double calc_deltaR(map[int, double])
        double calc_deltaR(cpp_mass_stream.MassStream)

        FluencePoint fluence_at_BU(double)
        double batch_average(double, std.string)
        double batch_average_k(double)
        void BUd_bisection_method()
        void run_P_NL(double)
        void calibrate_P_NL_to_BUd()

        cpp_mass_stream.MassStream calc()
        cpp_mass_stream.MassStream calc(map[int, double])
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream)

        void lattice_E_planar(double, double)
        void lattice_F_planar(double, double)
        void lattice_E_spherical(double, double)
        void lattice_F_spherical(double, double)
        void lattice_E_cylindrical(double, double)
        void lattice_F_cylindrical(double, double)

        void calc_zeta()
        void calc_zeta_planar()
        void calc_zeta_spherical()
        void calc_zeta_cylindrical()




cdef extern from "../LightWaterReactor1G.h":

    ReactorParameters filllwr_defaults()

    cdef cppclass LightWaterReactor1G(Reactor1G):
        # Constructors        
        LightWaterReactor1G()
        LightWaterReactor1G(std.string, std.string)
        LightWaterReactor1G(ReactorParameters, std.string)
        LightWaterReactor1G(std.string, ReactorParameters, std.string)

        # Methods
        void calc_params()




cdef extern from "../FastReactor1G.h":

    ReactorParameters fillfr_defaults()

    cdef cppclass FastReactor1G(Reactor1G):
        # Constructors        
        FastReactor1G()
        FastReactor1G(std.string, std.string)
        FastReactor1G(ReactorParameters, std.string)
        FastReactor1G(std.string, ReactorParameters, std.string)

        # Methods
        void calc_params()



cdef extern from "../FuelFabrication.h":

    cdef cppclass FuelFabrication(FCComp):
        # Constructors        
        FuelFabrication()
        FuelFabrication(std.string)
        FuelFabrication(set[std.string], std.string)
        FuelFabrication(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G, std.string)
        FuelFabrication(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G, set[std.string], std.string)

        # Attributes
        map[std.string, mass_stream.msp] mass_streams
        map[std.string, double] mass_weights_in
        map[std.string, double] mass_weights_out
        map[std.string, double] deltaRs

        Reactor1G reactor

        # Methods
        void initialize(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G)
        void calc_params()

        void calc_deltaRs()
        cpp_mass_stream.MassStream calc_core_input()
        void calc_mass_ratios()

        cpp_mass_stream.MassStream calc()
        cpp_mass_stream.MassStream calc(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G)
