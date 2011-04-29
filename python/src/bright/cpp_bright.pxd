"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

cimport std
cimport cpp_mass_stream
cimport mass_stream

cdef extern from "../../../cpp/bright.h" namespace "bright":
    std.string BRIGHT_DATA

    void bright_start() except +


cdef extern from "../../../cpp/FCComp.h" namespace "FCComps":
    set[int] track_isos
    void load_track_isos_hdf5(std.string, std.string, bint) except +
    void load_track_isos_text(std.string, bint) except +

    int verbosity
    bint write_hdf5
    bint write_text

    std.string output_filename


cdef extern from "../../../cpp/FCComp.h":
    cdef cppclass FCComp:
        # Constructors
        FCComp() except +
        FCComp(std.string) except +
        FCComp(set[std.string], std.string) except +

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
        void calc_params() except +
        void write_ms_pass() except +
        void write_params_pass() except +
        void write_text() except +
        void write_hdf5() except +
        void write() except +
        cpp_mass_stream.MassStream calc() except +


cdef extern from "../../../cpp/Enrichment.h":

    cdef cppclass EnrichmentParameters:
        # Constructors
        EnrichmentParameters() except +

        # Attributes
        double alpha_0
        double Mstar_0

        int j
        int k

        double N0
        double M0

        double xP_j
        double xW_j

    EnrichmentParameters fillUraniumEnrichmentDefaults() except +

    cdef cppclass Enrichment(FCComp): 
        # Constructors
        Enrichment() except +
        Enrichment(std.string) except +
        Enrichment(EnrichmentParameters, std.string) except +

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
        void initialize(EnrichmentParameters) except +
        void calc_params () except +
        cpp_mass_stream.MassStream calc () except +
        cpp_mass_stream.MassStream calc (map[int, double]) except +
        cpp_mass_stream.MassStream calc (cpp_mass_stream.MassStream) except +

        double PoverF (double, double, double) except +
        double WoverF (double, double, double) except +

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



cdef extern from "../../../cpp/Reprocess.h":

    cdef cppclass Reprocess(FCComp):
        # Constructors
        Reprocess() except +
        Reprocess(map[int, double], std.string) except +

        # Attributes
        map[int, double] sepeff

        # Methods
        void initialize(map[int, double]) except +
        void calc_params() except +
        cpp_mass_stream.MassStream calc() except +
        cpp_mass_stream.MassStream calc(map[int, double]) except +
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream) except +


cdef extern from "../../../cpp/Storage.h":

    cdef cppclass Storage(FCComp):
        # Constructors
        Storage() except +
        Storage(std.string) except +

        # Attributes
        double decay_time

        # Methods
        void calc_params() except +
        cpp_mass_stream.MassStream calc() except +
        cpp_mass_stream.MassStream calc(map[int, double]) except +
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream) except +
        cpp_mass_stream.MassStream calc(double) except +
        cpp_mass_stream.MassStream calc(map[int, double], double) except +
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream, double) except +



cdef extern from "../../../cpp/Reactor1G.h":

    cdef cppclass FluencePoint:
        # Constructors        
        FluencePoint() except +

        # Attributes
        int f
        double F
        double m


    cdef cppclass ReactorParameters:
        # Constructors        
        ReactorParameters() except +

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
        void initialize(ReactorParameters) except +
        void loadlib(std.string) except +
        void fold_mass_weights() except +

        void calc_Mj_F_() except +
        void calc_Mj_Fd_() except +

        void calc_ms_prod() except +
        void calcSubStreams() except +
        double calc_tru_cr() except +

        double calc_deltaR() except +
        double calc_deltaR(map[int, double]) except +
        double calc_deltaR(cpp_mass_stream.MassStream) except +

        FluencePoint fluence_at_BU(double) except +
        double batch_average(double, std.string) except +
        double batch_average_k(double) except +
        void BUd_bisection_method() except +
        void run_P_NL(double) except +
        void calibrate_P_NL_to_BUd() except +

        cpp_mass_stream.MassStream calc() except +
        cpp_mass_stream.MassStream calc(map[int, double]) except +
        cpp_mass_stream.MassStream calc(cpp_mass_stream.MassStream) except +

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




cdef extern from "../../../cpp/LightWaterReactor1G.h":

    ReactorParameters filllwr_defaults() except +

    cdef cppclass LightWaterReactor1G(Reactor1G):
        # Constructors        
        LightWaterReactor1G() except +
        LightWaterReactor1G(std.string, std.string) except +
        LightWaterReactor1G(ReactorParameters, std.string) except +
        LightWaterReactor1G(std.string, ReactorParameters, std.string) except +

        # Methods
        void calc_params() except +




cdef extern from "../../../cpp/FastReactor1G.h":

    ReactorParameters fillfr_defaults() except +

    cdef cppclass FastReactor1G(Reactor1G):
        # Constructors        
        FastReactor1G() except +
        FastReactor1G(std.string, std.string) except +
        FastReactor1G(ReactorParameters, std.string) except +
        FastReactor1G(std.string, ReactorParameters, std.string) except +

        # Methods
        void calc_params() except +



cdef extern from "../../../cpp/FuelFabrication.h":

    cdef cppclass FuelFabrication(FCComp):
        # Constructors        
        FuelFabrication() except +
        FuelFabrication(std.string) except +
        FuelFabrication(set[std.string], std.string) except +
        FuelFabrication(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G, std.string) except +
        FuelFabrication(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G, set[std.string], std.string) except +

        # Attributes
        map[std.string, mass_stream.msp] mass_streams
        map[std.string, double] mass_weights_in
        map[std.string, double] mass_weights_out
        map[std.string, double] deltaRs

        Reactor1G reactor

        # Methods
        void initialize(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G) except +
        void calc_params() except +

        void calc_deltaRs() except +
        cpp_mass_stream.MassStream calc_core_input() except +
        void calc_mass_ratios() except +

        cpp_mass_stream.MassStream calc() except +
        cpp_mass_stream.MassStream calc(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G) except +
