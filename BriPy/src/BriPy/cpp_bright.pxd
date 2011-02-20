"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

cimport std
cimport cpp_mass_stream
cimport mass_stream

cdef extern from "../bright.h" namespace "bright":
    std.string BRIGHT_DATA

    void BrightStart()


cdef extern from "../FCComp.h" namespace "FCComps":
    set[int] isos2track
    void load_isos2track_hdf5(std.string, std.string, bint)
    void load_isos2track_text(std.string, bint)

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

        cpp_mass_stream.MassStream IsosIn
        cpp_mass_stream.MassStream IsosOut

        map[std.string, double] ParamsIn
        map[std.string, double] ParamsOut

        int PassNum
        set[std.string] params2track

        # Methods
        void setParams()
        void writeIsoPass()
        void writeParamPass()
        void writeText()
        void writeHDF5()
        void writeout()
        cpp_mass_stream.MassStream doCalc()


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
        cpp_mass_stream.MassStream IsosTail

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
        void setParams ()
        cpp_mass_stream.MassStream doCalc ()
        cpp_mass_stream.MassStream doCalc (map[int, double])
        cpp_mass_stream.MassStream doCalc (cpp_mass_stream.MassStream)

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
        void setParams()
        cpp_mass_stream.MassStream doCalc()
        cpp_mass_stream.MassStream doCalc(map[int, double])
        cpp_mass_stream.MassStream doCalc(cpp_mass_stream.MassStream)


cdef extern from "../Storage.h":

    cdef cppclass Storage(FCComp):
        # Constructors
        Storage()
        Storage(std.string)

        # Attributes
        double decay_time

        # Methods
        void setParams()
        cpp_mass_stream.MassStream doCalc()
        cpp_mass_stream.MassStream doCalc(map[int, double])
        cpp_mass_stream.MassStream doCalc(cpp_mass_stream.MassStream)
        cpp_mass_stream.MassStream doCalc(double)
        cpp_mass_stream.MassStream doCalc(map[int, double], double)
        cpp_mass_stream.MassStream doCalc(cpp_mass_stream.MassStream, double)



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
        map[std.string, double] FuelForm
        map[std.string, double] CoolantForm
        double FuelDensity
        double CoolantDensity
        double pnl
        double BUt
        bint useDisadvantage
        std.string LatticeType
        bint HydrogenRescale
        double Radius
        double Length
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
        map[std.string, double] FuelChemicalForm
        map[std.string, double] CoolantChemicalForm
        double rhoF
        double rhoC
        double P_NL
        double TargetBU
        bint useZeta
        std.string Lattice
        bint H_XS_Rescale

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

        cpp_mass_stream.MassStream InU
        cpp_mass_stream.MassStream InTRU
        cpp_mass_stream.MassStream InLAN
        cpp_mass_stream.MassStream InACT
        cpp_mass_stream.MassStream OutU
        cpp_mass_stream.MassStream OutTRU
        cpp_mass_stream.MassStream OutLAN
        cpp_mass_stream.MassStream OutACT

        double deltaR
        double TruCR

        vector[double] SigmaFa_F_
        vector[double] SigmaFtr_F_
        vector[double] kappaF_F_

        vector[double] SigmaCa_F_
        vector[double] SigmaCtr_F_
        vector[double] kappaC_F_

        vector[double] LatticeE_F_
        vector[double] LatticeF_F_

        # Methods
        void initialize(ReactorParameters)
        void loadLib(std.string)
        void foldMassWeights()

        void mkMj_F_()
        void mkMj_Fd_()

        void calcOutIso()
        void calcSubStreams()
        double calcTruCR()

        double calc_deltaR()
        double calc_deltaR(map[int, double])
        double calc_deltaR(cpp_mass_stream.MassStream)

        FluencePoint FluenceAtBU(double)
        double batchAve(double, std.string)
        double batchAveK(double)
        void BUd_BisectionMethod()
        void Run_PNL(double)
        void Calibrate_PNL_2_BUd()

        cpp_mass_stream.MassStream doCalc()
        cpp_mass_stream.MassStream doCalc(map[int, double])
        cpp_mass_stream.MassStream doCalc(cpp_mass_stream.MassStream)

        void LatticeEPlanar(double, double)
        void LatticeFPlanar(double, double)
        void LatticeESpherical(double, double)
        void LatticeFSpherical(double, double)
        void LatticeECylindrical(double, double)
        void LatticeFCylindrical(double, double)

        void calcZeta()
        void calcZetaPlanar()
        void calcZetaSpherical()
        void calcZetaCylindrical()




cdef extern from "../LightWaterReactor1G.h":

    ReactorParameters fillLWRDefaults()

    cdef cppclass LightWaterReactor1G(Reactor1G):
        # Constructors        
        LightWaterReactor1G()
        LightWaterReactor1G(std.string, std.string)
        LightWaterReactor1G(ReactorParameters, std.string)
        LightWaterReactor1G(std.string, ReactorParameters, std.string)

        # Methods
        void setParams()




cdef extern from "../FastReactor1G.h":

    ReactorParameters fillFRDefaults()

    cdef cppclass FastReactor1G(Reactor1G):
        # Constructors        
        FastReactor1G()
        FastReactor1G(std.string, std.string)
        FastReactor1G(ReactorParameters, std.string)
        FastReactor1G(std.string, ReactorParameters, std.string)

        # Methods
        void setParams()



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
        void setParams()

        void calc_deltaRs()
        cpp_mass_stream.MassStream calc_core_input()
        void calc_mass_ratios()

        cpp_mass_stream.MassStream doCalc()
        cpp_mass_stream.MassStream doCalc(map[std.string, mass_stream.msp], map[std.string, double], Reactor1G)
