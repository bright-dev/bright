// One Group Reactor Component Class

#include "Reactor1G.h"

/**************************/
/*** FulencePoint Class ***/
/**************************/

FluencePoint::FluencePoint()
{
    f = 0;
    F = 0.0;
    m = 0.0;
};


FluencePoint::~FluencePoint()
{
};


/*******************************/
/*** ReactorParameters Class ***/
/*******************************/

ReactorParameters::ReactorParameters()
{
    batches = 0;
    flux = 0.0;
    FuelForm = std::map<std::string, double>();
    CoolantForm = std::map<std::string, double>();
    FuelDensity = 0.0;
    CoolantDensity = 0.0;
    pnl = 0.0;
    BUt = 0.0;
    useDisadvantage = false;
    LatticeType = std::string();
    HydrogenRescale = false;
    Radius = 0.0;
    Length = 0.0;
    open_slots = 0.0;
    total_slots = 0.0;
};


ReactorParameters::~ReactorParameters()
{
};

/***********************************************/
/*** Reactor1G Component Class and Functions ***/
/***********************************************/

Reactor1G::Reactor1G()
{
};

Reactor1G::Reactor1G(std::string n) : FCComp(n)
{
};

Reactor1G::Reactor1G(std::set<std::string> paramtrack, std::string n) : FCComp(paramtrack, n)
{
};

Reactor1G::Reactor1G(ReactorParameters rp, std::string n) : FCComp(n)
{
    initialize(rp);
};

Reactor1G::Reactor1G(ReactorParameters rp, std::set<std::string> paramtrack, std::string n) : FCComp(paramtrack, n)
{
    initialize(rp);
};

Reactor1G::~Reactor1G()
{
};

void Reactor1G::initialize(ReactorParameters rp)
{
    /** Sets reactor specific parameters.
     *  Must be done once at the beginning of reactor object life.
     */

    B = rp.batches;				//Total number of fuel loading batches
    phi = rp.flux;				//Flux used for Fluence
    FuelChemicalForm = rp.FuelForm;		//Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
    CoolantChemicalForm = rp.CoolantForm;	//Same a fuel chemical form but for coolant.  Should not have "IHM"
    rhoF = rp.FuelDensity;			//Fuel Density
    rhoC = rp.CoolantDensity;		//Coolant Density
    P_NL = rp.pnl;				//Non-Leakage Probability
    TargetBU = rp.BUt;			//Target Discharge Burnup, only used for graphing inside of this component
    useZeta = rp.useDisadvantage;		//Boolean value on whether or not the disadvantage factor should be used
    Lattice = rp.LatticeType;		//Lattice Type (Planar || Spherical || Cylindrical)
    H_XS_Rescale = rp.HydrogenRescale;	//Rescale the Hydrogen-1 XS?

    //Calculates Volumes
    r = rp.Radius;			//Fuel region radius
    l = rp.Length;			//Unit cell side length
    S_O = rp.open_slots;		//Number of open slots in fuel assembly
    S_T = rp.total_slots;		//Total number of Fuel assembly slots.
    //Fuel Volume
    VF = ((bright::pi*r*r)/(l*l)) * (1.0 - S_O/S_T); 
    //Coolant Volume
    VC = ((l*l - bright::pi*r*r)/(l*l)) * (1.0 - S_O/S_T) + (S_O/S_T);
};

void Reactor1G::loadLib(std::string libfile)
{
    //Loads Apporiate Libraries for Reactor and makes them into Burnup Parameters [F, pi(F), di(F), BUi(F), Tij(F)].

    //HDF5 types
    hid_t  rlib;
    herr_t rstat;

    rlib = H5Fopen (libfile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);		//Recator Library

    //Initializes Burnup Parameters...
    hsize_t dimFromIso[1];
    hsize_t dimToIso[1];

    rstat = H5LTget_dataset_info(rlib, "/FromIso_zz", dimFromIso, NULL, NULL);
    rstat = H5LTget_dataset_info(rlib, "/ToIso_zz",   dimToIso,   NULL, NULL);

    #ifdef _WIN32
        int * FromIso;
        int * ToIso;

        FromIso = new int [dimFromIso[0]];
        ToIso   = new int [dimToIso[0]];
    #else
        int FromIso [dimFromIso[0]];
        int ToIso   [dimToIso[0]];
    #endif

    rstat = H5LTread_dataset_int(rlib, "/FromIso_zz", FromIso);		
    rstat = H5LTread_dataset_int(rlib, "/ToIso_zz",   ToIso);		

    I.clear();
    I.insert(&FromIso[0], &FromIso[dimFromIso[0]]);
    J.clear();
    J.insert(&ToIso[0],   &ToIso[dimToIso[0]]);
    
    //Get Fluence Vector
    hsize_t dimsF[1];							//Read in number of data points
    rstat = H5LTget_dataset_info(rlib, "/Fluence", dimsF, NULL, NULL);
    int lenF = dimsF[0];

    //Make temp array
    #ifdef _WIN32
        float * tempF;
        tempF = new float [lenF];
    #else
        float tempF [lenF];
    #endif

    rstat = H5LTread_dataset_float(rlib, "/Fluence", tempF);		
    F.assign(&tempF[0], &tempF[lenF]);					//Fluence in [n/kb]

    for (IsoIter i = I.begin(); i != I.end(); i++ )
    {
        std::string iso = isoname::zzaaam_2_LLAAAM(*i);

        //Build BUi_F_
        #ifdef _WIN32
            float * tempBUi;
            tempBUi = new float [lenF];
        #else
            float tempBUi [lenF];
        #endif
        rstat = H5LTread_dataset_float(rlib, ("/Burnup/" + iso).c_str(), tempBUi);		
        BUi_F_[*i].assign(&tempBUi[0], &tempBUi[lenF]);

        //Build pi_F_
        #ifdef _WIN32
            float * temppi;
            temppi = new float [lenF];
        #else
            float temppi [lenF];
        #endif
        rstat = H5LTread_dataset_float(rlib, ("/Production/" + iso).c_str(), temppi);		
        pi_F_[*i].assign(&temppi[0], &temppi[lenF]);
        pi_F_[*i][0] = bright::SolveLine(0.0, F[2], pi_F_[*i][2], F[1], pi_F_[*i][1]);

        //Build di_F_
        #ifdef _WIN32
            float * tempdi;
            tempdi = new float [lenF];
        #else
            float tempdi [lenF];
        #endif
        rstat = H5LTread_dataset_float(rlib, ("/Destruction/" + iso).c_str(), tempdi);		
        di_F_[*i].assign(&tempdi[0], &tempdi[lenF]);
        di_F_[*i][0] = bright::SolveLine(0.0, F[2], di_F_[*i][2], F[1], di_F_[*i][1]);
        
        //Build Tij_F_
        for (int jn = 0; jn < dimToIso[0] ; jn++)
        {
            int j = ToIso[jn];
            std::string jso = isoname::zzaaam_2_LLAAAM(j);

            #ifdef _WIN32
                float * tempTij;
                tempTij = new float [lenF];
            #else
                float tempTij [lenF];
            #endif
            rstat = H5LTread_dataset_float(rlib, ("/Transmutation/" + iso + "/" + jso).c_str(), tempTij);
            Tij_F_[*i][j].assign(&tempTij[0], &tempTij[lenF]);
        };
    };
    rstat = H5Fclose(rlib);

    //Now get microscopic XS data from KAERI...
    //...But only if the disadvantage factor is used.
    if (!useZeta)
        return;

    //HDF5 types
    hid_t  kdblib;			//KaeriData.h5 file reference
    herr_t kdbstat;			//File status

    hsize_t xs_nfields, xs_nrows; 	//Number of rows and fields (named columns) in XS table

    //open the file
    kdblib = H5Fopen ( (bright::BRIGHT_DATA + "/KaeriData.h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);	//KAERI Data Library

    //Get Thermal Mawell Average Table & Field Data Dimensions 
    kdbstat = H5TBget_table_info(kdblib, "/XS/ThermalMaxwellAve", &xs_nfields, &xs_nrows);

    //Creating an empy array of character strings is tricky,
    //because character strings are arrays themselves!
    char ** xs_field_names = new char * [xs_nfields];
    for (int n = 0; n < xs_nfields; n++)
        xs_field_names[n] = new char [50]; 

    #ifdef _WIN32
        size_t * xs_field_sizes;
        size_t * xs_field_offsets;

        xs_field_sizes   = new size_t [xs_nfields];
        xs_field_offsets = new size_t [xs_nfields];
    #else
        size_t xs_field_sizes   [xs_nfields];
        size_t xs_field_offsets [xs_nfields];
    #endif

    size_t xs_type_size;
    kdbstat = H5TBget_field_info(kdblib, "/XS/ThermalMaxwellAve", xs_field_names, xs_field_sizes, xs_field_offsets, &xs_type_size);

    //Read the "isozz" column so that we can inteligently pick out our data
    int isozz_n = bright::find_index_char( (char *) "isozz", xs_field_names, xs_nfields);
    int * isozz = new int [xs_nrows];
    #ifdef _WIN32
        const size_t temp_xs_field_sizes_isozz_n [1] = {xs_field_sizes[isozz_n]};
        kdbstat = H5TBread_fields_name(kdblib, "/XS/ThermalMaxwellAve", xs_field_names[isozz_n], 0, xs_nrows, sizeof(int), 0, temp_xs_field_sizes_isozz_n, isozz);
    #else
        kdbstat = H5TBread_fields_name(kdblib, "/XS/ThermalMaxwellAve", xs_field_names[isozz_n], 0, xs_nrows, sizeof(int), 0, (const size_t [1]) {xs_field_sizes[isozz_n]}, isozz);
    #endif

    //Now, load the XS that we need.
    //NOTE: This maps metastable isotopes to their stable versions if they can't be found!
    int sigma_a_n = bright::find_index_char( (char *) "sigma_a", xs_field_names, xs_nfields);
    int sigma_s_n = bright::find_index_char( (char *) "sigma_s", xs_field_names, xs_nfields);

    for (std::set<int>::iterator i = FCComps::isos2track.begin(); i != FCComps::isos2track.end(); i++)
    {
        int iso_n = bright::find_index<int>(*i, isozz, xs_nrows);
        if (iso_n < 0)		
            iso_n = bright::find_index<int>(10*((*i)/10), isozz);
        if (iso_n < 0)
        {
            sigma_a_therm[*i] = 0.0;
            sigma_s_therm[*i] = 0.0;
            continue;
        };

        double * iso_sig_a = new double [1];
        double * iso_sig_s = new double [1];
        #ifdef _WIN32
            const size_t temp_xs_field_sizes_sigma_a_n [1] = {xs_field_sizes[sigma_a_n]};
            const size_t temp_xs_field_sizes_sigma_s_n [1] = {xs_field_sizes[sigma_s_n]};

            kdbstat = H5TBread_fields_name(kdblib, "/XS/ThermalMaxwellAve", xs_field_names[sigma_a_n], iso_n, 1, sizeof(double), 0, temp_xs_field_sizes_sigma_a_n, iso_sig_a);
            kdbstat = H5TBread_fields_name(kdblib, "/XS/ThermalMaxwellAve", xs_field_names[sigma_s_n], iso_n, 1, sizeof(double), 0, temp_xs_field_sizes_sigma_s_n, iso_sig_s);
        #else
            kdbstat = H5TBread_fields_name(kdblib, "/XS/ThermalMaxwellAve", xs_field_names[sigma_a_n], iso_n, 1, sizeof(double), 0, (const size_t [1]) {xs_field_sizes[sigma_a_n]}, iso_sig_a);
            kdbstat = H5TBread_fields_name(kdblib, "/XS/ThermalMaxwellAve", xs_field_names[sigma_s_n], iso_n, 1, sizeof(double), 0, (const size_t [1]) {xs_field_sizes[sigma_s_n]}, iso_sig_s);
        #endif
        sigma_a_therm[*i] = iso_sig_a[0];
        sigma_s_therm[*i] = iso_sig_s[0];
    };

    kdbstat = H5Fclose(kdblib);

    return;
};

void Reactor1G::foldMassWeights()
{
    /** Multiplies the burnup parameters by the mass weights."
     *  Calculates BU(F), P(F), D(F), and k(F)"
     */

    //First Things First, Let's calculate the atomic weight of the IHM
    double inverseA_IHM = 0.0;
    for (CompIter iso = IsosIn.comp.begin(); iso != IsosIn.comp.end(); iso++)
    {
        //Ensure that the isotope is officially allowed.
        if (0 == I.count(iso->first))
            continue;

        inverseA_IHM = inverseA_IHM + (iso->second)/ isoname::nuc_weight(iso->first);
    };
    A_IHM = 1.0 / inverseA_IHM;

    //Build the number density dictionaries
    niF.clear();
    niC.clear();
    //now for the ni in the Fuel
    for (std::map<std::string, double>::iterator key = FuelChemicalForm.begin(); key != FuelChemicalForm.end(); key++)
    {
        if ( (key->first) == "IHM")
        {
            for (CompIter iso = IsosIn.comp.begin(); iso != IsosIn.comp.end(); iso++)
            {
                //Ensure that the isotope is officially allowed.
                if (0 == I.count(iso->first))
                    continue;

                niF[iso->first] = FuelChemicalForm[key->first] * IsosIn.comp[iso->first];
            };
        }
        else
        {
            int key_zz = isoname::mixed_2_zzaaam(key->first);
            niF[key_zz] = FuelChemicalForm[key->first];
        }
    };
    //Note that the ni in the coolant is just CoolantChemicalForm
    for (std::map<std::string, double>::iterator key = CoolantChemicalForm.begin(); key != CoolantChemicalForm.end(); key++)
    {
        int key_zz = isoname::mixed_2_zzaaam(key->first);
        niC[key_zz] = CoolantChemicalForm[key->first];		
    };

    //Fuel mass weight
    miF.clear(); 
    for (CompIter iso = niF.begin(); iso != niF.end(); iso++)
    {
        if (niF[iso->first] == 0.0)
        {
            continue;
        }
        else
        {
            miF[iso->first] = niF[iso->first] * isoname::nuc_weight(iso->first) / A_IHM;
        }
    };

    //Coolant mass weight Calculation...requires MWF
    MWF = 0.0; 	//Fuel Molecular Weight
    for (std::map<std::string, double>::iterator key = FuelChemicalForm.begin(); key != FuelChemicalForm.end(); key++)
    {
        if ( (key->first) == "IHM")
            MWF = MWF + (FuelChemicalForm[key->first] * A_IHM);
        else
        {
            int key_zz = isoname::mixed_2_zzaaam(key->first);
            MWF = MWF + (FuelChemicalForm[key->first] * isoname::nuc_weight(key_zz));
        }
    };
    MWC = 0.0;	//Coolant Molecular Weight
    for (std::map<std::string, double>::iterator key = CoolantChemicalForm.begin(); key != CoolantChemicalForm.end(); key++)
    {
        int key_zz = isoname::mixed_2_zzaaam(key->first);
        MWC = MWC + (CoolantChemicalForm[key->first] * isoname::nuc_weight(key_zz));
    };
    miC.clear();
    double rel_Vol_coef = (rhoC * MWF * VC) / (rhoF * MWC * VF);
    for (CompIter iso = niC.begin(); iso != niC.end(); iso++)
    {
        if (niC[iso->first] == 0.0)
            continue;
        else
            miC[iso->first] = (niC[iso->first] * isoname::nuc_weight(iso->first) / A_IHM) * rel_Vol_coef;
    };

    //Fuel Number Density
    NiF.clear();
    for (CompIter iso = niF.begin(); iso != niF.end(); iso++)
    {
        NiF[iso->first] = niF[iso->first] * rhoF * (bright::N_A) / MWF;
    };

    //Coolant Number Density
    NiC.clear();
    for (CompIter iso = niC.begin(); iso != niC.end(); iso++)	{
        NiC[iso->first] = niC[iso->first] * rhoC * (bright::N_A) / MWC;
    };

    //BU(F)
    BU_F_.clear();
    BU_F_.assign( F.size(), 0.0 ); //re-initialize BU(F)
    for (CompIter i = miF.begin(); i != miF.end(); i++)
    {
        for (int f = 0; f < BU_F_.size(); f++)
        {
            BU_F_[f] = BU_F_[f] + (miF[i->first] * BUi_F_[i->first][f]);
        };
    };

    //P(F)
    P_F_.clear();
    P_F_.assign( F.size(), 0.0 );
    for (CompIter i = miF.begin(); i != miF.end(); i++)
    {
        for (int f = 0; f < P_F_.size(); f++)
        {
            P_F_[f] = P_F_[f] + (P_NL * miF[i->first] * pi_F_[i->first][f]);
        };
    };

    //d^F(F)
    dF_F_.clear();
    dF_F_.assign( F.size(), 0.0 );
    for (CompIter i = miF.begin(); i != miF.end(); i++)
    {
        for (int f = 0; f < dF_F_.size(); f++)
        {
            dF_F_[f] = dF_F_[f] + (miF[i->first] * di_F_[i->first][f]);
        };
    };

    //d^C(F)
    dC_F_.clear();
    dC_F_.assign( F.size(), 0.0 );
    for (CompIter i = miC.begin(); i != miC.end(); i++)
    {
        if (H_XS_Rescale && (i->first) == 10010)
        {
            for (int f = 0; f < dC_F_.size(); f++)
            {
                dC_F_[f] = dC_F_[f] + (miC[i->first] * di_F_[i->first][f] * (1.36927 - (0.01119 * BU_F_[f])));
            };
        }
        else
        {
            for (int f = 0; f < dC_F_.size(); f++)
            {
                dC_F_[f] = dC_F_[f] + (miC[i->first] * di_F_[i->first][f]);
            };
        };
    };

    //Implement the disadvantage factor, if needed.
    if (useZeta)
    {
        calcZeta();
        for (int f = 0; f < F.size(); f++)
        {
            dC_F_[f] = zeta_F_[f] * dC_F_[f];
        };
    }
    else
    {
        zeta_F_.clear();
        zeta_F_.assign( F.size(), 0.0 );
    };

    //D(F)
    D_F_.clear();
    D_F_.assign( F.size(), 0.0 );
    for (int f = 0; f < D_F_.size(); f++)
    {
        D_F_[f] = dF_F_[f] + dC_F_[f];
    };

    //k(F) -- almost meaningless
    k_F_.clear();
    k_F_.assign( F.size(), 0.0 );
    for (int f = 0; f < k_F_.size(); f++)
    {
        k_F_[f] = P_F_[f] / D_F_[f];
    };
};

void Reactor1G::mkMj_F_()
{
    //Generates the Mj(F) table.
    Mj_F_.clear();
    for (IsoIter j = J.begin(); j != J.end(); j++ )
    {
        Mj_F_[*j].assign( F.size(), 0.0 );
        for(CompIter i = IsosIn.comp.begin(); i != IsosIn.comp.end(); i++)
        {
            if (0 < I.count(i->first))
            {
                for (int f = 0; f < Mj_F_[*j].size(); f++)
                {
                    Mj_F_[*j][f] = Mj_F_[*j][f] + (miF[i->first] * Tij_F_[i->first][*j][f]);
                };
            };
        };
    };
};

void Reactor1G::mkMj_Fd_()
{
    /** Calculates the output isotopics of Mj(Fd).
     *  NOTE: Mj(Fd) is effectively the same variable as IsosOut before normalization!
     */

    CompDict tempOut;

    //Checks to see if the discharge index in at the end of the fluence table
    bool fdLast = false;
    if ( (fd+1) == F.size() )
        fdLast = true;			

    for (IsoIter j = J.begin(); j != J.end(); j++ )
    {
        if (fdLast)
            tempOut[*j] = bright::SolveLine(Fd, F[fd], Mj_F_[*j][fd], F[fd-1], Mj_F_[*j][fd-1]);
        else
            tempOut[*j] = bright::SolveLine(Fd, F[fd+1], Mj_F_[*j][fd+1], F[fd], Mj_F_[*j][fd]);
    };

    IsosOut = MassStream(tempOut);	
};

void Reactor1G::calcOutIso()
{
    //Wrapper to calculate discharge isotopics.
    mkMj_F_();
    mkMj_Fd_();
};

void Reactor1G::calcSubStreams()
{
    //Sets possibly relevant reactor input and output substreams.

    //Uranium
    InU  = IsosIn.get_u();
    OutU = IsosOut.get_u();

    //TRU
    try 
    {
        InTRU = IsosIn.get_tru();
    }
    catch (...)
    {
        CompDict cd;
        cd[942390] = 1.0;
        InTRU = MassStream(cd, 1.0);
        InTRU.mass = 0.0;
    };
    OutTRU = IsosOut.get_tru();

    //Lanthanides
    try
    {
        InLAN = IsosIn.get_lan();
    }
    catch (...)
    {
        CompDict cd;
        cd[581440] = 1.0;
        InLAN  = MassStream(cd, 1.0);
        InLAN.mass = 0.0;
    };
    OutLAN = IsosOut.get_lan();

    //Actinides
    InACT  = IsosIn.get_act();
    OutACT = IsosOut.get_act();
};


double Reactor1G::calc_deltaR()
{
    //Calculates the deltaR of the reactor with the current IsosIn
    foldMassWeights();
    deltaR = batchAve(TargetBU, "P") - batchAve(TargetBU, "D");
    return deltaR;
};

double Reactor1G::calc_deltaR(CompDict cd)
{
    //Calculates the deltaR of the reactor with the current IsosIn
    IsosIn = MassStream (cd);
    return calc_deltaR();
};

double Reactor1G::calc_deltaR(MassStream ms)
{
    //Calculates the deltaR of the reactor with the current IsosIn
    IsosIn = ms;
    return calc_deltaR();
};


double Reactor1G::calcTruCR()
{
    //Calculates the reactor's transuranic conversion ratio.
    TruCR = 1.0 - ((InTRU.mass - OutTRU.mass) / (BUd/931.46));
    return TruCR;
};



FluencePoint Reactor1G::FluenceAtBU(double BU)
{
    /** Gives the fluence at which the burnup BU occurs.
     *  Data returned as FluenceIndex stucture:
     *  	FI.f: index imeadiately lower than where BU achieved (int),
     *  	FI.F: flunece itself (double),
     *  	FI.m: slope dBU/dF bewteen points f and f+1 (double)
     */

    FluencePoint fp;

    //Finds the lower index
    fp.f = 0;
    while ( (fp.f < BU_F_.size()) && (BU_F_[fp.f] < BU) )
    {
        fp.f = fp.f + 1;
    };
    fp.f = fp.f - 1;

    if (fp.f < 0)
        fp.f = 0;
    
    if ( (fp.f + 1) == F.size() )
        fp.m = bright::slope(F[fp.f], BU_F_[fp.f], F[fp.f-1], BU_F_[fp.f-1]);
    else
        fp.m = bright::slope(F[fp.f+1], BU_F_[fp.f+1], F[fp.f], BU_F_[fp.f]);

    fp.F = ((BU - BU_F_[fp.f])/fp.m) + F[fp.f];

    return fp;
};

double Reactor1G::batchAve(double BUd, std::string PDk_flag)
{
    //Finds the batch-averaged P(F), D(F), or k(F) when at discharge burnup BUd.
    std::string PDk = bright::ToUpper(PDk_flag); 

    //For B batches (indexed from 0 -> B-1) calculate BU of each batch.
    std::vector<double> bu (B, 0.0);
    for (int b = 0; b < B; b++)
    {
        bu[b] = ((double) b + 1) * BUd / ((double) B);
    };

    //Calculate Fluence points for each BU
    std::vector<FluencePoint> fps; 
    fps.resize(B);
    for (int b = 0; b < B; b++)
    {
        fps[b] = FluenceAtBU(bu[b]);
    };

    std::vector<double> PDks (B, 0.0);
    for (int b = 0; b < B; b++)
    {
        double p;
        double d;

        if ( (fps[b].f + 1) == F.size() )
        {
            p = bright::SolveLine(fps[b].F, F[fps[b].f], P_F_[fps[b].f], F[fps[b].f-1], P_F_[fps[b].f-1]);
            d = bright::SolveLine(fps[b].F, F[fps[b].f], D_F_[fps[b].f], F[fps[b].f-1], D_F_[fps[b].f-1]);
        }
        else
        {
            p = bright::SolveLine(fps[b].F, F[fps[b].f+1], P_F_[fps[b].f+1], F[fps[b].f], P_F_[fps[b].f]);
            d = bright::SolveLine(fps[b].F, F[fps[b].f+1], D_F_[fps[b].f+1], F[fps[b].f], D_F_[fps[b].f]);
        };

        if (PDk == "K")
            PDks[b] = (p/d);
        else if (PDk == "P")
            PDks[b] = p;
        else if (PDk == "D")
            PDks[b] = d;
        else
        {
            if (1 < FCComps::verbosity) 
            {
                std::cout << "PDk flag is wrong: " << PDk << "\n";
                std::cout << "Using default of k.\n";
            };
            PDks[b] = (p/d);
        };
    };

    double numerator = 0.0;
    double denominator = 0.0;
    for (int b = 0; b < B; b++)
    {
        numerator   = numerator   + (PDks[b] / fps[b].m);
        denominator = denominator + (1.0 / fps[b].m);
    };

    return numerator/denominator;
};

double Reactor1G::batchAveK(double BUd) 
{
        return batchAve(BUd, "K");
};

void Reactor1G::BUd_BisectionMethod()
{
    //Calculates the maximum discharge burnup via the Bisection Method.
    int tempk = 1;
    double BUd_a, k_a, sign_a;
    double BUd_b, k_b, sign_b;
    double BUd_c, k_c, sign_c;
    
    //First Find a BUd that serves as an initial first guess
    while ( (tempk < k_F_.size()) && (1.0 < k_F_[tempk]) )
        tempk = tempk + 1;
    tempk = tempk - 1;
    if (tempk < 0)
        tempk = 0;

    BUd_a = BU_F_[tempk] * 2.0 * ((double) B) / ((double) B + 1.0);
    k_a = batchAveK( BUd_a );
    if (k_a == 1.0)
    {
        BUd = BUd_a;
        return;
    }
    else
        sign_a = (k_a - 1.0) / fabs(k_a - 1.0);
    
    //Find a BUd that serves as a valid second guess.  Remember, this is the bisection method here.
    BUd_b = BUd_a + sign_a * 5.0;
    k_b = batchAveK( BUd_b );
    if (k_b == 1.0)
    {
        BUd = BUd_b;
        return;
    }
    else
        sign_b = (k_b - 1.0) / fabs(k_b - 1.0);

    if (1 < FCComps::verbosity)
    {
        std::cout << "BUd_a = " << BUd_a << "\tk_a = " << k_a << "\tsign_a = " << sign_a << "\n";
        std::cout << "BUd_b = " << BUd_b << "\tk_b = " << k_b << "\tsign_b = " << sign_b << "\n";
        std::cout << "\n";
    };
    while (sign_a == sign_b)
    {
        BUd_b = BUd_b + sign_a * 5.0;
        k_b = batchAveK( BUd_b );
        if (k_b == 1.0)
        {
            BUd = BUd_b;
            return;
        }
        else
            sign_b = (k_b - 1.0) / fabs(k_b - 1.0);

        if (1 < FCComps::verbosity)
        {
            std::cout << "BUd_a = " << BUd_a << "\tk_a = " << k_a << "\tsign_a = " << sign_a << "\n";
            std::cout << "BUd_b = " << BUd_b << "\tk_b = " << k_b << "\tsign_b = " << sign_b << "\n";
            std::cout << "\n";
        };

        if ( (BUd_b < 0.0) || (1000.0 < BUd_b) )
            throw BadFuelForm ();
    };

    BUd_c = 0.0;
    k_c = 0.0;

    //Ok now that we have valid and hopefully close initial conditions, let's do it!
    double DoA = pow(10.0, -7);	//Degree of accuracy to carry out calculations to.
    int q = 0;			//index for number of iterations
    while ( ((DoA < fabs(1.0 - k_a)) || (DoA < fabs(1.0 - k_b))) && (0.0 < fabs(BUd_a - BUd_b)) && (q < 100) )
    {
        BUd_c = (BUd_a + BUd_b) / 2.0;
        k_c = batchAveK( BUd_c );
        if (k_c == 1.0)
        {
            sign_c = 0.0;
            break;
        }
        else
            sign_c = (k_c - 1.0) / fabs(k_c - 1.0);

        q = q + 1;

        if ( (sign_a == sign_c) && (sign_b != sign_c) )
        {
            BUd_a = BUd_c;
            k_a = k_c;
            sign_a = sign_c;
        }
        else if ( (sign_b == sign_c) && (sign_a != sign_c) )
        {
            BUd_b = BUd_c;
            k_b = k_c;
            sign_b = sign_c;
        }
        else
        {
            if (0 < FCComps::verbosity)
            {
                std::cout << "\n";
                std::cout << "SOMEWHERE WHILE FINDING k SOMETHING WENT WRONG!!!\n";
                std::cout << "Here is some information that might help you debug ^_^\n";
                std::cout << "BUd_a = " << BUd_a << "\tk_a = " << k_a << "\tsign_a = " << sign_a << "\n";
                std::cout << "BUd_b = " << BUd_b << "\tk_b = " << k_b << "\tsign_b = " << sign_b << "\n";
                std::cout << "BUd_c = " << BUd_c << "\tk_c = " << k_c << "\tsign_c = " << sign_c << "\n";
                std::cout << "\n";
            };
        };
    };

    //If c-set of variables wasn't altered, raise an exception.
    if ( (BUd_c == 0.0) && (k_c == 0.0) )
    {
        if (0 < FCComps::verbosity)
        {
            std::cout << "\n";
            std::cout << "SOMEWHERE WHILE FINDING k SOMETHING WENT WRONG!!!\n";
            std::cout << "Here is some information that might help you debug ^_^\n";
            std::cout << "BUd_a = " << BUd_a << "\tk_a = " << k_a << "\tsign_a = " << sign_a << "\n";
            std::cout << "BUd_b = " << BUd_b << "\tk_b = " << k_b << "\tsign_b = " << sign_b << "\n";
            std::cout << "BUd_c = " << BUd_c << "\tk_c = " << k_c << "\tsign_c = " << sign_c << "\n";
            std::cout << "\n";
        };
        throw BisectionMethodNotPerformed ("Burnup");
    };

    //print results, if desired.
    if (0 < FCComps::verbosity)
    {
        std::cout << "Final Result of Burnup Bisection Method Calculation:\n";
        std::cout << "BUd_a = " << BUd_a << "\tk_a = " << k_a << "\tsign_a = " << sign_a << "\n";
        std::cout << "BUd_b = " << BUd_b << "\tk_b = " << k_b << "\tsign_b = " << sign_b << "\n";
        std::cout << "BUd_c = " << BUd_c << "\tk_c = " << k_c << "\tsign_c = " << sign_c << "\n";
        std::cout << "\n";
    };

    //Sort to find closest value among results.
    if (sign_c == 0.0)
    {
        BUd = BUd_c;
        k = k_c;
    }
    else if ( fabs(k_a - 1.0) < DoA )
    {
        BUd = BUd_a;
        k = k_a;
    }
    else if ( fabs(k_b - 1.0) < DoA )
    {
        BUd = BUd_b;
        k = k_b;
    }
    else
    {
        if (0 < FCComps::verbosity)
            std::cout << "k did not converge with the Bisection Method to an accuracy of " << DoA << " in " << q << " iterations.\n";

        if ( (fabs(k_a - 1.0) < 0.01) && (fabs(k_a - 1.0) < fabs(k_b -1.0)) )
        {
            BUd = BUd_a;
            k = k_a;
            if (0 < FCComps::verbosity)
                std::cout << "However, k_a is within 1% of 1 and closer to 1 than k_b; using these values.\n";
        }
        else if ( (fabs(k_b - 1.0) < 0.01) && (fabs(k_b - 1.0) < fabs(k_a -1.0)) )
        {
            BUd = BUd_b;
            k = k_b;
            if (0 < FCComps::verbosity)
                std::cout << "However, k_b is within 1% of 1 and closer to 1 than k_a; using these values.\n";
        }
        else
        {
            if (0 < FCComps::verbosity)
                std::cout << "Alright.  It really didn't converge. Neither k_a nor k_b is within 1% of 1. Program will likely fail!\n";
        };
    };

    FluencePoint fp = FluenceAtBU(BUd);
    fd = fp.f;				//lower index of fluence at discharge
    Fd = fp.F;				//Fluence at discharge
    return;
};

void Reactor1G::Run_PNL(double temp_pnl)
{
    /** Does a reactor run for a specific P_NL.
     *  Requires that IsosIn be (meaningfully) set.
     *  For use with Calibrate_PNL_2_BUd
     */

    P_NL = temp_pnl;
    foldMassWeights();
    BUd_BisectionMethod();
};

void Reactor1G::Calibrate_PNL_2_BUd()
{
    /** Calibrates the non-leakage probability of a reactors to hit a target burnup.
     *  Calibration proceeds by bisection method...
     */

    double pnl_a, bud_a, sign_a;
    double pnl_b, bud_b, sign_b;
    double pnl_c, bud_c, sign_c;

    //Find an acceptable lower bound
    pnl_a = 0.05;
    bool FoundA = false;
    while (!FoundA)
    {
        try
        {
            Run_PNL(pnl_a);
            bud_a = BUd;
            sign_a = (bud_a - TargetBU) / fabs(bud_a - TargetBU);
            FoundA = true;
        }
        catch (BadFuelForm e)
        {
            pnl_a = pnl_a + 0.05;
        }
        catch (BisectionMethodNotPerformed e)
        {
            pnl_a = pnl_a + 0.05;
        };
    };

    //Find an acceptable upper bound
    pnl_b = 2.0;
    bool FoundB = false;
    while (!FoundB)
    {
        try
        {
            Run_PNL(pnl_b);
            bud_b = BUd;
            sign_b = (bud_b - TargetBU) / fabs(bud_b - TargetBU);
            FoundB = true;
        }
        catch (BadFuelForm e)
        {
            pnl_b = pnl_b - 0.05;
        }
        catch (BisectionMethodNotPerformed e)
        {
            pnl_b = pnl_b - 0.05;
        };
    };
        

    //Perform bisection method.
    double DoA = pow(10.0, -15);	//Degree of accuracy to carry out calculations to.
    int q = 0;
    while ( (DoA < fabs(pnl_a - pnl_b)) && (DoA < fabs(bud_a - bud_b)) && (q < 100) )
    {
        pnl_c = (pnl_a + pnl_b) / 2.0;
        Run_PNL(pnl_c);
        bud_c = BUd;
        sign_c = (bud_c - TargetBU) / fabs(bud_c - TargetBU);

        q = q + 1;

        if ( (sign_a == sign_c) && (sign_b != sign_c) )
        {
            pnl_a = pnl_c;
            bud_a = bud_c;
            sign_a = sign_c;
        }
        else if ( (sign_b == sign_c) && (sign_a != sign_c) )
        {
            pnl_b = pnl_c;
            bud_b = bud_c;
            sign_b = sign_c;
        }
        else
        {
            if (0 < FCComps::verbosity)
            {
                std::cout << "\n";
                std::cout << "SOMEWHERE WHILE FINDING k SOMETHING WENT WRONG!!!\n";
                std::cout << "Here is some information that might help you debug ^_^\n";
                std::cout << "pnl_a = " << pnl_a << "\tBUd_a = " << bud_a << "\tsign_a = " << sign_a << "\n";
                std::cout << "pnl_b = " << pnl_b << "\tBUd_b = " << bud_b << "\tsign_b = " << sign_b << "\n";
                std::cout << "pnl_c = " << pnl_c << "\tBUd_c = " << bud_c << "\tsign_c = " << sign_c << "\n";
                std::cout << "\n";
            };
        };
    };

    if (0 < FCComps::verbosity)
    {
        std::cout << "\n";
        std::cout << "Final Result P_NL Calibration to Burnup via Bisection Method Calculation:\n";
        std::cout << "Number of iterations q = " << q << "\n";
        std::cout << "pnl_a = " << pnl_a << "\tBUd_a = " << bud_a << "\tsign_a = " << sign_a << "\n";
        std::cout << "pnl_b = " << pnl_b << "\tBUd_b = " << bud_b << "\tsign_b = " << sign_b << "\n";
        std::cout << "pnl_c = " << pnl_c << "\tBUd_c = " << bud_c << "\tsign_c = " << sign_c << "\n";
    };

    return;	
};



MassStream Reactor1G::doCalc ()
{
    //Finds BUd and output isotopics.
    foldMassWeights();

    BUd_BisectionMethod();

    calcOutIso();

    return IsosOut;
};

MassStream Reactor1G::doCalc (CompDict incomp)
{
    //Finds BUd and output isotopics.
    IsosIn = MassStream (incomp);

    return doCalc();
};

MassStream Reactor1G::doCalc (MassStream instream)
{
    //Finds BUd and output isotopics.
    IsosIn = instream;

    return doCalc();
};


void Reactor1G::LatticeEPlanar(double a, double b)
{
    LatticeE_F_.clear();
        LatticeE_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        if (0.0 == kappaC_F_[f])
            LatticeE_F_[f] = 0.0;
        else
            LatticeE_F_[f] = kappaC_F_[f] * (b - a) * bright::COTH(kappaC_F_[f]*(b-a));
    };
    return;
};

void Reactor1G::LatticeFPlanar(double a, double b)
{
    LatticeF_F_.clear();
        LatticeF_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        if (0.0 == kappaF_F_[f])
            LatticeF_F_[f] = 0.0;
        else
            LatticeF_F_[f] = LatticeF_F_[f] * a * bright::COTH(LatticeF_F_[f]*a) ;
    };
    return; 
};

void Reactor1G::LatticeESpherical(double a, double b)
{
    LatticeE_F_.clear();
        LatticeE_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        if (0.0 == kappaC_F_[f])
            LatticeE_F_[f] = 0.0;
        else
        {
            double coef = pow(kappaC_F_[f], 3) * (pow(b,3) - pow(a,3)) / (3*kappaC_F_[f]*a);
            double num = 1.0 - ( kappaC_F_[f] * b * bright::COTH(kappaC_F_[f]*(b-a)) );
            double den = 1.0 - (pow(kappaC_F_[f], 2)*a*b) - ( kappaC_F_[f]*(b-a) * bright::COTH(kappaC_F_[f]*(b-a)) );
            LatticeE_F_[f] = coef * num / den;
        };
    };
    return;
};
    
void Reactor1G::LatticeFSpherical(double a, double b)
{
    LatticeF_F_.clear();
        LatticeF_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        double coef = pow(kappaF_F_[f], 2) * pow(a, 2) / 3.0;
        double num = bright::TANH(kappaF_F_[f]*a); 
        double den = (kappaF_F_[f]*a) - bright::TANH(kappaF_F_[f]*a); 
        LatticeF_F_[f] = coef * num / den;
    };
    return; 
};

void Reactor1G::LatticeECylindrical(double a, double b)
{
    namespace bm = boost::math;

    LatticeE_F_.clear();
        LatticeE_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        if (0.0 == kappaC_F_[f])
            LatticeE_F_[f] = 0.0;
        else
        {
            double coef = kappaC_F_[f] * (pow(b,2) - pow(a,2)) / (2.0*a);
            double num = ( bm::cyl_bessel_i(0, kappaC_F_[f]*a) * bm::cyl_bessel_k(1, kappaC_F_[f]*b) ) + \
                ( bm::cyl_bessel_k(0, kappaC_F_[f]*a) * bm::cyl_bessel_i(1, kappaC_F_[f]*b) );
            double den = ( bm::cyl_bessel_i(1, kappaC_F_[f]*b) * bm::cyl_bessel_k(1, kappaC_F_[f]*a) ) - \
                ( bm::cyl_bessel_k(1, kappaC_F_[f]*b) * bm::cyl_bessel_i(1, kappaC_F_[f]*a) );
            LatticeE_F_[f] = coef * num / den;
        };
    };
    return;
};
    
void Reactor1G::LatticeFCylindrical(double a, double b)
{
    namespace bm = boost::math;

    LatticeF_F_.clear();
        LatticeF_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        if (0.0 == kappaF_F_[f])
            LatticeE_F_[f] = 0.0;
        else
        {
            double num =  kappaF_F_[f] * a * bm::cyl_bessel_i(0, kappaF_F_[f]*a);
            double den = 2.0 * bm::cyl_bessel_i(1, kappaF_F_[f]*a);
            LatticeF_F_[f] = num / den;
        };
    };
    return;
};

void Reactor1G::calcZeta()
{
    // Computes the thermal disadvantage factor

    //Prepare data sets to calculate the disadvantage factor
    //For Fuel...
    SigmaFa_F_.clear();
    SigmaFtr_F_.clear();
    kappaF_F_.clear();

    SigmaFa_F_.assign( F.size(), 0.0 );
    SigmaFtr_F_.assign( F.size(), 0.0 );
    kappaF_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        //Calculate the following...
        //	* Macroscopic abspobtion XS in Fuel
        //	* Macroscopic transport XS in Fuel
        for(CompIter iso = niF.begin(); iso != niF.end(); iso++)
        {
            //If Lanthanide or Actinide, use ORIGEN Data as sigma_a
            //Else use KAERI Data for sigma_a
            if (570000 < iso->first < 720000 || 890000 < iso->first)
            {
                SigmaFa_F_[f]  = SigmaFa_F_[f]  + (NiF[iso->first] * di_F_[iso->first][f] * bright::bpcm2);

                SigmaFtr_F_[f] = SigmaFtr_F_[f] + (NiF[iso->first] * bright::bpcm2 * (di_F_[iso->first][f] + \
                    sigma_s_therm[iso->first]*(1.0 - 2.0/(3.0*isoname::nuc_weight(iso->first))) ) );
            }
            else
            {
                //renormalize sigma_a for this fluenece
                double sig_a = sigma_a_therm[iso->first] * di_F_[iso->first][f] / di_F_[iso->first][0];

                SigmaFa_F_[f]  = SigmaFa_F_[f]  + (NiF[iso->first] * sig_a * bright::bpcm2);

                SigmaFtr_F_[f] = SigmaFtr_F_[f] + (NiF[iso->first] * bright::bpcm2 * (sig_a + \
                    sigma_s_therm[iso->first]*(1.0 - 2.0/(3.0*isoname::nuc_weight(iso->first))) ) );
            };
        };

        //Calculate kappa
        kappaF_F_[f]   = sqrt( 3.0 * SigmaFtr_F_[f] * SigmaFa_F_[f] );
    };


    //For Coolant...
    SigmaCa_F_.clear();
    SigmaCtr_F_.clear();
    kappaC_F_.clear();

    SigmaCa_F_.assign( F.size(), 0.0 );
    SigmaCtr_F_.assign( F.size(), 0.0 );
    kappaC_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        //Calculate the following...
        //	* Macroscopic abspobtion XS in Coolant
        //	* Macroscopic transport XS in Coolant
        for(CompIter iso = niC.begin(); iso != niC.end(); iso++)
        {
            //If Lanthanide or Actinide, use ORIGEN Data as sigma_a
            //Else use KAERI Data for sigma_a
            if (570000 < iso->first < 720000 || 890000 < iso->first)
            {
                SigmaCa_F_[f]  = SigmaCa_F_[f]  + (NiC[iso->first] * di_F_[iso->first][f] * bright::bpcm2);

                SigmaCtr_F_[f] = SigmaCtr_F_[f] + (NiC[iso->first] * bright::bpcm2 * (di_F_[iso->first][f] + \
                    sigma_s_therm[iso->first]*(1.0 - 2.0/(3.0*isoname::nuc_weight(iso->first))) ) );
            }
            else
            {
                //renormalize sigma_a for this fluenece
                double sig_a = sigma_a_therm[iso->first] * di_F_[iso->first][f] / di_F_[iso->first][0];

                SigmaCa_F_[f]  = SigmaCa_F_[f]  + (NiC[iso->first] * sig_a * bright::bpcm2);

                SigmaCtr_F_[f] = SigmaCtr_F_[f] + (NiC[iso->first] * bright::bpcm2 * (sig_a + \
                    sigma_s_therm[iso->first]*(1.0 - 2.0/(3.0*isoname::nuc_weight(iso->first))) ) );
            };
        };

        //Calculate kappa
        kappaC_F_[f]   = sqrt( 3.0 * SigmaCtr_F_[f] * SigmaCa_F_[f] );
    };

    //Calculate the Lattice Functions
    double a, b;
    if (Lattice == "Planar")
    {
        a = r;
        b = l / 2.0;
    
        LatticeEPlanar(a, b);
        LatticeFPlanar(a, b);
    }
    else if (Lattice == "Spherical")
    {
        a = r;
        b = l / 2.0;
    
        LatticeESpherical(a, b);
        LatticeFSpherical(a, b);
    }
    else if (Lattice == "Cylindrical")
    {
        a = r;
        b = l / sqrt(bright::pi); //radius of cylinder with an equivilent cell volume

        LatticeECylindrical(a, b);
        LatticeFCylindrical(a, b);
    }
    else
    {
        if (0 < FCComps::verbosity)
            std::cout << "Did not specify use of planar or spheical or cylindrical lattice functions! Assuming cylindrical...\n";
        
        a = r;
        b = l / sqrt(bright::pi); //radius of cylinder with an equivilent cell volume

        LatticeECylindrical(a, b);
        LatticeFCylindrical(a, b);
    };

    //Finally, Calculate Zeta
    zeta_F_.clear();
    zeta_F_.assign( F.size(), 0.0);
    for (int f = 0; f < F.size(); f++) 
    {
        if (0.0 == SigmaCa_F_[f])
            zeta_F_[f] = 1.0;
        else
            zeta_F_[f] = LatticeF_F_[f] + ( SigmaFa_F_[f] * VF * (LatticeE_F_[f] - 1.0) / (SigmaCa_F_[f] * VC) );
    };


    //Unfortunately, the above formulation for the disadvantage factor is ONLY valid for a << b!!!
    //Often times in modern (thermal) reactors, this is not the case.
    //We have a 'thin moderator' situation.
    //
    //To fix this problem correctly requires going to a multi-region diffusion/transport calculation.
    //Doing so is beyond the scope of this code.
    //What is more in-line with current practice is to use the results of a more sophisticated method,
    //interpolate over them, and use them here.
    //
    //That is what is done here when 0.1 < VF / VC, (ie the fuel is greater than 10% of the coolant)
    //A baseline zeta is determined from data presented in "Thermal disadvantage factor calculation by 
    //the multiregion collision probability method" by B. Ozgener,  and H. A. Ozgener, Institute of 
    //Nuclear Energy, Istanbul Technical University 80626 Maslak, Istanbul, Turkey, Received 
    //28 January 2003;  accepted 20 May 2003.  Available online 6 March 2004.
    //This baseline is a function of (VF/VC).
    // 
    //The above calculation of zeta is then used as a scaling factor on the baseline function to 
    //account for variations in fuel composition and fluenece.

    //Check if we are in the proper Fuel-to-Coolant Regime

    double f2c = VF / VC;
    if (f2c < 0.1)
        return;

    double zetabase  = 1.30857959 - (0.10656299 * f2c);
    double zetaratio = zetabase / zeta_F_[0];	

    for (int f = 0; f < F.size(); f++) 
    {
        zeta_F_[f] = zeta_F_[f] * zetaratio;

        if (zeta_F_[f] < 1.0)
            zeta_F_[f] = 1.0;
    };

    return;
};

void Reactor1G::calcZetaPlanar()
{
    Lattice = "Planar";
    calcZeta();
    return;
};

void Reactor1G::calcZetaSpherical()
{
    Lattice = "Spherical";
    calcZeta();
    return;
};

void Reactor1G::calcZetaCylindrical()
{
    Lattice = "Cylindrical";
    calcZeta();
    return;
};
