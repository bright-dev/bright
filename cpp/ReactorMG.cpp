// One Group Reactor Component Class

#include "ReactorMG.h"

/***********************************************/
/*** ReactorMG Component Class and Functions ***/
/***********************************************/

ReactorMG::ReactorMG()
{
};


ReactorMG::ReactorMG(std::string n) : FCComp(n)
{
};


ReactorMG::ReactorMG(std::set<std::string> paramtrack, std::string n) : FCComp(paramtrack, n)
{
};


ReactorMG::ReactorMG(ReactorParameters rp, std::string n) : FCComp(n)
{
    initialize(rp);
};


ReactorMG::ReactorMG(ReactorParameters rp, std::set<std::string> paramtrack, std::string n) : FCComp(paramtrack, n)
{
    initialize(rp);
};


ReactorMG::~ReactorMG()
{
};


void ReactorMG::initialize(ReactorParameters rp)
{
    /** Sets reactor specific parameters.
     *  Must be done once at the beginning of reactor object life.
     */

    B = rp.batches;				//Total number of fuel loading batches
    flux = rp.flux;				//Flux used for Fluence
    fuel_chemical_form = rp.fuel_form;		//Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
    coolant_chemical_form = rp.coolant_form;	//Same a fuel chemical form but for coolant.  Should not have "IHM"

    rho_fuel = rp.fuel_density;			//Fuel Density
    rho_clad = rp.cladding_density;			// Cladding Density
    rho_cool = rp.coolant_density;		//Coolant Density

    P_NL = rp.pnl;				//Non-Leakage Probability
    target_BU = rp.BUt;			//Target Discharge Burnup, only used for graphing inside of this component
    burn_regions = rp.burn_regions;
    specific_power = rp.specific_power;
    burn_time = 0.0;
    burn_times = rp.burn_times;

    // Flags
    use_zeta = rp.use_disadvantage_factor;		//Boolean value on whether or not the disadvantage factor should be used
    lattice_flag = rp.lattice_type;		//lattice_flagType (Planar || Spherical || Cylindrical)
    rescale_hydrogen_xs = rp.rescale_hydrogen;	//Rescale the Hydrogen-1 XS?

    //Calculates Volumes
    r_fuel = rp.fuel_radius;			//Fuel region radius
    r_void = rp.void_radius;	//Void region radius
    r_clad = rp.clad_radius;	//clad region radius
    pitch = rp.unit_cell_pitch;			//Unit cell side length

    S_O = rp.open_slots;		//Number of open slots in fuel assembly
    S_T = rp.total_slots;		//Total number of Fuel assembly slots.
    //Fuel Volume
    VF = ((bright::pi * r_fuel * r_fuel)/(pitch*pitch)) * (1.0 - S_O/S_T); 
    //Coolant Volume
    VC = ((pitch*pitch - bright::pi * r_fuel * r_fuel)/(pitch*pitch)) * (1.0 - S_O/S_T) + (S_O/S_T);
};



void ReactorMG::loadlib(std::string libfile)
{
    //Loads Apporiate Libraries for ReactorMG

    // Check that the file is there
    if (!bright::FileExists(libfile))
        throw bright::FileNotFound(libfile);

    //Check to see if the file is in HDF5 format.
    bool isH5 = H5::H5File::isHdf5(libfile);
    if (!isH5)
    {
        std::cout << "!!!Warning!!! " << libfile << " is not a valid HDF5 file!\n";
        return;
    };

    // Turn off the exceptions
    H5::Exception::dontPrint();

    // Open file
    H5::H5File rmglib(libfile, H5F_ACC_RDONLY);

    // Load isos
    I = h5wrap::h5_array_to_cpp_set<int>(&rmglib, "/load_isos_zz", H5::PredType::NATIVE_INT);
    J = h5wrap::h5_array_to_cpp_set<int>(&rmglib, "/transmute_isos_zz", H5::PredType::NATIVE_INT);

    // Load perturbation table
    perturbations = h5wrap::HomogenousTypeTable<double>(&rmglib, "/perturbations");
    nperturbations = perturbations.shape[0];

    // Load in energy structure
    pert_data_g full_E_g = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/energy");
    E_g = full_E_g[0];
    G = E_g.size() - 1;

    // Load fluxes and fluence
    phi_g = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/phi_g");
    phi = h5wrap::h5_array_to_cpp_vector_1d<double>(&rmglib, "/phi");
    Phi = h5wrap::h5_array_to_cpp_vector_1d<double>(&rmglib, "/Phi");

    // Load time and burnup
    time0 = h5wrap::h5_array_to_cpp_vector_1d<double>(&rmglib, "/time0");
    BU0 = h5wrap::h5_array_to_cpp_vector_1d<double>(&rmglib, "/BU0");


    // Clear transmutation vectors and cross sections before reading in
    Ti0.clear();
    sigma_a.clear();
    sigma_s.clear();
    sigma_f.clear();
    nubar_sigma_f.clear();
    nubar.clear();

    // Load transmutation vectors and cross sections that are based off of isotope
    int iso_zz;
    std::string iso_LL;
    for(std::set<int>::iterator iso_iter = J.begin(); iso_iter != J.end(); iso_iter++)
    {
        iso_zz = *iso_iter;
        iso_LL = isoname::zzaaam_2_LLAAAM(iso_zz);

        // Add transmutation vector
        Ti0[iso_zz] = h5wrap::h5_array_to_cpp_vector_1d<double>(&rmglib, "/Ti0/" + iso_LL);

        // Add cross sections
        sigma_a[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_a/" + iso_LL);
        sigma_s[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_s/" + iso_LL);
        sigma_f[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_f/" + iso_LL);
        nubar_sigma_f[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/nubar_sigma_f/" + iso_LL);
        nubar[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/nubar/" + iso_LL);
        sigma_s_gh[iso_zz] = h5wrap::h5_array_to_cpp_vector_3d<double>(&rmglib, "/sigma_s_gh/" + iso_LL);
    };

    // close the reactor library
    rmglib.close();

    //Now get microscopic XS data from KAERI...
    //...But only if the disadvantage factor is used.
    if (!use_zeta)
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

    for (std::set<int>::iterator i = FCComps::track_isos.begin(); i != FCComps::track_isos.end(); i++)
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





void ReactorMG::calc_nearest_neighbors()
{
    /**
     * Returns a vector of the indices sorted, sorted by nearest neighbor
     */

    // Initialize
    std::map<std::string, std::vector<double> > norm_deltas;
    std::vector<double> rss (nperturbations, 0.0);
    nearest_neighbors.clear();

    // Calcumlate normailized deltas
    norm_deltas["fuel_density"] = bright::normalized_delta(rho_fuel, perturbations["fuel_density"]);
    norm_deltas["clad_density"] = bright::normalized_delta(rho_clad, perturbations["clad_density"]);
    norm_deltas["cool_density"] = bright::normalized_delta(rho_cool, perturbations["cool_density"]);

    norm_deltas["fuel_cell_radius"] = bright::normalized_delta(r_fuel, perturbations["fuel_cell_radius"]);
    norm_deltas["void_cell_radius"] = bright::normalized_delta(r_void, perturbations["void_cell_radius"]);
    norm_deltas["clad_cell_radius"] = bright::normalized_delta(r_clad, perturbations["clad_cell_radius"]);

    norm_deltas["unit_cell_pitch"] = bright::normalized_delta(pitch, perturbations["unit_cell_pitch"]);
    norm_deltas["burn_regions"] = bright::normalized_delta(burn_regions, perturbations["burn_regions"]);
    norm_deltas["fuel_specific_power"] = bright::normalized_delta(specific_power, perturbations["fuel_specific_power"]);

    // Calc pertubations for initial mass streams
    if (10 < perturbations.shape[1])
    {
        int iso_zz;
        std::string iso_LL;

        for (int p = 9; p < perturbations.shape[1] - 1; p++)
        {
            iso_LL = perturbations.cols[p];
            iso_zz = isoname::LLAAAM_2_zzaaam(iso_LL);

            if (0 < ms_feed.comp.count(iso_zz))
                norm_deltas[iso_LL] = bright::normalized_delta(ms_feed.comp[iso_zz], perturbations[iso_LL]);
            else
                norm_deltas[iso_LL] = bright::normalized_delta(0.0, perturbations[iso_LL]);                
        };
    };

    norm_deltas["burn_times"] = bright::normalized_delta(burn_time, perturbations["burn_times"]);


    // Now that we have the normalized deltas for each index
    // We want to calculate the 'distance' from the curent point
    // This is done by taking the root of the sum of squares for 
    // each index.
    std::string key;
    for (int p = 0; p < perturbations.shape[1]; p++)
    {
        key = perturbations.cols[p];

        // Take the sum of the squares
        for (int q = 0; q < perturbations.shape[0]; q++)
        {
            rss[q] = rss[q] + (norm_deltas[key][q] * norm_deltas[key][q]); 
        };
    };

    // Takes the root of the sum of squares
    for (int q = 0; q < perturbations.shape[0]; q++)
    {
        rss[q] = sqrt(rss[q]);
    };


    // Now that we have the root of the sum of squares, 
    // we need to sort the indices by rss from lowest 
    // to highest.
    nearest_neighbors = bright::sorted_index(rss);
};





void ReactorMG::fold_mass_weights()
{
    /** Multiplies the burnup parameters by the mass weights."
     *  Calculates BU(F), P(F), D(F), and k(F)"
     */

    //First Things First, Let's calculate the atomic weight of the IHM
    double inverseA_IHM = 0.0;
    for (CompIter iso = ms_feed.comp.begin(); iso != ms_feed.comp.end(); iso++)
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
    for (std::map<std::string, double>::iterator key = fuel_chemical_form.begin(); key != fuel_chemical_form.end(); key++)
    {
        if ( (key->first) == "IHM")
        {
            for (CompIter iso = ms_feed.comp.begin(); iso != ms_feed.comp.end(); iso++)
            {
                //Ensure that the isotope is officially allowed.
                if (0 == I.count(iso->first))
                    continue;

                niF[iso->first] = fuel_chemical_form[key->first] * ms_feed.comp[iso->first];
            };
        }
        else
        {
            int key_zz = isoname::mixed_2_zzaaam(key->first);
            niF[key_zz] = fuel_chemical_form[key->first];
        }
    };
    //Note that the ni in the coolant is just coolant_chemical_form
    for (std::map<std::string, double>::iterator key = coolant_chemical_form.begin(); key != coolant_chemical_form.end(); key++)
    {
        int key_zz = isoname::mixed_2_zzaaam(key->first);
        niC[key_zz] = coolant_chemical_form[key->first];		
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
    for (std::map<std::string, double>::iterator key = fuel_chemical_form.begin(); key != fuel_chemical_form.end(); key++)
    {
        if ( (key->first) == "IHM")
            MWF = MWF + (fuel_chemical_form[key->first] * A_IHM);
        else
        {
            int key_zz = isoname::mixed_2_zzaaam(key->first);
            MWF = MWF + (fuel_chemical_form[key->first] * isoname::nuc_weight(key_zz));
        }
    };
    MWC = 0.0;	//Coolant Molecular Weight
    for (std::map<std::string, double>::iterator key = coolant_chemical_form.begin(); key != coolant_chemical_form.end(); key++)
    {
        int key_zz = isoname::mixed_2_zzaaam(key->first);
        MWC = MWC + (coolant_chemical_form[key->first] * isoname::nuc_weight(key_zz));
    };
    miC.clear();
    double rel_Vol_coef = (rho_cool * MWF * VC) / (rho_fuel * MWC * VF);
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
        NiF[iso->first] = niF[iso->first] * rho_fuel * (bright::N_A) / MWF;
    };

    //Coolant Number Density
    NiC.clear();
    for (CompIter iso = niC.begin(); iso != niC.end(); iso++)	{
        NiC[iso->first] = niC[iso->first] * rho_cool * (bright::N_A) / MWC;
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
        if (rescale_hydrogen_xs && (i->first) == 10010)
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
    if (use_zeta)
    {
        calc_zeta();
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

void ReactorMG::calc_Mj_F_()
{
    //Generates the Mj(F) table.
    Mj_F_.clear();
    for (IsoIter j = J.begin(); j != J.end(); j++ )
    {
        Mj_F_[*j].assign( F.size(), 0.0 );
        for(CompIter i = ms_feed.comp.begin(); i != ms_feed.comp.end(); i++)
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

void ReactorMG::calc_Mj_Fd_()
{
    /** Calculates the output isotopics of Mj(Fd).
     *  NOTE: Mj(Fd) is effectively the same variable as ms_prod before normalization!
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

    ms_prod = MassStream(tempOut);	
};

void ReactorMG::calc_ms_prod()
{
    //Wrapper to calculate discharge isotopics.
    calc_Mj_F_();
    calc_Mj_Fd_();
};

void ReactorMG::calcSubStreams()
{
    //Sets possibly relevant reactor input and output substreams.

    //Uranium
    ms_feed_u  = ms_feed.get_u();
    ms_prod_u = ms_prod.get_u();

    //TRU
    try 
    {
        ms_feed_tru = ms_feed.get_tru();
    }
    catch (...)
    {
        CompDict cd;
        cd[942390] = 1.0;
        ms_feed_tru = MassStream(cd, 1.0);
        ms_feed_tru.mass = 0.0;
    };
    ms_prod_tru = ms_prod.get_tru();

    //Lanthanides
    try
    {
        ms_feed_lan = ms_feed.get_lan();
    }
    catch (...)
    {
        CompDict cd;
        cd[581440] = 1.0;
        ms_feed_lan  = MassStream(cd, 1.0);
        ms_feed_lan.mass = 0.0;
    };
    ms_prod_lan = ms_prod.get_lan();

    //Actinides
    ms_feed_act  = ms_feed.get_act();
    ms_prod_act = ms_prod.get_act();
};


double ReactorMG::calc_deltaR()
{
    //Calculates the deltaR of the reactor with the current ms_feed
    fold_mass_weights();
    deltaR = batch_average(target_BU, "P") - batch_average(target_BU, "D");
    return deltaR;
};

double ReactorMG::calc_deltaR(CompDict cd)
{
    //Calculates the deltaR of the reactor with the current ms_feed
    ms_feed = MassStream (cd);
    return calc_deltaR();
};

double ReactorMG::calc_deltaR(MassStream ms)
{
    //Calculates the deltaR of the reactor with the current ms_feed
    ms_feed = ms;
    return calc_deltaR();
};


double ReactorMG::calc_tru_cr()
{
    //Calculates the reactor's transuranic conversion ratio.
    tru_cr = 1.0 - ((ms_feed_tru.mass - ms_prod_tru.mass) / (BUd/931.46));
    return tru_cr;
};



FluencePoint ReactorMG::fluence_at_BU(double BU)
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

double ReactorMG::batch_average(double BUd, std::string PDk_flag)
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
        fps[b] = fluence_at_BU(bu[b]);
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

double ReactorMG::batch_average_k(double BUd) 
{
        return batch_average(BUd, "K");
};

void ReactorMG::BUd_bisection_method()
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
    k_a = batch_average_k( BUd_a );
    if (k_a == 1.0)
    {
        BUd = BUd_a;
        return;
    }
    else
        sign_a = (k_a - 1.0) / fabs(k_a - 1.0);
    
    //Find a BUd that serves as a valid second guess.  Remember, this is the bisection method here.
    BUd_b = BUd_a + sign_a * 5.0;
    k_b = batch_average_k( BUd_b );
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
        k_b = batch_average_k( BUd_b );
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
            throw bright::BadFuelForm ();
    };

    BUd_c = 0.0;
    k_c = 0.0;

    //Ok now that we have valid and hopefully close initial conditions, let's do it!
    double DoA = pow(10.0, -7);	//Degree of accuracy to carry out calculations to.
    int q = 0;			//index for number of iterations
    while ( ((DoA < fabs(1.0 - k_a)) || (DoA < fabs(1.0 - k_b))) && (0.0 < fabs(BUd_a - BUd_b)) && (q < 100) )
    {
        BUd_c = (BUd_a + BUd_b) / 2.0;
        k_c = batch_average_k( BUd_c );
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
        throw bright::BisectionMethodNotPerformed ("Burnup");
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

    FluencePoint fp = fluence_at_BU(BUd);
    fd = fp.f;				//lower index of fluence at discharge
    Fd = fp.F;				//Fluence at discharge
    return;
};

void ReactorMG::run_P_NL(double temp_pnl)
{
    /** Does a reactor run for a specific P_NL.
     *  Requires that ms_feed be (meaningfully) set.
     *  For use with calibrate_P_NL_to_BUd
     */

    P_NL = temp_pnl;
    fold_mass_weights();
    BUd_bisection_method();
};

void ReactorMG::calibrate_P_NL_to_BUd()
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
            run_P_NL(pnl_a);
            bud_a = BUd;
            sign_a = (bud_a - target_BU) / fabs(bud_a - target_BU);
            FoundA = true;
        }
        catch (bright::BadFuelForm e)
        {
            pnl_a = pnl_a + 0.05;
        }
        catch (bright::BisectionMethodNotPerformed e)
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
            run_P_NL(pnl_b);
            bud_b = BUd;
            sign_b = (bud_b - target_BU) / fabs(bud_b - target_BU);
            FoundB = true;
        }
        catch (bright::BadFuelForm e)
        {
            pnl_b = pnl_b - 0.05;
        }
        catch (bright::BisectionMethodNotPerformed e)
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
        run_P_NL(pnl_c);
        bud_c = BUd;
        sign_c = (bud_c - target_BU) / fabs(bud_c - target_BU);

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



MassStream ReactorMG::calc ()
{
    //Finds BUd and output isotopics.
    fold_mass_weights();

    BUd_bisection_method();

    calc_ms_prod();

    return ms_prod;
};

MassStream ReactorMG::calc (CompDict incomp)
{
    //Finds BUd and output isotopics.
    ms_feed = MassStream (incomp);

    return calc();
};

MassStream ReactorMG::calc (MassStream instream)
{
    //Finds BUd and output isotopics.
    ms_feed = instream;

    return calc();
};


void ReactorMG::lattice_E_planar(double a, double b)
{
    lattice_E_F_.clear();
        lattice_E_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        if (0.0 == kappaC_F_[f])
            lattice_E_F_[f] = 0.0;
        else
            lattice_E_F_[f] = kappaC_F_[f] * (b - a) * bright::COTH(kappaC_F_[f]*(b-a));
    };
    return;
};

void ReactorMG::lattice_F_planar(double a, double b)
{
    lattice_F_F_.clear();
        lattice_F_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        if (0.0 == kappaF_F_[f])
            lattice_F_F_[f] = 0.0;
        else
            lattice_F_F_[f] = lattice_F_F_[f] * a * bright::COTH(lattice_F_F_[f]*a) ;
    };
    return; 
};

void ReactorMG::lattice_E_spherical(double a, double b)
{
    lattice_E_F_.clear();
        lattice_E_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        if (0.0 == kappaC_F_[f])
            lattice_E_F_[f] = 0.0;
        else
        {
            double coef = pow(kappaC_F_[f], 3) * (pow(b,3) - pow(a,3)) / (3*kappaC_F_[f]*a);
            double num = 1.0 - ( kappaC_F_[f] * b * bright::COTH(kappaC_F_[f]*(b-a)) );
            double den = 1.0 - (pow(kappaC_F_[f], 2)*a*b) - ( kappaC_F_[f]*(b-a) * bright::COTH(kappaC_F_[f]*(b-a)) );
            lattice_E_F_[f] = coef * num / den;
        };
    };
    return;
};
    
void ReactorMG::lattice_F_spherical(double a, double b)
{
    lattice_F_F_.clear();
        lattice_F_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        double coef = pow(kappaF_F_[f], 2) * pow(a, 2) / 3.0;
        double num = bright::TANH(kappaF_F_[f]*a); 
        double den = (kappaF_F_[f]*a) - bright::TANH(kappaF_F_[f]*a); 
        lattice_F_F_[f] = coef * num / den;
    };
    return; 
};

void ReactorMG::lattice_E_cylindrical(double a, double b)
{
    namespace bm = boost::math;

    lattice_E_F_.clear();
        lattice_E_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        if (0.0 == kappaC_F_[f])
            lattice_E_F_[f] = 0.0;
        else
        {
            double coef = kappaC_F_[f] * (pow(b,2) - pow(a,2)) / (2.0*a);
            double num = ( bm::cyl_bessel_i(0, kappaC_F_[f]*a) * bm::cyl_bessel_k(1, kappaC_F_[f]*b) ) + \
                ( bm::cyl_bessel_k(0, kappaC_F_[f]*a) * bm::cyl_bessel_i(1, kappaC_F_[f]*b) );
            double den = ( bm::cyl_bessel_i(1, kappaC_F_[f]*b) * bm::cyl_bessel_k(1, kappaC_F_[f]*a) ) - \
                ( bm::cyl_bessel_k(1, kappaC_F_[f]*b) * bm::cyl_bessel_i(1, kappaC_F_[f]*a) );
            lattice_E_F_[f] = coef * num / den;
        };
    };
    return;
};
    
void ReactorMG::lattice_F_cylindrical(double a, double b)
{
    namespace bm = boost::math;

    lattice_F_F_.clear();
        lattice_F_F_.assign( F.size(), 0.0 );

    for (int f = 0; f < F.size(); f++)
    {
        if (0.0 == kappaF_F_[f])
            lattice_E_F_[f] = 0.0;
        else
        {
            double num =  kappaF_F_[f] * a * bm::cyl_bessel_i(0, kappaF_F_[f]*a);
            double den = 2.0 * bm::cyl_bessel_i(1, kappaF_F_[f]*a);
            lattice_F_F_[f] = num / den;
        };
    };
    return;
};

void ReactorMG::calc_zeta()
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

    //Calculate the lattice_flag Functions
    double a, b;
    if (lattice_flag == "Planar")
    {
        a = r_fuel;
        b = pitch / 2.0;
    
        lattice_E_planar(a, b);
        lattice_F_planar(a, b);
    }
    else if (lattice_flag == "Spherical")
    {
        a = r_fuel;
        b = pitch / 2.0;
    
        lattice_E_spherical(a, b);
        lattice_F_spherical(a, b);
    }
    else if (lattice_flag == "Cylindrical")
    {
        a = r_fuel;
        b = pitch / sqrt(bright::pi); //radius of cylinder with an equivilent cell volume

        lattice_E_cylindrical(a, b);
        lattice_F_cylindrical(a, b);
    }
    else
    {
        if (0 < FCComps::verbosity)
            std::cout << "Did not specify use of planar or spheical or cylindrical lattice functions! Assuming cylindrical...\n";
        
        a = r_fuel;
        b = pitch / sqrt(bright::pi); //radius of cylinder with an equivilent cell volume

        lattice_E_cylindrical(a, b);
        lattice_F_cylindrical(a, b);
    };

    //Finally, Calculate Zeta
    zeta_F_.clear();
    zeta_F_.assign( F.size(), 0.0);
    for (int f = 0; f < F.size(); f++) 
    {
        if (0.0 == SigmaCa_F_[f])
            zeta_F_[f] = 1.0;
        else
            zeta_F_[f] = lattice_F_F_[f] + ( SigmaFa_F_[f] * VF * (lattice_E_F_[f] - 1.0) / (SigmaCa_F_[f] * VC) );
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

void ReactorMG::calc_zeta_planar()
{
    lattice_flag = "Planar";
    calc_zeta();
    return;
};

void ReactorMG::calc_zeta_spherical()
{
    lattice_flag = "Spherical";
    calc_zeta();
    return;
};

void ReactorMG::calc_zeta_cylindrical()
{
    lattice_flag = "Cylindrical";
    calc_zeta();
    return;
};
