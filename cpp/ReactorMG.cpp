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
    chemical_form_fuel = rp.fuel_form;		//Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
    chemical_form_clad = rp.cladding_form;		//Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
    chemical_form_cool = rp.coolant_form;	//Same a fuel chemical form but for coolant.  Should not have "IHM"

    rho_fuel = rp.fuel_density;			//Fuel Density
    rho_clad = rp.cladding_density;			// Cladding Density
    rho_cool = rp.coolant_density;		//Coolant Density

    P_NL = rp.pnl;				//Non-Leakage Probability
    target_BU = rp.BUt;			//Target Discharge Burnup, only used for graphing inside of this component
    burn_regions = rp.burn_regions;
    specific_power = rp.specific_power;
    burn_time = 0.0;
    burn_times = rp.burn_times;
    S = burn_times.size();

    // Flags
    use_zeta = rp.use_disadvantage_factor;		//Boolean value on whether or not the disadvantage factor should be used
    lattice_flag = rp.lattice_type;		//lattice_flagType (Planar || Spherical || Cylindrical)
    rescale_hydrogen_xs = rp.rescale_hydrogen;	//Rescale the Hydrogen-1 XS?

    // Calculates Volumes
    r_fuel = rp.fuel_radius;    // Fuel region radius
    r_void = rp.void_radius;    // Void region radius
    r_clad = rp.clad_radius;    // Clad region radius
    pitch = rp.unit_cell_pitch; // Unit cell side length

    S_O = rp.open_slots;		//Number of open slots in fuel assembly
    S_T = rp.total_slots;		//Total number of Fuel assembly slots.

    // Fuel Volume Fraction
    V_fuel = ((bright::pi * r_fuel * r_fuel)/(pitch*pitch)) * (1.0 - S_O/S_T); 

    // Cladding Volume Fraction
    V_clad = ((bright::pi * (r_clad * r_clad - r_void * r_void))/(pitch*pitch)) * (1.0 - S_O/S_T); 

    // Coolant Volume Fraction
    V_cool = ((pitch*pitch - bright::pi * r_clad * r_clad)/(pitch*pitch)) * (1.0 - S_O/S_T) + (S_O/S_T);
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

    J_size = J.size();
    J_order = std::vector<int> (J.begin(), J.end());
    std::sort(J_order.begin(), J_order.end());
    for (int j = 0; j < J_size; j++)
        J_index[J_order[j]] = j;

    // Load perturbation table
    perturbations = h5wrap::HomogenousTypeTable<double>(&rmglib, "/perturbations");
    nperturbations = perturbations.shape[0];

    // Calculate perturbed fields
    std::vector<double> col_vec;
    perturbed_fields.clear();
    for (std::vector<std::string>::iterator col = perturbations.cols.begin(); col != perturbations.cols.end(); col++)
    {
        col_vec = perturbations[*col];
        perturbed_fields[*col] = std::vector<double> (3, -1.0);

        perturbed_fields[*col][0] = *std::min_element(col_vec.begin(), col_vec.end());
        perturbed_fields[*col][1] = *std::max_element(col_vec.begin(), col_vec.end());
        perturbed_fields[*col][2] = perturbed_fields[*col][1] - perturbed_fields[*col][0];
    };

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
    sigma_t_pg.clear();
    nubar_sigma_f_pg.clear();
    chi_pg.clear();
    sigma_s_pgh.clear();
    //sigma_a_pg.clear();
    //sigma_s_pg.clear();
    //sigma_f_pg.clear();
    //nubar_pg.clear();

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
        sigma_t_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_a/" + iso_LL);
        nubar_sigma_f_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/nubar_sigma_f/" + iso_LL);
        chi_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/chi/" + iso_LL);
        sigma_s_pgh[iso_zz] = h5wrap::h5_array_to_cpp_vector_3d<double>(&rmglib, "/sigma_s_gh/" + iso_LL);
        //sigma_a_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_a/" + iso_LL);
        //sigma_s_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_s/" + iso_LL);
        //sigma_f_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_f/" + iso_LL);
        //nubar_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/nubar/" + iso_LL);
    };

    // close the reactor library
    rmglib.close();


    //
    // Create a decay matrix from a file based off of the J isotopes
    //
    std::string decay_data_file = bright::BRIGHT_DATA + "/nuc_data.h5";
    decay_matrix = std::vector< std::vector<double> > (J_size, std::vector<double>(J_size, 0.0) );

    //Check to see if the file is in HDF5 format.
    if (!bright::FileExists(decay_data_file))
        throw bright::FileNotFound(decay_data_file);

    isH5 = H5::H5File::isHdf5(decay_data_file);
    if (!isH5)
    {
        std::cout << "!!!Warning!!! " << decay_data_file << " is not a valid HDF5 file!\n";
        return;
    };

    // Read in the decay data table as an array of FCComps::decay_iso_desc
    H5::H5File decay_data_h5 (decay_data_file.c_str(), H5F_ACC_RDONLY );
    H5::DataSet decay_data_set = decay_data_h5.openDataSet("/decay");
    H5::DataSpace decay_data_space = decay_data_set.getSpace();
    int decay_data_length = decay_data_space.getSimpleExtentNpoints(); 

    FCComps::decay_iso_stuct * decay_data_array = new FCComps::decay_iso_stuct [decay_data_length];
    decay_data_set.read(decay_data_array, FCComps::decay_iso_desc);

    // Make decay_martrix from this data.
    int i, j, ind, jnd;
    for (int l = 0; l < decay_data_length; l++)
    {
        i = decay_data_array[l].from_iso_zz;
        j = decay_data_array[l].to_iso_zz;

        // skip non-element from-isos
        if (J.count(i) < 1)
            continue;

        // skip non-element to-isos
        if (J.count(j) < 1)
            continue;

        // Get the indexes for these nulcides into the matrix
        ind = J_index[i];
        jnd = J_index[j];

        // Add diagonal elements
        decay_matrix[ind][ind] = -decay_data_array[l].decay_const;

        // Add i,j element to matrix
        if (i != j)
            decay_matrix[ind][jnd] = decay_data_array[l].branch_ratio * decay_data_array[l].decay_const;
    };

    return;
};





void ReactorMG::calc_nearest_neighbors()
{
    /**
     * Returns a vector of the indices sorted, sorted by nearest neighbor
     */

    // Initialize
    std::map<std::string, std::vector<double> > deltas;
    std::vector<double> rss (nperturbations, 0.0);
    nearest_neighbors.clear();

    // Calcumlate normailized deltas
    if (perturbed_fields["fuel_density"][2] != 0.0)
        deltas["fuel_density"] = bright::delta_vector(rho_fuel, perturbations["fuel_density"]);

    if (perturbed_fields["clad_density"][2] != 0.0)
        deltas["clad_density"] = bright::delta_vector(rho_clad, perturbations["clad_density"]);

    if (perturbed_fields["cool_density"][2] != 0.0)
        deltas["cool_density"] = bright::delta_vector(rho_cool, perturbations["cool_density"]);


    if (perturbed_fields["fuel_cell_radius"][2] != 0.0)
        deltas["fuel_cell_radius"] = bright::delta_vector(r_fuel, perturbations["fuel_cell_radius"]);

    if (perturbed_fields["void_cell_radius"][2] != 0.0)
        deltas["void_cell_radius"] = bright::delta_vector(r_void, perturbations["void_cell_radius"]);

    if (perturbed_fields["clad_cell_radius"][2] != 0.0)
        deltas["clad_cell_radius"] = bright::delta_vector(r_clad, perturbations["clad_cell_radius"]);


    if (perturbed_fields["unit_cell_pitch"][2] != 0.0)
        deltas["unit_cell_pitch"] = bright::delta_vector(pitch, perturbations["unit_cell_pitch"]);

    if (perturbed_fields["burn_regions"][2] != 0.0)
        deltas["burn_regions"] = bright::delta_vector(burn_regions, perturbations["burn_regions"]);

    if (perturbed_fields["fuel_specific_power"][2] != 0.0)
        deltas["fuel_specific_power"] = bright::delta_vector(specific_power, perturbations["fuel_specific_power"]);


    // Calc pertubations for initial mass streams
    if (10 < perturbations.shape[1])
    {
        int iso_zz;
        std::string iso_LL;
        double iso_mass;

        for (int p = 9; p < perturbations.shape[1] - 1; p++)
        {
            // Grab some names
            iso_LL = perturbations.cols[p];
            iso_zz = isoname::LLAAAM_2_zzaaam(iso_LL);

            // Determine the mass of the isotope in the feed
            if (0 < ms_feed.comp.count(iso_zz))
                iso_mass = ms_feed.comp[iso_zz];
            else
                iso_mass = 0.0;

            // Calculate the delta if appropriate.
            if (perturbed_fields[iso_LL][2] != 0.0)
                deltas[iso_LL] = bright::delta_vector(iso_mass, perturbations[iso_LL]);                
        };
    };

    if (perturbed_fields["burn_times"][2] != 0.0)
        deltas["burn_times"] = bright::delta_vector(burn_time, perturbations["burn_times"]);


    // Now that we have the normalized deltas for each index
    // We want to calculate the 'distance' from the curent point
    // This is done by taking the root of the sum of squares for 
    // each index.
    std::string key;
    double norm_factor_sqrd;
    for (std::map<std::string, std::vector<double> >::iterator d = deltas.begin(); d != deltas.end(); d++)
    {
        key = (*d).first;
        norm_factor_sqrd = (perturbed_fields[key][2] * perturbed_fields[key][2]);

        // Take the sum of the squares
        for (int q = 0; q < perturbations.shape[0]; q++)
        {
            rss[q] = rss[q] + (deltas[key][q] * deltas[key][q] / norm_factor_sqrd);
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












void ReactorMG::interpolate_cross_sections()
{
    // Grab the nearest and next nearest neighbor maps
    int a0 = nearest_neighbors[0]; 
    int a1 = nearest_neighbors[1]; 
    std::map<std::string, double> nn0 = perturbations[a0];
    std::map<std::string, double> nn1 = perturbations[a1];

    // Calculate the x-factors to interpolate against.
    // For every variable that is perturbed, 
    // this factor is the same for all y variables
    // and is equal to:
    //
    //     (xa - xa1)     (xb - xb1)
    //     ----------  x  ----------  ...
    //     (xa2 - xa1)    (xb2 - xb2)
    //
    // Then the value of the interpolation is 
    // calculated via
    //
    // y = x_factor * (y2 - y1) + y1
    //
    double x_factor = 0.0;    

    if (nn0["fuel_density"] != nn1["fuel_density"])
        x_factor = x_factor + ((rho_fuel - nn0["fuel_density"])/(nn1["fuel_density"] - nn0["fuel_density"]));

    if (nn0["clad_density"] != nn1["clad_density"])
        x_factor = x_factor + ((rho_clad - nn0["clad_density"])/(nn1["clad_density"] - nn0["clad_density"]));

    if (nn0["cool_density"] != nn1["cool_density"])
        x_factor = x_factor + ((rho_cool - nn0["cool_density"])/(nn1["cool_density"] - nn0["cool_density"]));


    if (nn0["fuel_cell_radius"] != nn1["fuel_cell_radius"])
        x_factor = x_factor + ((r_fuel - nn0["fuel_cell_radius"])/(nn1["fuel_cell_radius"] - nn0["fuel_cell_radius"]));

    if (nn0["void_cell_radius"] != nn1["void_cell_radius"])
        x_factor = x_factor + ((r_void - nn0["void_cell_radius"])/(nn1["void_cell_radius"] - nn0["void_cell_radius"]));

    if (nn0["clad_cell_radius"] != nn1["clad_cell_radius"])
        x_factor = x_factor + ((r_clad - nn0["clad_cell_radius"])/(nn1["clad_cell_radius"] - nn0["clad_cell_radius"]));


    if (nn0["unit_cell_pitch"] != nn1["unit_cell_pitch"])
        x_factor = x_factor + ((pitch - nn0["unit_cell_pitch"])/(nn1["unit_cell_pitch"] - nn0["unit_cell_pitch"]));

    if (nn0["burn_regions"] != nn1["burn_regions"])
        x_factor = x_factor + ((burn_regions - nn0["burn_regions"])/(nn1["burn_regions"] - nn0["burn_regions"]));

    if (nn0["fuel_specific_power"] != nn1["fuel_specific_power"])
        x_factor = x_factor + ((specific_power - nn0["fuel_specific_power"])/(nn1["fuel_specific_power"] - nn0["fuel_specific_power"]));


    // Calc x-factor for initial mass streams
    if (10 < perturbations.shape[1])
    {
        int iso_zz;
        std::string iso_LL;
        double iso_mass;

        for (int p = 9; p < perturbations.shape[1] - 1; p++)
        {
            // Grab some names
            iso_LL = perturbations.cols[p];
            iso_zz = isoname::LLAAAM_2_zzaaam(iso_LL);

            // Determine the mass of the isotope in the feed
            if (0 < ms_feed.comp.count(iso_zz))
                iso_mass = ms_feed.comp[iso_zz];
            else
                iso_mass = 0.0;

            // Calculate the x-factor if appropriate.
            if (nn0[iso_LL] != nn1[iso_LL])
                x_factor = x_factor + ((iso_mass - nn0[iso_LL])/(nn1[iso_LL] - nn0[iso_LL]));
        };
    };

    if (nn0["burn_times"] != nn1["burn_times"])
        x_factor = x_factor + ((burn_time - nn0["burn_times"])/(nn1["burn_times"] - nn0["burn_times"]));



    // Now that we have found the x-factor, we get to do the actual interpolations. Oh Joy!
    for (iso_iter iso = J.begin(); iso != J.end(); iso++)
    {
        // Interpolate the cross-sections
        sigma_t_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_t_pg[*iso][a1], sigma_t_pg[*iso][a0]);
        nubar_sigma_f_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, nubar_sigma_f_pg[*iso][a1], nubar_sigma_f_pg[*iso][a0]);
        chi_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, chi_pg[*iso][a1], chi_pg[*iso][a0]);
        //sigma_a_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_a_pg[*iso][a1], sigma_a_pg[*iso][a0]);
        //sigma_s_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_s_pg[*iso][a1], sigma_s_pg[*iso][a0]);
        //sigma_f_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_f_pg[*iso][a1], sigma_f_pg[*iso][a0]);
        //nubar_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, nubar_pg[*iso][a1], nubar_pg[*iso][a0]);

        for (int g = 0; g < G; g++)
            sigma_s_itgh[*iso][bt_s][g] = bright::y_x_factor_interpolation(x_factor, sigma_s_pgh[*iso][a1][g], sigma_s_pgh[*iso][a0][g]);
    };
    
};









void ReactorMG::calc_mass_weights()
{
    /** 
     *  Calculates the appropriate mass fractions for this time step
     */

    // First things first, let's calculate the atomic weight of the HM
    double inverse_A_HM = 0.0;
    double mass_HM = 0.0;
    for (iso_iter iso = J.begin(); iso != J.end(); iso++)
    {
        mass_HM += T_it[*iso][bt_s];
        inverse_A_HM += (T_it[*iso][bt_s] / isoname::nuc_weight(*iso));
    };
    A_HM_t[bt_s] = mass_HM / inverse_A_HM;


    //
    // Calculate the molecular weight value
    //
    int key_zz;

    // Fuel Molecular Weight
    for (std::map<std::string, double>::iterator key = chemical_form_fuel.begin(); key != chemical_form_fuel.end(); key++)
    {
        if ( (key->first) == "IHM")
            MW_fuel_t[bt_s] += chemical_form_fuel[key->first] * A_HM_t[bt_s];
        else
        {
            key_zz = isoname::mixed_2_zzaaam(key->first);
            MW_fuel_t[bt_s] += chemical_form_fuel[key->first] * isoname::nuc_weight(key_zz);
        };
    };

    // Cladding Molecular Weight
    for (std::map<std::string, double>::iterator key = chemical_form_clad.begin(); key != chemical_form_clad.end(); key++)
    {
        key_zz = isoname::mixed_2_zzaaam(key->first);
        MW_clad_t[bt_s] += chemical_form_clad[key->first] * isoname::nuc_weight(key_zz);
    };

    // Coolant Molecular Weight
    for (std::map<std::string, double>::iterator key = chemical_form_cool.begin(); key != chemical_form_cool.end(); key++)
    {
        key_zz = isoname::mixed_2_zzaaam(key->first);
        MW_cool_t[bt_s] += chemical_form_cool[key->first] * isoname::nuc_weight(key_zz);
    };



    //
    // Build the atom number density dictionaries
    //

    // now for the n_it in the Fuel
    for (std::map<std::string, double>::iterator key = chemical_form_fuel.begin(); key != chemical_form_fuel.end(); key++)
    {
        if ( (key->first) == "IHM")
        {
            for (iso_iter iso = J.begin(); iso != J.end(); iso++)
                n_fuel_it[*iso][bt_s] += chemical_form_fuel[key->first] * T_it[*iso][bt_s];
        }
        else
        {
            key_zz = isoname::mixed_2_zzaaam(key->first);
            n_fuel_it[key_zz][bt_s] += chemical_form_fuel[key->first];
        }
    };

    // Note that the n_it in the cladding is just chemical_form_clad
    for (std::map<std::string, double>::iterator key = chemical_form_clad.begin(); key != chemical_form_clad.end(); key++)
    {
        key_zz = isoname::mixed_2_zzaaam(key->first);
        n_clad_it[key_zz][bt_s] = chemical_form_clad[key->first];
    };

    // Note that the n_it in the coolant is just chemical_form_cool
    for (std::map<std::string, double>::iterator key = chemical_form_cool.begin(); key != chemical_form_cool.end(); key++)
    {
        key_zz = isoname::mixed_2_zzaaam(key->first);
        n_cool_it[key_zz][bt_s] = chemical_form_cool[key->first];
    };



    //
    // Build the number density dictionaries
    //
    for (iso_iter iso = J.begin(); iso != J.end(); iso++)
    {
        // Fuel Number Density
        N_fuel_it[*iso][bt_s] = n_fuel_it[*iso][bt_s] * rho_fuel * (bright::N_A / MW_fuel_t[bt_s]);

        // Cladding Number Density
        N_clad_it[*iso][bt_s] = n_clad_it[*iso][bt_s] * rho_clad * (bright::N_A / MW_clad_t[bt_s]);

        // Coolant Number Density
        N_cool_it[*iso][bt_s] = n_cool_it[*iso][bt_s] * rho_cool * (bright::N_A / MW_cool_t[bt_s]);
    };



    //
    // Build the mass weight dictionaries
    //
    double iso_weight;
    double relative_volume_clad = (rho_clad * MW_fuel_t[bt_s] * V_clad) / (rho_fuel * MW_clad_t[bt_s] * V_fuel);
    double relative_volume_cool = (rho_cool * MW_fuel_t[bt_s] * V_cool) / (rho_fuel * MW_cool_t[bt_s] * V_fuel);

    for (iso_iter iso = J.begin(); iso != J.end(); iso++)
    {
        iso_weight = isoname::nuc_weight(*iso);

        // Fuel mass weight
        m_fuel_it[*iso][bt_s] =  n_fuel_it[*iso][bt_s] * iso_weight / A_HM_t[bt_s];

        // Cladding mass weight 
        m_clad_it[*iso][bt_s] = (n_clad_it[*iso][bt_s] * iso_weight / A_HM_t[bt_s]) * relative_volume_clad;

        // Coolant mass weight 
        m_cool_it[*iso][bt_s] = (n_cool_it[*iso][bt_s] * iso_weight / A_HM_t[bt_s]) * relative_volume_cool;
    };
};







void ReactorMG::fold_mass_weights()
{
    // Folds mass weight in with cross-sections for current time step

    double V_total = V_fuel + V_clad + V_cool;
    double V_frac_fuel = V_fuel / V_total;
    double V_frac_clad = V_clad / V_total;
    double V_frac_cool = V_cool / V_total;

    double n_i = 0.0;
    double N_i_cm2pb = 0.0;

    for (iso_iter iso = J.begin(); iso != J.end(); iso++)
    {
        // Calculate the core-average number density for this isotope.
        n_i = ((n_fuel_it[*iso][bt_s] * V_frac_fuel) \
            +  (n_clad_it[*iso][bt_s] * V_frac_clad) \
            +  (n_cool_it[*iso][bt_s] * V_frac_cool));

        N_i_cm2pb = bright::cm2_per_barn * ((N_fuel_it[*iso][bt_s] * V_frac_fuel) \
                                         +  (N_clad_it[*iso][bt_s] * V_frac_clad) \
                                         +  (N_cool_it[*iso][bt_s] * V_frac_cool));
        // Loop over all groups
        for (int g = 0; g < G; g++)
        {
            Sigma_t_tg[bt_s][g] += N_i_cm2pb * sigma_t_itg[*iso][bt_s][g];
            nubar_Sigma_f_tg[bt_s][g] += N_i_cm2pb * nubar_sigma_f_itg[*iso][bt_s][g];
// TESTME FIXME  I possibly neede to think more about the 
// Normalization factor here for chi
            chi_tg[bt_s][g] += n_i * chi_itg[*iso][bt_s][g];
            //Sigma_a_tg[bt_s][g] += N_i_cm2pb * sigma_a_itg[*iso][bt_s][g];
            //Sigma_s_tg[bt_s][g] += N_i_cm2pb * sigma_s_itg[*iso][bt_s][g];
            //Sigma_f_tg[bt_s][g] += N_i_cm2pb * sigma_f_itg[*iso][bt_s][g];
            //nubar_tg[bt_s][g] += N_i_cm2pb * nubar_itg[*iso][bt_s][g];

            for (int h =0; h < G; h++)
                Sigma_s_tgh[bt_s][g][h] += N_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];
        };
    };
};





void ReactorMG::assemble_multigroup_matrices()
{
    // Assembles the cross section matrices needed for multigroup 
    // Burnup-criticality calculations.

    // Assemble the A matrix 
    for (int g = 0; g < G; g++)
    {
        // Add the total cross section
        A_tgh[bt_s][g][g] += Sigma_t_tg[bt_s][g];

        // Subtract the scattering kernel
        for (int h = 0; h < G; h++)
            A_tgh[bt_s][g][h] -= Sigma_s_tgh[bt_s][g][h];
    };


    // Assemble the F matrix
    F_tgh[bt_s] = bright::vector_outer_product(chi_tg[bt_s], nubar_Sigma_f_tg[bt_s]);

    // Grab the inverse of the A matrix
    A_inv_tgh[bt_s] = bright::matrix_inverse(A_tgh[bt_s]);

    // Multiply the inverse of A by F
    A_inv_F_tgh[bt_s] = bright::matrix_multiplication(A_inv_tgh[bt_s], F_tgh[bt_s]);
};






void ReactorMG::calc_criticality()
{
    // Init values
    int n = 0;
    int N = 100;

    float epsik = 1.0;
    float tmp_epsiphi; 
    float epsiphi = 1.0;
    float epsilon = 0.005;

    double k0 = 1.0;
    std::vector<double> phi0 (1.0, G);
    double k1;
    std::vector<double> phi1;

    int g = 0;
    double invPk;
    double nu_Sigma_f_phi0;
    double nu_Sigma_f_phi1;

    // Solve for k and phi simeltaneoulsy
    while ((n < N) && ((epsilon < epsik) || (epsilon < epsiphi)))
    {
        // Calculate the next eigen-flux
        invPk = 1.0 / (P_NL * k0);
        phi1 = bright::scalar_matrix_vector_product(invPk, A_inv_F_tgh[bt_s], phi0);

        // Calculate the next eigen-k
        nu_Sigma_f_phi0 = 0.0;
        nu_Sigma_f_phi1 = 0.0;
        for (g = 0; g < G; g++)
        {
            nu_Sigma_f_phi0 += nubar_Sigma_f_tg[bt_s][g] * phi0[g]; 
            nu_Sigma_f_phi1 += nubar_Sigma_f_tg[bt_s][g] * phi1[g]; 
        };
        k1 = k0 * nu_Sigma_f_phi1 / nu_Sigma_f_phi0;

        // Calculate the epsilon value of k 
        epsik = fabs(1.0 - (k0/k1));

        // Calulate the maximum epsilon of phi over all groups
        epsiphi = fabs(1.0 - (phi0[0]/phi1[0]));
        for (g = 1; g < G; g++)
        {
            tmp_epsiphi = fabs(1.0 - (phi0[g]/phi1[g]));
            if (epsiphi < tmp_epsiphi)
                epsiphi = tmp_epsiphi;
        };

        // Set the next eigens to the previous values befor looping
        k0 = k1;
        phi0 = phi1;
    };

    // Set the final values to the class members
    k_t[bt_s] = k1;
    phi_tg[bt_s] = phi1;

    phi_t[bt_s] = 0.0;
    for (g = 0; g < G; g++)
        phi_t[bt_s] += phi1[0];

    Phi_t[bt_s] = phi_t[bt_s] * (burn_times[bt_s] - burn_times[bt_s]) * bright::sec_per_day;
};








void ReactorMG::burnup_core()
{
    // Burns up the core and fills in parameter values as we go.


    // Initialize the transmutation matrix with values from ms_feed
    T_it.clear();

    // Also initialize the cross-section matrices as a function of time.
    sigma_t_itg.clear();
    nubar_sigma_f_itg.clear();
    chi_itg.clear();
    sigma_s_itgh.clear();
    //sigma_a_itg.clear();
    //sigma_s_itg.clear();
    //sigma_f_itg.clear();
    //nubar_itg.clear();

    // Also initilaize the mass weights
    A_HM_t = std::vector<double>(S, 0.0);
    MW_fuel_t = std::vector<double>(S, 0.0);
    MW_clad_t = std::vector<double>(S, 0.0);
    MW_cool_t = std::vector<double>(S, 0.0);
    n_fuel_it.clear();
    n_clad_it.clear();
    n_cool_it.clear();
    m_fuel_it.clear();
    m_clad_it.clear();
    m_cool_it.clear();
    N_fuel_it.clear();
    N_clad_it.clear();
    N_cool_it.clear();

    for (iso_iter iso = J.begin(); iso != J.end(); iso++)
    {
        // Init the transmutation matrix
        T_it[*iso] = time_data(S, -1.0);

        if (0 < ms_feed.comp.count(*iso))
            T_it[*iso][0] = ms_feed.comp[*iso];
        else
            T_it[*iso][0] = 0.0;

        // Init the cross-sections
        sigma_t_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        nubar_sigma_f_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        chi_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        sigma_s_itgh[*iso] = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, -1.0)));
        //sigma_a_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        //sigma_s_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        //sigma_f_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        //nubar_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));

        // Init the mass weights
        n_fuel_it[*iso] = std::vector<double>(S, 0.0);
        n_clad_it[*iso] = std::vector<double>(S, 0.0);
        n_cool_it[*iso] = std::vector<double>(S, 0.0);
        m_fuel_it[*iso] = std::vector<double>(S, 0.0);
        m_clad_it[*iso] = std::vector<double>(S, 0.0);
        m_cool_it[*iso] = std::vector<double>(S, 0.0);
        N_fuel_it[*iso] = std::vector<double>(S, 0.0);
        N_clad_it[*iso] = std::vector<double>(S, 0.0);
        N_cool_it[*iso] = std::vector<double>(S, 0.0);
    };

    // Init the macroscopic cross sections
    Sigma_t_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    nubar_Sigma_f_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    chi_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_s_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    //Sigma_a_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    //Sigma_s_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    //Sigma_f_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    //nubar_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));


    // Init the criticality matrices
    A_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    F_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    A_inv_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    A_inv_F_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));


    // Loop through all time steps
    for (int s = 0; s < S; s++)
    {
        // Set the current time
        bt_s = s;
        burn_time = burn_times[s];

        // Find the nearest neightbors for this time.
        calc_nearest_neighbors();

        // Interpolate cross section in preparation for 
        // criticality calculation.
        interpolate_cross_sections();

        // Fold the mass weights for this time step
        calc_mass_weights();
        fold_mass_weights();

        // Make the cross-section matrices before the criticality calulation
        assemble_multigroup_matrices();
        calc_criticality();
    };

};






void ReactorMG::old_burnup_core()
{
};
/*
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
*/





void ReactorMG::calc_T_itd()
{
};
    /** Calculates the output isotopics of Mj(Fd).
     *  NOTE: Mj(Fd) is effectively the same variable as ms_prod before normalization!
     */

/*
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
*/

void ReactorMG::calc_ms_prod()
{
    //Wrapper to calculate discharge isotopics.
    calc_T_itd();
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




double ReactorMG::calc_tru_cr()
{
    //Calculates the reactor's transuranic conversion ratio.
    tru_cr = 1.0 - ((ms_feed_tru.mass - ms_prod_tru.mass) / (BUd/931.46));
    return tru_cr;
};





double ReactorMG::calc_deltaR()
{
};
/*
    //Calculates the deltaR of the reactor with the current ms_feed
    fold_mass_weights();
    deltaR = batch_average(target_BU, "P") - batch_average(target_BU, "D");
    return deltaR;
};
*/

double ReactorMG::calc_deltaR(CompDict cd)
{
};
/*
    //Calculates the deltaR of the reactor with the current ms_feed
    ms_feed = MassStream (cd);
    return calc_deltaR();
};
*/

double ReactorMG::calc_deltaR(MassStream ms)
{
};
/*
    //Calculates the deltaR of the reactor with the current ms_feed
    ms_feed = ms;
    return calc_deltaR();
};
*/


FluencePoint ReactorMG::fluence_at_BU(double BU)
{
};
    /** Gives the fluence at which the burnup BU occurs.
     *  Data returned as FluenceIndex stucture:
     *  	FI.f: index imeadiately lower than where BU achieved (int),
     *  	FI.F: flunece itself (double),
     *  	FI.m: slope dBU/dF bewteen points f and f+1 (double)
     */

/*
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
*/


double ReactorMG::batch_average(double BUd, std::string PDk_flag)
{
};
/*
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
*/

double ReactorMG::batch_average_k(double BUd) 
{
};
/*
        return batch_average(BUd, "K");
};
*/


void ReactorMG::BUd_bisection_method()
{
};
/*
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
*/

void ReactorMG::run_P_NL(double temp_pnl)
{
};
    /** Does a reactor run for a specific P_NL.
     *  Requires that ms_feed be (meaningfully) set.
     *  For use with calibrate_P_NL_to_BUd
     */

/*
    P_NL = temp_pnl;
    fold_mass_weights();
    BUd_bisection_method();
};
*/


void ReactorMG::calibrate_P_NL_to_BUd()
{
};
    /** Calibrates the non-leakage probability of a reactors to hit a target burnup.
     *  Calibration proceeds by bisection method...
     */
/*
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
*/



MassStream ReactorMG::calc()
{
};
/*
    //Finds BUd and output isotopics.
    fold_mass_weights();

    BUd_bisection_method();

    calc_ms_prod();

    return ms_prod;
};
*/


MassStream ReactorMG::calc (CompDict incomp)
{
};
/*
    //Finds BUd and output isotopics.
    ms_feed = MassStream (incomp);

    return calc();
};
*/


MassStream ReactorMG::calc (MassStream instream)
{
};
/*
    //Finds BUd and output isotopics.
    ms_feed = instream;

    return calc();
};
*/


void ReactorMG::lattice_E_planar(double a, double b)
{
};
/*
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
*/



void ReactorMG::lattice_F_planar(double a, double b)
{
};
/*
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
*/


void ReactorMG::lattice_E_spherical(double a, double b)
{
};
/*
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
*/
  
  
void ReactorMG::lattice_F_spherical(double a, double b)
{
};
/*
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
*/

void ReactorMG::lattice_E_cylindrical(double a, double b)
{
};
/*
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
*/  
  
void ReactorMG::lattice_F_cylindrical(double a, double b)
{
};
/*
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
*/

void ReactorMG::calc_zeta()
{
};
/*
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
                SigmaFa_F_[f]  = SigmaFa_F_[f]  + (NiF[iso->first] * di_F_[iso->first][f] * bright::cm2_per_barn);

                SigmaFtr_F_[f] = SigmaFtr_F_[f] + (NiF[iso->first] * bright::cm2_per_barn * (di_F_[iso->first][f] + \
                    sigma_s_therm[iso->first]*(1.0 - 2.0/(3.0*isoname::nuc_weight(iso->first))) ) );
            }
            else
            {
                //renormalize sigma_a for this fluenece
                double sig_a = sigma_a_therm[iso->first] * di_F_[iso->first][f] / di_F_[iso->first][0];

                SigmaFa_F_[f]  = SigmaFa_F_[f]  + (NiF[iso->first] * sig_a * bright::cm2_per_barn);

                SigmaFtr_F_[f] = SigmaFtr_F_[f] + (NiF[iso->first] * bright::cm2_per_barn * (sig_a + \
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
                SigmaCa_F_[f]  = SigmaCa_F_[f]  + (NiC[iso->first] * di_F_[iso->first][f] * bright::cm2_per_barn);

                SigmaCtr_F_[f] = SigmaCtr_F_[f] + (NiC[iso->first] * bright::cm2_per_barn * (di_F_[iso->first][f] + \
                    sigma_s_therm[iso->first]*(1.0 - 2.0/(3.0*isoname::nuc_weight(iso->first))) ) );
            }
            else
            {
                //renormalize sigma_a for this fluenece
                double sig_a = sigma_a_therm[iso->first] * di_F_[iso->first][f] / di_F_[iso->first][0];

                SigmaCa_F_[f]  = SigmaCa_F_[f]  + (NiC[iso->first] * sig_a * bright::cm2_per_barn);

                SigmaCtr_F_[f] = SigmaCtr_F_[f] + (NiC[iso->first] * bright::cm2_per_barn * (sig_a + \
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
*/

void ReactorMG::calc_zeta_planar()
{
};
/*
    lattice_flag = "Planar";
    calc_zeta();
    return;
};
*/

void ReactorMG::calc_zeta_spherical()
{
};
/*
    lattice_flag = "Spherical";
    calc_zeta();
    return;
};
*/

void ReactorMG::calc_zeta_cylindrical()
{
};
/*
    lattice_flag = "Cylindrical";
    calc_zeta();
    return;
};
*/
