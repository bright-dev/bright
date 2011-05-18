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
    sigma_f_pg.clear();
    sigma_gamma_pg.clear();
    sigma_2n_pg.clear();
    sigma_3n_pg.clear();
    sigma_alpha_pg.clear();
    sigma_proton_pg.clear();
    sigma_gamma_x_pg.clear();
    sigma_2n_x_pg.clear();

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
        sigma_t_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_t/" + iso_LL);
        nubar_sigma_f_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/nubar_sigma_f/" + iso_LL);
        chi_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/chi/" + iso_LL);
        sigma_s_pgh[iso_zz] = h5wrap::h5_array_to_cpp_vector_3d<double>(&rmglib, "/sigma_s_gh/" + iso_LL);
        sigma_f_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_f/" + iso_LL);
        sigma_gamma_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_gamma/" + iso_LL);
        sigma_2n_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_2n/" + iso_LL);
        sigma_3n_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_3n/" + iso_LL);
        sigma_alpha_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_alpha/" + iso_LL);
        sigma_proton_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_proton/" + iso_LL);
        sigma_gamma_x_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_gamma_x/" + iso_LL);
        sigma_2n_x_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(&rmglib, "/sigma_2n_x/" + iso_LL);
    };

    // close the reactor library
    rmglib.close();


    //
    // Create a decay matrix from a file based off of the J isotopes
    //
    std::string nuc_data_file = bright::BRIGHT_DATA + "/nuc_data.h5";

    //Check to see if the file is in HDF5 format.
    if (!bright::FileExists(nuc_data_file))
        throw bright::FileNotFound(nuc_data_file);

    isH5 = H5::H5File::isHdf5(nuc_data_file);
    if (!isH5)
    {
        std::cout << "!!!Warning!!! " << nuc_data_file << " is not a valid HDF5 file!\n";
        return;
    };

    // Open the HDF5 file
    H5::H5File nuc_data_h5 (nuc_data_file.c_str(), H5F_ACC_RDONLY );


    //
    // Read in the decay data table as an array of FCComps::decay_iso_desc
    //
    H5::DataSet decay_data_set = nuc_data_h5.openDataSet("/decay");
    H5::DataSpace decay_data_space = decay_data_set.getSpace();
    int decay_data_length = decay_data_space.getSimpleExtentNpoints(); 

    FCComps::decay_iso_stuct * decay_data_array = new FCComps::decay_iso_stuct [decay_data_length];
    decay_data_set.read(decay_data_array, FCComps::decay_iso_desc);


    // Make decay_martrix from this data.
    decay_matrix = std::vector< std::vector<double> > (J_size, std::vector<double>(J_size, 0.0) );

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


    //
    // Read in the fission table
    //
    H5::DataSet fission_set = nuc_data_h5.openDataSet("/neutron/xs_mg/fission");
    H5::DataSpace fission_space = fission_set.getSpace();
    int fission_length = fission_space.getSimpleExtentNpoints(); 

    FCComps::fission_struct * fission_array = new FCComps::fission_struct [fission_length];
    fission_set.read(fission_array, FCComps::fission_desc);

    // Run through the array and make join maps
    //  key = fission index
    //  value = vector of J_indexs
    std::map<int, std::vector<int> > thermal_join;
    std::map<int, std::vector<int> > fast_join;

    int ty, fy;
    for (int l = 0; l < fission_length; l++)
    {
        i = fission_array[l].iso_zz;

        // skip non-element from-isos
        if (J.count(i) < 1)
            continue;


        // make thermal join
        ty = fission_array[l].thermal_yield;

        if (thermal_join.count(ty) < 1)
            thermal_join[ty] = std::vector<int>();

        thermal_join[ty].push_back(J_index[i]);


        // make fast join
        fy = fission_array[l].fast_yield;

        if (fast_join.count(fy) < 1)
            fast_join[fy] = std::vector<int>();

        fast_join[ty].push_back(J_index[i]);
    };


    // Read in fission product yeilds
    H5::DataSet fp_yields_set = nuc_data_h5.openDataSet("/neutron/fission_products/yields");
    H5::DataSpace fp_yields_space = fp_yields_set.getSpace();
    int fp_yields_length = fp_yields_space.getSimpleExtentNpoints(); 

    FCComps::fission_product_yields_struct * fp_yields_array = new FCComps::fission_product_yields_struct [fp_yields_length];
    fp_yields_set.read(fp_yields_array, FCComps::fission_product_yields_desc);


    // Run through the array and make yield matrices
    thermal_yield_matrix = std::vector< std::vector<double> > (J_size, std::vector<double>(J_size, 0.0) );
    fast_yield_matrix = std::vector< std::vector<double> > (J_size, std::vector<double>(J_size, 0.0) );

    int index, tj, fj, TJ, FJ;
    double mf;
    for (int l = 0; l < fp_yields_length; l++)
    {
        // Get important data from struct
        index = fp_yields_array[l].index;
        j = fp_yields_array[l].to_iso_zz;
        jnd = J_index[j];
        mf = fp_yields_array[l].mass_frac;

        // Add to thermal yields
        if (0 < thermal_join.count(index))
        {
            TJ = thermal_join[index].size();
            for (tj = 0; tj < TJ; tj++)
            {
                ind = thermal_join[index][tj];
                thermal_yield_matrix[ind][jnd] = mf;
            };
        };

        // Add to fast yields.
        if (0 < fast_join.count(index))
        {
            FJ = fast_join[index].size();
            for (fj = 0; fj < FJ; fj++)
            {
                ind = fast_join[index][fj];
                fast_yield_matrix[ind][jnd] = mf;
            };
        };
    };


    // Make fission product yield matrix
    fission_product_yield_matrix = std::vector< std::vector< std::vector<double> > > (J_size, std::vector< std::vector<double> >(J_size, std::vector<double>(G, 0.0) ) );
    for (ind = 0; ind < J_size; ind++)
    {
        for (jnd = 0; jnd < J_size; jnd++)
        {
            for (int g = 0; g < G; g++)
            {
                // Test which regime we are in
                if (0.1 < E_g[g])
                    fission_product_yield_matrix[ind][jnd][g] = fast_yield_matrix[ind][jnd];
                else
                    fission_product_yield_matrix[ind][jnd][g] = thermal_yield_matrix[ind][jnd];
            };
        };
    };


    // close the nuc_data library
    nuc_data_h5.close();

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
        std::string iso_col;
        double iso_mass;

        for (int p = 9; p < perturbations.shape[1] - 1; p++)
        {
            // Grab some names
            iso_col = perturbations.cols[p];
            iso_LL = perturbations.cols[p];
            iso_LL.replace(0, 8, "");
            iso_zz = isoname::LLAAAM_2_zzaaam(iso_LL);

            // Determine the mass of the isotope in the feed
            if (0 < ms_feed.comp.count(iso_zz))
                iso_mass = ms_feed.comp[iso_zz];
            else
                iso_mass = 0.0;

            // Calculate the delta if appropriate.
            if (perturbed_fields[iso_col][2] != 0.0)
                deltas[iso_col] = bright::delta_vector(iso_mass, perturbations[iso_col]);                
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
        std::string iso_col;
        double iso_mass;

        for (int p = 9; p < perturbations.shape[1] - 1; p++)
        {
            // Grab some names
            iso_col = perturbations.cols[p];
            iso_LL = perturbations.cols[p];
            iso_LL.replace(0, 8, "");
            iso_zz = isoname::LLAAAM_2_zzaaam(iso_LL);

            // Determine the mass of the isotope in the feed
            if (0 < ms_feed.comp.count(iso_zz))
                iso_mass = ms_feed.comp[iso_zz];
            else
                iso_mass = 0.0;

            // Calculate the x-factor if appropriate.
            if (nn0[iso_col] != nn1[iso_col])
                x_factor = x_factor + ((iso_mass - nn0[iso_col])/(nn1[iso_col] - nn0[iso_col]));
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
        sigma_f_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_f_pg[*iso][a1], sigma_f_pg[*iso][a0]);
        sigma_gamma_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_gamma_pg[*iso][a1], sigma_gamma_pg[*iso][a0]);
        sigma_2n_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_2n_pg[*iso][a1], sigma_2n_pg[*iso][a0]);
        sigma_3n_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_3n_pg[*iso][a1], sigma_3n_pg[*iso][a0]);
        sigma_alpha_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_alpha_pg[*iso][a1], sigma_alpha_pg[*iso][a0]);
        sigma_proton_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_proton_pg[*iso][a1], sigma_proton_pg[*iso][a0]);
        sigma_gamma_x_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_gamma_x_pg[*iso][a1], sigma_gamma_x_pg[*iso][a0]);
        sigma_2n_x_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_2n_x_pg[*iso][a1], sigma_2n_x_pg[*iso][a0]);

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
    int g;
    double mu;
    double nuc_weight;

    double V_total = V_fuel + V_clad + V_cool;
    double V_frac_fuel = V_fuel / V_total;
    double V_frac_clad = V_clad / V_total;
    double V_frac_cool = V_cool / V_total;

    double N_fuel_i_cm2pb = 0.0;
    double N_clad_i_cm2pb = 0.0;
    double N_cool_i_cm2pb = 0.0;
    double N_i_cm2pb = 0.0;

    for (iso_iter iso = J.begin(); iso != J.end(); iso++)
    {
        nuc_weight = isoname::nuc_weight(*iso);

        N_fuel_i_cm2pb = bright::cm2_per_barn * N_fuel_it[*iso][bt_s];

        N_clad_i_cm2pb = bright::cm2_per_barn * N_clad_it[*iso][bt_s];

        N_cool_i_cm2pb = bright::cm2_per_barn * N_cool_it[*iso][bt_s];

        N_i_cm2pb = bright::cm2_per_barn * ((N_fuel_it[*iso][bt_s] * V_frac_fuel) \
                                         +  (N_clad_it[*iso][bt_s] * V_frac_clad) \
                                         +  (N_cool_it[*iso][bt_s] * V_frac_cool));



        // Loop over all groups
        for (g = 0; g < G; g++)
        {
            Sigma_t_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_t_itg[*iso][bt_s][g];
            nubar_Sigma_f_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * nubar_sigma_f_itg[*iso][bt_s][g];
            chi_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * chi_itg[*iso][bt_s][g];
            Sigma_f_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_f_itg[*iso][bt_s][g];
            Sigma_gamma_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_gamma_itg[*iso][bt_s][g];
            Sigma_2n_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_2n_itg[*iso][bt_s][g];
            Sigma_3n_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_3n_itg[*iso][bt_s][g];
            Sigma_alpha_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_alpha_itg[*iso][bt_s][g];
            Sigma_proton_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_proton_itg[*iso][bt_s][g];
            Sigma_gamma_x_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_gamma_x_itg[*iso][bt_s][g];
            Sigma_2n_x_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_2n_x_itg[*iso][bt_s][g];

            Sigma_t_clad_tg[bt_s][g] += N_clad_i_cm2pb * sigma_t_itg[*iso][bt_s][g];

            Sigma_t_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_t_itg[*iso][bt_s][g];

            Sigma_t_tg[bt_s][g] += N_i_cm2pb * sigma_t_itg[*iso][bt_s][g];
            nubar_Sigma_f_tg[bt_s][g] += N_i_cm2pb * nubar_sigma_f_itg[*iso][bt_s][g];
            chi_tg[bt_s][g] += N_i_cm2pb * chi_itg[*iso][bt_s][g];
            Sigma_f_tg[bt_s][g] += N_i_cm2pb * sigma_f_itg[*iso][bt_s][g];
            Sigma_gamma_tg[bt_s][g] += N_i_cm2pb * sigma_gamma_itg[*iso][bt_s][g];
            Sigma_2n_tg[bt_s][g] += N_i_cm2pb * sigma_2n_itg[*iso][bt_s][g];
            Sigma_3n_tg[bt_s][g] += N_i_cm2pb * sigma_3n_itg[*iso][bt_s][g];
            Sigma_alpha_tg[bt_s][g] += N_i_cm2pb * sigma_alpha_itg[*iso][bt_s][g];
            Sigma_proton_tg[bt_s][g] += N_i_cm2pb * sigma_proton_itg[*iso][bt_s][g];
            Sigma_gamma_x_tg[bt_s][g] += N_i_cm2pb * sigma_gamma_x_itg[*iso][bt_s][g];
            Sigma_2n_x_tg[bt_s][g] += N_i_cm2pb * sigma_2n_x_itg[*iso][bt_s][g];


            for (int h =0; h < G; h++)
            {
                Sigma_s_fuel_tgh[bt_s][g][h] += N_fuel_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];
                Sigma_s_clad_tgh[bt_s][g][h] += N_clad_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];
                Sigma_s_cool_tgh[bt_s][g][h] += N_cool_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];
                Sigma_s_tgh[bt_s][g][h] += N_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];
            };
        };
    };

    // Re-Normalize chi
    double chi_fuel_tot = 0.0;
    for (g = 0; g < G; g++)
        chi_fuel_tot += chi_fuel_tg[bt_s][g]; 
    for (g = 0; g < G; g++)
        chi_fuel_tg[bt_s][g] = chi_fuel_tg[bt_s][g] / chi_fuel_tot; 

    double chi_tot = 0.0;
    for (g = 0; g < G; g++)
        chi_tot += chi_tg[bt_s][g]; 
    for (g = 0; g < G; g++)
        chi_tg[bt_s][g] = chi_tg[bt_s][g] / chi_tot; 

    // Calculate kappas as the inverse of the diffusion length
    // k = 1 / D = 3 * Sigma_tr = 3 (Sigma_t - mu * Sigma_s_gg)
    for (g = 0; g < G; g++)
    {
        kappa_fuel_tg[bt_s][g] = (3.0 * Sigma_t_fuel_tg[bt_s][g]) - (2.0 * Sigma_s_fuel_tgh[bt_s][g][g] / MW_fuel_t[bt_s]); 
        kappa_clad_tg[bt_s][g] = (3.0 * Sigma_t_clad_tg[bt_s][g]) - (2.0 * Sigma_s_clad_tgh[bt_s][g][g] / MW_clad_t[bt_s]); 
        kappa_cool_tg[bt_s][g] = (3.0 * Sigma_t_cool_tg[bt_s][g]) - (2.0 * Sigma_s_cool_tgh[bt_s][g][g] / MW_cool_t[bt_s]); 
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
        A_fuel_tgh[bt_s][g][g] += Sigma_t_fuel_tg[bt_s][g];
        A_clad_tgh[bt_s][g][g] += Sigma_t_clad_tg[bt_s][g];
        A_cool_tgh[bt_s][g][g] += Sigma_t_cool_tg[bt_s][g];
        A_tgh[bt_s][g][g] += Sigma_t_tg[bt_s][g];

        // Subtract the scattering kernel
        for (int h = 0; h < G; h++)
        {
            A_fuel_tgh[bt_s][g][h] -= Sigma_s_fuel_tgh[bt_s][g][h];
            A_clad_tgh[bt_s][g][h] -= Sigma_s_clad_tgh[bt_s][g][h];
            A_cool_tgh[bt_s][g][h] -= Sigma_s_cool_tgh[bt_s][g][h];
            A_tgh[bt_s][g][h] -= Sigma_s_tgh[bt_s][g][h];
        };
    };


    // Assemble the F matrix
    F_fuel_tgh[bt_s] = bright::vector_outer_product(nubar_Sigma_f_fuel_tg[bt_s], chi_fuel_tg[bt_s]);
    F_tgh[bt_s] = bright::vector_outer_product(nubar_Sigma_f_tg[bt_s], chi_tg[bt_s]);

    // Grab the inverse of the A matrix
    A_inv_fuel_tgh[bt_s] = bright::matrix_inverse(A_fuel_tgh[bt_s]);
    A_inv_clad_tgh[bt_s] = bright::matrix_inverse(A_clad_tgh[bt_s]);
    A_inv_cool_tgh[bt_s] = bright::matrix_inverse(A_cool_tgh[bt_s]);
    A_inv_tgh[bt_s] = bright::matrix_inverse(A_tgh[bt_s]);

    // Multiply the inverse of A by F
    A_inv_F_fuel_tgh[bt_s] = bright::matrix_multiplication(A_inv_fuel_tgh[bt_s], F_fuel_tgh[bt_s]);
    A_inv_F_tgh[bt_s] = bright::matrix_multiplication(A_inv_tgh[bt_s], F_tgh[bt_s]);

    //
    // Assemble the energy integral of transmutation matrix
    //
    int g, i, j, ind, jnd;
    std::vector< std::vector< std::vector<double> > > T_matrix = std::vector< std::vector< std::vector<double> > > (J_size, std::vector< std::vector<double> >(J_size, std::vector<double>(G, 0.0) ) );
    
    // Add the fission yields first
    for (ind = 0; ind < J_size; ind++)
    {
        for (jnd = 0; jnd < J_size; jnd++)
        {
            for (g = 0; g < G; g++)
                T_matrix[ind][jnd][g] = fission_product_yield_matrix[ind][jnd][g] * sigma_f_itg[J_order[ind]][bt_s][g];
        };
    };
    // Add the other cross sections
    int j_gamma, j_2n, j_3n, j_alpha, j_proton, j_gamma_x, j_2n_x; 
    for (ind = 0; ind < J_size; ind++)
    {
        i = J_order[ind];

        // Get from isos
        j_gamma = ((i/10) + 1) * 10;
        j_2n = j_gamma - 20;
        j_3n = j_2n - 10;
        j_alpha = j_3n - 20010;
        j_proton = j_gamma - 10010;
        j_gamma_x = j_gamma + 1;
        j_2n_x = j_2n + 1;

        // Add the capture cross-section
        if (0 < J.count(j_gamma))
        {
            jnd = J_index[j_gamma];
            for (g = 0; g < G; g++)
                T_matrix[ind][jnd][g] += sigma_gamma_itg[i][bt_s][g];
        };

        // Add the (n, 2n) cross-section
        if (0 < J.count(j_2n))
        {
            jnd = J_index[j_2n];
            for (g = 0; g < G; g++)
                T_matrix[ind][jnd][g] += sigma_2n_itg[i][bt_s][g];
        };

        // Add the (n, 3n) cross-section
        if (0 < J.count(j_3n))
        {
            jnd = J_index[j_3n];
            for (g = 0; g < G; g++)
                T_matrix[ind][jnd][g] += sigma_3n_itg[i][bt_s][g];
        };

        // Add the (n, alpha) cross-section
        if (0 < J.count(j_2n))
        {
            jnd = J_index[j_alpha];
            for (g = 0; g < G; g++)
                T_matrix[ind][jnd][g] += sigma_alpha_itg[i][bt_s][g];
        };

        // Add the (n, proton) cross-section
        if (0 < J.count(j_proton))
        {
            jnd = J_index[j_2n];
            for (g = 0; g < G; g++)
                T_matrix[ind][jnd][g] += sigma_proton_itg[i][bt_s][g];
        };

        // Add the capture (excited) cross-section
        if (0 < J.count(j_gamma_x))
        {
            jnd = J_index[j_gamma_x];
            for (g = 0; g < G; g++)
                T_matrix[ind][jnd][g] += sigma_gamma_x_itg[i][bt_s][g];
        };

        // Add the (n, 2n *) cross-section
        if (0 < J.count(j_2n_x))
        {
            jnd = J_index[j_2n_x];
            for (g = 0; g < G; g++)
                T_matrix[ind][jnd][g] += sigma_2n_x_itg[i][bt_s][g];
        };
    };

    // Multiply by flux and integrate over energy
    for (ind = 0; ind < J_size; ind++)
    {
        for (jnd = 0; jnd < J_size; jnd++)
        {
            for (g = 0; g < G; g++)
                T_int_tij[bt_s][ind][jnd] += T_matrix[ind][jnd][g] * phi_tg[bt_s][g];
        };
    };

    // Make the transmutation matrix for this time step
    M_tij[bt_s] = bright::matrix_addition(T_int_tij[bt_s], decay_matrix);
};









void ReactorMG::calc_criticality()
{
    // Init values
    int n = 0;
    int N = 5;

    float epsik = 1.0;
    float tmp_epsiphi; 
    float epsiphi = 1.0;
    float epsilon = 0.005;

    double k0 = 1.0;
    std::vector<double> phi0 (G, 1.0);
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

/*
        for (g = 0; g < G; g++)
            std::cout << phi0[g] << "   ";
        std::cout << "\n";
        for (g = 0; g < G; g++)
            std::cout << phi1[g] << "   ";
        std::cout << "\n----------\n\n\n";
*/

        // Calculate the next eigen-k
        nu_Sigma_f_phi0 = 0.0;
        nu_Sigma_f_phi1 = 0.0;
        for (g = 0; g < G; g++)
        {
            nu_Sigma_f_phi0 += nubar_Sigma_f_tg[bt_s][g] * phi0[g]; 
            nu_Sigma_f_phi1 += nubar_Sigma_f_tg[bt_s][g] * phi1[g]; 
        };
        k1 = k0 * nu_Sigma_f_phi1 / nu_Sigma_f_phi0;

//        std::cout << k0 << "   " << k1 << "\n";

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
        n++;
    };

//    std::cout << "N = " << N << "\n";

    // Set the final values to the class members
    k_t[bt_s] = k1;
    phi_tg[bt_s] = phi1;

    phi_t[bt_s] = 0.0;
    for (g = 0; g < G; g++)
        phi_t[bt_s] += phi1[0];

    Phi_t[bt_s] = phi_t[bt_s] * (burn_times[bt_s] - burn_times[bt_s]) * bright::sec_per_day;
};










void ReactorMG::calc_transmutation()
{
    // Calculates a tranmutation step via the Pade method

    // Get the transmutation matrix for this time delta
    double dt = burn_times[bt_s] - burn_times[bt_s - 1];
    std::vector< std::vector<double> > Mt = bright::scalar_matrix_product(dt, M_tij[bt_s]);

    // Set Pade coefficients, for p = q = 6
    double N_coef_n [7] = {1.00000000e+00,   5.00000000e-01,   1.13636364e-01, \
                           1.51515152e-02,   1.26262626e-03,   6.31313131e-05, \
                           1.50312650e-06};
    double D_coef_n [7] = {1.00000000e+00,   5.00000000e-01,   1.13636364e-01, \
                           1.51515152e-02,   1.26262626e-03,   6.31313131e-05, \
                           1.50312650e-06};

    // Init the new matrices
    std::vector< std::vector< std::vector<double> > > Mt_n = std::vector< std::vector< std::vector<double> > > (7, std::vector< std::vector<double> >(J_size, std::vector<double>(J_size, 0.0) ) );
    std::vector< std::vector< std::vector<double> > > neg_Mt_n = std::vector< std::vector< std::vector<double> > > (7, std::vector< std::vector<double> >(J_size, std::vector<double>(J_size, 0.0) ) );
    std::vector< std::vector<double> > N_pq = std::vector< std::vector<double> >(J_size, std::vector<double>(J_size, 0.0) );
    std::vector< std::vector<double> > D_pq = std::vector< std::vector<double> >(J_size, std::vector<double>(J_size, 0.0) );

    int n, i, j, ind, jnd;
    for (ind = 0; ind < J_size; ind++)
    {
        Mt_n[0][ind][ind] = 1.0;
        neg_Mt_n[0][ind][ind] = 1.0;
    };

    // Calculate the exponential matrices
    for (n = 1; n < 7; n++)
    {
        Mt_n[n] = bright::matrix_multiplication(Mt_n[n-1], Mt);
        neg_Mt_n[n] = bright::matrix_multiplication(Mt_n[n-1], bright::scalar_matrix_product(-1.0, Mt));
    };

    // Calculate the pade numerator and denom matrices
    for (n = 0; n < 7; n++)
    {
        N_pq = bright::matrix_addition(N_pq, bright::scalar_matrix_product(N_coef_n[n], Mt_n[n]));
        D_pq = bright::matrix_addition(D_pq, bright::scalar_matrix_product(D_coef_n[n], neg_Mt_n[n]));
    };

    // Invert the denominator
    std::vector< std::vector<double> > inv_D_pq = bright::matrix_inverse(D_pq);

    // Approximate the exponential e^(Mt) = D^-1 * N
    std::vector< std::vector<double> > exp_Mt = bright::matrix_multiplication(inv_D_pq, N_pq);

    // Make mass vectors
    std::vector<double> comp_prev (J_size, 0.0);
    for (ind = 0; ind < J_size; ind++)
    {
        i = J_order[ind];
        comp_prev[ind] = T_it[i][bt_s-1];
    };

    // Get the composition for the next time step
    std::vector<double> comp_next = bright::scalar_matrix_vector_product(1.0, exp_Mt, comp_prev);

    // Copy this composition back to the tranmutuation matrix
    for (ind = 0; ind < J_size; ind++)
    {
        i = J_order[ind];
        T_it[i][bt_s] = comp_next[ind];
    };

    // Calculate the burnup 
    CompDict cd_prev, cd_next;
    for (ind = 0; ind < J_size; ind++)
    {
        i = J_order[ind];
        cd_prev[i] = comp_prev[ind];
        cd_next[i] = comp_next[ind];
    };

    MassStream ms_prev (cd_prev);
    MassStream act_prev = ms_prev.get_act();

    MassStream ms_next (cd_next);
    MassStream act_next = ms_next.get_act();

    double delta_BU = (act_prev.mass - act_next.mass) * 931.46;

    BU_t[bt_s] = delta_BU + BU_t[bt_s - 1];
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
    sigma_f_itg.clear();
    sigma_gamma_itg.clear();
    sigma_2n_itg.clear();
    sigma_3n_itg.clear();
    sigma_alpha_itg.clear();
    sigma_proton_itg.clear();
    sigma_gamma_x_itg.clear();
    sigma_2n_x_itg.clear();

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

    // Also init attributes caluclated from burnup
    phi_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0) );
    phi_t = std::vector<double>(S, -1.0);
    Phi_t = std::vector<double>(S, -1.0);
    BU_t = std::vector<double>(S, -1.0);

    zeta_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 1.0) );
    lattice_E_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0) );
    lattice_F_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0) );


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
        sigma_f_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        sigma_gamma_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        sigma_2n_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        sigma_3n_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        sigma_alpha_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        sigma_proton_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        sigma_gamma_x_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
        sigma_2n_x_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));

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
    Sigma_t_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    nubar_Sigma_f_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    chi_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_s_fuel_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    Sigma_f_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_gamma_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_2n_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_3n_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_alpha_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_proton_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_gamma_x_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_2n_x_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));

    Sigma_t_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_s_clad_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));

    Sigma_t_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_s_cool_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));

    Sigma_t_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    nubar_Sigma_f_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    chi_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_s_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    Sigma_f_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_gamma_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_2n_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_3n_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_alpha_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_proton_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_gamma_x_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
    Sigma_2n_x_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));


    // Init the criticality matrices
    A_fuel_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    F_fuel_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    A_inv_fuel_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    A_inv_F_fuel_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));

    A_clad_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    A_inv_clad_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));

    A_cool_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    A_inv_cool_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));

    A_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    F_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    A_inv_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
    A_inv_F_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));

    T_int_tij = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(J_size, std::vector<double>(J_size, 0.0)));
    M_tij = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(J_size, std::vector<double>(J_size, 0.0)));

    // Init the multilplcation factpr
    k_t = std::vector<double>(S, -1.0);

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

        // Preform the criticality and burnup calulations
        assemble_multigroup_matrices();
        calc_criticality();

        if (s == 0)
            BU_t[0] = 0.0;
        else        
            calc_transmutation();
    };

};









void ReactorMG::calc_T_itd()
{
    /** Calculates the output isotopics of Mj(Fd).
     *  NOTE: Mj(Fd) is effectively the same variable as ms_prod before normalization!
     */

    CompDict tempOut;

    //Checks to see if the discharge index in at the end of the fluence table
    bool td_nLast = false;
    if ( (td_n+1) == S )
        td_nLast = true;			

    for (iso_iter j = J.begin(); j != J.end(); j++ )
    {
        if (td_nLast)
            tempOut[*j] = bright::SolveLine(Phid, Phi_t[td_n], T_it[*j][td_n], Phi_t[td_n-1], T_it[*j][td_n-1]);
        else
            tempOut[*j] = bright::SolveLine(Phid, Phi_t[td_n+1], T_it[*j][td_n+1], Phi_t[td_n], T_it[*j][td_n]);
    };

    ms_prod = MassStream(tempOut);	
};



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
    while ( (fp.f < S) && (BU_t[fp.f] < BU) )
    {
        fp.f = fp.f + 1;
    };
    fp.f = fp.f - 1;

    if (fp.f < 0)
        fp.f = 0;
    
    if ( (fp.f + 1) == S )
        fp.m = bright::slope(Phi_t[fp.f], BU_t[fp.f], Phi_t[fp.f-1], BU_t[fp.f-1]);
    else
        fp.m = bright::slope(Phi_t[fp.f+1], BU_t[fp.f+1], Phi_t[fp.f], BU_t[fp.f]);

    fp.F = ((BU - BU_t[fp.f])/fp.m) + Phi_t[fp.f];

    return fp;
};






double ReactorMG::batch_average_k(double BUd)
{
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

    // Get the multiplication factot and flux for each batch.
    std::vector<double> k_b (B, 0.0);
    std::vector<double> phi_b (B, 0.0);
    for (int b = 0; b < B; b++)
    {
        if ( (fps[b].f + 1) == S )
        {
            k_b[b] = bright::SolveLine(fps[b].F, Phi_t[fps[b].f], k_t[fps[b].f], Phi_t[fps[b].f-1], k_t[fps[b].f-1]);
            phi_b[b] = bright::SolveLine(fps[b].F, Phi_t[fps[b].f], phi_t[fps[b].f], Phi_t[fps[b].f-1], phi_t[fps[b].f-1]);
        }
        else
        {
            k_b[b] = bright::SolveLine(fps[b].F, Phi_t[fps[b].f+1], k_t[fps[b].f+1], Phi_t[fps[b].f], k_t[fps[b].f]);
            phi_b[b] = bright::SolveLine(fps[b].F, Phi_t[fps[b].f+1], phi_t[fps[b].f+1], Phi_t[fps[b].f], phi_t[fps[b].f]);
        };
    };

    // Calculate the flux weighted avearge
    double numerator = 0.0;
    double denominator = 0.0;
    for (int b = 0; b < B; b++)
    {
        numerator   += (k_b[b] * phi_b[b]);
        denominator += phi_b[b];
    };

    double batch_ave_k = numerator/denominator;

    return batch_ave_k;
};




void ReactorMG::BUd_bisection_method()
{
    //Calculates the maximum discharge burnup via the Bisection Method.
    int tempk = 1;
    double BUd_a, k_a, sign_a;
    double BUd_b, k_b, sign_b;
    double BUd_c, k_c, sign_c;
    
    //First Find a BUd that serves as an initial first guess
    while ( (tempk < S) && (1.0 < k_t[tempk]) )
        tempk = tempk + 1;
    tempk = tempk - 1;
    if (tempk < 0)
        tempk = 0;

    BUd_a = BU_t[tempk] * 2.0 * ((double) B) / ((double) B + 1.0);
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
    td_n = fp.f;    // lower index of discharge time
    Phid = fp.F;    // Discharge fluence

    // time at discharge
    if ( (fp.f + 1) == S )
        td = bright::SolveLine(fp.F, Phi_t[fp.f], burn_times[fp.f], Phi_t[fp.f-1], burn_times[fp.f-1]);
    else
        td = bright::SolveLine(fp.F, Phi_t[fp.f+1], burn_times[fp.f+1], Phi_t[fp.f], burn_times[fp.f]);

    return;
};




void ReactorMG::run_P_NL(double temp_pnl)
{
    /** Does a reactor run for a specific P_NL.
     *  Requires that ms_feed be (meaningfully) set.
     *  For use with calibrate_P_NL_to_BUd
     */

    P_NL = temp_pnl;
    burnup_core();
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



MassStream ReactorMG::calc()
{
    // Finds BUd and output isotopics.
    burnup_core();

    BUd_bisection_method();

    calc_ms_prod();

    return ms_prod;
};


MassStream ReactorMG::calc (CompDict incomp)
{
    // Finds BUd and output isotopics.
    ms_feed = MassStream (incomp);

    return calc();
};


MassStream ReactorMG::calc (MassStream instream)
{
    // Finds BUd and output isotopics.
    ms_feed = instream;

    return calc();
};










// 
// Lattice and Zeta functions below
//



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
        a = r;
        b = l / 2.0;
    
        lattice_E_planar(a, b);
        lattice_F_planar(a, b);
    }
    else if (lattice_flag == "Spherical")
    {
        a = r;
        b = l / 2.0;
    
        lattice_E_spherical(a, b);
        lattice_F_spherical(a, b);
    }
    else if (lattice_flag == "Cylindrical")
    {
        a = r;
        b = l / sqrt(bright::pi); //radius of cylinder with an equivilent cell volume

        lattice_E_cylindrical(a, b);
        lattice_F_cylindrical(a, b);
    }
    else
    {
        if (0 < FCComps::verbosity)
            std::cout << "Did not specify use of planar or spheical or cylindrical lattice functions! Assuming cylindrical...\n";
        
        a = r;
        b = l / sqrt(bright::pi); //radius of cylinder with an equivilent cell volume

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

