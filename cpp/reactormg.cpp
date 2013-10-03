// Multi-Group Reactor Component Class
#include "reactormg.h"



  /***********************************************/
  /*** ReactorMG Component Class and Functions ***/
  /***********************************************/

bright::ReactorMG::ReactorMG(std::string n) : FCComp(n)
{
};


bright::ReactorMG::ReactorMG(std::set<std::string> paramtrack, std::string n) : bright::FCComp(paramtrack, n)
{
};


bright::ReactorMG::ReactorMG(ReactorParameters rp, std::string n) : bright::FCComp(n)
{
  initialize(rp);
};


bright::ReactorMG::ReactorMG(ReactorParameters rp, std::set<std::string> paramtrack, std::string n) : bright::FCComp(paramtrack, n)
{
  initialize(rp);
};


bright::ReactorMG::~ReactorMG()
{
};




void bright::ReactorMG::initialize(ReactorParameters rp)
{
  /** Sets reactor specific parameters.
   *  Must be done once at the beginning of reactor object life.
   */

  B = rp.batches;				// Total number of fuel loading batches
  flux = rp.flux;				// Flux used for Fluence
  chemical_form_fuel = rp.fuel_form;		// Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
  chemical_form_clad = rp.cladding_form;  // Chemical form of Fuel as Dictionary.  Keys are elements or isotopes while values represent mass weights.  Denote heavy metal by key "IHM".
  chemical_form_cool = rp.coolant_form;	// Same a fuel chemical form but for coolant.  Should not have "IHM"

  rho_fuel = rp.fuel_density;     // Fuel Density
  rho_clad = rp.cladding_density; // Cladding Density
  rho_cool = rp.coolant_density;  // Coolant Density

  P_NL = rp.pnl;				// Non-Leakage Probability
  target_BU = rp.BUt;			// Target Discharge Burnup, only used for graphing inside of this component
  burn_regions = rp.burn_regions;
  specific_power = rp.specific_power;
  burn_time = 0.0;
  burn_times = rp.burn_times;
  S = burn_times.size();

  // Flags
  use_zeta = rp.use_disadvantage_factor;		//Boolean value on whether or not the disadvantage factor should be used
  lattice_flag = rp.lattice_type;		//lattice_flagType (Planar || Spherical || Cylindrical)
  rescale_hydrogen_xs = rp.rescale_hydrogen;	//Rescale the Hydrogen-1 XS?
  burnup_via_constant = rp.burnup_via_constant;  // power or flux
  branch_ratio_cutoff = rp.branch_ratio_cutoff; // Cut-off for bateman chains

  // Calculates Volumes
  r_fuel = rp.fuel_radius;    // Fuel region radius
  r_void = rp.void_radius;    // Void region radius
  r_clad = rp.clad_radius;    // Clad region radius
  pitch = rp.unit_cell_pitch; // Unit cell side length

  S_O = rp.open_slots;		// Number of open slots in fuel assembly
  S_T = rp.total_slots;		// Total number of Fuel assembly slots.

  // Fuel Volume Fraction
  V_fuel = ((pyne::pi * r_fuel * r_fuel)/(pitch*pitch)) * (1.0 - S_O/S_T); 

  // Cladding Volume Fraction
  V_clad = ((pyne::pi * (r_clad * r_clad - r_void * r_void))/(pitch*pitch)) * (1.0 - S_O/S_T); 

  // Coolant Volume Fraction
  V_cool = ((pitch*pitch - pyne::pi * r_clad * r_clad)/(pitch*pitch)) * (1.0 - S_O/S_T) + (S_O/S_T);
};



void bright::ReactorMG::loadlib(std::string lib)
{
  // Loads Apporiate Libraries for ReactorMG

  // Check that the file is there
  if (!pyne::file_exists(lib))
    throw pyne::FileNotFound(lib);

  //Check to see if the file is in HDF5 format.
  bool isH5 = H5::H5File::isHdf5(lib);
  if (!isH5)
  {
    std::cout << "!!!Warning!!! " << lib << " is not a valid HDF5 file!\n";
    return;
  };

  // Turn off the exceptions
  H5::Exception::dontPrint();

  // Open file
  H5::H5File rmglib(lib, H5F_ACC_RDONLY);
  hid_t rmglibid = rmglib.getId();

  // Load isos
  std::string load_zz, transmute_zz;
  if (h5wrap::path_exists(rmglibid, "/load_isos_zz"))
  {
    load_zz = "/load_isos_zz";
    transmute_zz = "/transmute_isos_zz";
  }
  else if (h5wrap::path_exists(rmglibid, "/load_nucs_zz"))
  {
    load_zz = "/load_nucs_zz";
    transmute_zz = "/transmute_nucs_zz";
  }
  else
    throw h5wrap::PathNotFound(lib, "/load_nucs_zz or /load_nucs_zz");

  I = h5wrap::h5_array_to_cpp_set<int>(rmglibid, load_zz,      H5T_NATIVE_INT);
  J = h5wrap::h5_array_to_cpp_set<int>(rmglibid, transmute_zz, H5T_NATIVE_INT);
  K = h5wrap::h5_array_to_cpp_set<int>(rmglibid, transmute_zz, H5T_NATIVE_INT);

  // Load perturbation table
  perturbations = h5wrap::HomogenousTypeTable<double>(rmglibid, "/perturbations");
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
  pert_data_g full_E_g = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/energy");
  E_g = full_E_g[0];
  G = E_g.size() - 1;

  // Load fluxes and fluence
  phi_g = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/phi_g");
  phi = h5wrap::h5_array_to_cpp_vector_1d<double>(rmglibid, "/phi");
  Phi = h5wrap::h5_array_to_cpp_vector_1d<double>(rmglibid, "/Phi");

  // Load time and burnup
  time0 = h5wrap::h5_array_to_cpp_vector_1d<double>(rmglibid, "/time0");
  BU0 = h5wrap::h5_array_to_cpp_vector_1d<double>(rmglibid, "/BU0");

  // Clear transmutation vectors and cross sections before reading in
  Ti0.clear();
  sigma_t_pg.clear();
  sigma_a_pg.clear();
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
  for(nuc_iter nuciter = J.begin(); nuciter != J.end(); nuciter++)
  {
    iso_zz = *nuciter;
    iso_LL = pyne::nucname::name(iso_zz);

    // Add transmutation vector
    Ti0[iso_zz] = h5wrap::h5_array_to_cpp_vector_1d<double>(rmglibid, "/Ti0/" + iso_LL);

    // Add cross sections
    sigma_t_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/sigma_t/" + iso_LL);
    sigma_a_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/sigma_a/" + iso_LL);
    nubar_sigma_f_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/nubar_sigma_f/" + iso_LL);
    chi_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/chi/" + iso_LL);
    sigma_s_pgh[iso_zz] = h5wrap::h5_array_to_cpp_vector_3d<double>(rmglibid, "/sigma_s_gh/" + iso_LL);
    sigma_f_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/sigma_f/" + iso_LL);
    sigma_gamma_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/sigma_gamma/" + iso_LL);
    sigma_2n_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/sigma_2n/" + iso_LL);
    sigma_3n_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/sigma_3n/" + iso_LL);
    sigma_alpha_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/sigma_alpha/" + iso_LL);
    sigma_proton_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/sigma_proton/" + iso_LL);
    sigma_gamma_x_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/sigma_gamma_x/" + iso_LL);
    sigma_2n_x_pg[iso_zz] = h5wrap::h5_array_to_cpp_vector_2d<double>(rmglibid, "/sigma_2n_x/" + iso_LL);
  };

  // close the reactor library
  rmglib.close();


  //
  // Create a decay matrix from a file based off of the J isotopes
  //
  std::string nuc_data_file; 
  std::string nuc_data_file_bright = bright::BRIGHT_DATA + "/nuc_data.h5";

  //Check to see if the file is in HDF5 format.
  if (pyne::file_exists(pyne::NUC_DATA_PATH))
    nuc_data_file = pyne::NUC_DATA_PATH;
  else if (pyne::file_exists(bright::BRIGHT_DATA + "/nuc_data.h5"))
    nuc_data_file = bright::BRIGHT_DATA + "/nuc_data.h5";
  else
    throw pyne::FileNotFound("nuc_data.h5");

  isH5 = H5::H5File::isHdf5(nuc_data_file);
  if (!isH5)
    throw h5wrap::FileNotHDF5(nuc_data_file);

  // Open the HDF5 file
  H5::H5File nuc_data_h5 (nuc_data_file.c_str(), H5F_ACC_RDONLY );


  //
  // Read in the decay data table as an array of pyne::atomic_decay_struct
  //
  int i, j, k, ind, jnd, knd, l, g;
  // Get the HDF5 compound type (table) description
  H5::CompType atom_dec_desc(sizeof(pyne::atomic_decay_struct));
  atom_dec_desc.insertMember("from_nuc",   HOFFSET(pyne::atomic_decay_struct, from_nuc),   H5::PredType::NATIVE_INT);
  atom_dec_desc.insertMember("level", HOFFSET(pyne::atomic_decay_struct, level), H5::PredType::NATIVE_DOUBLE);
  atom_dec_desc.insertMember("to_nuc",   HOFFSET(pyne::atomic_decay_struct, to_nuc),   H5::PredType::NATIVE_INT);
  atom_dec_desc.insertMember("half_life", HOFFSET(pyne::atomic_decay_struct, half_life), H5::PredType::NATIVE_DOUBLE);
  atom_dec_desc.insertMember("decay_const", HOFFSET(pyne::atomic_decay_struct, decay_const), H5::PredType::NATIVE_DOUBLE);
  atom_dec_desc.insertMember("branch_ratio", HOFFSET(pyne::atomic_decay_struct, branch_ratio), H5::PredType::NATIVE_DOUBLE);

  H5::DataSet decay_data_set = nuc_data_h5.openDataSet("/atomic_decay");
  H5::DataSpace decay_data_space = decay_data_set.getSpace();
  int decay_data_length = decay_data_space.getSimpleExtentNpoints(); 

  pyne::atomic_decay_struct * decay_data_array = new pyne::atomic_decay_struct [decay_data_length];
  decay_data_set.read(decay_data_array, atom_dec_desc);

  // Finish initializing K, based on decay info    
  for (l = 0; l < decay_data_length; l++)
  {
    K.insert(decay_data_array[l].from_nuc);
    K.insert(decay_data_array[l].to_nuc);
  };

  K_num = K.size();
  K_ord = std::vector<int> (K.begin(), K.end());
  std::sort(K_ord.begin(), K_ord.end());
  for (k = 0; k < K_num; k++)
    K_ind[K_ord[k]] = k;

  // Make decay_martrix from this data.
  decay_matrix = bright::SparseMatrix<double>(2*decay_data_length, K_num, K_num);

  for (l = 0; l < decay_data_length; l++)
  {
    i = decay_data_array[l].from_nuc;
    j = decay_data_array[l].to_nuc;

    if (i == j)
      continue;

    // Get the indexes for these nulcides into the matrix
    ind = K_ind[i];
    jnd = K_ind[j];

    // Add diagonal elements
    decay_matrix.push_back(ind, ind, -decay_data_array[l].decay_const);

    // Add i,j element to matrix
    if (i != j)
      decay_matrix.push_back(ind, jnd, decay_data_array[l].branch_ratio * decay_data_array[l].decay_const);
  };

  decay_matrix.clean_up();


  //
  // Read in the fission table
  //
  H5::DataSet fission_set = nuc_data_h5.openDataSet("/neutron/cinder_xs/fission");
  H5::DataSpace fission_space = fission_set.getSpace();
  int fission_length = fission_space.getSimpleExtentNpoints(); 

  bright::fission_struct * fission_array = new bright::fission_struct [fission_length];
  fission_set.read(fission_array, bright::fission_desc);

  // Run through the array and make join maps
  //  key = fission index
  //  value = vector of J_indexs
  std::map<int, std::vector<int> > thermal_join;
  std::map<int, std::vector<int> > fast_join;

  int ty, fy;
  for (l = 0; l < fission_length; l++)
  {
    i = fission_array[l].nuc_zz;

    // skip non-element from-isos
    if (K.count(i) < 1)
      continue;

    // make thermal join
    ty = fission_array[l].thermal_yield;

    if (thermal_join.count(ty) < 1)
      thermal_join[ty] = std::vector<int>();

    thermal_join[ty].push_back(K_ind[i]);

    // make fast join
    fy = fission_array[l].fast_yield;

    if (fast_join.count(fy) < 1)
      fast_join[fy] = std::vector<int>();

    fast_join[fy].push_back(K_ind[i]);
  };


  // Read in fission product yeilds
  H5::DataSet fp_yields_set = nuc_data_h5.openDataSet("/neutron/cinder_fission_products/yields");
  H5::DataSpace fp_yields_space = fp_yields_set.getSpace();
  int fp_yields_length = fp_yields_space.getSimpleExtentNpoints(); 

  bright::fission_product_yields_struct * fp_yields_array = new bright::fission_product_yields_struct [fp_yields_length];
  fp_yields_set.read(fp_yields_array, bright::fission_product_yields_desc);


  // Run through the array and make yield matrices
  thermal_yield_matrix = bright::SparseMatrix<double>(fp_yields_length, K_num, K_num);
  fast_yield_matrix = bright::SparseMatrix<double>(fp_yields_length, K_num, K_num);

  int index, tj, fj, TJ, FJ;
  double mf;
  for (l = 0; l < fp_yields_length; l++)
  {
    // Get important data from struct
    index = fp_yields_array[l].index;
    j = fp_yields_array[l].to_nuc_zz;
    jnd = K_ind[j];
    mf = fp_yields_array[l].mass_frac;

    // Add to thermal yields
    if (0 < thermal_join.count(index))
    {
      TJ = thermal_join[index].size();
      for (tj = 0; tj < TJ; tj++)
      {
        ind = thermal_join[index][tj];
        thermal_yield_matrix.push_back(ind, jnd, mf);
      };
    };

    // Add to fast yields.
    if (0 < fast_join.count(index))
    {
      FJ = fast_join[index].size();
      for (fj = 0; fj < FJ; fj++)
      {
        ind = fast_join[index][fj];
        fast_yield_matrix.push_back(ind, jnd, mf);
      };
    };
  };

  thermal_yield_matrix.clean_up();
  fast_yield_matrix.clean_up();


  // Make fission product yield matrix
  fission_product_yield_matrix = std::vector< bright::SparseMatrix<double> > (G);

  // Set the mass fraction between thermal and fast data.
  // Do not interpolate here, you'll get negative masses...
  for (g = 0; g < G; g++)
  {
    // fission_product_yield_matrix[g].push_back(0, 0, 0.0);
/*
*/
    if (0.001 < E_g[g])
      fission_product_yield_matrix[g] = fast_yield_matrix;
    else
      fission_product_yield_matrix[g] = thermal_yield_matrix;
  };

  //
  // Read in the one group cross sections
  //
  // Thermal
  H5::DataSet xs_1g_thermal_set = nuc_data_h5.openDataSet("/neutron/simple_xs/thermal");
  H5::DataSpace xs_1g_thermal_space = xs_1g_thermal_set.getSpace();
  int xs_1g_thermal_length = xs_1g_thermal_space.getSimpleExtentNpoints(); 

  bright::xs_1g_struct * xs_1g_thermal_array = new bright::xs_1g_struct [xs_1g_thermal_length];
  xs_1g_thermal_set.read(xs_1g_thermal_array, bright::xs_1g_desc);

  // Fast
  H5::DataSet xs_1g_fast_set = nuc_data_h5.openDataSet("/neutron/simple_xs/fission_spectrum_ave");
  H5::DataSpace xs_1g_fast_space = xs_1g_fast_set.getSpace();
  int xs_1g_fast_length = xs_1g_fast_space.getSimpleExtentNpoints(); 

  bright::xs_1g_struct * xs_1g_fast_array = new bright::xs_1g_struct [xs_1g_fast_length];
  xs_1g_fast_set.read(xs_1g_fast_array, bright::xs_1g_desc);

  // Copy the data over
  double Eng_g;
  nuc_set xs_isos;
  std::vector<double> sig_t, sig_a, sig_f, nu_sig_f, sig_gamma, sig_2n, sig_3n, sig_alpha, sig_proton;
  std::vector< std::vector<double> > zeros_pg;
  std::vector< std::vector< std::vector<double> > > zeros_pgh;

  zeros_pg = std::vector< std::vector<double> >(nperturbations, std::vector<double> (G,  0.0));
  zeros_pgh =  std::vector< std::vector< std::vector<double> > >(nperturbations, std::vector< std::vector<double> > (G,  std::vector<double> (G, 0.0)));

  for (l = 0; l < xs_1g_fast_length; l++)
  {
    i = xs_1g_thermal_array[l].nuc_zz;

    if (J.count(i) == 1)
      continue;

    if (K.count(i) == 0)
      continue;

    xs_isos.insert(i);

    // Init the interpolation arrays
    sig_t = std::vector<double>(G);
    sig_a = std::vector<double>(G);
    sig_f = std::vector<double>(G);
    nu_sig_f = std::vector<double>(G);
    sig_gamma = std::vector<double>(G);
    sig_2n = std::vector<double>(G);
    sig_3n = std::vector<double>(G);
    sig_alpha = std::vector<double>(G);
    sig_proton = std::vector<double>(G);

    // Fill the ineterpolation array
    for (g = 0; g < G; g++)
    {
      Eng_g = E_g[g];

      sig_t[g] = pyne::solve_line(Eng_g, 1.0, xs_1g_fast_array[l].sigma_t, 2.53e-08, xs_1g_thermal_array[l].sigma_t);
      if (sig_t[g] < 0.0)
        sig_t[g] = 0.0;

      sig_a[g] = pyne::solve_line(Eng_g, 1.0, xs_1g_fast_array[l].sigma_a, 2.53e-08, xs_1g_thermal_array[l].sigma_a);
      if (sig_a[g] < 0.0)
        sig_a[g] = 0.0;

      sig_f[g] = pyne::solve_line(Eng_g, 1.0, xs_1g_fast_array[l].sigma_f, 2.53e-08, xs_1g_thermal_array[l].sigma_f);
      if (sig_f[g] < 0.0)
        sig_f[g] = 0.0;

      nu_sig_f[g] = 2.5 * sig_f[g];

      sig_gamma[g] = pyne::solve_line(Eng_g, 1.0, xs_1g_fast_array[l].sigma_gamma, 2.53e-08, xs_1g_thermal_array[l].sigma_gamma);
      if (sig_gamma[g] < 0.0)
        sig_gamma[g] = 0.0;

      sig_2n[g] = pyne::solve_line(Eng_g, 1.0, xs_1g_fast_array[l].sigma_2n, 2.53e-08, xs_1g_thermal_array[l].sigma_2n);
      if (sig_2n[g] < 0.0)
        sig_2n[g] = 0.0;

      sig_3n[g] = pyne::solve_line(Eng_g, 1.0, xs_1g_fast_array[l].sigma_3n, 2.53e-08, xs_1g_thermal_array[l].sigma_3n);
      if (sig_3n[g] < 0.0)
        sig_3n[g] = 0.0;

      sig_alpha[g] = pyne::solve_line(Eng_g, 1.0, xs_1g_fast_array[l].sigma_alpha, 2.53e-08, xs_1g_thermal_array[l].sigma_alpha);
      if (sig_alpha[g] < 0.0)
        sig_alpha[g] = 0.0;

      sig_proton[g] = pyne::solve_line(Eng_g, 1.0, xs_1g_fast_array[l].sigma_proton, 2.53e-08, xs_1g_thermal_array[l].sigma_proton);
      if (sig_proton[g] < 0.0)
        sig_proton[g] = 0.0;
    };

    // Copy back the data to the XS library
    sigma_t_pg[i] = std::vector< std::vector<double> >(nperturbations, sig_t);
    sigma_a_pg[i] = std::vector< std::vector<double> >(nperturbations, sig_a);
    sigma_f_pg[i] = std::vector< std::vector<double> >(nperturbations, sig_f);
    nubar_sigma_f_pg[i] = std::vector< std::vector<double> >(nperturbations, nu_sig_f);
    sigma_gamma_pg[i] = std::vector< std::vector<double> >(nperturbations, sig_gamma);
    sigma_2n_pg[i] = std::vector< std::vector<double> >(nperturbations, sig_2n);
    sigma_3n_pg[i] = std::vector< std::vector<double> >(nperturbations, sig_3n);
    sigma_alpha_pg[i] = std::vector< std::vector<double> >(nperturbations, sig_alpha);
    sigma_proton_pg[i] = std::vector< std::vector<double> >(nperturbations, sig_proton);

    // Fill in zeros is places where data is not avilable
    chi_pg[i] = zeros_pg;
    sigma_s_pgh[i] = zeros_pgh;
    sigma_gamma_x_pg[i] = zeros_pg;
    sigma_2n_x_pg[i] = zeros_pg;
  };

  // Zero out XS for isos present in decay but not XS data
  for (nuc_iter iso = K.begin(); iso != K.end(); iso++)
  {
    i = *iso;
    if ((xs_isos.count(i) == 1) || (J.count(i) == 1))
      continue;

    sigma_t_pg[i] = zeros_pg;
    sigma_a_pg[i] = zeros_pg;
    nubar_sigma_f_pg[i] = zeros_pg;
    chi_pg[i] = zeros_pg;
    sigma_f_pg[i] = zeros_pg;
    sigma_s_pgh[i] = zeros_pgh;
    sigma_gamma_pg[i] = zeros_pg;
    sigma_2n_pg[i] = zeros_pg;
    sigma_3n_pg[i] = zeros_pg;
    sigma_alpha_pg[i] = zeros_pg;
    sigma_proton_pg[i] = zeros_pg;
    sigma_gamma_x_pg[i] = zeros_pg;
    sigma_2n_x_pg[i] = zeros_pg;
  };

  // close the nuc_data library
  nuc_data_h5.close();

  return;
};





void bright::ReactorMG::calc_nearest_neighbors()
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
      iso_zz = pyne::nucname::id(iso_LL);

      // Determine the mass of the isotope in the feed
      if (0 < mat_feed.comp.count(iso_zz))
        iso_mass = mat_feed.comp[iso_zz];
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












void bright::ReactorMG::interpolate_cross_sections()
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
      iso_zz = pyne::nucname::id(iso_LL);

      // Determine the mass of the isotope in the feed
      if (0 < mat_feed.comp.count(iso_zz))
        iso_mass = mat_feed.comp[iso_zz];
      else
        iso_mass = 0.0;

      // Calculate the x-factor if appropriate.
      if (nn0[iso_col] != nn1[iso_col])
        x_factor = x_factor + ((iso_mass - nn0[iso_col])/(nn1[iso_col] - nn0[iso_col]));
    };
  };

  if (nn0["burn_times"] != nn1["burn_times"])
    x_factor = x_factor + ((burn_time - nn0["burn_times"])/(nn1["burn_times"] - nn0["burn_times"]));

  // Let's flesh out this time step a bit
  // for the nuclides not in the data library
  for (nuc_iter iso = K.begin(); iso != K.end(); iso++)
  {
    if (J.count(*iso) == 1)
      continue;

    sigma_t_itg[*iso][bt_s] = sigma_t_pg[*iso][a0];
    sigma_a_itg[*iso][bt_s] = sigma_a_pg[*iso][a0];
    nubar_sigma_f_itg[*iso][bt_s] = nubar_sigma_f_pg[*iso][a0];
    chi_itg[*iso][bt_s] = chi_pg[*iso][a0];
    sigma_f_itg[*iso][bt_s] = sigma_f_pg[*iso][a0];
    sigma_gamma_itg[*iso][bt_s] = sigma_gamma_pg[*iso][a0];
    sigma_2n_itg[*iso][bt_s] = sigma_2n_pg[*iso][a0];
    sigma_3n_itg[*iso][bt_s] = sigma_3n_pg[*iso][a0];
    sigma_alpha_itg[*iso][bt_s] = sigma_alpha_pg[*iso][a0];
    sigma_proton_itg[*iso][bt_s] = sigma_proton_pg[*iso][a0];
    sigma_gamma_x_itg[*iso][bt_s] = sigma_gamma_x_pg[*iso][a0];
    sigma_2n_x_itg[*iso][bt_s] = sigma_2n_x_pg[*iso][a0];

    for (int g = 0; g < G; g++)
      sigma_s_itgh[*iso][bt_s][g] = sigma_s_pgh[*iso][a0][g];
  };

  // Now that we have found the x-factor, we get to do the actual interpolations. Oh Joy!
  for (nuc_iter iso = J.begin(); iso != J.end(); iso++)
  {
    // Interpolate the cross-sections
    sigma_t_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_t_pg[*iso][a1], sigma_t_pg[*iso][a0]);
    sigma_a_itg[*iso][bt_s] = bright::y_x_factor_interpolation(x_factor, sigma_a_pg[*iso][a1], sigma_a_pg[*iso][a0]);
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









void bright::ReactorMG::calc_mass_weights()
{
  /** 
   *  Calculates the appropriate mass fractions for this time step
   */

  // First things first, let's calculate the atomic weight of the HM
  double inverse_A_HM = 0.0;
  double mass_HM = 0.0;
  for (nuc_iter iso = K.begin(); iso != K.end(); iso++)
  {
    mass_HM += T_it[*iso][bt_s];
    inverse_A_HM += (T_it[*iso][bt_s] / pyne::atomic_mass(*iso));
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
      key_zz = pyne::nucname::id(key->first);
      MW_fuel_t[bt_s] += chemical_form_fuel[key->first] * pyne::atomic_mass(key_zz);
    };
  };

  // Cladding Molecular Weight
  for (std::map<std::string, double>::iterator key = chemical_form_clad.begin(); key != chemical_form_clad.end(); key++)
  {
    key_zz = pyne::nucname::id(key->first);
    MW_clad_t[bt_s] += chemical_form_clad[key->first] * pyne::atomic_mass(key_zz);
  };

  // Coolant Molecular Weight
  for (std::map<std::string, double>::iterator key = chemical_form_cool.begin(); key != chemical_form_cool.end(); key++)
  {
    key_zz = pyne::nucname::id(key->first);
    MW_cool_t[bt_s] += chemical_form_cool[key->first] * pyne::atomic_mass(key_zz);
  };



  //
  // Build the atom number density dictionaries
  //

  // now for the n_it in the Fuel
  for (std::map<std::string, double>::iterator key = chemical_form_fuel.begin(); key != chemical_form_fuel.end(); key++)
  {
    if ( (key->first) == "IHM")
    {
      for (nuc_iter iso = K.begin(); iso != K.end(); iso++)
        n_fuel_it[*iso][bt_s] += chemical_form_fuel[key->first] * T_it[*iso][bt_s];
    }
    else
    {
      key_zz = pyne::nucname::id(key->first);
      n_fuel_it[key_zz][bt_s] += chemical_form_fuel[key->first];
    }
  };

  // Note that the n_it in the cladding is just chemical_form_clad
  for (std::map<std::string, double>::iterator key = chemical_form_clad.begin(); key != chemical_form_clad.end(); key++)
  {
    key_zz = pyne::nucname::id(key->first);
    n_clad_it[key_zz][bt_s] = chemical_form_clad[key->first];
  };

  // Note that the n_it in the coolant is just chemical_form_cool
  for (std::map<std::string, double>::iterator key = chemical_form_cool.begin(); key != chemical_form_cool.end(); key++)
  {
    key_zz = pyne::nucname::id(key->first);
    n_cool_it[key_zz][bt_s] = chemical_form_cool[key->first];
  };



  //
  // Build the number density dictionaries
  //
  for (nuc_iter iso = K.begin(); iso != K.end(); iso++)
  {
    // Fuel Number Density
    N_fuel_it[*iso][bt_s] = n_fuel_it[*iso][bt_s] * rho_fuel * (pyne::N_A / MW_fuel_t[bt_s]);

    // Cladding Number Density
    N_clad_it[*iso][bt_s] = n_clad_it[*iso][bt_s] * rho_clad * (pyne::N_A / MW_clad_t[bt_s]);

    // Coolant Number Density
    N_cool_it[*iso][bt_s] = n_cool_it[*iso][bt_s] * rho_cool * (pyne::N_A / MW_cool_t[bt_s]);
  };



  //
  // Build the mass weight dictionaries
  //
  double iso_weight;
  double relative_volume_clad = (rho_clad * MW_fuel_t[bt_s] * V_clad) / (rho_fuel * MW_clad_t[bt_s] * V_fuel);
  double relative_volume_cool = (rho_cool * MW_fuel_t[bt_s] * V_cool) / (rho_fuel * MW_cool_t[bt_s] * V_fuel);

  for (nuc_iter iso = K.begin(); iso != K.end(); iso++)
  {
    iso_weight = pyne::atomic_mass(*iso);

    // Fuel mass weight
    m_fuel_it[*iso][bt_s] =  n_fuel_it[*iso][bt_s] * iso_weight / A_HM_t[bt_s];

    // Cladding mass weight 
    m_clad_it[*iso][bt_s] = (n_clad_it[*iso][bt_s] * iso_weight / A_HM_t[bt_s]) * relative_volume_clad;

    // Coolant mass weight 
    m_cool_it[*iso][bt_s] = (n_cool_it[*iso][bt_s] * iso_weight / A_HM_t[bt_s]) * relative_volume_cool;
  };
};







void bright::ReactorMG::fold_mass_weights()
{
  // Folds mass weight in with cross-sections for current time step
  int g, h;
  double mu;

  double N_fuel_i_cm2pb = 0.0;
  double N_clad_i_cm2pb = 0.0;
  double N_cool_i_cm2pb = 0.0;

  double Sig_s_fuel_ig, Sig_s_clad_ig, Sig_s_cool_ig;
  double AW_ig;

  for (nuc_iter iso = K.begin(); iso != K.end(); iso++)
  {
    N_fuel_i_cm2pb = pyne::cm2_per_barn * N_fuel_it[*iso][bt_s];
    N_clad_i_cm2pb = pyne::cm2_per_barn * N_clad_it[*iso][bt_s];
    N_cool_i_cm2pb = pyne::cm2_per_barn * N_cool_it[*iso][bt_s];

    AW_ig = pyne::atomic_mass(*iso);

    // Loop over all groups
    for (g = 0; g < G; g++)
    {
      Sigma_t_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_t_itg[*iso][bt_s][g];
      Sigma_a_fuel_tg[bt_s][g] += N_fuel_i_cm2pb * sigma_a_itg[*iso][bt_s][g];
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
      Sigma_a_clad_tg[bt_s][g] += N_clad_i_cm2pb * sigma_a_itg[*iso][bt_s][g];
      nubar_Sigma_f_clad_tg[bt_s][g] += N_clad_i_cm2pb * nubar_sigma_f_itg[*iso][bt_s][g];
      chi_clad_tg[bt_s][g] += N_clad_i_cm2pb * chi_itg[*iso][bt_s][g];
      Sigma_f_clad_tg[bt_s][g] += N_clad_i_cm2pb * sigma_f_itg[*iso][bt_s][g];
      Sigma_gamma_clad_tg[bt_s][g] += N_clad_i_cm2pb * sigma_gamma_itg[*iso][bt_s][g];
      Sigma_2n_clad_tg[bt_s][g] += N_clad_i_cm2pb * sigma_2n_itg[*iso][bt_s][g];
      Sigma_3n_clad_tg[bt_s][g] += N_clad_i_cm2pb * sigma_3n_itg[*iso][bt_s][g];
      Sigma_alpha_clad_tg[bt_s][g] += N_clad_i_cm2pb * sigma_alpha_itg[*iso][bt_s][g];
      Sigma_proton_clad_tg[bt_s][g] += N_clad_i_cm2pb * sigma_proton_itg[*iso][bt_s][g];
      Sigma_gamma_x_clad_tg[bt_s][g] += N_clad_i_cm2pb * sigma_gamma_x_itg[*iso][bt_s][g];
      Sigma_2n_x_clad_tg[bt_s][g] += N_clad_i_cm2pb * sigma_2n_x_itg[*iso][bt_s][g];

      Sigma_t_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_t_itg[*iso][bt_s][g];
      Sigma_a_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_a_itg[*iso][bt_s][g];
      nubar_Sigma_f_cool_tg[bt_s][g] += N_cool_i_cm2pb * nubar_sigma_f_itg[*iso][bt_s][g];
      chi_cool_tg[bt_s][g] += N_cool_i_cm2pb * chi_itg[*iso][bt_s][g];
      Sigma_f_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_f_itg[*iso][bt_s][g];
      Sigma_gamma_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_gamma_itg[*iso][bt_s][g];
      Sigma_2n_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_2n_itg[*iso][bt_s][g];
      Sigma_3n_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_3n_itg[*iso][bt_s][g];
      Sigma_alpha_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_alpha_itg[*iso][bt_s][g];
      Sigma_proton_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_proton_itg[*iso][bt_s][g];
      Sigma_gamma_x_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_gamma_x_itg[*iso][bt_s][g];
      Sigma_2n_x_cool_tg[bt_s][g] += N_cool_i_cm2pb * sigma_2n_x_itg[*iso][bt_s][g];


      Sig_s_fuel_ig = 0.0;
      Sig_s_clad_ig = 0.0;
      Sig_s_cool_ig = 0.0;

      for (h =0; h < G; h++)
      {
        Sigma_s_fuel_tgh[bt_s][g][h] += N_fuel_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];
        Sigma_s_clad_tgh[bt_s][g][h] += N_clad_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];
        Sigma_s_cool_tgh[bt_s][g][h] += N_cool_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];

        Sig_s_fuel_ig += N_fuel_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];
        Sig_s_clad_ig += N_clad_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];
        Sig_s_cool_ig += N_cool_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][h];

/*
        Sig_s_fuel_ig += N_fuel_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][g];
        Sig_s_clad_ig += N_clad_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][g];
        Sig_s_cool_ig += N_cool_i_cm2pb * sigma_s_itgh[*iso][bt_s][g][g];

        Sig_s_fuel_ig += N_fuel_i_cm2pb * sigma_s_itgh[*iso][bt_s][h][h];
        Sig_s_clad_ig += N_clad_i_cm2pb * sigma_s_itgh[*iso][bt_s][h][h];
        Sig_s_cool_ig += N_cool_i_cm2pb * sigma_s_itgh[*iso][bt_s][h][h];

        Sig_s_fuel_ig += N_fuel_i_cm2pb * sigma_s_itgh[*iso][bt_s][h][g];
        Sig_s_clad_ig += N_clad_i_cm2pb * sigma_s_itgh[*iso][bt_s][h][g];
        Sig_s_cool_ig += N_cool_i_cm2pb * sigma_s_itgh[*iso][bt_s][h][g];
*/
      };

      // Calculate kappas as the inverse of the diffusion length
      // k = 1 / D = 3 * Sigma_tr = 3 (Sigma_t - mu * Sigma_s_gg)
/*
      kappa_fuel_tg[bt_s][g] += (3.0 * N_fuel_i_cm2pb * sigma_t_itg[*iso][bt_s][g]) - (2.0 * Sig_s_fuel_ig / AW_ig); 
      kappa_clad_tg[bt_s][g] += (3.0 * N_fuel_i_cm2pb * sigma_t_itg[*iso][bt_s][g]) - (2.0 * Sig_s_clad_ig / AW_ig); 
      kappa_cool_tg[bt_s][g] += (3.0 * N_fuel_i_cm2pb * sigma_t_itg[*iso][bt_s][g]) - (2.0 * Sig_s_cool_ig / AW_ig); 
*/

      kappa_fuel_tg[bt_s][g] += (N_fuel_i_cm2pb * sigma_a_itg[*iso][bt_s][g]) * Sig_s_fuel_ig * (3.0 - (2.0 / AW_ig));
      kappa_clad_tg[bt_s][g] += (N_clad_i_cm2pb * sigma_a_itg[*iso][bt_s][g]) * Sig_s_clad_ig * (3.0 - (2.0 / AW_ig));
      kappa_cool_tg[bt_s][g] += (N_cool_i_cm2pb * sigma_a_itg[*iso][bt_s][g]) * Sig_s_cool_ig * (3.0 - (2.0 / AW_ig));
    };
  };


  // sqrt(kappa)
  for (g=0; g < G; g++)
  {
    kappa_fuel_tg[bt_s][g] = sqrt(kappa_fuel_tg[bt_s][g]);
    kappa_clad_tg[bt_s][g] = sqrt(kappa_clad_tg[bt_s][g]);
    kappa_cool_tg[bt_s][g] = sqrt(kappa_cool_tg[bt_s][g]);
  };

  // Re-Normalize chi
  double chi_fuel_tot = 0.0;
  double chi_clad_tot = 0.0;
  double chi_cool_tot = 0.0;

  for (g = 0; g < G; g++)
  {
    chi_fuel_tot += chi_fuel_tg[bt_s][g]; 
    chi_clad_tot += chi_clad_tg[bt_s][g]; 
    chi_cool_tot += chi_cool_tg[bt_s][g]; 
  };
  for (g = 0; g < G; g++)
  {
    // Normailze Fuel
    if (chi_fuel_tot == 0.0)
      chi_fuel_tg[bt_s][g] = 0.0;
    else
      chi_fuel_tg[bt_s][g] = chi_fuel_tg[bt_s][g] / chi_fuel_tot;

    // Normailze Cladding
    if (chi_clad_tot == 0.0)
      chi_clad_tg[bt_s][g] = 0.0;
    else
      chi_clad_tg[bt_s][g] = chi_clad_tg[bt_s][g] / chi_clad_tot;

    // Normailze Coolant
    if (chi_cool_tot == 0.0)
      chi_cool_tg[bt_s][g] = 0.0;
    else
      chi_cool_tg[bt_s][g] = chi_cool_tg[bt_s][g] / chi_cool_tot;
  };


  // Calculate the disadvantage factor, if required.
  if (use_zeta)
    calc_zeta();

  // Calculate Core averaged XS
  double zeta_V_cool_g, denom_g;
  for (g = 0; g < G; g++)
  {
    zeta_V_cool_g = zeta_tg[bt_s][g] * V_cool;
    denom_g = V_fuel + zeta_V_cool_g; 

    Sigma_t_tg[bt_s][g] = (V_fuel * Sigma_t_fuel_tg[bt_s][g] + zeta_V_cool_g * Sigma_t_cool_tg[bt_s][g]) / denom_g;
    Sigma_a_tg[bt_s][g] = (V_fuel * Sigma_a_fuel_tg[bt_s][g] + zeta_V_cool_g * Sigma_a_cool_tg[bt_s][g]) / denom_g;
    nubar_Sigma_f_tg[bt_s][g] = (V_fuel * nubar_Sigma_f_fuel_tg[bt_s][g] + zeta_V_cool_g * nubar_Sigma_f_cool_tg[bt_s][g]) / denom_g;
    chi_tg[bt_s][g] = (V_fuel * chi_fuel_tg[bt_s][g] + zeta_V_cool_g * chi_cool_tg[bt_s][g]) / denom_g;
    Sigma_f_tg[bt_s][g] = (V_fuel * Sigma_f_fuel_tg[bt_s][g] + zeta_V_cool_g * Sigma_f_cool_tg[bt_s][g]) / denom_g;
    Sigma_gamma_tg[bt_s][g] = (V_fuel * Sigma_gamma_fuel_tg[bt_s][g] + zeta_V_cool_g * Sigma_gamma_cool_tg[bt_s][g]) / denom_g;
    Sigma_2n_tg[bt_s][g] = (V_fuel * Sigma_2n_fuel_tg[bt_s][g] + zeta_V_cool_g * Sigma_2n_cool_tg[bt_s][g]) / denom_g;
    Sigma_3n_tg[bt_s][g] = (V_fuel * Sigma_3n_fuel_tg[bt_s][g] + zeta_V_cool_g * Sigma_3n_cool_tg[bt_s][g]) / denom_g;
    Sigma_alpha_tg[bt_s][g] = (V_fuel * Sigma_alpha_fuel_tg[bt_s][g] + zeta_V_cool_g * Sigma_alpha_cool_tg[bt_s][g]) / denom_g;
    Sigma_proton_tg[bt_s][g] = (V_fuel * Sigma_proton_fuel_tg[bt_s][g] + zeta_V_cool_g * Sigma_proton_cool_tg[bt_s][g]) / denom_g;
    Sigma_gamma_x_tg[bt_s][g] = (V_fuel * Sigma_gamma_x_fuel_tg[bt_s][g] + zeta_V_cool_g * Sigma_gamma_x_cool_tg[bt_s][g]) / denom_g;
    Sigma_2n_x_tg[bt_s][g] = (V_fuel * Sigma_2n_x_fuel_tg[bt_s][g] + zeta_V_cool_g * Sigma_2n_x_cool_tg[bt_s][g]) / denom_g;

    for (h = 0; h < G; h++)
    {
      Sigma_s_tgh[bt_s][g][h] = (V_fuel * Sigma_s_fuel_tgh[bt_s][g][h] + zeta_V_cool_g * Sigma_s_cool_tgh[bt_s][g][h]) / denom_g;
    };
  };

  // Re-Normalize chi
  double chi_tot = 0.0;
  for (g = 0; g < G; g++)
    chi_tot += chi_tg[bt_s][g]; 
  for (g = 0; g < G; g++)
    chi_tg[bt_s][g] = chi_tg[bt_s][g] / chi_tot;
  if (chi_tot == 0.0)
  {
    for (g = 0; g < G; g++)
      chi_tg[bt_s][g] = 0.0;
  };
};





void bright::ReactorMG::assemble_multigroup_matrices()
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
  F_fuel_tgh[bt_s] = bright::vector_outer_product(chi_fuel_tg[bt_s], nubar_Sigma_f_fuel_tg[bt_s]);

  //F_tgh[bt_s] = bright::vector_outer_product(chi_tg[bt_s], nubar_Sigma_f_tg[bt_s]);
  //F_tgh[bt_s] = bright::vector_outer_product(chi_tg[bt_s], nubar_Sigma_f_fuel_tg[bt_s]);
  F_tgh[bt_s] = bright::vector_outer_product(chi_fuel_tg[bt_s], nubar_Sigma_f_fuel_tg[bt_s]);

  //F_fuel_tgh[bt_s] = bright::vector_outer_product(nubar_Sigma_f_fuel_tg[bt_s], chi_fuel_tg[bt_s]);
  //F_tgh[bt_s] = bright::vector_outer_product(nubar_Sigma_f_tg[bt_s], chi_tg[bt_s]);

  // Grab the inverse of the A matrix
  A_inv_fuel_tgh[bt_s] = bright::matrix_inverse(A_fuel_tgh[bt_s]);
  A_inv_clad_tgh[bt_s] = bright::matrix_inverse(A_clad_tgh[bt_s]);
  A_inv_cool_tgh[bt_s] = bright::matrix_inverse(A_cool_tgh[bt_s]);
  A_inv_tgh[bt_s] = bright::matrix_inverse(A_tgh[bt_s]);

  // Multiply the inverse of A by F
  A_inv_F_fuel_tgh[bt_s] = bright::matrix_multiplication(A_inv_fuel_tgh[bt_s], F_fuel_tgh[bt_s]);
  A_inv_F_tgh[bt_s] = bright::matrix_multiplication(A_inv_tgh[bt_s], F_tgh[bt_s]);

};






void bright::ReactorMG::assemble_transmutation_matrices()
{
  //
  // Assemble the energy integral of transmutation matrix
  //
  int g, i, j, ind, jnd;
  std::vector< bright::SparseMatrix<double> > T_matrix = std::vector< bright::SparseMatrix<double> > (G,  bright::SparseMatrix<double>(fast_yield_matrix.size(), K_num, K_num));
  std::vector< bright::sparse_matrix_entry<double> >::iterator fpy_iter, fpy_end;
  
  
  // Add the cross sections
  double fpy, sig;
  int fpy_ind;
  int ind_PU239 = K_ind[942390];
  int j_gamma, j_2n, j_3n, j_alpha, j_proton, j_gamma_x, j_2n_x; 
  int jnd_gamma, jnd_2n, jnd_3n, jnd_alpha, jnd_proton, jnd_gamma_x, jnd_2n_x;
  for (ind = 0; ind < K_num; ind++)
  {
    i = K_ord[ind];

    // Add diag entries
    for (g = 0; g < G; g++)
      T_matrix[g].push_back(ind, ind, -(sigma_f_itg[i][bt_s][g] + \
                                        sigma_gamma_itg[i][bt_s][g] + \
                                        sigma_2n_itg[i][bt_s][g] + \
                                        sigma_3n_itg[i][bt_s][g] + \
                                        sigma_alpha_itg[i][bt_s][g] + \
                                        sigma_proton_itg[i][bt_s][g] + \
                                        sigma_gamma_x_itg[i][bt_s][g] + \
                                        sigma_2n_x_itg[i][bt_s][g]));

    // Add the fission source
    for (g = 0; g < G; g++)
    {
      sig = sigma_f_itg[i][bt_s][g];
      if (sig == 0.0)
        continue;

      fpy_end = fission_product_yield_matrix[g].sm.end();
      fpy_iter = bright::find_row(fission_product_yield_matrix[g].sm.begin(), fpy_end, ind);
      fpy_ind = ind;

      // Deafult to Pu239 FP if yields not available
      if (fpy_iter == fpy_end)
      {
        fpy_iter = bright::find_row(fission_product_yield_matrix[g].sm.begin(), fpy_end, ind_PU239);
        fpy_ind = ind_PU239;
      };

      while((*fpy_iter).row == fpy_ind)
      {
        jnd = (*fpy_iter).col;
        fpy = (*fpy_iter).val;
        T_matrix[g].push_back(ind, jnd, fpy * sig);
        fpy_iter++;
      };
    };

    // Add the capture cross-section
    j_gamma = ((i/10) + 1) * 10;
    if (K.count(j_gamma) == 1)
    {
      jnd_gamma = K_ind[j_gamma];
      for (g = 0; g < G; g++)
      {
        sig = sigma_gamma_itg[i][bt_s][g];
        if (sig == 0.0)
          continue;

        T_matrix[g].push_back(ind, jnd_gamma, sig);
      };
    };

    // Add the (n, 2n) cross-section
    j_2n = j_gamma - 20;
    if (K.count(j_2n) == 1)
    {
      jnd_2n = K_ind[j_2n];
      for (g = 0; g < G; g++)
      {
        sig = sigma_2n_itg[i][bt_s][g];
        if (sig == 0.0)
          continue;

        T_matrix[g].push_back(ind, jnd_2n, sig);
      };
    };

    // Add the (n, 3n) cross-section
    j_3n = j_2n - 10;
    if (K.count(j_3n) == 1)
    {
      jnd_3n = K_ind[j_3n];
      for (g = 0; g < G; g++)
      {
        sig = sigma_3n_itg[i][bt_s][g];
        if (sig == 0.0)
          continue;

        T_matrix[g].push_back(ind, jnd_3n, sig);
      };
    };

    // Add the (n, alpha) cross-section
    j_alpha = j_3n - 20010;
    if (K.count(j_alpha) == 1)
    {
      jnd_alpha = K_ind[j_alpha];
      for (g = 0; g < G; g++)
      {
        sig = sigma_alpha_itg[i][bt_s][g];
        if (sig == 0.0)
          continue;

        T_matrix[g].push_back(ind, jnd_alpha, sig);
      };
    };

    // Add the (n, proton) cross-section
    j_proton = j_gamma - 10010;
    if (K.count(j_proton) == 1)
    {
      jnd_proton = K_ind[j_proton];
      for (g = 0; g < G; g++)
      {
        sig = sigma_proton_itg[i][bt_s][g];
        if (sig == 0.0)
          continue;

        T_matrix[g].push_back(ind, jnd_proton, sig);
      };
    };

    // Add the capture (excited) cross-section
    j_gamma_x = j_gamma + 1;
    if (K.count(j_gamma_x) == 1)
    {
      jnd_gamma_x = K_ind[j_gamma_x];
      for (g = 0; g < G; g++)
      {
        sig = sigma_gamma_x_itg[i][bt_s][g];
        if (sig == 0.0)
          continue;

        T_matrix[g].push_back(ind, jnd_gamma_x, sig);
      };
    };

    // Add the (n, 2n *) cross-section
    j_2n_x = j_2n + 1;
    if (K.count(j_2n_x) == 1)
    {
      jnd_2n_x = K_ind[j_2n_x];
      for (g = 0; g < G; g++)
      {
        sig = sigma_2n_x_itg[i][bt_s][g];
        if (sig == 0.0)
          continue;

        T_matrix[g].push_back(ind, jnd_2n_x, sig);
      };
    };
  };

  for (g = 0; g < G; g++)
    T_matrix[g].clean_up();



  // Multiply by flux and integrate over energy
  double adj_phi;
  for (g = 0; g < G; g++)
  {
    // Adjust the flux value for burning
    adj_phi = pyne::cm2_per_barn * phi_tg[bt_s][g];

    // Skip this group if there is no flux 
    if (adj_phi == 0.0)
      continue;

    T_int_tij[bt_s] = T_int_tij[bt_s] + (T_matrix[g] * adj_phi);
  };
  

  // Make the transmutation matrix for this time step
  M_tij[bt_s] = (T_int_tij[bt_s] + decay_matrix);

  // Add initial transmutatio chains
  if (bt_s == 0)
  {
    std::vector< bright::sparse_matrix_entry<double> >::iterator M_iter, M_end;
    M_iter = M_tij[bt_s].sm.begin();
    M_end = M_tij[bt_s].sm.end();

    // Initialize the chains container
    for ( ; M_iter != M_end; M_iter++)
    {
      ind = (*M_iter).row;
      i = K_ord[ind];

      jnd = (*M_iter).col;
      j = K_ord[jnd];

      if (i == j)
        continue;

      if (transmutation_chains.count(i) == 0)
        transmutation_chains[i] = std::map<int, std::vector< std::vector<int> > > ();

      transmutation_chains[i][j] = std::vector< std::vector<int> >(1,  std::vector<int>(2));
      transmutation_chains[i][j][0][0] = i;
      transmutation_chains[i][j][0][1] = j;
    };
  };

};





void bright::ReactorMG::add_transmutation_chains(std::vector<int> tc)
{
  int n, nik, Nik;
  int chain_len = tc.size();
  int i = tc[0];
  int j = tc[chain_len - 1];
  int k, jnd, knd;
  bool k_in_chain = false;
  bool chain_present = false;
  bool chain_ind_same = false;
  double br_ij, br_ik;
  double branch_ratio_cutoff_point;

  jnd = K_ind[j];


  // I swear this 42 number is meaningful...
  if (21 < chain_len)
    return;

  if (j < 860000)
  {
    if (i < 860000)
      branch_ratio_cutoff_point = 5E-4;
    else
      branch_ratio_cutoff_point = 5E-2;
  }
  else
    branch_ratio_cutoff_point = 5E-3;

  br_ij = 1.0;
  for (n = 1; n < chain_len; n++)
    br_ij *= branch_ratios[K_ind[tc[n-1]]][K_ind[tc[n]]]; 

  if (br_ij < branch_ratio_cutoff_point)
    return;

  std::vector<int> next_chain;

  for (knd = 0; knd < K_num; knd++)
  {
    br_ik = br_ij * branch_ratios[jnd][knd];
    if (br_ik < branch_ratio_cutoff_point)
      continue;

    if (jnd == knd)
      continue;

    k = K_ord[knd];

    if (transmutation_chains[i].count(k) == 0)
      transmutation_chains[i][k] = std::vector< std::vector<int> >();

    // Don't allow cyclic chains
    k_in_chain = false;
    for (n = 0; n < chain_len; n++)
      if (k == tc[n])
        k_in_chain = true;

    if (k_in_chain)
      continue;

    // construct new chains
    next_chain = std::vector<int>(tc);
    next_chain.push_back(k);
    int next_chain_size = next_chain.size();
  
    // chack if the chain is already present
    chain_present = false;
    Nik = transmutation_chains[i][k].size();
    for (nik = 0; nik < Nik; nik++)
    {
      if (next_chain_size == transmutation_chains[i][k][nik].size())
      {
        int ncp  = 0;
        chain_ind_same = true;
        //while (chain_ind_same && ncp < next_chain_size)
        for (ncp = 0; chain_ind_same && ncp < next_chain_size; ncp++)
        {
          chain_ind_same = chain_ind_same && (next_chain[ncp] == transmutation_chains[i][k][nik][ncp]);
//          std::cout << "            chain = " << next_chain[ncp] << " --- " << transmutation_chains[i][k][nik][ncp] << "  " << chain_ind_same << "  " << ncp << "/" << next_chain_size << "\n";
          //ncp++;
        };

//        std::cout << "          ---------------\n";

        chain_present = (chain_present || chain_ind_same);
      };
    };

    if (chain_present)
    {
      std::cout << "        Present chains = " << i << " --> " << j << " --> " << k << "  " << chain_present << "  " << chain_ind_same << "  " << next_chain_size << "  " << Nik << "\n";
      continue;
    };

    std::cout << "      Adding chains = " << i << " --> " << j << " --> " << k << "  " << chain_present << "  " << chain_ind_same << "  " << next_chain.size() << "  " << Nik << "  " << branch_ratio_cutoff_point << "\n";

    // add new chains
    transmutation_chains[i][k].push_back(next_chain);
    add_transmutation_chains(next_chain);
  };

};






double bright::ReactorMG::bateman_chain(int i, int j, int c, double t)
{
  // Solves the Bateman Equations for a unit mass of isotope i -> 
  // decaying into isotope j.
  int ind = K_ind[i];
  int jnd = K_ind[j];
  
  std::vector<int> chain = transmutation_chains[i][j][c];

  int n;
  int N = chain.size();

  int qnd, rnd;
  qnd = ind;

  double B = 1.0;
  double alpha_num = 1.0;

  double trans_cutoff = 1E+3 / t;

  for (n = 0; n < N - 1; n++)
  {
    rnd = K_ind[chain[n+1]];
    B *= branch_ratios[qnd][rnd];

    if (trans_consts[qnd] < trans_cutoff)
      alpha_num *= trans_consts[qnd];
    //alpha_num *= trans_consts[qnd];

    qnd = rnd;
  };

  // Check trivial results
  if ((B < branch_ratio_cutoff) || (alpha_num == 0.0))
    return 0.0;

  int m;
  double alpha_den, sum_part;
  for (n = 0; n < N; n++)
  {
    qnd = K_ind[chain[n]];
    alpha_den = 1.0;

    if (trans_cutoff < trans_consts[qnd])
      continue;

    for (m = 0; m < N; m++)
    {
      rnd = K_ind[chain[m]];

      if (trans_cutoff < trans_consts[rnd])
        continue;

      // Debug for equal lambdas
      //if (trans_consts[rnd] == trans_consts[qnd] && rnd != qnd)
      //    std::cout << "    trans constants equal for " << chain[n] << ", " << chain[m] << " = " << trans_consts[qnd] << ", " << trans_consts[rnd] << "\n";

      if (n != m)
        alpha_den *= (trans_consts[rnd] - trans_consts[qnd]);
     };

    sum_part += (exp(-trans_consts[qnd] * t) / alpha_den);
  };

  double mass_frac = B * alpha_num * sum_part;

  // Sanity check
  if (!(0.0 < mass_frac) || 1.0 < mass_frac)
    mass_frac = 0.0;

  return mass_frac;
}





double bright::ReactorMG::bateman(int i, int j, double t)
{
  int c, C;
  double total_mass_frac = 0.0;

  C = transmutation_chains[i][j].size();
  for (c = 0; c < C; c++)
  {
    total_mass_frac += bateman_chain(i, j, c, t);
    //std::cout <<  i << " --> " << j << "  " << c << "/" << C << "   " << total_mass_frac << "\n";
  };

  return total_mass_frac;
};






void bright::ReactorMG::calc_criticality()
{
  // Init values
  int n = 0;
  int N = 100;

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
    phi1 = bright::scalar_matrix_vector_product(1.0 / k0, A_inv_F_tgh[bt_s], phi0);

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
    n++;
  };

  // Set the final flux values to the class members
  std::cout << "   k0 = " << k0 << "\n";

  // Normalize the flux
  double phi1_tot = 0.0;
  for (g = 0; g < G; g++)
    phi1_tot += phi1[g];
  for (g = 0; g < G; g++)
    phi_tg[bt_s][g] = phi1[g] / phi1_tot;

  // Rescale t flux dependingon the burnup method used.
  if (burnup_via_constant == "flux")
  {
    phi_t[bt_s] = flux;
    for (g = 0; g < G; g++)
      phi_tg[bt_s][g] *= flux;
  }
  else if (burnup_via_constant == "power")
  {
    double norm_fission_reaction_rate = 0.0;
    for (g = 0; g < G; g++)
      norm_fission_reaction_rate += (Sigma_f_fuel_tg[bt_s][g] * phi_tg[bt_s][g]);

    // Rescale the flux
    phi_t[bt_s] = specific_power * rho_fuel * 1E3 / (3.28446179835e-11 * norm_fission_reaction_rate);
    for (g = 0; g < G; g++)
      phi_tg[bt_s][g] *= phi_t[bt_s];

    std::cout << "   nfrr = " << norm_fission_reaction_rate << "\n";
    std::cout << "   flux = " << phi_t[bt_s] << "\n";
  }
  else
    std::cout << "burnup_via_constant is not setup properly\n";


  if (bt_s == 0)
    Phi_t[bt_s] = 0.0;
  else
    Phi_t[bt_s] = Phi_t[bt_s - 1] + (phi_t[bt_s] * (burn_times[bt_s] - burn_times[bt_s-1]) * pyne::sec_per_day * pyne::cm2_per_barn * 1E+3);

  // Calculate the multiplication factor, physically, and not from the eigenvalue
  double k_num = 0.0;
  double k_den = 0.0;

  for (g = 0; g < G; g++)
  {
    k_num += V_fuel * nubar_Sigma_f_fuel_tg[bt_s][g] * phi_tg[bt_s][g];
    k_den += phi_tg[bt_s][g] * ((V_fuel * Sigma_a_fuel_tg[bt_s][g]) + (zeta_tg[bt_s][g] * V_cool * Sigma_a_cool_tg[bt_s][g]));
  };

  k_t[bt_s] = P_NL * k_num / k_den;
};









void bright::ReactorMG::calc_transmutation()
{
  // Calculates a tranmutation step via the Pade method
  int i, j, ind, jnd;

  // Get the transmutation matrix for this time delta
  double dt = (burn_times[bt_s + 1] - burn_times[bt_s]) * pyne::sec_per_day;

  int knd, qnd;
  trans_consts = std::vector<double>(K_num, 0.0);
  branch_ratios = M_tij[bt_s].todense();
  for (knd = 0; knd < K_num; knd++)
  {
    trans_consts[knd] = -branch_ratios[knd][knd];

    if (trans_consts[knd] == 0.0)
      continue;

    branch_ratios[knd][knd] = 0.0;
    for (qnd = 0; qnd < K_num; qnd++)
      branch_ratios[knd][qnd] = branch_ratios[knd][qnd] / trans_consts[knd];
  };


  // Fill in the chains container with more than one-step values
  if (bt_s == 0)
  {
    for (nuc_iter iso = K.begin(); iso != K.end(); iso++)
    {
      i = (*iso);
      ind = K_ind[i];

      for (jnd = 0; jnd < K_num; jnd++)
      {
        if (branch_ratios[ind][jnd] == 0.0)
          continue;

        j = K_ord[jnd];

        std::cout << "    Adding chains for " << i << " --> " << j << "\n";
        for (int ncp = 0; ncp < transmutation_chains[i][j].size(); ncp++)
          add_transmutation_chains(transmutation_chains[i][j][ncp]);
      };
    };
  };



  // Make mass vectors
  std::vector<double> comp_prev (K_num, 0.0);
  for (ind = 0; ind < K_num; ind++)
  {
    i = K_ord[ind];
    comp_prev[ind] = T_it[i][bt_s];
  };


  std::vector<double> comp_next (K_num, 0.0);
  for (ind = 0; ind < K_num; ind++)
  {
    if (comp_prev[ind] == 0.0)
      continue;

    i = K_ord[ind];
    if (J.count(i) == 0)
      continue;

    for (jnd = 0; jnd < K_num; jnd++)
    {
      j = K_ord[jnd];

      //if (J.count(j) == 0)
      //    continue;


      if (i == j)
        comp_next[ind] += comp_prev[ind] * exp(-trans_consts[ind] * dt);
      else if (transmutation_chains[i].count(j) == 0)
        continue;
      else
        comp_next[jnd] += comp_prev[ind] * bateman(i, j, dt);
    };
  };

  // Copy this composition back to the tranmutuation matrix
  for (ind = 0; ind < K_num; ind++)
  {
    i = K_ord[ind];
    T_it[i][bt_s+1] = comp_next[ind];
  };

  // Calculate the burnup 
  pyne::comp_map cd_prev, cd_next;
  for (ind = 0; ind < K_num; ind++)
  {
    i = K_ord[ind];
    cd_prev[i] = comp_prev[ind];
    cd_next[i] = comp_next[ind];
  };

  double delta_BU;

  if (burnup_via_constant == "flux")
  {
    pyne::Material ms_prev (cd_prev);
    pyne::Material act_prev = ms_prev.sub_act();

    pyne::Material ms_next (cd_next);
    pyne::Material act_next = ms_next.sub_act();

    delta_BU = (act_prev.mass - act_next.mass) * 931.46;
  }
  else if (burnup_via_constant == "power")
  {
    delta_BU = specific_power * (burn_times[bt_s + 1] - burn_times[bt_s]);
  }
  else
    std::cout << "burnup_via_constant not set properly!\n";

  BU_t[bt_s+1] = delta_BU + BU_t[bt_s];

  std::cout << "   BU_t = " << BU_t[bt_s+1] << "\n";
};












void bright::ReactorMG::init_core()
{
  // Burns up the core and fills in parameter values as we go.


  // Initialize the transmutation matrix with values from mat_feed
  T_it.clear();

  // Also initialize the cross-section matrices as a function of time.
  sigma_t_itg.clear();
  sigma_a_itg.clear();
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
  BU_t = std::vector<double>(S, 0.0);

  zeta_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 1.0) );
  lattice_E_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0) );
  lattice_F_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0) );


  for (nuc_iter iso = K.begin(); iso != K.end(); iso++)
  {
    // Init the transmutation matrix
    T_it[*iso] = time_data(S, 0.0);

    if (0 < mat_feed.comp.count(*iso))
      T_it[*iso][0] = mat_feed.comp[*iso];

    // Init the cross-sections
    sigma_t_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
    sigma_a_itg[*iso] = std::vector< std::vector<double> >(S, std::vector<double>(G, -1.0));
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
  Sigma_a_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
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
  kappa_fuel_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));

  Sigma_t_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_a_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  nubar_Sigma_f_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  chi_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_s_clad_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
  Sigma_f_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_gamma_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_2n_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_3n_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_alpha_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_proton_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_gamma_x_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_2n_x_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  kappa_clad_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));

  Sigma_t_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_a_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  nubar_Sigma_f_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  chi_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_s_cool_tgh = std::vector< std::vector< std::vector<double> > >(S, std::vector< std::vector<double> >(G, std::vector<double>(G, 0.0)));
  Sigma_f_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_gamma_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_2n_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_3n_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_alpha_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_proton_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_gamma_x_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_2n_x_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  kappa_cool_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));


  Sigma_t_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
  Sigma_a_tg = std::vector< std::vector<double> >(S, std::vector<double>(G, 0.0));
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

  T_int_tij = std::vector< bright::SparseMatrix<double> > (S,  bright::SparseMatrix<double> (0, K_num, K_num));
  M_tij = std::vector< bright::SparseMatrix<double> > (S,  bright::SparseMatrix<double> (0, K_num, K_num));

  // Init the multilplcation factpr
  k_t = std::vector<double>(S, -1.0);
};

void bright::ReactorMG::burnup_core()
{
  // prep the core
  init_core();
  
  // Loop through all time steps
  for (int s = 0; s < S; s++)
  {
    // Set the current time
    bt_s = s;
    burn_time = burn_times[s];

    if (2 <= bright::verbosity)
      std::cout << "Time step " << s << " = " << burn_times[s] << " days\n";

    // Find the nearest neightbors for this time.
    calc_nearest_neighbors();

    // Interpolate cross section in preparation for 
    // criticality calculation.
    interpolate_cross_sections();

    // Fold the mass weights for this time step
    calc_mass_weights();
    fold_mass_weights();

    // Preform the criticality calulation
    assemble_multigroup_matrices();
    calc_criticality();

    // Preform the burnup calulation
    assemble_transmutation_matrices();

    if (s != S - 1)
      calc_transmutation();
  };

};









void bright::ReactorMG::calc_T_itd()
{
  /** Calculates the output isotopics of Mj(Fd).
   *  NOTE: Mj(Fd) is effectively the same variable as mat_prod before normalization!
   */

  pyne::comp_map tempOut;

  //Checks to see if the discharge index in at the end of the fluence table
  bool td_nLast = false;
  if ( (td_n+1) == S )
    td_nLast = true;			

  for (nuc_iter j = J.begin(); j != J.end(); j++ )
  {
    if (td_nLast)
      tempOut[*j] = pyne::solve_line(Phid, Phi_t[td_n], T_it[*j][td_n], Phi_t[td_n-1], T_it[*j][td_n-1]);
    else
      tempOut[*j] = pyne::solve_line(Phid, Phi_t[td_n+1], T_it[*j][td_n+1], Phi_t[td_n], T_it[*j][td_n]);
  };

  mat_prod = pyne::Material(tempOut);	
};



void bright::ReactorMG::calc_mat_prod()
{
  //Wrapper to calculate discharge isotopics.
  calc_T_itd();
};






void bright::ReactorMG::calc_sub_mats()
{
  //Sets possibly relevant reactor input and output substreams.

  //Uranium
  mat_feed_u  = mat_feed.sub_u();
  mat_prod_u = mat_prod.sub_u();

  //TRU
  try 
  {
    mat_feed_tru = mat_feed.sub_tru();
  }
  catch (...)
  {
    pyne::comp_map cd;
    cd[942390] = 1.0;
    mat_feed_tru = pyne::Material(cd, 1.0);
    mat_feed_tru.mass = 0.0;
  };
  mat_prod_tru = mat_prod.sub_tru();

  //Lanthanides
  try
  {
    mat_feed_lan = mat_feed.sub_lan();
  }
  catch (...)
  {
    pyne::comp_map cd;
    cd[581440] = 1.0;
    mat_feed_lan  = pyne::Material(cd, 1.0);
    mat_feed_lan.mass = 0.0;
  };
  mat_prod_lan = mat_prod.sub_lan();

  //Actinides
  mat_feed_act  = mat_feed.sub_act();
  mat_prod_act = mat_prod.sub_act();
};




double bright::ReactorMG::calc_tru_cr()
{
  //Calculates the reactor's transuranic conversion ratio.
  tru_cr = 1.0 - ((mat_feed_tru.mass - mat_prod_tru.mass) / (BUd/931.46));
  return tru_cr;
};






bright::FluencePoint bright::ReactorMG::fluence_at_BU(double BU)
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
    fp.m = pyne::slope(Phi_t[fp.f], BU_t[fp.f], Phi_t[fp.f-1], BU_t[fp.f-1]);
  else
    fp.m = pyne::slope(Phi_t[fp.f+1], BU_t[fp.f+1], Phi_t[fp.f], BU_t[fp.f]);

  fp.F = ((BU - BU_t[fp.f])/fp.m) + Phi_t[fp.f];

  return fp;
};






double bright::ReactorMG::batch_average_k(double BUd)
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
      k_b[b] = pyne::solve_line(fps[b].F, Phi_t[fps[b].f], k_t[fps[b].f], Phi_t[fps[b].f-1], k_t[fps[b].f-1]);
      phi_b[b] = pyne::solve_line(fps[b].F, Phi_t[fps[b].f], phi_t[fps[b].f], Phi_t[fps[b].f-1], phi_t[fps[b].f-1]);
    }
    else
    {
      k_b[b] = pyne::solve_line(fps[b].F, Phi_t[fps[b].f+1], k_t[fps[b].f+1], Phi_t[fps[b].f], k_t[fps[b].f]);
      phi_b[b] = pyne::solve_line(fps[b].F, Phi_t[fps[b].f+1], phi_t[fps[b].f+1], Phi_t[fps[b].f], phi_t[fps[b].f]);
    };
  };

  // Calculate the flux weighted avearge
  double numerator = 0.0;
  double denominator = 0.0;
  for (int b = 0; b < B; b++)
  {
    numerator   += (k_b[b] * phi_b[b]);
    denominator += phi_b[b];

    //numerator   += (k_b[b] * fps[b].F);
    //denominator += fps[b].F;
  };

  double batch_ave_k = numerator/denominator;

  return batch_ave_k;
};




void bright::ReactorMG::BUd_bisection_method()
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

  if (1 < bright::verbosity)
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

    if (1 < bright::verbosity)
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
      if (0 < bright::verbosity)
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
    if (0 < bright::verbosity)
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
  if (0 < bright::verbosity)
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
    if (0 < bright::verbosity)
      std::cout << "k did not converge with the Bisection Method to an accuracy of " << DoA << " in " << q << " iterations.\n";

    if ( (fabs(k_a - 1.0) < 0.01) && (fabs(k_a - 1.0) < fabs(k_b -1.0)) )
    {
      BUd = BUd_a;
      k = k_a;
      if (0 < bright::verbosity)
        std::cout << "However, k_a is within 1% of 1 and closer to 1 than k_b; using these values.\n";
    }
    else if ( (fabs(k_b - 1.0) < 0.01) && (fabs(k_b - 1.0) < fabs(k_a -1.0)) )
    {
      BUd = BUd_b;
      k = k_b;
      if (0 < bright::verbosity)
        std::cout << "However, k_b is within 1% of 1 and closer to 1 than k_a; using these values.\n";
    }
    else
    {
      if (0 < bright::verbosity)
        std::cout << "Alright.  It really didn't converge. Neither k_a nor k_b is within 1% of 1. Program will likely fail!\n";
    };
  };

  FluencePoint fp = fluence_at_BU(BUd);
  td_n = fp.f;    // lower index of discharge time
  Phid = fp.F;    // Discharge fluence

  // time at discharge
  if ( (fp.f + 1) == S )
    td = pyne::solve_line(fp.F, Phi_t[fp.f], burn_times[fp.f], Phi_t[fp.f-1], burn_times[fp.f-1]);
  else
    td = pyne::solve_line(fp.F, Phi_t[fp.f+1], burn_times[fp.f+1], Phi_t[fp.f], burn_times[fp.f]);

  return;
};




void bright::ReactorMG::run_P_NL(double temp_pnl)
{
  /** Does a reactor run for a specific P_NL.
   *  Requires that mat_feed be (meaningfully) set.
   *  For use with calibrate_P_NL_to_BUd
   */

  P_NL = temp_pnl;
  burnup_core();
  BUd_bisection_method();
};




void bright::ReactorMG::calibrate_P_NL_to_BUd()
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
      if (0 < bright::verbosity)
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

  if (0 < bright::verbosity)
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



pyne::Material bright::ReactorMG::calc()
{
  // Finds BUd and output isotopics.
  burnup_core();

  BUd_bisection_method();

  calc_mat_prod();

  return mat_prod;
};


pyne::Material bright::ReactorMG::calc (pyne::comp_map incomp)
{
  // Finds BUd and output isotopics.
  mat_feed = pyne::Material (incomp);

  return calc();
};


pyne::Material bright::ReactorMG::calc (pyne::Material mat)
{
  // Finds BUd and output isotopics.
  mat_feed = mat;

  return calc();
};










// 
// Lattice and Zeta functions below
//



void bright::ReactorMG::lattice_E_planar(double a, double b)
{
  for (int g = 0; g < G; g++)
  {
    if (0.0 == kappa_cool_tg[bt_s][g])
      lattice_E_tg[bt_s][g] = 0.0;
    else
      lattice_E_tg[bt_s][g] = kappa_cool_tg[bt_s][g] * (b - a) * pyne::coth(kappa_cool_tg[bt_s][g]*(b-a));
  };
};






void bright::ReactorMG::lattice_F_planar(double a, double b)
{
  for (int g = 0; g < G; g++)
  {
    if (0.0 == kappa_fuel_tg[bt_s][g])
      lattice_F_tg[bt_s][g] = 0.0;
    else
      lattice_F_tg[bt_s][g] = kappa_fuel_tg[bt_s][g] * a * pyne::coth(kappa_fuel_tg[bt_s][g]*a) ;
  };
};




void bright::ReactorMG::lattice_E_spherical(double a, double b)
{
  double coef, num, den;

  for (int g = 0; g < G; g++)
  {
    if (0.0 == kappa_cool_tg[bt_s][g])
      lattice_E_tg[bt_s][g] = 0.0;
    else
    {
      coef = pow(kappa_cool_tg[bt_s][g], 3) * (pow(b,3) - pow(a,3)) / (3.0 * kappa_cool_tg[bt_s][g] * a);
      num = 1.0 - (kappa_cool_tg[bt_s][g] * b * pyne::coth(kappa_cool_tg[bt_s][g] * (b - a)));
      den = 1.0 - (pow(kappa_cool_tg[bt_s][g], 2) * a * b) - (kappa_cool_tg[bt_s][g] * (b - a) * pyne::coth(kappa_cool_tg[bt_s][g] * (b - a)));
      lattice_E_tg[bt_s][g] = coef * num / den;
    };
  };
};
  




void bright::ReactorMG::lattice_F_spherical(double a, double b)
{
  double coef, num, den;

  for (int g = 0; g < G; g++)
  {
    coef = pow(kappa_fuel_tg[bt_s][g], 2) * pow(a, 2) / 3.0;
    num = pyne::tanh(kappa_fuel_tg[bt_s][g] * a); 
    den = (kappa_fuel_tg[bt_s][g] * a) - pyne::tanh(kappa_fuel_tg[bt_s][g] * a); 
    lattice_F_tg[bt_s][g] = coef * num / den;
  };
};





void bright::ReactorMG::lattice_E_cylindrical(double a, double b)
{
  namespace bm = boost::math;
  double coef, num, den;

  for (int g = 0; g < G; g++)
  {
    if (0.0 == kappa_cool_tg[bt_s][g])
      lattice_E_tg[bt_s][g] = 0.0;
    else
    {
      coef = kappa_cool_tg[bt_s][g] * (pow(b, 2) - pow(a, 2)) / (2.0 * a);

      num = (bm::cyl_bessel_i(0, kappa_cool_tg[bt_s][g] * a) * bm::cyl_bessel_k(1, kappa_cool_tg[bt_s][g] * b)) + \
          (bm::cyl_bessel_k(0, kappa_cool_tg[bt_s][g] * a) * bm::cyl_bessel_i(1, kappa_cool_tg[bt_s][g] * b));

      den = (bm::cyl_bessel_i(1, kappa_cool_tg[bt_s][g] * b) * bm::cyl_bessel_k(1, kappa_cool_tg[bt_s][g] * a)) - \
          (bm::cyl_bessel_k(1, kappa_cool_tg[bt_s][g] * b) * bm::cyl_bessel_i(1, kappa_cool_tg[bt_s][g] * a));

      lattice_E_tg[bt_s][g] = coef * num / den;
    };
  };
};



  
void bright::ReactorMG::lattice_F_cylindrical(double a, double b)
{
  namespace bm = boost::math;
  double num, den;

  for (int g = 0; g < G; g++)
  {
    if (0.0 == kappa_fuel_tg[bt_s][g])
      lattice_E_tg[bt_s][g] = 0.0;
    else
    {
      num =  kappa_fuel_tg[bt_s][g] * a * bm::cyl_bessel_i(0, kappa_fuel_tg[bt_s][g] * a);
      den = 2.0 * bm::cyl_bessel_i(1, kappa_fuel_tg[bt_s][g] * a);
      lattice_F_tg[bt_s][g] = num / den;
    };
  };
};



void bright::ReactorMG::calc_zeta()
{
  // Computes the thermal disadvantage factor

  // Calculate the lattice_flag Functions
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
    //a = r_fuel + ((r_clad - r_fuel) / 2.0);
    b = pitch / sqrt(pyne::pi); // radius of cylinder with an equivilent cell volume

    lattice_E_cylindrical(a, b);
    lattice_F_cylindrical(a, b);
  }
  else
  {
    if (0 < bright::verbosity)
      std::cout << "Did not specify use of planar or spheical or cylindrical lattice functions! Assuming cylindrical...\n";
    
    a = r_fuel;
    b = pitch / sqrt(pyne::pi); //radius of cylinder with an equivilent cell volume

    lattice_E_cylindrical(a, b);
    lattice_F_cylindrical(a, b);
  };

  // Finally, Calculate Zeta
  int g, h;
  for (g = 0; g < G; g++) 
  {
    if (0.0 == Sigma_a_cool_tg[bt_s][g])
      zeta_tg[bt_s][g] = 1.0;
    else
      //zeta_tg[bt_s][g] = lattice_F_tg[bt_s][g] + (Sigma_a_fuel_tg[bt_s][g] * (lattice_E_tg[bt_s][g] - 1.0) / (Sigma_a_cool_tg[bt_s][g]));
      zeta_tg[bt_s][g] = lattice_F_tg[bt_s][g] + (Sigma_a_fuel_tg[bt_s][g] * V_fuel * (lattice_E_tg[bt_s][g] - 1.0) / (Sigma_a_cool_tg[bt_s][g] * V_cool));
  };



  // Unfortunately, the above formulation for the disadvantage factor is ONLY valid for a << b!!!
  // Often times in modern (thermal) reactors, this is not the case.
  // We have a 'thin moderator' situation.
  //
  // To fix this problem correctly requires going to a multi-region diffusion/transport calculation.
  // Doing so is beyond the scope of this code.
  // What is more in-line with current practice is to use the results of a more sophisticated method,
  // interpolate over them, and use them here.
  //
  // That is what is done here when 0.1 < VF / VC, (ie the fuel is greater than 10% of the coolant)
  // A baseline zeta is determined from data presented in "Thermal disadvantage factor calculation by 
  // the multiregion collision probability method" by B. Ozgener,  and H. A. Ozgener, Institute of 
  // Nuclear Energy, Istanbul Technical University 80626 Maslak, Istanbul, Turkey, Received 
  // 28 January 2003;  accepted 20 May 2003.  Available online 6 March 2004.
  // This baseline is a function of (VF/VC).
  // 
  // The above calculation of zeta is then used as a scaling factor on the baseline function to 
  // account for variations in fuel composition and fluenece.

  // Check if we are in the proper Fuel-to-Coolant Regime
/*
  double f2c = V_fuel / V_cool;
  if (f2c < 0.1)
    return;

  double zetabase  = 1.30857959 - (0.10656299 * f2c);

  // Find an index that is hopefully thermal
  int g_therm = (2 * G) / 3;
  while (g_therm != G-1 && 0.0001 < E_g[g_therm])
  {
    std::cout << g_therm << "\n";
    g_therm = g_therm + 1;
  };

  double zetaratio = zetabase / zeta_tg[bt_s][g_therm];

  for (g = 0; g < G; g++) 
  {
    zeta_tg[bt_s][g] = zeta_tg[bt_s][g] * zetaratio;

    if (zeta_tg[bt_s][g] < 1.0)
      zeta_tg[bt_s][g] = 1.0;
  };

  // try something else
  for (g = 0; g < G; g++) 
  {
    if (g < g_therm)
      zeta_tg[bt_s][g] = 1.0;
//    else if (G <= g + 1)
//      zeta_tg[bt_s][g] = 1.0;
    else
      zeta_tg[bt_s][g] = 1.3;
  };
*/

};



//template class bright::sparse_matrix_entry<double>;

template class bright::SparseMatrix<double>;
