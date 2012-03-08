// General Library 

#include "bright.h"

//Bright Globals

std::string bright::BRIGHT_DATA = "";

void bright::bright_start()
{
  #ifdef _WIN32
    char * tmpBRIGHT_DATA;
    size_t lenBRIGHT_DATA;
    errno_t errBRIGHT_DATA = _dupenv_s(&tmpBRIGHT_DATA, &lenBRIGHT_DATA, "BRIGHT_DATA");
    if (errBRIGHT_DATA) std::cout << "BRIGHT_DATA Enviromental Variable could not be found\n";
      BRIGHT_DATA = (std::string) tmpBRIGHT_DATA;
  #else
    BRIGHT_DATA = getenv("BRIGHT_DATA");
  #endif
  return;
};



#ifdef _WIN32
  int null_set [1] = {922350};
  std::set<int> bright::track_nucs (null_set, null_set+1);
  std::vector<int> bright::track_nucs_order (null_set, null_set+1);
#else
  int null_set [0] = {};
  std::set<int> bright::track_nucs (null_set, null_set+0);
  std::vector<int> bright::track_nucs_order (null_set, null_set+0);
#endif

int bright::verbosity  = 0;
int bright::write_text = 1;
int bright::write_hdf5 = 0;

std::string bright::output_filename = "fuel_cycle.h5";


void bright::sort_track_nucs()
{
  track_nucs_order = std::vector<int> (track_nucs.begin(), track_nucs.end());
  std::sort(track_nucs_order.begin(), track_nucs_order.end());
};



void bright::load_track_nucs_hdf5(std::string filename, std::string datasetname, bool clear_prev)
{
  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  //Load values into track_nucs from an hdf5 file.
  //If the dataspace name is not given, try some defaults.
  int dslen = 14;
  std::string defaultsets [14] = {
    "/track_nucs",
    "/Isos2Track",
    "/isostrack",   
    "/IsosTrack",
    "/isotrack",   
    "/IsoTrack",    
    "/ToIso",
    "/ToIsos",
    "/ToIso_zz",
    "/ToIso_MCNP",  
    "/FromIso",  
    "/FromIsos",  
    "/FromIso_zz",
    "/FromIso_MCNP"
    };

  //Open file
  H5::H5File isofile(filename, H5F_ACC_RDONLY);

  //Open dataset
  H5::DataSet isoset;
  if (datasetname.length() != 0)
    isoset = isofile.openDataSet(datasetname);
  else
  {
    //Try to grab one of the default data sets
    //...while suppressing HDF5 errors.
    H5::FileIException not_found_error;
    not_found_error.dontPrint();

    int n = 0;
    while (n < dslen)
    {
      try {  
        isoset = isofile.openDataSet(defaultsets[n]);
        break;
      }
      catch(H5::FileIException not_found_error)
      {
      };
      n++;
    };
    if (n == dslen)
      throw H5::FileIException("load_track_nucs", "Dataset not found!");
  };

  //Read in isos from dataset.
  H5::DataSpace isospace = isoset.getSpace();
  hsize_t isolen[1];
  int isodim = isospace.getSimpleExtentDims(isolen, NULL);

  //Try native int data type first
  #ifdef _WIN32
    // Windows VC++ doesn't accept variable length arrays!
    // So let's make the read in array greater than the number of known isotopes...
    int         iso_out_int [5000];
  #else
    int         iso_out_int [isolen[0]];
  #endif
  isoset.read(iso_out_int, H5::PredType::NATIVE_INT);
  //Maybe add other data types in the future... 

  //Clear previous entries
  if (clear_prev)
    track_nucs.clear();

  //load into track_nucs
  for(int n = 0; n < isolen[0]; n++)
    track_nucs.insert(pyne::nucname::zzaaam(iso_out_int[n]));

  // Sort the results
  sort_track_nucs();
};

void bright::load_track_nucs_text(std::string filename, bool clear_prev)
{
  // Check that the file is there
  if (!pyne::file_exists(filename))
    throw pyne::FileNotFound(filename);

  //Clear previous entries
  if (clear_prev)
    track_nucs.clear();

  //open file
  std::fstream isofile;
  isofile.open(filename.c_str(), std::fstream::in);

  char isoraw [20];
  std::string isostr;
  while (!isofile.eof())
  {
    isofile.width(20);
    isofile >> isoraw;
    isostr.assign(isoraw);
    isostr = pyne::remove_characters(isostr, "()[],.;{}!#|");
    track_nucs.insert(pyne::nucname::zzaaam(isostr));
  };

  //close file
  isofile.close();

  // Sort the results
  sort_track_nucs();
};




//
// Some HDF5 helpers
//
H5::StrType bright::nuc_name_type = H5::StrType(0, 6);

hsize_t bright::cinder_g_dims [1] = {63};
H5::ArrayType bright::cinder_g_type = H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 1, bright::cinder_g_dims);


//
// Fission Product HDF5 interface
//

H5::CompType bright::make_fission_desc()
{
  //  Makes a fission compound datatype
  H5::CompType fdesc( sizeof(fission_struct) );

  fdesc.insertMember( "nuc_name", HOFFSET(fission_struct, nuc_name), bright::nuc_name_type);
  fdesc.insertMember( "nuc_zz", HOFFSET(fission_struct, nuc_zz), H5::PredType::NATIVE_INT);

  fdesc.insertMember( "thermal_yield", HOFFSET(fission_struct, thermal_yield), H5::PredType::NATIVE_INT8);
  fdesc.insertMember( "fast_yield", HOFFSET(fission_struct, fast_yield), H5::PredType::NATIVE_INT8);
  fdesc.insertMember( "high_energy_yield", HOFFSET(fission_struct, high_energy_yield), H5::PredType::NATIVE_INT8);

  fdesc.insertMember( "xs", HOFFSET(fission_struct, xs), bright::cinder_g_type);
 
  return fdesc;
};

H5::CompType bright::fission_desc = bright::make_fission_desc();



H5::CompType bright::make_fission_product_yields_desc()
{
  //  Makes a decay isotope compound datatype
  H5::CompType fpydesc( sizeof(fission_product_yields_struct) );

  fpydesc.insertMember( "index", HOFFSET(fission_product_yields_struct, index), H5::PredType::NATIVE_INT16);

  fpydesc.insertMember( "from_nuc_name", HOFFSET(fission_product_yields_struct, from_nuc_name), bright::nuc_name_type);
  fpydesc.insertMember( "from_nuc_zz", HOFFSET(fission_product_yields_struct, from_nuc_zz), H5::PredType::NATIVE_INT);

  fpydesc.insertMember( "to_nuc_name", HOFFSET(fission_product_yields_struct, to_nuc_name), bright::nuc_name_type);
  fpydesc.insertMember( "to_nuc_zz", HOFFSET(fission_product_yields_struct, to_nuc_zz), H5::PredType::NATIVE_INT);

  fpydesc.insertMember( "mass_frac", HOFFSET(fission_product_yields_struct, mass_frac), H5::PredType::NATIVE_DOUBLE);
 
  return fpydesc;
};

H5::CompType bright::fission_product_yields_desc = bright::make_fission_product_yields_desc();





H5::CompType bright::make_xs_1g_desc()
{
  //  Makes a decay isotope compound datatype
  H5::CompType xs1gdesc( sizeof(xs_1g_struct) );

  xs1gdesc.insertMember( "nuc_name", HOFFSET(xs_1g_struct, nuc_name), bright::nuc_name_type);
  xs1gdesc.insertMember( "nuc_zz", HOFFSET(xs_1g_struct, nuc_zz), H5::PredType::NATIVE_INT);

  xs1gdesc.insertMember( "sigma_t", HOFFSET(xs_1g_struct, sigma_t), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_s", HOFFSET(xs_1g_struct, sigma_s), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_e", HOFFSET(xs_1g_struct, sigma_e), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_i", HOFFSET(xs_1g_struct, sigma_i), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_a", HOFFSET(xs_1g_struct, sigma_a), H5::PredType::NATIVE_DOUBLE);

  xs1gdesc.insertMember( "sigma_gamma", HOFFSET(xs_1g_struct, sigma_gamma), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_f", HOFFSET(xs_1g_struct, sigma_f), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_alpha", HOFFSET(xs_1g_struct, sigma_alpha), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_proton", HOFFSET(xs_1g_struct, sigma_proton), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_duet", HOFFSET(xs_1g_struct, sigma_duet), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_trit", HOFFSET(xs_1g_struct, sigma_trit), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_2n", HOFFSET(xs_1g_struct, sigma_2n), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_3n", HOFFSET(xs_1g_struct, sigma_3n), H5::PredType::NATIVE_DOUBLE);
  xs1gdesc.insertMember( "sigma_4n", HOFFSET(xs_1g_struct, sigma_4n), H5::PredType::NATIVE_DOUBLE);
 
  return xs1gdesc;
};

H5::CompType bright::xs_1g_desc = bright::make_xs_1g_desc();




/*
 *  Vectorized functions
 */

std::vector<double> bright::delta_vector(double x, std::vector<double> vec)
{
  // This functions finds the 
  // value of (x - vec[i]) for all elements i 
  // in the vector.
  std::vector<double> d (vec.size(), 0.0);

  // Calculate the normalized delta for 
  // all i elements.
  for(int i = 0; i < vec.size(); i++)
    d[i] = (x - vec[i]);

  return d;
};





std::vector<double> bright::normalized_delta(double x, std::vector<double> vec)
{
  // This functions find the normalized 
  // value of (x - vec[i]) for all elements i 
  // in the vector.
  //
  // This is equivelent to the fraction:
  //     (x - vec[i])
  //     ------------
  //      norm_factor
  //
  // Where the normalization factor is 
  //   norm_factor = (vec_max - vec_min) 
  // if the min does not equal the max.
  // and norm_factor = vec_min = vec_max 
  // if it does.

  double norm_factor;
  std::vector<double> nd (vec.size(), 0.0);

  // Get the min and max out of the vector
  double vec_min = *std::min_element(vec.begin(), vec.end());
  double vec_max = *std::max_element(vec.begin(), vec.end());

  if (vec_min == vec_max)
    norm_factor = vec_min;
  else
    norm_factor = vec_max - vec_min;

  // Calculate the normalized delta for 
  // all i elements.
  for(int i = 0; i < vec.size(); i++)
    nd[i] = (x - vec[i]) / norm_factor;

  return nd;
};





bool bright::sorted_index_comparator(std::pair<int, double> i, std::pair<int, double> j)
{
  return i.second < j.second;
};


std::vector<int> bright::sorted_index(std::vector<double> vec)
{
  // Make an indexed vector
  int I = vec.size();
  std::vector< std::pair<int, double> > ind_vec (I);
  for (int i = 0; i < I; i++)
    ind_vec[i] = std::pair<int, double>(i, vec[i]);

  // Sort the indexed vector
  std::sort(ind_vec.begin(), ind_vec.end(), sorted_index_comparator);

  // Grab the indicies out of ind_vec
  std::vector<int> ind (I);
  for (int i = 0; i < I; i++)
    ind[i] = ind_vec[i].first;

  return ind;
};




std::vector<double> bright::y_x_factor_interpolation(double x_factor, std::vector<double> y2, std::vector<double> y1)
{
  // This function calculates the following equation in a vectorized way
  //
  //      y = x(y2 - y1) * x_factor + y1
  //
  // y1 must be of the same size as y2

  int N = y1.size();

  std::vector<double> y (N, -1.0);

  for (int n = 0; n < N; n++)
    y[n] = ((y2[n] - y1[n]) * x_factor) + y1[n];

  return y;
};






std::vector< std::vector<double> > bright::vector_outer_product(std::vector<double> a, std::vector<double> b)
{
  // Performs outer product operation on two vectors
  int I = a.size(); 

  if (I != b.size())
    throw VectorSizeError();

  std::vector< std::vector<double> > c (I, std::vector<double>(I, 0.0)); 

  for (int i = 0; i < I; i++)
  {
    for (int j = 0; j < I; j++)
      c[i][j] = a[i] * b[j];
  };

  return c;
};






std::vector< std::vector<double> > bright::matrix_inverse(std::vector< std::vector<double> > a)
{
  // Performs outer product operation on two vectors
  int I = a.size(); 

  std::vector< std::vector<double> > a_inv (I, std::vector<double>(I, 0.0)); 

  /* This function calculates the inverse of a square matrix
   *
   * Code is rewritten from c++ template code Mike Dinolfo
   * by D. Kroon which was rewritten by Anthony Scopatz
   * which was found at http://snippets.dzone.com/posts/show/7558
   *
   */
  /* Loop variables */
  int i, j, k;

  /* Sum variables */
  double sum, x;
    
  /*  Copy the input matrix to output matrix */
  for (i = 0; i < I; i++) 
  {
    for (j = 0; j < I; j++)
      a_inv[i][j] = a[i][j]; 
  };
    
  /* Add small value to diagonal if diagonal is zero */
  for(i = 0; i < I; i++)
  { 
    if((a_inv[i][i] < 1e-12) && (a_inv[i][i] > -1e-12))
      a_inv[i][i] = 1e-12; 
  }
    
  /* Matrix size of one is special cased */
  if (I == 1)
  {
    a_inv[0][0] = 1.0 / a_inv[0][0];
    return a_inv;
  };

  /* Matrix size must be larger than zero */
  if (I <= 0)
    throw VectorSizeError();

  /* normalize row 0 */
  for (i = 1; i < I; i++) 
    a_inv[0][i] /= a_inv[0][0];

  /* Do LU separation */    
  for (i = 1; i < I; i++)  
  {
    /* do a column of L */
    for (j = i; j < I; j++)  
    { 
      sum = 0.0;
      for (k = 0; k < i; k++) 
        sum += a_inv[j][k] * a_inv[k][i];

      a_inv[j][i] -= sum;
    };

    if (i == I-1)
      continue;

        
    /* do a row of U */
    for (j = i+1; j < I; j++)
    {
      sum = 0.0;
      for (k = 0; k < i; k++)
        sum += a_inv[i][k] * a_inv[k][j];

      a_inv[i][j] = (a_inv[i][j] - sum) / a_inv[i][i];
    };
  };

    /* invert L */ 
    for ( i = 0; i < I; i++ )  
    {
        for ( j = i; j < I; j++ )  
        {
            x = 1.0;

            if ( i != j ) 
            {
                x = 0.0;
                for ( k = i; k < j; k++ ) 
                    x -= a_inv[j][k] * a_inv[k][i];
            };

            a_inv[j][i] = x / a_inv[j][j];
        };
    };

  /* invert U */ 
  for ( i = 0; i < I; i++ ) 
  {
    for ( j = i; j < I; j++ )  
    {
      if ( i == j ) 
        continue;

      sum = 0.0;
      for ( k = i; k < j; k++ )
        sum += a_inv[k][j] * ( (i==k) ? 1.0 : a_inv[i][k] );

      a_inv[i][j] = -sum;
    };
  };

  /* final inversion */ 
  for ( i = 0; i < I; i++ ) 
  {
    for ( j = 0; j < I; j++ )  
    {
      sum = 0.0;

      for ( k = ((i>j)?i:j); k < I; k++ ) 
        sum += ((j==k)?1.0:a_inv[j][k]) * a_inv[k][i];

      a_inv[j][i] = sum;
    };
  };
 
  return a_inv;
};





std::vector< std::vector<double> > bright::matrix_addition(std::vector< std::vector<double> > a, std::vector< std::vector<double> > b)
{
  // Adds two matrices together

  int I = a.size();

  if ( I != a[0].size() || I != b.size() || I != b[0].size())
    throw VectorSizeError();
    
  std::vector< std::vector<double> > c (I, std::vector<double>(I, 0.0)); 

  int i, j;

  for (i = 0; i < I; i++)
  {
    for (j = 0; j < I; j++)
      c[i][j] = a[i][j] + b[i][j];
  };

  return c;
};




std::vector< std::vector<double> > bright::matrix_multiplication(std::vector< std::vector<double> > a, std::vector< std::vector<double> > b)
{
  // Multiplies two matrices together

  int I = a.size();

  if ( I != a[0].size() || I != b.size() || I != b[0].size())
    throw VectorSizeError();
    
  std::vector< std::vector<double> > c (I, std::vector<double>(I, 0.0)); 

  int i, j, k;

  for (i = 0; i < I; i++)
  {
    for (j = 0; j < I; j++)
    {
      for (k = 0; k < I; k++)        
        c[i][j] += a[i][k] * b[k][j];
    };
  };

  return c;
};





std::vector< std::vector<double> > bright::scalar_matrix_product(double a, std::vector< std::vector<double> > M)
{
  // Solves the equation r = aM for a scalar a and Matrix M.
  // Returns the resultant vector r.

  int I = M.size();

  if (I != M[0].size())
    throw VectorSizeError();

  std::vector< std::vector<double> > r (I, std::vector<double>(I, 0.0)); 

  for (int i = 0; i < I; i++)
  {
    for (int j = 0; j < I; j++)
      r[i][j] += (a * M[i][j]);
  };

  return r;
};





std::vector<double> bright::scalar_matrix_vector_product(double a, std::vector< std::vector<double> > M, std::vector<double> v)
{
  // Solves the equation r = aMv for a scalar a, Matrix M, and vector v.
  // Returns the resultant vector r.

  int I = M.size();

  if ( I != M[0].size() || I != v.size())
    throw VectorSizeError();

  std::vector<double> r (I, 0.0);

  for (int i = 0; i < I; i++)
  {
    for (int j = 0; j < I; j++)
      r[i] += (M[i][j] * v[j]);

    r[i] = (r[i] * a);
  };

  return r;
};





/* 
 * Array Helpers
 */

int bright::find_index_char(char * val, char ** arr, int arr_len)
{
  // Finds an element 'val' in array 'arr'
  // returns the index of val's first location
  // returns -1 if not found.
  // For Arrays of char strings

  if (arr_len < 0)
    arr_len = length_array(arr);

  for (int n = 0; n < arr_len; n++)
  {
    if (strcmp(arr[n], val) == 0)
       return n;
  };

  return -1;
};





