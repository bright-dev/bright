// Storage facility class

#include "storage.h"

/******************************/
/*** Storage Facility Class ***/
/******************************/

/***************************/
/*** Protected Functions ***/
/***************************/

void bright::Storage::initialize ()
{
  char decay_file[500];
  strcpy(decay_file, getenv("BRIGHT_DATA") );
  #ifdef _WIN32
    strcat(decay_file, "\\decay.h5");
  #else
    strcat(decay_file, "/decay.h5");
  #endif

  // NOTE: The 'decay.h5' librray probably should be rewritten more heirarchically! Changing the following...
  // open the 'decay.h5' file
  hid_t  decay_file_id, decay_dset_id, decay_dspc_id, decay_data_id;
  herr_t decay_status;

  decay_file_id  = H5Fopen(decay_file, H5F_ACC_RDONLY, H5P_DEFAULT);  // Opens the hdf5 file
  decay_dset_id  = H5Dopen2(decay_file_id, "/Decay", H5P_DEFAULT);    // Opens the Dataset
  decay_dspc_id  = H5Dget_space(decay_dset_id);	                      // Gets the filespace in order to...
  decay_data_len = H5Sget_simple_extent_npoints(decay_dspc_id);       // Calculate the number of data entries.

  decay_data_id = H5Tcreate(H5T_COMPOUND, sizeof(decay_nuc) );	// Maps the file entries to a the data structure.
  decay_status  = H5Tinsert(decay_data_id, "fromiso",     HOFFSET(decay_nuc, fromiso),     H5T_STD_I32LE );
  decay_status  = H5Tinsert(decay_data_id, "halflife",    HOFFSET(decay_nuc, halflife),    H5T_IEEE_F64LE);
  decay_status  = H5Tinsert(decay_data_id, "decayconst",  HOFFSET(decay_nuc, decayconst),  H5T_IEEE_F64LE);
  decay_status  = H5Tinsert(decay_data_id, "toiso",       HOFFSET(decay_nuc, toiso),       H5T_STD_I32LE );
  decay_status  = H5Tinsert(decay_data_id, "branchratio", HOFFSET(decay_nuc, branchratio), H5T_IEEE_F64LE);

  // Initializes an array of the data struct and fills it!
  decay_data   = new decay_nuc [decay_data_len];
  decay_status = H5Dread(decay_dset_id, decay_data_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, decay_data);

  // Closes the hdf5 file.
  decay_status = H5Dclose(decay_dset_id);
  decay_status = H5Fclose(decay_file_id);

  // Put the library in a sexier form...
  //					...for Katy.
  from_nuc_struct * fis;
  for (int i = 0; i < decay_data_len; i++)
  {
    if (0 < decay.count(decay_data[i].fromiso))
      decay[decay_data[i].fromiso].toiso[decay_data[i].toiso] = &(decay_data[i].branchratio);
    else
    {
      fis = new from_nuc_struct;
      (*fis).halflife = &(decay_data[i].halflife);
      (*fis).decayconst = &decay_data[i].decayconst;
      (*fis).toiso[decay_data[i].toiso] = &(decay_data[i].branchratio);
      decay[decay_data[i].fromiso] = (*fis);
    };
  };	
};


void bright::Storage::print_chain (nuc_chain nc)
{
  if (nc.empty())
    return;

  std::cout << "[";
  for (int n = 0; n < nc.size(); n++)
    std::cout << nc[n] << ", ";

  std::cout << "]\n";
  return;
};	


double bright::Storage::bateman (int iso, double mass, nuc_chain nucchain)
{
  // Solves the Bateman Equations for a isotope and what it decays into.
  double coef  = mass;
  double sumpart = 0.0;
  for (int n = 0; n < nucchain.size(); n++)
  {
    if (nucchain[n] != iso)
    {
      // Note: that decay[isochain[n]].decayconst = the decay constant, while...
      // ...decay[isochain[n]].toiso[isochain[n+1]] represents the branch ratio.
      coef = coef * (*decay[nucchain[n]].decayconst) * (*decay[nucchain[n]].toiso[nucchain[n+1]]);
    };

    double prodpart = 1.0;
    for (int m = 0; m < nucchain.size(); m++)
    {
      if (n != m)
        prodpart = prodpart * ((*decay[nucchain[m]].decayconst) - (*decay[nucchain[n]].decayconst));
    }
    sumpart = sumpart + ((exp(-(*decay[nucchain[n]].decayconst) * decay_time))/prodpart);
  };
  return coef * sumpart;
};


void bright::Storage::addchains(nuc_chain nc)
{
  int lastnuc = nc.back();
  if (decay[lastnuc].toiso.begin()->first == 0)
    return;

  // continue on with next IsoChain if the end of this chain has been reached.
  to_nuc_dict tonucs = decay[lastnuc].toiso;
  for (to_nuc_iter tni = tonucs.begin(); tni != tonucs.end(); tni++)
  {
    nuc_chain nucchain (nc);
    nucchain.push_back(tni->first);
    nucchains.insert(nucchain);
    addchains(nucchain);
  }
  return;
};


void bright::Storage::addchains(int i)
{
  nuc_chain nc (1, i);
  nucchains.insert(nc);
  addchains(nc);
};



/****************************/
/*** Storage Constructors ***/
/****************************/

bright::Storage::Storage () : bright::FCComp(stor_p2track, "")
{
  // Empty storage component
  initialize();
};


bright::Storage::Storage(std::string n) : bright::FCComp (stor_p2track, n)
{
  initialize();
}


bright::Storage::~Storage ()
{
  //Should close the 'decay.h5' file
}


/************************/
/*** Public Functions ***/
/************************/

void bright::Storage::calc_params()
{
  params_prior_calc["Mass"]  = mat_feed.mass;
  params_after_calc["Mass"] = mat_prod.mass;
}


pyne::Material bright::Storage::calc()
{
  // Main part of the cooling code.
  // mat is a mass stream of nuclides as the keys with the mass as a float as the value.
  // decay_time is a float value for the time in seconds.
  // bright::track_nucs throws out any values not in the list before returning vector

  // Initialize the components.
  pyne::comp_map cdin, cdout;
  cdin = mat_feed.mult_by_mass();

  // Adds decay chains to isochains set that aren't already there.
  for (pyne::comp_iter ci = cdin.begin(); ci != cdin.end(); ci++)
  {
    nuc_chain nc (1, ci->first);
    if (0 == nucchains.count(nc))
    {
      nucchains.insert(nc);
      addchains(nc);
    };
  };

  int mom, daughter;
  for (nuc_chain_set_iter ncsi = nucchains.begin(); ncsi != nucchains.end(); ncsi++)
  {
    mom = (*ncsi)[0];
    daughter = (*ncsi)[(*ncsi).size()-1];
    if ( (0 < cdin.count(mom)) && (0 < bright::track_nucs.count(daughter)) )
    {
      if (0 < cdout.count(daughter))
        cdout[daughter] = cdout[daughter] + bateman(daughter, cdin[mom], *ncsi);
      else
        cdout[daughter] = bateman(daughter, cdin[mom], *ncsi);
    };
  };

  mat_prod = pyne::Material (cdout);
  return mat_prod;
};


pyne::Material bright::Storage::calc(pyne::comp_map cd)
{
  mat_feed = pyne::Material (cd);
  return calc();
}

pyne::Material bright::Storage::calc(pyne::Material mat)
{
  mat_feed = mat;
  return calc();
}


pyne::Material bright::Storage::calc(double t)
{
    decay_time = t;
    return calc();
}


pyne::Material bright::Storage::calc(pyne::comp_map cd, double t)
{
    decay_time = t;
    mat_feed = pyne::Material (cd);
    return calc();
}


pyne::Material bright::Storage::calc(pyne::Material mat, double t)
{
    decay_time = t;
    mat_feed = mat;
    return calc();
}
