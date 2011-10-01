// Reprocessing Component Class

#include "Reprocess.h"

/**************************************************/
/*** Reprocessing Component Class and Functions ***/
/**************************************************/

/***************************/
/*** Protected Functions ***/
/***************************/

void Reprocess::initialize(SepEffDict sed)
{
    //Initializes the reprocessing component with specific separation efficiencies.
    //sepeff = dictioanry of separation efficiencies.  Of form {aazzzm: 0.99}, eg {922350, 0.999, 942390: 0.99}
    for (std::set<int>::iterator iso = bright::track_isos.begin(); iso != bright::track_isos.end(); iso++)
    {
        if (0 < sed.count(*iso))
            sepeff[*iso] = sed[*iso];
        else if (0 < sed.count((*iso)/10000))
            sepeff[*iso] = sed[(*iso)/10000];
        else
            sepeff[*iso] = 1.0;
    }
}

/******************************/
/*** Reprocess Constructors ***/
/******************************/

Reprocess::Reprocess () : FCComp (rep_p2track, "")
{
    //Reprocessing Fuel Cycle Component.  Applies Separation Efficiencies.
}

Reprocess::Reprocess (SepEffDict sed, std::string n) : FCComp (rep_p2track, n)
{
    //Reprocessing Fuel Cycle Component.  Applies Separation Efficiencies.
    //sed = dictioanry of separation efficiencies.  Of form {zz: 0.99}, eg {92: 0.999, 94: 0.99} or of form {aazzzm: 0.99}, eg {922350, 0.999, 942390: 0.99}
    initialize(sed);
}

Reprocess::Reprocess (std::map<std::string, double> ssed, std::string n) : FCComp (rep_p2track, n)
{
    //Reprocessing Fuel Cycle Component.  Applies Separation Efficiencies.
    //ssed = string dictioanry of separation efficiencies.  Of form {zz: 0.99}, eg {92: 0.999, 94: 0.99} or of form {aazzzm: 0.99}, eg {922350, 0.999, 942390: 0.99}

    SepEffDict sed;
    for (std::map<std::string, double>::iterator i = ssed.begin(); i != ssed.end(); i++)
    {
        if (0 < pyne::nucname::LLzz.count(i->first))
            sed[pyne::nucname::LLzz[i->first]] = i->second;
        else
        {
            try
            {
                sed[pyne::nucname::zzaaam(i->first)] = i->second;
            }
            catch (std::exception& e)
            {
                continue;
            }
        }
    }
    initialize(sed);
}

Reprocess::~Reprocess ()
{
}


/************************/
/*** Public Functions ***/
/************************/

void Reprocess::calc_params ()
{
    params_prior_calc["Mass"]  = mat_feed.mass;
    params_after_calc["Mass"] = mat_prod.mass;	
}

pyne::Material Reprocess::calc ()
{
    //Does the Reprocessing
    pyne::comp_map incomp  = mat_feed.mult_by_mass();
    pyne::comp_map outcomp;
    for (CompIter i = incomp.begin(); i != incomp.end(); i++)
    {
        outcomp[i->first] = (i->second) * sepeff[i->first];
    }
    mat_prod = pyne::Material (outcomp);
    return mat_prod;
}

pyne::Material Reprocess::calc (pyne::comp_map incomp)
{
    //Does the Reprocessing
    //incomp = input component dictionary of all nuclides. Standard pyne::comp_map object. Assigns this to mat_feed.
    mat_feed = pyne::Material (incomp);
    pyne::comp_map outcomp;
    for (CompIter i = incomp.begin(); i != incomp.end(); i++)
    {
        outcomp[i->first] = (i->second) * sepeff[i->first];
    }
    mat_prod = pyne::Material (outcomp);
    return mat_prod;
}

pyne::Material Reprocess::calc (pyne::Material instream)
{
    //Does the Reprocessing
    //instream = input stream of all nuclides. Standard pyne::Material object.
    mat_feed = instream;
    return calc();
}
