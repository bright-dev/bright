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
    for (std::set<int>::iterator iso = FCComps::track_isos.begin(); iso != FCComps::track_isos.end(); iso++)
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
        if (0 < isoname::LLzz.count(i->first))
            sed[isoname::LLzz[i->first]] = i->second;
        else
        {
            try
            {
                sed[isoname::mixed_2_zzaaam(i->first)] = i->second;
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
    params_prior_calc["Mass"]  = ms_feed.mass;
    params_after_calc["Mass"] = ms_prod.mass;	
}

MassStream Reprocess::calc ()
{
    //Does the Reprocessing
    CompDict incomp  = ms_feed.mult_by_mass();
    CompDict outcomp;
    for (CompIter i = incomp.begin(); i != incomp.end(); i++)
    {
        outcomp[i->first] = (i->second) * sepeff[i->first];
    }
    ms_prod = MassStream (outcomp);
    return ms_prod;
}

MassStream Reprocess::calc (CompDict incomp)
{
    //Does the Reprocessing
    //incomp = input component dictionary of all nuclides. Standard CompDict object. Assigns this to ms_feed.
    ms_feed = MassStream (incomp);
    CompDict outcomp;
    for (CompIter i = incomp.begin(); i != incomp.end(); i++)
    {
        outcomp[i->first] = (i->second) * sepeff[i->first];
    }
    ms_prod = MassStream (outcomp);
    return ms_prod;
}

MassStream Reprocess::calc (MassStream instream)
{
    //Does the Reprocessing
    //instream = input stream of all nuclides. Standard MassStream object.
    ms_feed = instream;
    return calc();
}
