// Storage facility class

#include "Storage.h"

/******************************/
/*** Storage Facility Class ***/
/******************************/

/***************************/
/*** Protected Functions ***/
/***************************/

void Storage::initialize ()
{
    char Decay_File[500];
    strcpy(Decay_File, getenv("BRIGHT_DATA") );
    strcat(Decay_File, "/decay.h5");

    //NOTE: The 'decay.h5' librray probably should be rewritten more heirarchically! Changing the following...
    //open the 'decay.h5' file
    hid_t  decay_file_id, decay_dset_id, decay_dspc_id, decay_data_id;
    herr_t decay_status;

    decay_file_id  = H5Fopen(Decay_File, H5F_ACC_RDONLY, H5P_DEFAULT);	//Opens the hdf5 file
    decay_dset_id  = H5Dopen2(decay_file_id, "/Decay", H5P_DEFAULT);		//Opens the Dataset
    decay_dspc_id  = H5Dget_space(decay_dset_id);				//Gets the filespace in order to...
    decay_data_len = H5Sget_simple_extent_npoints(decay_dspc_id);		//Calculate the number of data entries.

    decay_data_id = H5Tcreate(H5T_COMPOUND, sizeof(DecayIso) );	//Maps the file entries to a the data structure.
    decay_status  = H5Tinsert(decay_data_id, "fromiso",     HOFFSET(DecayIso, fromiso),     H5T_STD_I32LE );
    decay_status  = H5Tinsert(decay_data_id, "halflife",    HOFFSET(DecayIso, halflife),    H5T_IEEE_F64LE);
    decay_status  = H5Tinsert(decay_data_id, "decayconst",  HOFFSET(DecayIso, decayconst),  H5T_IEEE_F64LE);
    decay_status  = H5Tinsert(decay_data_id, "toiso",       HOFFSET(DecayIso, toiso),       H5T_STD_I32LE );
    decay_status  = H5Tinsert(decay_data_id, "branchratio", HOFFSET(DecayIso, branchratio), H5T_IEEE_F64LE);

    //Initializes an array of the data struct and fills it!
    decay_data   = new DecayIso [decay_data_len];
    decay_status = H5Dread(decay_dset_id, decay_data_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, decay_data);

    //Closes the hdf5 file.
    decay_status = H5Dclose(decay_dset_id);
    decay_status = H5Fclose(decay_file_id);

    //Put the library in a sexier form...
    //					...for Katy.
    FromIsoStruct * fis;
    for (int i = 0; i < decay_data_len; i++)
    {
        if (0 < decay.count(decay_data[i].fromiso))
        {
            decay[decay_data[i].fromiso].toiso[decay_data[i].toiso] = &(decay_data[i].branchratio);
        }
        else
        {
            fis = new FromIsoStruct;
            (*fis).halflife = &(decay_data[i].halflife);
            (*fis).decayconst = &decay_data[i].decayconst;
            (*fis).toiso[decay_data[i].toiso] = &(decay_data[i].branchratio);
            decay[decay_data[i].fromiso] = (*fis);
        }
    }	
}

void Storage::PrintChain (IsoChain ic)
{
    if (ic.empty())
        return;

    std::cout << "[";
    for (int n = 0; n < ic.size(); n++)
    {
        std::cout << ic[n] << ", ";
    }
    std::cout << "]\n";
    return;
}		

double Storage::bateman (int iso, double mass, IsoChain isochain)
{
    //Solves the Bateman Equations for a isotope and what it decays into.
    double coef  = mass;
    double sumpart = 0.0;
    for (int n = 0; n < isochain.size(); n++)
    {
        if (isochain[n] != iso)
        {
            //Note: that decay[isochain[n]].decayconst = the decay constant, while...
            //...decay[isochain[n]].toiso[isochain[n+1]] represents the branch ratio.
            coef = coef * (*decay[isochain[n]].decayconst) * (*decay[isochain[n]].toiso[isochain[n+1]]);
        }

        double prodpart = 1.0;
        for (int m = 0; m < isochain.size(); m++)
        {
            if (n != m)
                prodpart = prodpart * ((*decay[isochain[m]].decayconst) - (*decay[isochain[n]].decayconst));
        }
        sumpart = sumpart + ((exp(-(*decay[isochain[n]].decayconst) * decay_time))/prodpart);
    }
    return coef * sumpart;
}

void Storage::addchains(IsoChain ic)
{
    int lastiso = ic.back();
    if ( decay[lastiso].toiso.begin()->first == 0)
    {
        return;
    }

    //continue on with next IsoChain if the end of this chain has been reached.
    ToIsoDict toisos = decay[lastiso].toiso;
    for (ToIsoIter tii = toisos.begin(); tii != toisos.end(); tii++)
    {
        IsoChain isochain (ic);
        isochain.push_back(tii->first);
        isochains.insert(isochain);
        addchains(isochain);
    }
    return;
}

void Storage::addchains(int i)
{
    IsoChain ic (1, i);
    isochains.insert(ic);
    addchains(ic);
}

/****************************/
/*** Storage Constructors ***/
/****************************/

Storage::Storage () : FCComp (stor_p2track, "")
{
    //Empty storage component
    initialize();
}

Storage::Storage(std::string n) : FCComp (stor_p2track, n)
{
    initialize();
}

Storage::~Storage ()
{
    //Should close the 'decay.h5' file
}

/************************/
/*** Public Functions ***/
/************************/

void Storage::setParams()
{
    ParamsIn["Mass"]  = IsosIn.mass;
    ParamsOut["Mass"] = IsosOut.mass;
}

MassStream Storage::doCalc()
{
    //Main part of the cooling code.
    //	instream is a mass stream of nuclides as the keys with the mass as a float as the value.
    //	decay_time is a float value for the time in seconds.
    //	FCComps::isos2track throws out any values not in the list before returning vector

    //Initialize the components.
    CompDict cdin, cdout;
    cdin = IsosIn.multByMass();

    //Adds decay chains to isochains set that aren't already there.
    for (CompIter ci = cdin.begin(); ci != cdin.end(); ci++)
    {
        IsoChain ic (1, ci->first);
        if (0 == isochains.count(ic))
        {
            	isochains.insert(ic);
            addchains(ic);
        }
    }

    int mom, daughter;
    for (IsoChainSetIter icsi = isochains.begin(); icsi != isochains.end(); icsi++)
    {
        mom = (*icsi)[0];
        daughter = (*icsi)[(*icsi).size()-1];
        if ( (0 < cdin.count(mom)) && (0 < FCComps::isos2track.count(daughter)) )
        {
            if (0 < cdout.count(daughter))
                cdout[daughter] = cdout[daughter] + bateman(daughter, cdin[mom], *icsi);
            else
                cdout[daughter] = bateman(daughter, cdin[mom], *icsi);
        }
    }
    IsosOut = MassStream (cdout);
    return IsosOut;
}

MassStream Storage::doCalc(CompDict cd)
{
    IsosIn = MassStream (cd);
    return doCalc();
}

MassStream Storage::doCalc(MassStream ms)
{
    IsosIn = ms;
    return doCalc();
}

MassStream Storage::doCalc(double t)
{
    decay_time = t;
    return doCalc();
}

MassStream Storage::doCalc(CompDict cd, double t)
{
    decay_time = t;
    IsosIn = MassStream (cd);
    return doCalc();
}

MassStream Storage::doCalc(MassStream ms, double t)
{
    decay_time = t;
    IsosIn = ms;
    return doCalc();
}
