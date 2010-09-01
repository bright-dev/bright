########################
#### CHAR Libraries ####
########################
from char import *

def Write_ORIGEN_Libs():
    """Writes MCNP output to Coarse ORIGEN libraries."""
    global InFlux, InXS, InBurnSum, AtMult, AtBin, LibFileOpen, libfile, mat, fm, n, flux_g, flux, E_up, mult, burn
    global InFluxCinder, E_up_cinder, flux_g_cinder, flux_cinder
    global BR_NG, BR_N2N

    if mult == {}:
        Parse_MCNP()

    if not ( 'libs' in os.listdir('.') ):
        os.mkdir('libs/')
    os.chdir('libs/')
    if ( 'ORIGEN' in os.listdir('.') ):
        metasci.SafeRemove('ORIGEN', IsDir=True)
        os.removedirs('ORIGEN')
    os.mkdir('ORIGEN/')
    os.chdir('ORIGEN/')

    for t in CoarseTime:
        Make_Input_Deck.Make_TAPE9(t)

    os.chdir('../../')
    return

def Parse_TAPE6(p = ""):
    "Parses an ORIGEN TAPE6.OUT file that is in the current directory + path p."
    InTable5 = False
    outvec = {}

    tape6 = open("%sTAPE6.OUT"%p, 'r')
    for line in tape6:
        if "BURNUP,MWD" in line:
            ls = line.split()
            BU = float(ls[-1])
            continue
        elif "K INFINITY" in line:
            ls = line.split()
            k = float(ls[-1])
            continue
        elif "NEUT PRODN" in line:
            ls = line.split()
            Pro = float(ls[-1])
            continue
        elif "NEUT DESTN" in line:
            ls = line.split()
            Des = float(ls[-1])
            continue
        elif "5 SUMMARY TABLE:  CONCENTRATIONS, GRAMS" in line:
            InTable5 = True
            continue
        elif InTable5 and ("OUTPUT UNIT =  6" in line):
            InTable5 = False
            continue
        elif InTable5:
            ls = line.split()
            try:
                iso = isoname.LLAAAM_2_zzaaam(ls[0])
            except:
                try:
                    iso = isoname.LLAAAM_2_zzaaam(ls[0] + ls[1])
                except:
                    continue
            outvec[iso] = float(ls[-1])
        else:
            continue
    tape6.close()
    return BU, k, Pro, Des, outvec

def Write_TXT_Lib_ORIGEN(BU, k, Pro, Des, Tij):
    """Writes ORIGEN output to text libraries.
    Not Really helpful.
    Inputs are dictionaries of tables."""

    os.chdir('libs/ORIGEN/')

    G = len(BU[BU.keys()[0]][0])

    for iso in BU.keys():
        for n_g in range(G):
            os.chdir('E%s'%(G - n_g) )
            os.chdir('%s'%iso)

            #Write Burnup
            BUfile = open('Burnup.lib', 'w')
            BUfile.write("Time\tBurnup\n")
            for n_t in FineTimeIndex:
                BUfile.write('%s\t%.6E\n'%(FineTime[n_t], BU[iso][n_t][n_g]) )
            BUfile.close()

            #Write Multiplication factor
            kfile = open('k.lib', 'w')
            kfile.write("Time\tk\n")
            for n_t in FineTimeIndex:
                kfile.write('%s\t%.6E\n'%(FineTime[n_t], k[iso][n_t][n_g]) )
            kfile.close()

            #Write Neutron Production Rate 
            Profile = open('Pro.lib', 'w')
            Profile.write("Time\tPro\n")
            for n_t in FineTimeIndex:
                Profile.write('%s\t%.6E\n'%(FineTime[n_t], Pro[iso][n_t][n_g]) )
            Profile.close()

            #Write Neutron Destruction Rate 
            Desfile = open('Des.lib', 'w')
            Desfile.write("Time\tDes\n")
            for n_t in FineTimeIndex:
                Desfile.write('%s\t%.6E\n'%(FineTime[n_t], Des[iso][n_t][n_g]) )
            Desfile.close()

            #Write Transformation Matrix
            Tijfile = open('Tij.lib', 'w')
            Tijfile.write("Time")
            for n_t in FineTimeIndex:
                Tijfile.write('\t%s\t'%FineTime[n_t])
            Tijfile.write('\n')
            for j in CoreTran_zzaaam:
                if not (j in Tij[iso][n_t][n_g].keys() ):
                    continue
                Tijfile.write(CoreTran_LLAAAM[CoreTran_zzaaam.index(j)])
                for n_t in FineTimeIndex:
                    Tijfile.write('\t%.6E'%Tij[iso][n_t][n_g][j]) 
                Tijfile.write('\n')
            Tijfile.close()

            os.chdir('../') #Returns to the Energy Directory
            os.chdir('../') #Returns to the ORIGEN Directory

    os.chdir('../../') #Return to 'reactor' root
    return 


def Write_HDF5_Lib_ORIGEN(BU, k, Pro, Des, Tij):
    """Writes ORIGEN output to the HDF5 library."""
    os.chdir('libs/')
    libfile = tb.openFile(reactor + ".h5", mode = "r+")
    lfrf = libfile.root.Fine

    #Write Data
    libfile.createGroup("/Fine", "Burnup")
    libfile.createGroup("/Fine", "k")
    libfile.createGroup("/Fine", "Production")
    libfile.createGroup("/Fine", "Destruction")
    libfile.createGroup("/Fine", "Transmutation")
    for iso in BU.keys(): 
        isoLL = isoname.zzaaam_2_LLAAAM(iso)

        #Writes the easy data
        libfile.createArray("/Fine/Burnup",      isoLL, BU[iso],  "Burnup BU(F) for {0}".format(isoLL))                     #Burnup Data
        libfile.createArray("/Fine/k",           isoLL, k[iso],   "Multiplication Factor k(F) for {0}".format(isoLL))       #Multiplication Factor Data
        libfile.createArray("/Fine/Production",  isoLL, Pro[iso], "Neutron Production Rate p(F) for {0}".format(isoLL))     #Production Rate Data
        libfile.createArray("/Fine/Destruction", isoLL, Des[iso], "Neutron Destruction Rate d(F) for {0}".format(isoLL))    #Destruction Rate Data

        #Writes the transmutation matrices
        libfile.createGroup("/Fine/Transmutation", isoLL)
        for jso in CoreTran_zzaaam:
            if (jso not in Tij[iso][0].keys()):
                continue
            jsoLL = isoname.zzaaam_2_LLAAAM(jso)
    
            M = []
            for n_t in FineTimeIndex:
                M.append(Tij[iso][n_t][jso])
            #Transmutation Matrix
            libfile.createArray("/Fine/Transmutation/{0}".format(isoLL), jsoLL, M, 
                "Transmutation Matrix Tij(F) from {0} to {1}".format(isoLL, jsoLL) )	

    libfile.close()
    os.chdir('../') #Returns to 'reactor' root directory
    return
