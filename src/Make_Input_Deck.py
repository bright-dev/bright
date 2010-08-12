########################
### Custom Libraries ###
########################
import metasci.nuke.Origen as msno

######################
### CHAR Libraries ###
######################
from char import *


###################################
#### MCNP Input Deck Functions ####
###################################

def Make_Burn():
    """Generates the Burn card string for an MCNP inputfile."""

    #make list of time steps
    TimeSteps = []
    for n in CoarseTimeIndex:
        if (0 == n%2) and not (n == 0):
            	TimeSteps.append(CoarseTime[n] - CoarseTime[n-2])

    #Loads the burn card.
    burnline = "BURN TIME=" + str(TimeSteps)[1:-1].replace(',', '') 
    burnline = burnline + " MAT=1 MATVOL=" + "{0:G}".format(FuelCellVolume) 
    burnline = burnline + " POWER={0:G}".format(msn.CellPower(InitialFuelStream.comp, FuelSpecificPower, FuelCellVolume, FuelDensity)) 
    burnline = burnline + " PFRAC="
    for t in TimeSteps:
        burnline = burnline + "1.0 " 
    burnline = burnline + "OMIT=1 27 6014 7016 8018 9018 87223 88227 89228 90234 91229 91230 91232 92229 92230 92231 94245 "
    burnline = burnline + "95240 95642 95644 97245 97246 97247 97248 99249 99250 99251 99252 99253 "
    burnline = burnline + "BOPT=1.0 -04 -1"  

    burnline = msn.Line2MCNP(burnline)

    return burnline

def Make_Mat1():
    """Generates the Material 1 card string for an MCNP inputfile."""

    #writes m1, the initial fuel composition
    m1line = "m1 "
    for key in InitialFuelStream.comp.keys():
        m1line = m1line + '{0}{1} -{2:G} '.format(isoname.zzaaam_2_MCNP( key ), ENDF_FLAG, InitialFuelStream.comp[key])
    for iso in CoreLoad_MCNP:
        if not (iso in isoname.zzaaam_2_MCNP_List(InitialFuelStream.comp.keys()) ) and not (iso in metastableMCNP):
            m1line = m1line + '{0}{1} -{2:G} '.format(iso, ENDF_FLAG, 1.0E-36)
    m1line = msn.Line2MCNP(m1line)

    return m1line

def Make_MatOther():
    """Generates the Other Material card string for an MCNP inputfile."""

    #writes remaining materials
    moline = ""
    for mnum in range(moffset+1, moffset+1+len(mat_number)):
        for iso in mat_number.keys():
            if mat_number[iso] == mnum:
                moline = moline + "m{0} {1}{2} 1.0\n".format(mat_number[iso], iso, ENDF_FLAG) 
                break
    return moline[:-1]

def Make_Pert():
    """Generates the perturbation cards string for an MCNP inputfile."""

    #writes remaining materials
    pline = ""
    pline = pline + "PERT1:N CELL=1 RHO=-{0:G}\n".format(FuelDensity * 0.0)
    pline = pline + "PERT2:N CELL=1 RHO=-{0:G}\n".format(FuelDensity * 0.5)
    pline = pline + "PERT3:N CELL=1 RHO=-{0:G}\n".format(FuelDensity * 1.5)
    pline = pline + "PERT4:N CELL=1 RHO=-{0:G}\n".format(FuelDensity * 2.0)
    return pline[:-1]

def Make_Tally():
    """Generates the Tally card string for an MCNP inputfile."""

    tline = ""

    #Write Group Flux Tally
    tline = tline + "f14:n 1\n"
    tline = tline + "FC14 flux\n"
    tline = tline + "E14 {0} NT\n".format(GroupStructure)
    tline = tline + "sd14 {0:G}\n".format(FuelCellVolume)
    tline = tline + "c \n" 

    #Write Multiplier Tallies
    tline = tline + "f24:n 1\n"
    tline = tline + "FC24 XS\n"
    l2w = "fm24 "
    for mnum in sorted(mat_number.values()):
        for fm in FMdic.values():
            l2w = l2w + "(1 {0:d} {1:d}) ".format(mnum, fm)
    tline = tline + msn.Line2MCNP(l2w) 
    tline = tline + "\n" 
    tline = tline + "E24 {0} NT\n".format(GroupStructure) 
    tline = tline + "sd24 {0:G}\n".format(FuelCellVolume) 
    tline = tline + "c \n"

    #Write Group Flux Tally for CINDER data
    #This is to fix meta-stable read-in flags for ORIGEN
    tline = tline + "f34:n 1\n" 
    tline = tline + "FC34 cinderflux\n" 
    CinderGroupLine = "E34 1.00000E-11 5.00000E-09 1.00000E-08 1.50000E-08 2.00000E-08 2.50000E-08 " + \
                        "3.00000E-08 3.50000E-08 4.20000E-08 5.00000E-08 5.80000E-08 6.70000E-08 " + \
                        "8.00000E-08 1.00000E-07 1.52000E-07 2.51000E-07 4.14000E-07 6.83000E-07 " + \
                        "1.12500E-06 1.85500E-06 3.05900E-06 5.04300E-06 8.31500E-06 1.37100E-05 " + \
                        "2.26000E-05 3.72700E-05 6.14400E-05 1.01300E-04 1.67000E-04 2.75400E-04 " + \
                        "4.54000E-04 7.48500E-04 1.23400E-03 2.03500E-03 2.40400E-03 2.84000E-03 " + \
                        "3.35500E-03 5.53100E-03 9.11900E-03 1.50300E-02 1.98900E-02 2.55400E-02 " + \
                        "4.08700E-02 6.73800E-02 1.11100E-01 1.83200E-01 3.02000E-01 3.88700E-01 " + \
                        "4.97900E-01 6.39279E-01 8.20850E-01 1.10803E+00 1.35335E+00 1.73774E+00 " + \
                        "2.23130E+00 2.86505E+00 3.67879E+00 4.96585E+00 6.06500E+00 1.00000E+01 " + \
                        "1.49182E+01 1.69046E+01 2.00000E+01 2.50000E+01 NT"
    tline = tline + msn.Line2MCNP(CinderGroupLine)
    tline = tline + "\n" 
    tline = tline + "sd34 {0:G}\n".format(FuelCellVolume)
    tline = tline + "c " 

    return tline


def Make_ScatTally():
    """Generates Group-to-Group Scattering XS Tallies."""

    tline = ""

    tline = tline + "F44:n 1\n"
    tline = tline + "FC44 Group-to-Group Scattering XS\n"
    l2w = "FM44 "
    for mnum in sorted(mat_number.values()):
        l2w = l2w + "(1 {0:d} 2:4) ".format(mnum)
    tline = tline + msn.Line2MCNP(l2w) 
    tline = tline + "\n" 
    tline = tline + "E44 {0} NT\n".format(GroupStructure) 
    tline = tline + "SD44 {0:G}\n".format(FuelCellVolume) 
    tline = tline + "FT44 SCX 1 \n"
    tline = tline + "c \n"

    tline = tline + "c " 

    return tline


###################################
### ORIGEN Input Deck Functions ###
###################################

def Make_TAPE4(isovec, n = "TAPE4.INP"):
    """Writes a TAPE4.INP file."""
    msno.writeTAPE4(isovec, n)
    return

def Make_TAPE5(t, p = ""):
    """Writes a TAPE5 files to the current directory + p.
    Only use this function from within the libs/ORIGEN/ directory."""

    #grab important data from library file.
    libfile = tb.openFile('../{0}.h5'.format(reactor), 'r')
    lfr = libfile.root
    n_t = FineTime.index(t)
    IRFtime = FineTime[n_t] - FineTime[n_t-1]
    IRFflux = lfr.Fine.flux[n_t]
    libfile.close()

    #Grab the library numbers used in the original TAPE9 file...
    NLB = []
    tape9_orig = open("../../{0}.tape9".format(reactor), 'r')
    for line in tape9_orig:
        ls = line.split()
        if ls == []:
            continue
        elif (3 < int(ls[0])) and not (ls[0] in NLB):
            NLB.append(ls[0])
    tape9_orig.close()

    #Make template file fill-value dictionary
    tape5_kw = {
        'ORIGEN_Cut_Off': ORIGEN_Cut_Off,
        'IRFtime': '{0:.10E}'.format(IRFtime),
        'IRFflux': '{0:.10E}'.format(IRFflux),
        'NLB1': NLB[0],
        'NLB2': NLB[1],
        'NLB3': NLB[2],
        }

    #Fill the template
    with open("../../{0}.tape5".format(reactor), 'r') as f:
        template_file = f.read()

    with open("{0}{1}_T{2}.tape5".format(p, reactor, t), 'w') as f:
        f.write(template_file.format(**tape5_kw))

    return

def Make_TAPE9(t, p = ""):
    """Makes a new ORIGEN TAPE9.INP file based on original TAPE9 file and stored XS.
    Places the new TAPE9.INP library in the current directory + path p.

    Args:
        * `t` (float): Time at which to evaluate the new library.

    Keyword Args:
    	* `p` (str):  Path to append to the new tape9 file.

    Only use this function from within the libs/ORIGEN/ directory.
    """

    libfile = tb.openFile('../{0}.h5'.format(reactor), 'r')
    lfr = libfile.root
    lfrf = lfr.Fine

    n_t    = FineTime.index(t)
    phi_t  = lfrf.flux[n_t]
    phi_gt = lfrf.flux_g[n_t]

    tape9_kw = {
        'name_new': "{0}{1}_T{2}.tape9".format(p, reactor, t),
        'name_org': "../../{0}.tape9".format(reactor),
        'iso_list': list(lfr.CoreLoad_zzaaam),
        'SNG':   {},
        'SN2N':  {},
        'SN3N':  {},
        'SNF':   {},
        'SNA':   {},
        'SNP':   {},
        'SNGX':  {},
        'SN2NX': {},
        }

    #Get Nomral XS
    for iso in tape9_kw['iso_list']:
        iso_LLAAAM = isoname.zzaaam_2_LLAAAM(iso)
        tape9_kw['SNG'][iso]  = msn.GroupCollapse(getattr(lfrf.sigma_gamma,  iso_LLAAAM)[n_t], phi_gt, phi_t)
        tape9_kw['SN2N'][iso] = msn.GroupCollapse(getattr(lfrf.sigma_2n,     iso_LLAAAM)[n_t], phi_gt, phi_t)
        tape9_kw['SN3N'][iso] = msn.GroupCollapse(getattr(lfrf.sigma_3n,     iso_LLAAAM)[n_t], phi_gt, phi_t)
        tape9_kw['SNF'][iso]  = msn.GroupCollapse(getattr(lfrf.sigma_f,      iso_LLAAAM)[n_t], phi_gt, phi_t)
        tape9_kw['SNA'][iso]  = msn.GroupCollapse(getattr(lfrf.sigma_alpha,  iso_LLAAAM)[n_t], phi_gt, phi_t)
        tape9_kw['SNP'][iso]  = msn.GroupCollapse(getattr(lfrf.sigma_proton, iso_LLAAAM)[n_t], phi_gt, phi_t)

    #(n, g*) XS
    for metaiso_table in lfrf.BranchRatio.NG:
        branch_ratio = metaiso_table[n_t]
        metaiso = isoname.mixed_2_zzaaam(metaiso_table.name)
        initiso = int(metaiso/10)*10 - 10

        if initiso in tape9_kw['iso_list']:
            tape9_kw['SNGX'][initiso] = tape9_kw['SNG'][initiso] * branch_ratio
        
    #(n, 2n*) XS
    for metaiso_table in lfrf.BranchRatio.N2N:
        branch_ratio = metaiso_table[n_t]
        metaiso = isoname.mixed_2_zzaaam(metaiso_table.name)
        initiso = int(metaiso/10)*10 + 10

        if initiso in tape9_kw['iso_list']:
            tape9_kw['SN2NX'][initiso] = tape9_kw['SN2N'][initiso] * branch_ratio

    #Close out the library file
    libfile.close()

    #Finally, make the tape9 file
    msno.writeTAPE9(**tape9_kw)

    return
