########################
### Custom Libraries ###
########################
import metasci.nuke.Origen as msno

######################
### CHAR Libraries ###
######################
from char import *


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
