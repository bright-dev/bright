########################
#### CHAR Libraries ####
########################
from char import *
import Make_Input_Deck 


############################
### MCNP Parse Functions ###
############################
InFlux = False
InXS = False
InBurnSum = False
AtMult = False
AtBin = False
LibFileOpen = False
libfile = None
mat = None
fm = None
n = 0
E_up = []
flux_g = []
flux = []
mult = {}
burn = {"step": [], "duration": [], "time": [], "power": [], "keff": [], "flux": [], "ave. nu": [], "ave. q": [], "burnup": [], "source": []}

#Cinder Flux Data
InFluxCinder = False
E_up_cinder = []
flux_g_cinder = []
flux_cinder = []

#Meta-Stable to Ground State Branch Ratios
BR_NG  = {}
BR_N2N = {}

def Init_MCNP():
    """Initializes variables for MCNP parsing"""
    global InFlux, InXS, InBurnSum, AtMult, AtBin, LibFileOpen, libfile, mat, fm, n, flux_g, flux, E_up, mult, burn
    global InFluxCinder, E_up_cinder, flux_g_cinder, flux_cinder
    global BR_NG, BR_N2N

    InFlux = False
    InXS = False
    InBurnSum = False
    AtMult = False
    AtBin = False
    LibFileOpen = False
    libfile = None
    mat = None
    fm = None
    n = 0
    flux_g = []
    E_up = []
    mult = {}
    burn = {"step": [], "duration": [], "time": [], "power": [], "keff": [], "flux": [], "ave. nu": [], "ave. q": [], "burnup": [], "source": []}

    InFluxCinder = False
    E_up_cinder = []
    flux_g_cinder = []
    flux_cinder = []

    BR_NG  = {}
    BR_N2N = {}

    return

def Parse_MCNP():
    """Parse MCNPX Output"""
    global InFlux, InXS, InBurnSum, AtMult, AtBin, LibFileOpen, libfile, mat, fm, n, flux_g, flux, E_up, mult, burn
    global InFluxCinder, E_up_cinder, flux_g_cinder, flux_cinder
    global BR_NG, BR_N2N

    if not (mult == {}):
        Init_MCNP()

    outfile = open(reactor + '.o', 'r')

    for line in outfile:
        ls = line.split()
        lslen = len(ls)

        if lslen == 0:
            if (InFlux or InFluxCinder) and AtBin:
                InFlux = False
                InFluxCinder = False
                AtBin = False
                continue
            elif AtMult and AtBin:
                AtMult = False
                AtBin = False
                mat = None
                fm = None
                n = 0
                continue
            elif InBurnSum:
                InBurnSum = False
            else:
                continue

        #Get Flux Data
        elif ls[0] == 'energy' and (InFlux or InFluxCinder):
            AtBin = True
            continue
        elif (1 < lslen) and (ls[1] == 'flux'):
            InFlux = True
            flux_g.append([])
            continue
        elif (1 < lslen) and (ls[1] == 'cinderflux'):
            InFluxCinder = True
            flux_g_cinder.append([])
            continue
        elif InFlux and AtBin:
            if len(flux_g) == 1:
                E_up.append(float(ls[0]))
            flux_g[-1].append(float(ls[1]))
            continue
        elif InFluxCinder and AtBin:
            if len(flux_g_cinder) == 1:
                E_up_cinder.append(float(ls[0]))
            flux_g_cinder[-1].append(float(ls[1]))
            continue

        #Get Multiplier Data
        elif (1 < lslen) and (ls[1] == 'XS'):
            InXS = True
            continue
        elif InXS and (ls[0] == 'multiplier'):
            AtMult = True
            mat = isoname.MCNP_2_LLAAAM( number_mat[ ls[-2] ] )
            fm = dicFM[ ls[-1] ]
            lib = '%s_%s.lib'%(mat, fm)
            if lib in mult.keys():
                mult[lib].append([])
            else:
                mult[lib] = [[]]
            continue
        elif ls[0] == 'energy' and AtMult:
            AtBin = True
            continue
        elif AtMult and AtBin:
            mult[lib][-1].append( float(ls[1]) / flux_g[-1][n] )
            n = n + 1

        #Get Burnup Summary Table
        elif ls == ["(days)", "(days)", "(MW)", "(GWd/MTU)", "(nts/sec)"]:
            InBurnSum = True
            continue
        elif InBurnSum and (0 == len(burn['step'])):
            burn["step"].append( float(ls[0]) )
            burn["duration"].append( float(ls[1]) )
            burn["time"].append( float(ls[2]) )
            burn["power"].append( float(ls[3]) )
            burn["keff"].append( float(ls[4]) )
            burn["flux"].append( float(ls[5]) )
            burn["ave. nu"].append( float(ls[6]) )
            burn["ave. q"].append( float(ls[7]) )
            burn["burnup"].append( float(ls[8]) )
            burn["source"].append( float(ls[9]) )
        elif InBurnSum and (float(ls[0]) <= burn['step'][-1]):
            continue
        elif InBurnSum and (burn['step'][-1] < float(ls[0])):
            #Half Step Data Interpolation
            burn["step"].append( (burn["step"][-1] + float(ls[0])) / 2.0)
            burn["duration"].append( (burn["duration"][-1] + float(ls[1])) / 2.0)
            burn["time"].append( (burn["time"][-1] + float(ls[2])) / 2.0)
            burn["power"].append( (burn["power"][-1] + float(ls[3])) /2.0)
            burn["keff"].append( (burn["keff"][-1] + float(ls[4])) / 2.0)
            burn["flux"].append( (burn["flux"][-1] + float(ls[5])) / 2.0)
            burn["ave. nu"].append( (burn["ave. nu"][-1] + float(ls[6])) / 2.0)
            burn["ave. q"].append( (burn["ave. q"][-1] + float(ls[7])) / 2.0)
            burn["burnup"].append( (burn["burnup"][-1] + float(ls[8])) / 2.0)
            burn["source"].append( (burn["source"][-1] + float(ls[9])) / 2.0)

            #Full Step Data
            burn["step"].append( float(ls[0]) )
            burn["duration"].append( float(ls[1]) )
            burn["time"].append( float(ls[2]) )
            burn["power"].append( float(ls[3]) )
            burn["keff"].append( float(ls[4]) )
            burn["flux"].append( float(ls[5]) )
            burn["ave. nu"].append( float(ls[6]) )
            burn["ave. q"].append( float(ls[7]) )
            burn["burnup"].append( float(ls[8]) )
            burn["source"].append( float(ls[9]) )
            continue

    outfile.close()

    #Calculate Total Fluxes	
    #User-specified groups
    for gth_flux_spectrum in flux_g:
        f = 0.0
        for f_g in gth_flux_spectrum:
            f = f + f_g
        flux.append(f)
    #Cinder Energy Groups
    for gth_flux_spectrum in flux_g_cinder:
        f = 0.0
        for f_g in gth_flux_spectrum:
            f = f + f_g
        flux_cinder.append(f)

    #Renormalize Fluxes if burned...
    if 0 < len( burn["step"] ):
        #User-specified groups
        for t in range(len(flux_g)):
            for g in range(len(flux_g[t])):
                flux_g[t][g] = flux_g[t][g] * burn["flux"][t] / flux[t]
        flux = burn["flux"]
        #Cinder Energy Groups
        for t in range(len(flux_g_cinder)):
            for g in range(len(flux_g_cinder[t])):
                flux_g_cinder[t][g] = flux_g_cinder[t][g] * burn["flux"][t] / flux_cinder[t]
        flux_cinder = burn["flux"]

    #Calculate Meta-Stable to Ground State branch ratios, if metastables are suppossed to be tracked...
    for iso in metastabletrak:
        anum = (iso%10000)/10
        znum = iso/10000
        NG_grd_str  = "{0}{1:03d}0 produced by the following C-X  (g   )".format(anum, znum)
        NG_mss_str  = "{0}{1:03d}1 produced by the following C-X  (g  *)".format(anum, znum)
        N2N_grd_str = "{0}{1:03d}0 produced by the following C-X  (2n  )".format(anum, znum)
        N2N_mss_str = "{0}{1:03d}1 produced by the following C-X  (2n *)".format(anum, znum)

        InGrndNG  = False
        InMetaNG  = False
        InGrndN2N = False
        InMetaN2N = False

        GrndXS_NG  = []
        MetaXS_NG  = []
        GrndXS_N2N = []
        MetaXS_N2N = []

        cinderfile = open(CINDER_DAT, 'r')
        for line in cinderfile:
            if len(line) <= 4:
                continue
            elif NG_grd_str in line:
                InGrndNG = True
                InMetaNG = False
                continue
            elif NG_mss_str in line:
                InGrndNG = False
                InMetaNG = True
                continue
            elif N2N_grd_str in line:
                InGrndN2N = True
                InMetaN2N = False
                continue
            elif N2N_mss_str in line:
                InGrndN2N = False
                InMetaN2N = True
                continue
            elif "_______________________" in line:
                InGrndNG  = False
                InMetaNG  = False
                InGrndN2N = False
                InMetaN2N = False
                continue
            elif "   #" == line[:4]:
                InGrndNG  = False
                InMetaNG  = False
                InGrndN2N = False
                InMetaN2N = False
                continue
            elif InGrndNG:
                for xs in line.split():
                    GrndXS_NG.append( float(xs) )
            elif InMetaNG:
                for xs in line.split():
                    MetaXS_NG.append( float(xs) )
            elif InGrndN2N:
                for xs in line.split():
                    GrndXS_N2N.append( float(xs) )
            elif InMetaN2N:
                for xs in line.split():
                    MetaXS_N2N.append( float(xs) )
        cinderfile.close()

        if not (GrndXS_NG == []) and not (MetaXS_NG == []):
            BR_NG[iso] = []
            for t in range( len(flux_g_cinder) ):
                gXS_NG = msn.GroupCollapse(GrndXS_NG, flux_g_cinder[t], flux_cinder[t])
                mXS_NG = msn.GroupCollapse(MetaXS_NG, flux_g_cinder[t], flux_cinder[t])
                BR_NG[iso].append(mXS_NG / gXS_NG)

        if not (GrndXS_N2N == []) and not (MetaXS_N2N == []):
            BR_N2N[iso] = []
            for t in range( len(flux_g_cinder) ):
                gXS_N2N = msn.GroupCollapse(GrndXS_N2N, flux_g_cinder[t], flux_cinder[t])
                mXS_N2N = msn.GroupCollapse(MetaXS_N2N, flux_g_cinder[t], flux_cinder[t])
                BR_N2N[iso].append(mXS_N2N / gXS_N2N)
    return

def Write_TXT_Lib_MCNP():
    """Writes MCNP output to text libraries."""
    global InFlux, InXS, InBurnSum, AtMult, AtBin, LibFileOpen, libfile, mat, fm, n, flux_g, flux, E_up, mult, burn
    global InFluxCinder, E_up_cinder, flux_g_cinder, flux_cinder
    global BR_NG, BR_N2N
    
    if mult == {}:
        Parse_MCNP()

    if not ( 'libs' in os.listdir('.') ):
        os.mkdir('libs/')
    os.chdir('libs/')
    for lib in os.listdir('.'):
        if lib[-4:] == ".lib":
            metasci.SafeRemove(lib)

    headline = "E_up"
    for t in CoarseTime:
        headline = headline + "\t\t" + str(t)
    headline = headline + "\n"

    ### Write Flux Files
    #Group Fluxes
    libfile = open('flux_g.lib', 'w')
    libfile.write(headline)
    for e in range(len(E_up)):
        line = '%.5E'%E_up[e]
        for f in range(len(flux_g)):
            line = line + "\t%.6E"%flux_g[f][e]
        libfile.write('%s\n'%line)
    libfile.close()	 

    #Total Flux
    libfile = open('flux.lib', 'w')
    libfile.write('Time\t' + headline[5:])
    line = "Flux"
    for f in range(len(flux)):
        line = line + "\t%.6E"%flux[f]
    libfile.write('%s\n'%line)
    libfile.close()	 

    #Cinder Group Fluxes
    libfile = open('flux_g_cinder.lib', 'w')
    libfile.write(headline)
    for e in range(len(E_up_cinder)):
        line = '%.5E'%E_up_cinder[e]
        for f in range(len(flux_g_cinder)):
            line = line + "\t%.6E"%flux_g_cinder[f][e]
        libfile.write('%s\n'%line)
    libfile.close()	 

    #Cinder Total Flux
    libfile = open('flux_cinder.lib', 'w')
    libfile.write('Time\t' + headline[5:])
    line = "Flux"
    for f in range(len(flux_cinder)):
        line = line + "\t%.6E"%flux_cinder[f]
    libfile.write('%s\n'%line)
    libfile.close()	 

    #Write Multiplier Files
    for lib in mult.keys():
        libfile = open(lib, 'w')
        libfile.write(headline)
        for e in range(len(E_up)):
            line = '%.5E'%E_up[e]
            for ml in range(len(mult[lib])):
                line = line + "\t%.6E"%mult[lib][ml][e]
            libfile.write('%s\n'%line)
        libfile.close()

    #Meta-Stable Branch Ratio Files
    for iso in BR_NG.keys():
        libfile = open(isoname.zzaaam_2_LLAAAM(iso) + '_NGamma_Branch_Ratio.lib', 'w')
        libfile.write('Time\t' + headline[5:])
        line = "Ratio"
        for br in range(len(BR_NG[iso])):
            line = line + "\t%.6E"%BR_NG[iso][br]
        libfile.write('%s\n'%line)
        libfile.close()	 
    for iso in BR_N2N.keys():
        libfile = open(isoname.zzaaam_2_LLAAAM(iso) + '_N2N_Branch_Ratio.lib', 'w')
        libfile.write('Time\t' + headline[5:])
        line = "Ratio"
        for br in range(len(BR_N2N[iso])):
            line = line + "\t%.6E"%BR_N2N[iso][br]
        libfile.write('%s\n'%line)
        libfile.close()	 

    os.chdir('..')
    return 

def Write_HDF5_Lib_MCNP():
    """Writes MCNP output to an HDF5 library."""
    global InFlux, InXS, InBurnSum, AtMult, AtBin, LibFileOpen, libfile, mat, fm, n, flux_g, flux, E_up, mult, burn
    global InFluxCinder, E_up_cinder, flux_g_cinder, flux_cinder
    global BR_NG, BR_N2N
    
    if mult == {}:
        Parse_MCNP()

    if not ( 'libs' in os.listdir('.') ):
        os.mkdir('libs/')
    os.chdir('libs/')
    for h5 in os.listdir('.'):
        if h5[-3:] == ".h5":
            metasci.SafeRemove(h5)

    # Initialize the HDF5 file
    libfile = tb.openFile(reactor + ".h5", mode = "w", title = '[CHAR] %s'%reactor)
    root = libfile.root

    ##########################################
    ### Make arrays that apply to all data ###
    ##########################################
    #Add Energy Group Array
    libfile.createArray(root, "E_up", E_up, "Upper Energy Limit [MeV]") 					#Upper Energy Limit
    libfile.createArray(root, "CoreLoad_zzaaam", CoreLoad_zzaaam, "Core Loading Isotopes [zzaaam]")		#Core load isotopes in zzaaam form
    libfile.createArray(root, "CoreLoad_LLAAAM", CoreLoad_LLAAAM, "Core Loading Isotopes [LLAAAM]")		#Core load isotopes in LLAAAM form
    libfile.createArray(root, "CoreLoad_MCNP",   CoreLoad_MCNP,   "Core Loading Isotopes [MCNP]")		#Core load isotopes in MCNP form
    libfile.createArray(root, "CoreTran_zzaaam", CoreTran_zzaaam, "Core Transformation Isotopes [zzaaam]")	#Core transformation isotopes in zzaaam form
    libfile.createArray(root, "CoreTran_LLAAAM", CoreTran_LLAAAM, "Core Transformation Isotopes [LLAAAM]")	#Core transformation isotopes in LLAAAM form
    libfile.createArray(root, "CoreTran_MCNP",   CoreTran_MCNP,   "Core Transformation Isotopes [MCNP]")	#Core transformation isotopes in MCNP form

    #Add Coarsely binned time and flux data
    libfile.createGroup(root, "Coarse")
    libfile.createArray("/Coarse", "time",   CoarseTime, "Burn Time, Coarse [days]")		#Burn up times
    libfile.createArray("/Coarse", "flux",   flux,       "Neutron Flux, Coarse [n/s/cm2]")	 	#Flux
    libfile.createArray("/Coarse", "flux_g", flux_g,     "Neutron Group Flux, Coarse [n/s/cm2]")	#Group fluxes

    #Add Finely Binned time and flux data
    libfile.createGroup(root, "Fine")
    flux_fine = []
    flux_g_fine = []
    for t in FineTime:
        n = 0
        while (CoarseTime[n] <= t) and not (t <= CoarseTime[n+1]):
                n = n + 1
        flux_fine.append( metasci.SolveLine(t, CoarseTime[n+1], flux[n+1], CoarseTime[n], flux[n]) )
        flux_g_fine.append([])
        for e in range( len(E_up) ):
            flux_g_fine[-1].append( metasci.SolveLine(t, CoarseTime[n+1], flux_g[n+1][e], CoarseTime[n], flux_g[n][e]) )
    libfile.createArray("/Fine", "time",   FineTime,    "Burn Time, Fine [days]")			#Burn up times
    libfile.createArray("/Fine", "flux",   flux_fine,   "Neutron Flux, Fine [n/s/cm2]") 		#Flux
    libfile.createArray("/Fine", "flux_g", flux_g_fine, "Neutron Group Flux, Fine [n/s/cm2]")	#Group fluxes


    #############################################
    ### Make arrays that apply to Cinder data ###
    #############################################
    libfile.createGroup(root, "CINDER")
    libfile.createArray("/CINDER", "time",   CoarseTime,    "Cinder Burn Time [days]")			#Burn up times
    libfile.createArray("/CINDER", "E_up",   E_up_cinder,   "Upper Energy Limit - Cinder [MeV]") 		#Upper Energy Limit
    libfile.createArray("/CINDER", "flux_g", flux_g_cinder, "Neutron Group Flux - Cinder [n/s/cm2]")	#Group fluxes
    libfile.createArray("/CINDER", "flux",   flux_cinder,   "Neutron Flux - Cinder [n/s/cm2]")		#Flux

    ##########################################
    ### Make Meta-Stable Branch Ratio Data ###
    ##########################################
    libfile.createGroup("/Coarse", "BranchRatio")
    libfile.createGroup("/Coarse/BranchRatio", "NG")
    libfile.createGroup("/Coarse/BranchRatio", "N2N")
    libfile.createGroup("/Fine", "BranchRatio")
    libfile.createGroup("/Fine/BranchRatio", "NG")
    libfile.createGroup("/Fine/BranchRatio", "N2N")

    #Coarse Data
    for iso in BR_NG.keys():
        iname = isoname.zzaaam_2_LLAAAM(iso)
        libfile.createArray("/Coarse/BranchRatio/NG", iname, BR_NG[iso], "(n, gamma) Meta-Stable to Ground %s Branch Ratio"%iname)	#Meta-Stable Branch Ratio Data
    for iso in BR_N2N.keys():
        iname = isoname.zzaaam_2_LLAAAM(iso)
        libfile.createArray("/Coarse/BranchRatio/N2N", iname, BR_N2N[iso], "(n, 2n) Meta-Stable to Ground %s Branch Ratio"%iname)	#Meta-Stable Branch Ratio Data

    #Fine Data
    BR_NG_Fine  = {}
    for iso in BR_NG.keys():
        finetemp = []
        for t in FineTime:
            n = 0
            while (CoarseTime[n] <= t) and not (t <= CoarseTime[n+1]):
                    n = n + 1
            finetemp.append( metasci.SolveLine(t, CoarseTime[n+1], BR_NG[iso][n+1], CoarseTime[n], BR_NG[iso][n]) )
        iname = isoname.zzaaam_2_LLAAAM(iso)
        libfile.createArray("/Fine/BranchRatio/NG", iname, finetemp, "(n, gamma) Meta-Stable to Ground %s Branch Ratio"%iname)	#Meta-Stable Branch Ratio Data
    BR_N2N_Fine = {}
    for iso in BR_N2N.keys():
        finetemp = []
        for t in FineTime:
            n = 0
            while (CoarseTime[n] <= t) and not (t <= CoarseTime[n+1]):
                    n = n + 1
            finetemp.append( metasci.SolveLine(t, CoarseTime[n+1], BR_N2N[iso][n+1], CoarseTime[n], BR_N2N[iso][n]) )
        iname = isoname.zzaaam_2_LLAAAM(iso)
        libfile.createArray("/Fine/BranchRatio/N2N", iname, finetemp, "(n, 2n) Meta-Stable to Ground %s Branch Ratio"%iname)	#Meta-Stable Branch Ratio Data

    ##############################
    ### Make Multiplier Arrays ###
    ##############################
    for key in FMdic.keys():
        libfile.createGroup("/Coarse", key)
    for key in FMdic.keys():
        libfile.createGroup("/Fine", key)

    #Coarse Data
    for key in mult.keys():
        iso, nada, mult_type = key[:-4].partition("_")
        if mult_type[:5] == "sigma":
            mult_str = mult_type + " [barns]"
        else:
            mult_str = mult_type + " [unitless]"
        libfile.createArray("/Coarse/"+mult_type, iso, mult[key], mult_str)

        #Fine Data
        finetemp = []
        for t in FineTime:
            n = 0
            while (CoarseTime[n] <= t) and not (t <= CoarseTime[n+1]):
                    n = n + 1
            finetemp.append([])
            for e in range( len(E_up) ):
                finetemp[-1].append( metasci.SolveLine(t, CoarseTime[n+1], mult[key][n+1][e], CoarseTime[n], mult[key][n][e]) )
        libfile.createArray("/Fine/"+mult_type, iso, finetemp, mult_str)
    
    libfile.close()
    os.chdir('..')
    return

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
