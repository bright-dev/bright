##########################
### Standard Libraries ###
##########################
import os
import subprocess
from math import pi

##########################
#### Custom Libraries ####
##########################
import isoname
import MassStream
import metasci

######################
### CHAR Libraries ###
######################
from defchar import *

########################
### Global Functions ###
########################

def QPrint(s, quiet=False):
    """Quick, Quiet Print override."""
    global Quiet

    if not (quiet or Quiet):
        if type(s) == str:
            print s
        elif type(s) == list:
            for el in s:
                print el[:-1]
        else:
            pass
    return s

##########################
#### Global Variables ####
##########################
dicFM = metasci.ReverseDic(FMdic)

inputfile = reactor + ".i"
runscript = 'run_{0}.sh'.format(reactor)
RemoteDir = RemoteDir + "runchar/{0}/".format(reactor)

class RC:
    def __init__(self):
        self.RemoteURL  = RemoteURL
        self.RemoteUser = RemoteUser
        self.RemoteDir  = RemoteDir
    def run(self, cmd):
        return subprocess.call("ssh {rc.RemoteUser}@{rc.RemoteURL} \"{remcmd}\"".format(rc=self, remcmd=cmd), shell=True)
    def put(self, loc_file, rem_file):
        return subprocess.call("scp {lf} {rc.RemoteUser}@{rc.RemoteURL}:{rf}".format(rc=self, lf=loc_file, rf=rem_file), shell=True)
    def get(self, rem_file, loc_file):
        return subprocess.call("scp {rc.RemoteUser}@{rc.RemoteURL}:{rf} {lf}".format(rc=self, lf=loc_file, rf=rem_file), shell=True)
RemoteConnection = RC()

ORIGEN_FASTorTHERM = ORIGEN_FASTorTHERM.lower()
if not (ORIGEN_FASTorTHERM in ['fast', 'therm']):
    print "ORIGEN_FASTorTHERM not properly set:", ORIGEN_FASTorTHERM
    print "Setting to \'therm\'."

#########################
### Make Cell Volumes ###
#########################
#Division by 4 because template only looks at a quarter of the unit cell
FuelCellVolume = 0.25 * UnitCellHeight * pi * (FuelCellRadius**2) 
CladCellVolume = (0.25 * UnitCellHeight * pi * (CladCellRadius**2)) - FuelCellVolume
UnitCellVolume = 0.25 * UnitCellHeight * (UnitCellPitch**2)
CoolCellVolume = UnitCellVolume - FuelCellVolume - CladCellVolume

UnitCellHalfPitch = UnitCellPitch / 2.0 

#######################
### Make Time Steps ###
#######################
CoarseTime = range(0, BurnTime, CoarseStep/2)
CoarseTime.append(BurnTime)
CoarseTimeIndex = range(len(CoarseTime))

FineTime = range(0, BurnTime, FineStep)
FineTime.append(BurnTime)
FineTimeIndex = range(len(FineTime))

######################################################
### Makes the Core Loading isotopic tracking lists ###
######################################################
coreload = []
f = open(coreloadtrackfile, 'r')
for line in f:
    coreload.append(line.split()[0])
f.close()
coreload = sorted( isoname.mixed_2_zzaaam_List(coreload) )
metastabletrak = []
for iso in coreload:
    if not ( (iso%10) == 0):
        metastabletrak.append(iso)
        NGammaParent = ((iso/10) - 1) * 10
        if not (NGammaParent in coreload):
            coreload.append(NGammaParent)
        N2NParent = ((iso/10) + 1) * 10 
        if not (N2NParent in coreload):
            coreload.append(N2NParent)
coreload = isoname.zzaaam_2_MCNP_List(coreload)
metastableMCNP = isoname.zzaaam_2_MCNP_List(metastabletrak)
InXSDIR = {}
for iso in coreload:
    InXSDIR[iso] = False
xsdir = open(os.getenv("DATAPATH") + "/xsdir", 'r' )
for line in xsdir:
    ls = line.split()
    if ls == []:
        continue
    elif not ('.' in ls[0]):
        continue
    else:
        i, p, l = ls[0].partition('.')
        try:
            xs_i = int(i)
        except:
            continue
        for iso in InXSDIR.keys():
            if xs_i == int(iso):
                InXSDIR[iso] = True
xsdir.close()
coreload = []
for iso in InXSDIR.keys():
    if InXSDIR[iso]:
        coreload.append(iso)
    else:
        if Verbose:
            QPrint("The following nuclide could not be found in $DATAPATH/xsdir: %s."%isoname.MCNP_2_LLAAAM(iso) )
CoreLoad_zzaaam = isoname.MCNP_2_zzaaam_List(coreload)
CoreLoad_LLAAAM = isoname.MCNP_2_LLAAAM_List(coreload)
CoreLoad_MCNP   = coreload

#############################################################
### Makes the Core Transformation isotopic tracking lists ###
#############################################################
coretran = []
f = open(coretrantrackfile, 'r')
for line in f:
    coretran.append(line.split()[0])
f.close()
coretran = sorted( isoname.mixed_2_zzaaam_List(coretran) )
for iso in coretran:
    if not (iso%10 == 0):
        NGammaParent = ((iso/10) - 1) * 10 
        if not (NGammaParent in coretran):
            coretran.append(NGammaParent)
        N2NParent = ((iso/10) + 1) * 10
        if not (N2NParent in coretran):
            coretran.append(N2NParent)
coretran = isoname.zzaaam_2_MCNP_List(coretran)
#InXSDIR = {}
#for iso in coretran:
#	InXSDIR[iso] = False
#xsdir = open(os.getenv("DATAPATH") + "/xsdir", 'r' )
#for line in xsdir:
#	ls = line.split()
#	if ls == []:
#		continue
#	elif not ('.' in ls[0]):
#		continue
#	else:
#		i, p, l = ls[0].partition('.')
#		try:
#			xs_i = int(i)
#		except:
#			continue
#		for iso in InXSDIR.keys():
#			if xs_i == int(iso):
#				InXSDIR[iso] = True
#xsdir.close()
#coretran = []
#for iso in InXSDIR.keys():
#	if InXSDIR[iso]:
#		coretran.append(iso)
#	else:
#		libcom.QPrint("The following nuclide could not be found in $DATAPATH/xsdir: %s."%isoname.MCNP_2_LLAAAM(iso) )
CoreTran_zzaaam = isoname.MCNP_2_zzaaam_List(coretran)
CoreTran_LLAAAM = isoname.MCNP_2_LLAAAM_List(coretran)
CoreTran_MCNP   = coretran

############################################
### Map from isotopes to material number ###
############################################
mat_number = {}		
mnum = moffset
for iso in CoreLoad_MCNP:
    mnum = mnum + 1
    mat_number[iso] = mnum
del mnum 
number_mat = metasci.ReverseDic(mat_number)

InitialFuelStream = MassStream.MassStream(InitialFuelForm)


######################
### ORIGEN Options ###
######################

#Format/Set Concentration cut off
try:
    ORIGEN_Cut_Off = "{0:.3E}".format(ORIGEN_Concentration_Cut_Off)
except NameError:
    ORIGEN_Cut_Off = '1.-10'

