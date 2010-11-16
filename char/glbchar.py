##########################
### Standard Libraries ###
##########################
from __future__ import print_function
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
from char import defchar

########################
### Global Functions ###
########################


##########################
#### Global Variables ####
##########################

class RemoteConnection(object):
    def __init__(self, url='', user='', dir=''):
        self.url  = url
        self.user = user
        self.dir  = dir

    def run(self, cmd):
        return subprocess.call("ssh {user}@{url} \"{remcmd}\"".format(remcmd=cmd, **self), shell=True)

    def put(self, loc_file, rem_file):
        return subprocess.call("scp {lf} {user}@{url}:{rf}".format(lf=loc_file, rf=rem_file, **self), shell=True)

    def get(self, rem_file, loc_file):
        return subprocess.call("scp {user}@{url}:{rf} {lf}".format(lf=loc_file, rf=rem_file, **self), shell=True)


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

try:
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
except:
    pass

coreload = []
for iso in InXSDIR.keys():
    if InXSDIR[iso]:
        coreload.append(iso)
    else:
        if 0 < verbosity:
            print("The following nuclide could not be found in $DATAPATH/xsdir: {0}.".format(isoname.MCNP_2_LLAAAM(iso)))

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


def defchar_update(defchar):
    """Takes the defchar namespace, updates it, and returns it."""
    defchar.tallies_reversed = metasci.ReverseDic(defchar.tallies)

    # Name some files and directories
    defchar.input_file = defchar.reactor + ".i"
    defchar.run_script = 'run_{0}.sh'.format(defchar.reactor)

    if hasattr(defchar, 'remote_dir'):
        defchar.remote_dir = defchar.remote_dir + "runchar/{0}/".format(defchar.reactor)

    # Setup a remote connection instance
    rckw = {}
    if hasattr(defchar, 'remote_url'):
        rckw['url'] = defcahr.remote_url
    if hasattr(defchar, 'remote_user'):
        rckw['user'] = defcahr.remote_user
    if hasattr(defchar, 'remote_dir'):
        rckw['dir'] = defcahr.remote_dir
    defchar.remote_connection = RemoteConnection(**rckw)

    defchar.initial_fuel_stream = MassStream.MassStream(defchar.initial_fuel_form)
    return defchar
