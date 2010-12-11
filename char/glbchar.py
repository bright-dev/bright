##########################
### Standard Libraries ###
##########################
from __future__ import print_function
import subprocess

##########################
#### Custom Libraries ####
##########################
import numpy as np

import isoname
import MassStream

import metasci
import metasci.nuke as msn

######################
### CHAR Libraries ###
######################
from char import defchar

########################
### Global Functions ###
########################

def iso_list_conversions(iso_list):
    """Converts an isotopic list from a mixed from to zzaaam, LLAAAM, MCNP form as well as doing the
    having a separate lists fo just the metastable isotopes.  Returns a dictionary."""

    zzaaam = sorted( isoname.mixed_2_zzaaam_List(iso_list) )
    metastable = []

    for iso in zzaaam:
        if not ( (iso%10) == 0):
            metastable.append(iso)

            NGammaParent = ((iso/10) - 1) * 10
            if not (NGammaParent in zzaaam):
                zzaaam.append(NGammaParent)

            N2NParent = ((iso/10) + 1) * 10 
            if not (N2NParent in zzaaam):
                zzaaam.append(N2NParent)

    zzaaam = sorted(zzaaam)
    metastable = sorted(metastable)

    iso_dict = {'zzaaam': zzaaam, 
                'LLAAAM': isoname.zzaaam_2_LLAAAM_List(zzaaam),
                'MCNP':    isoname.zzaaam_2_MCNP_List(zzaaam),

                'metastable_zzaaam': metastable, 
                'metastable_LLAAAM': isoname.zzaaam_2_LLAAAM_List(metastable),
                'metastable_MCNP':   isoname.zzaaam_2_MCNP_List(metastable),
                }

    return iso_dict

def iso_file_conversions(filename):
    """Takes a file that contains whitespace separated isotope names and runs iso_list_conversions on it."""
    with open(filename, 'r') as f:
        s = f.read()

    iso_list = s.split()
    iso_dict = iso_list_conversions(iso_list)
    return iso_dict

# Old code that masked isotopes that are not present in the MCNP xsdir file.
# I longer think that this is the right design pattern.  If an isotope should 
# not be included, the user simply shouldn't include it!
""" 
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
"""


##########################
#### Global Variables ####
##########################

class RemoteConnection(object):
    def __init__(self, url='', user='', dir=''):
        self.url  = url
        self.user = user
        self.dir  = dir

    def run(self, cmd):
        return subprocess.call("ssh {user}@{url} \"{remcmd}\"".format(remcmd=cmd, **self.__dict__), shell=True)

    def put(self, loc_file, rem_file):
        return subprocess.call("rsync -r --partial --progress --rsh=ssh {lf} {user}@{url}:{rf}".format(
                                lf=loc_file, rf=rem_file, **self.__dict__), shell=True)

    def get(self, rem_file, loc_file):
        return subprocess.call("rsync -r --partial --progress --rsh=ssh {user}@{url}:{rf} {lf}".format(
                                lf=loc_file, rf=rem_file, **self.__dict__), shell=True)


def defchar_update(defchar):
    """Takes the defchar namespace, updates it, and returns it."""
    if hasattr(defchar, 'tallies'):
        defchar.tallies_reversed = metasci.ReverseDic(defchar.tallies)

    if not hasattr(defchar, 'scheduler'):
        defchar.scheduler = ''

    # Name some files and directories
    defchar.input_file = defchar.reactor + ".i"
    defchar.run_script = 'run_{0}.sh'.format(defchar.reactor)

    if hasattr(defchar, 'remote_dir'):
        defchar.remote_dir = defchar.remote_dir + "runchar/{0}/".format(defchar.reactor)

    # Setup a remote connection instance
    rckw = {}
    if hasattr(defchar, 'remote_url'):
        rckw['url'] = defchar.remote_url
    if hasattr(defchar, 'remote_user'):
        rckw['user'] = defchar.remote_user
    if hasattr(defchar, 'remote_dir'):
        rckw['dir'] = defchar.remote_dir
    defchar.remote_connection = RemoteConnection(**rckw)

    # Make Time Steps 
    defchar.coarse_time = np.arange(0, defchar.burn_time + defchar.coarse_step/10.0, defchar.coarse_step)
    defchar.coarse_time_index = range(len(defchar.coarse_time))

    defchar.fine_time = np.arange(0, defchar.burn_time + defchar.fine_step/10.0, defchar.fine_step)
    defchar.fine_time_index = range(len(defchar.fine_time))

    # Make isotopic lists
    if isinstance(defchar.core_load_isos, basestring):
        defchar.core_load = iso_file_conversions(defchar.core_load_isos)
    elif isinstance(defchar.core_load_isos, list):
        defchar.core_load = iso_list_conversions(defchar.core_load_isos)
    else:
        raise TypeError("The core_load_isos type was not correct.")

    if isinstance(defchar.core_transmute_isos, basestring):
        defchar.core_transmute = iso_file_conversions(defchar.core_transmute_isos)
    elif isinstance(defchar.core_transmute_isos, list):
        defchar.core_transmute = iso_list_conversions(defchar.core_transmute_isos)
    else:
        raise TypeError("The core_transmute_isos type was not correct.")

    # Make fuel stream
    defchar.IHM_stream = MassStream.MassStream(defchar.initial_heavy_metal)
    isovec, AW, MW = msn.convolve_initial_fuel_form(defchar.IHM_stream, defchar.fuel_chemical_form)

    defchar.initial_fuel_stream = MassStream.MassStream(isovec)
    defchar.IHM_weight = AW
    defchar.fuel_weight = MW

    # Make arrays out of quatities that are allowed to vary.
    defchar.fuel_density = np.atleast_1d(defchar.fuel_density)
    defchar.clad_density = np.atleast_1d(defchar.clad_density)
    defchar.cool_density = np.atleast_1d(defchar.cool_density)
    
    defchar.fuel_cell_radius = np.atleast_1d(defchar.fuel_cell_radius)
    defchar.void_cell_radius = np.atleast_1d(defchar.void_cell_radius)
    defchar.clad_cell_radius = np.atleast_1d(defchar.clad_cell_radius)

    defchar.unit_cell_pitch = np.atleast_1d(defchar.unit_cell_pitch)

    defchar.fuel_specific_power = np.atleast_1d(defchar.fuel_specific_power)

    # Set up tuple of parameters to perform a burnup step for
    defchar.perturbation_params = ('fuel_density', 
                                   'clad_density', 
                                   'cool_density',

                                   'fuel_cell_radius', 
                                   'void_cell_radius', 
                                   'clad_cell_radius', 

                                   'unit_cell_pitch', 
                                   'fuel_specific_power', 
                                   'coarse_time')   # coarse_time needs to be the last element

    return defchar
