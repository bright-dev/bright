##########################
### Standard Libraries ###
##########################
from __future__ import print_function
import re
import subprocess

##########################
#### Custom Libraries ####
##########################
import numpy as np

import isoname
import MassStream

import metasci
import metasci.nuke as msn
from metasci.nuke import ace
from metasci.colortext import message, failure

######################
### CHAR Libraries ###
######################
from char import defchar
from n_code_serpent import zzaaam_2_serpent

########################
### Global Functions ###
########################
initial_iso_pattern = 'initial_([A-Za-z]{1,2}\d{1,3}[Mm]?)'


def iso_list_conversions(iso_list):
    """Converts an isotopic list from a mixed from to zzaaam, LLAAAM, MCNP form as well as doing the
    having a separate lists fo just the metastable isotopes.  Returns a dictionary."""

    zzaaam = sorted( isoname.mixed_2_zzaaam_List(iso_list) )
    metastable = []

    for iso in zzaaam:
        if not ( (iso%10) == 0):
            continue
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


def serpent_xs_isos_available(xsdata):
    """Finds the isotopes available to serpent for cross-section generation.

    Args:
        * xsdata (str): path to serpent *.xsdata file that will be used.

    Returns:
        * serpent_isos (set): Set of isotopes serpent has available.
    """
    # Grabs the ZAID and metastable flag separately
    xsdata_pattern = "\s*[\dA-Za-z-]+\.\d{2}[a-z]\s+\d{4,6}\.\d{2}[a-z]  \d\s+(\d{4,6})  (\d)\s+.*"

    with open(xsdata, 'r') as f:
        raw_xsdata = f.read()

    serpent_iso = set(int(''.join(m.groups())) for m in re.finditer(xsdata_pattern, raw_xsdata))
    return serpent_iso


serpent_mt_always = set(range(-9, 3))
"""A set of MT numbers that is always available in Serpent."""

def serpent_mt_avaliable(xsdata, isos, temp_flag, verbosity=100):
    """Finds the MT numbers available for each isotope.

    Args:
        * xsdata (str): path to serpent *.xsdata file that will be used.
        * isos (list of zzaaam): List of isotopes to find MT numbers for. 
          isotopes must be valid for serpent.
        * temp_flag (3-character string): Flag for the temperature.

    Returns:
        * iso_mt (dict of sets): A dictionary whose keys are isotopes (zzaaam) and whose 
          keys are sets of MT numbers that serpent has available.
    """
    if 0 < verbosity:
        print(message("Grabbing valid MT numbers for available isotopes:"))

    # First, read in the xsdata file
    xsdata_dict = {}
    with open(xsdata, 'r') as f:
        for line in f:
            ls = line.split()
            xsdata_dict[ls[0]] = (ls[1], ls[-1])

    # Now, find the MTs for each iso
    iso_mts = {}
    for iso_zz in isos:
        # Convert iso 
        iso_serp = zzaaam_2_serpent(iso_zz)
        iso_serp_flag = "{0}.{1}".format(iso_serp, temp_flag)

        if 0 < verbosity:
            print("  Isotope {0:>7} {1:>11}".format(iso_zz, iso_serp_flag))

        # Get the MT numbers
        mts = ace.mt(*xsdata_dict[iso_serp_flag])
        iso_mt = mts | serpent_mt_always

        # Add this iso to the dict
        iso_mts[iso_zz] = iso_mt

    if 0 < verbosity:
        print(message("Done!"))
        print()

    return iso_mts

def temperature_flag(t):
    """Converts a temperature into a the proper continuous energy flag used in ACE files.

    Args:
        * t (int): Temperature, multiple of 300 K.

    Returns: 
        * temp_flag (3-character string)
    """

    t = int(t)

    # Check temperature value validity
    if t%300 != 0:
        raise ValueError("The temperature value must be a multiple of 300 K!")
    elif t <= 0:
        raise ValueError("The temperature value must be positive!")
    elif 9999 < t:
        raise ValueError("The temperature value must less than 10000 K!")

    # Make the temperature flag
    temp_flag = "{0:02}c".format(t/100)

    return temp_flag


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
        return subprocess.call("rsync -rh --partial --progress --rsh=ssh {lf} {user}@{url}:{rf}".format(
                                lf=loc_file, rf=rem_file, **self.__dict__), shell=True)

    def get(self, rem_file, loc_file):
        return subprocess.call("rsync -rh --partial --progress --rsh=ssh {user}@{url}:{rf} {lf}".format(
                                lf=loc_file, rf=rem_file, **self.__dict__), shell=True)



def defchar_update_for_execution(defchar):
    """Updates the defchar namespace for runs where an execution is going to occur."""
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

    # Find which isotopes are available in serpent
    # and which ones must be handled manually.
    core_transmute_set = set(defchar.core_transmute['zzaaam'])
    serpent_xs_isos_set = serpent_xs_isos_available(defchar.serpent_xsdata)

    core_transmute_in_serpent = list(core_transmute_set & serpent_xs_isos_set)
    core_transmute_not_in_serpent = list(core_transmute_set - serpent_xs_isos_set)

    defchar.core_transmute_in_serpent = iso_list_conversions(core_transmute_in_serpent)
    defchar.core_transmute_not_in_serpent = iso_list_conversions(core_transmute_not_in_serpent)

    # Make temperature flag
    if not hasattr(defchar, 'temperature'):
        defchar.temperature = 600

    defchar.temp_flag = temperature_flag(defchar.temperature)

    # Grab the MT numbers that are available for valid isotopes.
    defchar.iso_mts = serpent_mt_avaliable(defchar.serpent_xsdata, 
                                           defchar.core_transmute_in_serpent['zzaaam'], 
                                           defchar.temp_flag, 
                                           defchar.verbosity)

    # Make fuel stream
    defchar.IHM_stream = MassStream.MassStream(defchar.initial_heavy_metal)

    if hasattr(defchar, 'sensitivity_mass_fractions'):
        defchar.deltam = np.atleast_1d(defchar.sensitivity_mass_fractions)
        defchar.deltam.sort()


    # Make arrays out of quatities that are allowed to vary.
    defchar.fuel_density = np.atleast_1d(defchar.fuel_density)
    defchar.clad_density = np.atleast_1d(defchar.clad_density)
    defchar.cool_density = np.atleast_1d(defchar.cool_density)
    
    defchar.fuel_cell_radius = np.atleast_1d(defchar.fuel_cell_radius)
    defchar.void_cell_radius = np.atleast_1d(defchar.void_cell_radius)
    defchar.clad_cell_radius = np.atleast_1d(defchar.clad_cell_radius)

    defchar.unit_cell_pitch = np.atleast_1d(defchar.unit_cell_pitch)

    defchar.burn_regions = np.atleast_1d(defchar.burn_regions)
    defchar.fuel_specific_power = np.atleast_1d(defchar.fuel_specific_power)

    # Grab initial iso perturbation
    max_mass = 0.0
    initial_iso_vars = []
    for var in defchar.__dict__:
        m = re.match(initial_iso_pattern, var)
        if m is None:
            continue

        defchar_initial_iso = getattr(defchar, var)
        defchar_initial_iso = np.atleast_1d(defchar_initial_iso)

        initial_iso_vars.append(var)
        max_mass += np.max(defchar_initial_iso)

    initial_iso_vars.sort()
    defchar.initial_iso_vars = initial_iso_vars

    if 1.0 < max_mass:
        print(failure("The maxium mass of initial heavy metal perturbations exceeds 1.0 kg!"))
        raise SystemExit

    # Set up tuple of parameters to perform a burnup step for
    defchar.perturbation_params = ['fuel_density', 
                                   'clad_density', 
                                   'cool_density',

                                   'fuel_cell_radius', 
                                   'void_cell_radius', 
                                   'clad_cell_radius', 

                                   'unit_cell_pitch', 
                                   'burn_regions', 
                                   'fuel_specific_power',]

    defchar.perturbation_params.extend(initial_iso_vars)
 
    defchar.perturbation_params.append('burn_times')   # burn_times needs to be the last element

    return defchar


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
    defchar.burn_times = np.arange(0, defchar.burn_time + defchar.time_step/10.0, defchar.time_step)
    defchar.burn_times_index = range(len(defchar.burn_times))

    # Upddate the namespace if we are going to execute a neutronics code
    if (defchar.options.MAKE_INPUT or defchar.options.RUN_TRANSPORT or 
        defchar.options.RUN_BURNUP or defchar.options.RUN_XS_GEN or 
        defchar.options.RUN_DELTAM):

        defchar = defchar_update_for_execution(defchar)

    return defchar


