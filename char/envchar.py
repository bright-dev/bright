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

import nucname
from pyne.material import Material

import metasci
import metasci.nuke as msn
from metasci.nuke import ace
from metasci.colortext import message, failure

######################
### CHAR Libraries ###
######################
from n_code_serpent import zzaaam_2_serpent
from tally_types import restricted_tallies
from char import utils

########################
### Global Functions ###
########################
initial_iso_pattern = 'initial_([A-Za-z]{1,2}\d{1,3}[Mm]?)'


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

serpent_mt_fission = set([-6, 18, 19, 20, 21, 38])
"""A set of MT numbers for fission cross-sections in Serpent."""

serpent_mt_nubar = set([-7])
"""A set of MT numbers for the number of neutrons per fission times the fission cross-sections in Serpent."""

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
        iso_mt = (mts | serpent_mt_always)

        if (iso_zz, temp_flag) in restricted_tallies:
            iso_mt = iso_mt - restricted_tallies[(iso_zz, temp_flag)]

        if 0 == len(iso_mt & serpent_mt_fission - set([-6])):
            # if isotopic fission not avilable, remove material 
            # fission and nubar from avilable tallies
            iso_mt = iso_mt - serpent_mt_fission
            iso_mt = iso_mt - serpent_mt_nubar

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




def update_env_for_execution(env):
    """Updates the env namespace for runs where an execution is going to occur."""
    # Make isotopic lists
    if isinstance(env['core_load_isos'], basestring):
        env['core_load'] = iso_file_conversions(env['core_load_isos'])
    elif isinstance(env['core_load_isos'], list):
        env['core_load'] = iso_list_conversions(env['core_load_isos'])
    else:
        raise TypeError("The core_load_isos type was not correct.")

    if isinstance(env['core_transmute_isos'], basestring):
        env['core_transmute'] = iso_file_conversions(env['core_transmute_isos'])
    elif isinstance(env['core_transmute_isos'], list):
        env['core_transmute'] = iso_list_conversions(env['core_transmute_isos'])
    else:
        raise TypeError("The core_transmute_isos type was not correct.")

    # Find which isotopes are available in serpent
    # and which ones must be handled manually.
    core_transmute_set = set(env['core_transmute']['zzaaam'])
    serpent_xs_isos_set = serpent_xs_isos_available(env['serpent_xsdata'])

    core_transmute_in_serpent = list(core_transmute_set & serpent_xs_isos_set)
    core_transmute_not_in_serpent = list(core_transmute_set - serpent_xs_isos_set)

    env['core_transmute_in_serpent'] = iso_list_conversions(core_transmute_in_serpent)
    env['core_transmute_not_in_serpent'] = iso_list_conversions(core_transmute_not_in_serpent)

    #env['xs_models_needed'] = (0 < len(env['core_transmute_not_in_serpent']))

    # Make temperature flag
    if 'temperature' not in env:
        env['temperature'] = 600

    env['temp_flag'] = temperature_flag(env['temperature'])

    # Grab the MT numbers that are available for valid isotopes.
    env['iso_mts'] = serpent_mt_avaliable(env['serpent_xsdata'], 
                                          env['core_transmute_in_serpent']['zzaaam'], 
                                          env['temp_flag'], 
                                          env['verbosity'])
    #print(env['iso_mts'][501250]); raise SystemExit

    # Make fuel stream
    env['IHM_stream'] = Material(env['initial_heavy_metal'])

    if 'sensitivity_mass_fractions' in env:
        env['deltam'] = np.atleast_1d(env['sensitivity_mass_fractions'])
        env['deltam'].sort()


    # Make arrays out of quatities that are allowed to vary.
    env['fuel_density'] = np.atleast_1d(env['fuel_density'])
    env['clad_density'] = np.atleast_1d(env['clad_density'])
    env['cool_density'] = np.atleast_1d(env['cool_density'])
    
    env['fuel_cell_radius'] = np.atleast_1d(env['fuel_cell_radius'])
    env['void_cell_radius'] = np.atleast_1d(env['void_cell_radius'])
    env['clad_cell_radius'] = np.atleast_1d(env['clad_cell_radius'])

    env['unit_cell_pitch'] = np.atleast_1d(env['unit_cell_pitch'])

    env['burn_regions'] = np.atleast_1d(env['burn_regions'])
    env['fuel_specific_power'] = np.atleast_1d(env['fuel_specific_power'])

    # Grab initial iso perturbation
    max_mass = 0.0
    initial_iso_keys = []
    for key in env:
        m = re.match(initial_iso_pattern, key)
        if m is None:
            continue

        env_initial_iso = env[key]
        env_initial_iso = np.atleast_1d(env_initial_iso)

        initial_iso_keys.append(key)
        max_mass += np.max(env_initial_iso)

    initial_iso_keys.sort()
    env['initial_iso_keys'] = initial_iso_keys

    if 1.0 < max_mass:
        print(failure("The maxium mass of initial heavy metal perturbations exceeds 1.0 kg!"))
        raise SystemExit

    # Set up tuple of parameters to perform a burnup step for
    env['perturbation_params'] = ['fuel_density', 
                                  'clad_density', 
                                  'cool_density',

                                  'fuel_cell_radius', 
                                  'void_cell_radius', 
                                  'clad_cell_radius', 

                                  'unit_cell_pitch', 
                                  'burn_regions', 
                                  'fuel_specific_power',]

    env['perturbation_params'].extend(initial_iso_keys)
 
    # burn_times needs to be the last element
    env['perturbation_params'].append('burn_times')   

    return env


def update_env(env):
    """Takes the env[' namespace, updates it, and returns it."""
    if 'transporter' not in env:
        env['transporter'] = ''

    if 'scheduler' not in env:
        env['scheduler'] = ''

    # Name some files and directories
    env['input_file'] = env['reactor'] + ".i"
    env['run_script'] = 'run_{0}.sh'.format(env['reactor'])

    if 'remote_dir' in env:
        env['remote_dir'] = env['remote_dir'] + "runchar/{0}/".format(env['reactor'])

    # Setup a remote connection instance
    rckw = {}
    if 'remote_url' in env:
        rckw['url'] = env['remote_url']
    if 'remote_user' in env:
        rckw['user'] = env['remote_user']
    if 'remote_dir' in env:
        rckw['dir'] = env['remote_dir']
    env['remote_connection'] = RemoteConnection(**rckw)

    # Make Time Steps 
    env['burn_times'] = np.arange(0, env['burn_time'] + env['time_step']/10.0, env['time_step'])
    env['burn_times_index'] = range(len(env['burn_times']))

    # Upddate the namespace if we are going to execute a neutronics code
    if (env['options'].MAKE_INPUT or env['options'].RUN_TRANSPORT or 
        env['options'].RUN_BURNUP or env['options'].RUN_XS_GEN or 
        env['options'].RUN_DELTAM):

        env = update_env_for_execution(env)

    return env


