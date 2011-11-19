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

#import metasci
#import metasci.nuke as msn
#from metasci.nuke import ace
#from metasci.colortext import message, failure

######################
### CHAR Libraries ###
######################
from tally_types import restricted_tallies
from char import utils

##########################
#### Global Variables ####
##########################
initial_iso_pattern = 'initial_([A-Za-z]{1,2}\d{1,3}[Mm]?)'
    
class RemoteConnection(object):
    def __init__(self, url='', user='', dir=''):
        self.url  = url
        self.user = user
        self.dir  = dir

    def run(self, cmd):
        return subprocess.call('ssh {user}@{url} \"{remcmd}\"'.format(remcmd=cmd, **self.__dict__), shell=True)

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
        env['core_load'] = load_nuc_file(env['core_load_isos'])
    else:
        env['core_load'] = sorted(nucname.zzaaam(nuc) for nuc in env['core_load_isos'])

    if isinstance(env['core_transmute_isos'], basestring):
        env['core_transmute'] = load_nuc_file(env['core_transmute_isos'])
    else:
        env['core_transmute'] = sorted(nucname.zzaaam(nuc) for nuc in env['core_transmute_isos'])

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


