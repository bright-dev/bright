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

from pyne import nucname
from pyne.material import Material

######################
### CHAR Libraries ###
######################
import bright.xsgen.utils as utils
from char.tally_types import restricted_tallies

##########################
#### Global Variables ####
##########################
initial_nuc_pattern = 'initial_([A-Za-z]{0,2}\d{1,7}[Mm]?)'


def update_env_for_execution(env):
    """Updates the env namespace for runs where an execution is going to occur."""
    # Make nuclide lists
    if isinstance(env.core_load_nucs, basestring):
        core_load = load_nuc_file(env.core_load_nucs)
    else:
        core_load = [nucname.zzaaam(nuc) for nuc in env.core_load_nucs]
    env.core_load = sorted(set(core_load))

    if isinstance(env.core_transmute_nucs, basestring):
        core_transmute = load_nuc_file(env.core_transmute_nucs)
    else:
        core_transmute = [nucname.zzaaam(nuc) for nuc in env.core_transmute_nucs]
    env.core_transmute = sorted(set(core_transmute))

    # Make temperature
    env.temperature = env.get('temperature', 600)

    # Make fuel stream
    env.ihm_mat = Material(env.initial_heavy_metal)

    # make sensitivity mass fractions
    if 'sensitivity_mass_fractions' in env:
        env.deltam = np.atleast_1d(env.sensitivity_mass_fractions)
        env.deltam.sort()


    # Make arrays out of quatities that are allowed to vary.
    env.fuel_density = np.atleast_1d(env.fuel_density)
    env.clad_density = np.atleast_1d(env.clad_density)
    env.cool_density = np.atleast_1d(env.cool_density)

    env.fuel_cell_radius = np.atleast_1d(env.fuel_cell_radius)
    env.void_cell_radius = np.atleast_1d(env.void_cell_radius)
    env.clad_cell_radius = np.atleast_1d(env.clad_cell_radius)

    env.unit_cell_pitch = np.atleast_1d(env.unit_cell_pitch)

    env.burn_regions = np.atleast_1d(env.burn_regions)
    env.fuel_specific_power = np.atleast_1d(env.fuel_specific_power)

    # Grab initial nuc perturbation
    max_mass = 0.0
    initial_nuc_keys = []
    for key in env:
        m = re.match(initial_nuc_pattern, key)
        if m is None:
            continue

        env_initial_nuc = env[key]
        env_initial_nuc = np.atleast_1d(env_initial_nuc)
        env[key] = env_initial_nuc

        initial_nuc_keys.append(key)
        max_mass += np.max(env_initial_nuc)

    initial_nuc_keys.sort()
    env.initial_nuc_keys = initial_nuc_keys

    if 1.0 < max_mass:
        msg = "The maxium mass of initial heavy metal perturbations exceeds 1.0 kg!"
        raise SystemExit(utils.failure(msg))

    # Set up tuple of parameters to perform a burnup step for
    env.perturbation_params = ['fuel_density',
                                  'clad_density',
                                  'cool_density',

                                  'fuel_cell_radius',
                                  'void_cell_radius',
                                  'clad_cell_radius',

                                  'unit_cell_pitch',
                                  'burn_regions',
                                  'fuel_specific_power',]

    env.perturbation_params.extend(initial_nuc_keys)

    # burn_times needs to be the last element
    env.perturbation_params.append('burn_times')

    return env


def update_env(env):
    """Takes the env namespace, updates it, and returns it."""
    env.transporter = env.get('transporter', '')
    env.scheduler = env.get('scheduler', '')

    # Name some files and directories
    env.run_script = 'run_{0}.sh'.format(env.reactor)

    # Setup a remote connection instance
    rckw = {}
    if 'remote_url' in env:
        rckw['url'] = env.remote_url
    if 'remote_user' in env:
        rckw['user'] = env.remote_user
    if 'remote_dir' in env:
        env.remote_dir = env.remote_dir + "runchar/{0}/".format(env.reactor)
        rckw['dir'] = env.remote_dir
    env.remote_connection = utils.RemoteConnection(**rckw)

    # Make Time Steps
    if 'burn_times' in env:
        env.burn_times = np.asarray(env.burn_times, dtype=float)
    else:
        bt_upper_lim = env.burn_time + env.time_step/10.0
        env.burn_times = np.arange(0, bt_upper_lim, env.time_step)
    env.burn_times_index = range(len(env.burn_times))

    # Upddate the namespace if we are going to execute a neutronics code
    if (env.MAKE_INPUT or env.RUN_TRANSPORT or
        env.RUN_BURNUP or env.RUN_XS_GEN or
        env.RUN_DELTAM):

        env = update_env_for_execution(env)

    return env


