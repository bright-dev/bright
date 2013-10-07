"""A class to setup, run, and parse Serpent."""
from __future__ import print_function
import os
import re
import sys
import subprocess
from itertools import product

import numpy as np
import tables as tb

from pyne import nucname
from pyne import ace
from pyne.material import Material, from_atom_frac
from pyne.pyne_config import pyne_conf
from pyne.xs.cache import xs_cache
import pyne.xs.channels

from scipy.integrate import cumtrapz

import tally_types
from pyne.utils import message, failure
import pyne.serpent

from char import utils
from char.tally_types import restricted_tallies

# Hide warnings from numpy
np.seterr(divide='ignore')

# Setup serpent running
try:
    import serpent
except ImportError:
    pass

def run_serpent_subprocess(argstr):
    rtn = subprocess.check_call('sss-dev ' + argstr, shell=True)
    return rtn

run_serpent_in_switch = {
    '': serpent.main, 
    'sh': run_serpent_subprocess,
    'bash': run_serpent_subprocess,
    'subprocess': run_serpent_subprocess,
    'python': serpent.main, 
    'module': serpent.main, 
    }

partial_fission_mts = set([19, 20, 21, 38])


#
# Helper Functions
#


def serpent_xs_nucs_available(xsdata):
    """Finds the nuclides available to serpent for cross-section generation.

    Args:
        * xsdata (str): path to serpent *.xsdata file that will be used.

    Returns:
        * serpent_nucs (set): Set of nuclides serpent has available.
    """
    # Grabs the ZAID and metastable flag separately
    xsdata_pattern = "\s*[\dA-Za-z-]+\.\d{2}[a-z]\s+\d{4,6}\.\d{2}[a-z]  \d\s+(\d{4,6})  (\d)\s+.*"

    with open(xsdata, 'r') as f:
        raw_xsdata = f.read()

    serpent_nuc = set(int(''.join(m.groups())) for m in re.finditer(xsdata_pattern, raw_xsdata))
    return serpent_nuc


serpent_mt_always = set(range(-9, 3))
"""A set of MT numbers that is always available in Serpent."""

serpent_mt_fission = set([-6, 18, 19, 20, 21, 38])
"""A set of MT numbers for fission cross-sections in Serpent."""

serpent_mt_nubar = set([-7])
"""A set of MT numbers for the number of neutrons per fission times the fission cross-sections in Serpent."""

def serpent_mt_avaliable(xsdata, nucs, temp_flag, verbosity=100):
    """Finds the MT numbers available for each nuclide.

    Args:
        * xsdata (str): path to serpent *.xsdata file that will be used.
        * nucs (list of zzaaam): List of nuclides to find MT numbers for. 
          nuclides must be valid for serpent.
        * temp_flag (3-character string): Flag for the temperature.

    Returns:
        * nuc_mt (dict of sets): A dictionary whose keys are nuclides (zzaaam) and whose 
          keys are sets of MT numbers that serpent has available."""
    if 0 < verbosity:
        print(message("Grabbing valid MT numbers for available nuclides:"))

    # First, read in the xsdata file
    xsdata_dict = {}
    with open(xsdata, 'r') as f:
        for line in f:
            ls = line.split()
            xsdata_dict[ls[0]] = (ls[1], ls[-1])

    # Now, find the MTs for each nuc
    nuc_mts = {}
    for nuc_zz in nucs:
        # Convert nuc 
        nuc_serp = nucname.serpent(nuc_zz)
        nuc_serp_flag = "{0}.{1}".format(nuc_serp, temp_flag)

        if 0 < verbosity:
            print("  Nuclide {0:>7} {1:>11}".format(nuc_zz, nuc_serp_flag))

        # Get the MT numbers
        ace_lib = ace.Library(xsdata_dict[nuc_serp_flag][1])
        ace_lib.read()
        nuc_tab = ace_lib.find_table(xsdata_dict[nuc_serp_flag][0])
        nuc_tab._read_mtr()
        mts = set(rx.MT for rx in nuc_tab.reactions)
        nuc_mt = (mts | serpent_mt_always)

        if (nuc_zz, temp_flag) in restricted_tallies:
            nuc_mt = nuc_mt - restricted_tallies[(nuc_zz, temp_flag)]

        if 0 == len(nuc_mt & serpent_mt_fission - set([-6])):
            # if nuclide fission not avilable, remove material 
            # fission and nubar from avilable tallies
            nuc_mt = nuc_mt - serpent_mt_fission
            nuc_mt = nuc_mt - serpent_mt_nubar

        # Add this nuc to the dict
        nuc_mts[nuc_zz] = nuc_mt

    if 0 < verbosity:
        print(message("Done!"))
        print()

    return nuc_mts



#
# Neutron Code Class
#


class NCodeSerpent(object):
    """A Serpent neutronics code wrapper class."""

    def __init__(self, env):
        self.env = self.update_env(env)

        self.name = "Serpent"
        self.run_str = "sss-dev"

        # Get run serpent function.
        if 'run_serpent_in' not in self.env:
            self.env['run_serpent_in'] = 'python'
        self.run_serpent = run_serpent_in_switch[self.env['run_serpent_in']]

        # Remote file lists
        self.place_remote_files = ['.']
        self.fetch_remote_files = ['.']

        # Set some helpful attributes 
        self.mpi_flag = self.get_mpi_flag()
        self.ntimes = len(self.env['burn_times'])
        self.G = len(self.env['group_structure']) - 1

        # Grab stock high-resolution group structure
        with tb.openFile(pyne_conf.NUC_DATA_PATH, 'r') as f:
            self.hi_res_group_structure = np.array(f.root.neutron.cinder_xs.E_g)

        # Set Stock tallies for serpent
        if 'tallies' not in env:
            self.env['tallies'] = tally_types.serpent_default
        self.env['tallies'] = {t: tally_types.serpent_tallies[t] for t in self.env['tallies']}

        # Make perturbation table
        data = [self.env[a] for a in self.env['perturbation_params']]
        self.perturbations = [p for p in product(*data)]
        self.nperturbations = len(self.perturbations)
        self.pert_cols = {p: np.array([row[i] for row in self.perturbations])
                                              for i, p in enumerate(self.env['perturbation_params'])}


    #
    # Environment updater
    #
    def update_env(self, env):
        """Update environment to work with serpent."""

        # Find which nuclides are available in serpent
        # and which ones must be handled manually.
        core_transmute_set = set(env['core_transmute'])
        serpent_xs_nucs_set = serpent_xs_nucs_available(env['serpent_xsdata'])

        core_transmute_in_serpent = sorted(core_transmute_set & serpent_xs_nucs_set)
        core_transmute_not_in_serpent = sorted(core_transmute_set - serpent_xs_nucs_set)

        env['core_transmute_in_serpent'] = core_transmute_in_serpent
        env['core_transmute_not_in_serpent'] = core_transmute_not_in_serpent

        # Make temperature
        env['temp_flag'] = utils.temperature_flag(env['temperature'])

        # Grab the MT numbers that are available for all valid nuclides
        env['nuc_mts'] = serpent_mt_avaliable(env['serpent_xsdata'], 
                                              env['core_transmute_in_serpent'], 
                                              env['temp_flag'], 
                                              env['verbosity'])
        return env


    #
    # Serpent input file generation methods
    #

    def make_input_material_weights(self, comp, mass_weighted=True):
        """This function takes an nuclide vector, comp, and returns a serpent string representation.
        Note that by default these are mass weights, rather than atomic weights."""
        comp_str = ''

        for nuc in comp.keys():
            nuc_serp = nucname.serpent(nuc)

            comp_str += "{0:>11}".format("{0}.{1}".format(nuc_serp, self.env['temp_flag']))

            if mass_weighted:
                comp_str += "  -{0:.5G}\n".format(comp[nuc])
            else:
                comp_str += "   {0:.5G}\n".format(comp[nuc])
        return comp_str
            

    def update_initial_fuel(self, n):
        """Required to allow for initial fuel nuclide concentration perturbations."""
        if len(self.env['initial_nuc_keys']) == 0:
            # Don't bother readin the perturbation table if there 
            # is nothing there to read!
            ihm_mat = 1.0 * self.env['ihm_mat']
        else:
            # Read the perturbations from the file
            pert_nucs = set()
            init_nuc_conc = {} 

            for iiv in self.env['initial_nuc_keys']:
                iiv_col = self.pert_cols[iiv]

                nuc_zz = nucname.zzaaam(iiv.partition('_')[2])
                pert_nucs.add(nuc_zz)

                init_nuc_conc[nuc_zz] = iiv_col[n]

            # generate a pertubed stream
            pert_stream = Material(init_nuc_conc)

            # generate a non-pertubed stream
            all_nucs = set(self.env['ihm_mat'].comp.keys())
            non_pert_nucs = all_nucs - pert_nucs
            non_pert_stream = self.env['ihm_mat'][non_pert_nucs]
            non_pert_stream.mass = 1.0 - pert_stream.mass

            # generate an initial heavy metal stream
            ihm_mat = pert_stream + non_pert_stream

        self.ihm_mat = ihm_mat

        # Convolve the streams
        atom_frac_fuel = {k: v for k, v in self.env['fuel_chemical_form'].items() if k != "IHM"}
        atom_frac_fuel[ihm_mat] = self.env['fuel_chemical_form'].get("IHM", 0.0)
        self.initial_fuel_stream = from_atom_frac(atom_frac_fuel)
        self.IHM_weight = self.ihm_mat.molecular_weight()
        self.fuel_weight = self.initial_fuel_stream.molecular_weight()
    

    def make_input_fuel(self, ms=None):
        if 'fuel_form_mass_weighted' in self.env:
            mass_weighted = self.env['fuel_form_mass_weighted']
        else:
            mass_weighted = True

        if ms == None:
            comp_str = self.make_input_material_weights(self.initial_fuel_stream.comp, mass_weighted)
        else:
            comp_str = self.make_input_material_weights(ms.mult_by_mass(), True)

        return comp_str


    def make_input_cladding(self):
        # Try to load cladding stream
        if 'clad_form' in self.env:
            clad_stream = Material(self.env['clad_form'])
        else:
            if 0 < self.env['verbosity']:
                print(message("Cladding not found.  Proceeding with standard zircaloy mixture."))
                print()
            clad_form = {
                # Natural Zirconium
                400900: 0.98135 * 0.5145,
                400910: 0.98135 * 0.1122,
                400920: 0.98135 * 0.1715,
                400940: 0.98135 * 0.1738,
                400960: 0.98135 * 0.0280,
                # The plastic is all melted and the natural Chromium too..
                240500: 0.00100 * 0.04345,
                240520: 0.00100 * 0.83789,
                240530: 0.00100 * 0.09501,
                240540: 0.00100 * 0.02365,
                # Natural Iron
                260540: 0.00135 * 0.05845,
                260560: 0.00135 * 0.91754,
                260570: 0.00135 * 0.02119,
                260580: 0.00135 * 0.00282,
                # Natural Nickel
                280580: 0.00055 * 0.68077,
                280600: 0.00055 * 0.26223,
                280610: 0.00055 * 0.01140,
                280620: 0.00055 * 0.03634,
                280640: 0.00055 * 0.00926,
                # Natural Tin
                501120: 0.01450 * 0.0097,
                501140: 0.01450 * 0.0065,
                501150: 0.01450 * 0.0034,
                501160: 0.01450 * 0.1454,
                501170: 0.01450 * 0.0768,
                501180: 0.01450 * 0.2422,
                501190: 0.01450 * 0.0858,
                501200: 0.01450 * 0.3259,
                501220: 0.01450 * 0.0463,
                501240: 0.01450 * 0.0579,
                # We Need Oxygen!
                80160:  0.00125,
                }
            clad_stream = Material(clad_form)
            self.env['clad_form'] = clad_form

        # Try to find if the cladding form is mass or atomic weighted
        if 'clad_form_mass_weighted' in self.env:
            mass_weighted = self.env['clad_form_mass_weighted']
        else:
            mass_weighted = True
        
        return self.make_input_material_weights(clad_stream.comp, mass_weighted)


    def make_input_coolant(self):
        # Try to load coolant stream
        if 'cool_form' in self.env:
            cool_stream = Material(self.env['cool_form'])
        else:
            if 0 < self.env['verbosity']:
                print(message("Coolant not found.  Proceeding with borated light water.\n"))

            MW = (2 * 1.0) + (1 * 16.0) + (0.199 * 550 * 10.0**-6 * 10.0) + (0.801 * 550 * 10.0**-6 * 11.0)

            cool_form = {
                10010: (2 * 1.0) / MW,
                80160: (1 * 16.0) / MW,
                50100: (0.199 * 550 * 10.0**-6 * 10.0) / MW,
                50110: (0.801 * 550 * 10.0**-6 * 11.0) / MW,
                }
            cool_stream = Material(cool_form)
            self.env['cool_form'] = cool_form

        # Try to find if the coolant form is mass or atomic weighted
        if 'cool_form_mass_weighted' in self.env:
            mass_weighted = self.env['cool_form_mass_weighted']
        else:
            mass_weighted = True
        
        return self.make_input_material_weights(cool_stream.comp, mass_weighted)


    def make_input_geometry(self, n):
        # Require
        geom = {
            'fuel_radius': self.pert_cols['fuel_cell_radius'][n],
            'clad_radius': self.pert_cols['clad_cell_radius'][n],
            'void_radius': self.pert_cols['void_cell_radius'][n],
            'cell_pitch':  self.pert_cols['unit_cell_pitch'][n],
            }

        # Tries to get the lattice specification
        # If it isn't present, use a default 17x17 PWR lattice
        if ('lattice' in self.env) and ('lattice_xy' in self.env):
            geom['lattice']    = self.env['lattice']
            geom['lattice_xy'] = self.env['lattice_xy']
        else:
            if 0 < self.env['verbosity']:
                print(message("Lattice specification not found."))
                print(message("Using a default 17x17 PWR lattice."))
                print()
            geom['lattice_xy'] = 17
            geom['lattice']    = ("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" 
                                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" 
                                  "1 1 1 1 1 2 1 1 2 1 1 2 1 1 1 1 1 \n" 
                                  "1 1 1 2 1 1 1 1 1 1 1 1 1 2 1 1 1 \n" 
                                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" 
                                  "1 1 2 1 1 2 1 1 2 1 1 2 1 1 2 1 1 \n" 
                                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" 
                                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" 
                                  "1 1 2 1 1 2 1 1 2 1 1 2 1 1 2 1 1 \n" 
                                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" 
                                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" 
                                  "1 1 2 1 1 2 1 1 2 1 1 2 1 1 2 1 1 \n" 
                                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" 
                                  "1 1 1 2 1 1 1 1 1 1 1 1 1 2 1 1 1 \n" 
                                  "1 1 1 1 1 2 1 1 2 1 1 2 1 1 1 1 1 \n" 
                                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" 
                                  "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")

        # Determine if lattice is symetric
        lat = np.array(geom['lattice'].split(), dtype='int')
        lat = lat.reshape((geom['lattice_xy'], geom['lattice_xy']))
        lat_trans = lat.transpose()
        if (lat == lat_trans).all():    # Condition for symmetry; A = A^T
            geom['sym_flag'] = ''
        else:
            geom['sym_flag'] = '% The lattice is not symmetric! Forced to use whole geometry...\n%'

        # Set half of the lattice pitch.
        half_lat_pitch = (float(geom['lattice_xy']) * geom['cell_pitch']) / 2.0
        geom['half_lattice_pitch'] = "{0}".format(half_lat_pitch)

        return geom


    def make_input_energy_groups(self, group_structure=None):
        """Makes the energy group structure.

        CHAR and most other neutronics codes sepecify this using 
        upper energy bounds.  That way the number of bounds equals the number
        of groups G. Serpent, however, uses only the internal boundaries, so there
        are G-1 energies given for G groups.  Additionally, the highest energy group 
        has the range Bound[G-1] <= Group 1 < inifinity.  The lowest energy group 
        thus covers 0.0 MeV <= Group G < Bound[1].
        """
        if group_structure is None:
            group_structure = self.env['group_structure']

        e = {}
        gs = ["{0:.6G}".format(float(gb)) for gb in group_structure]

        # Set number of (serpent) groups
        e['n_groups'] = len(group_structure) - 1
        e['group_lower_bound'] = gs[0]
        e['group_upper_bound'] = gs[-1]
        e['group_inner_structure'] = "  " + "\n  ".join(gs[1:-1])
        e['group_structure'] = "  " + "\n  ".join(gs)

        return e        


    def make_burnup(self, n):
        """Generates a dictionary of values that fill the burnup portion of the serpent template."""
        # make burnup dictionary
        bu = {'decay_lib': self.env['serpent_decay_lib'],
              'fission_yield_lib': self.env['serpent_fission_yield_lib'],
              'num_burn_regions':  int(self.pert_cols['burn_regions'][n]), 
              'fuel_specific_power': self.pert_cols['fuel_specific_power'][n],
              }

        bu['depletion_times'] = ''
        for bt in self.env['burn_times'][1:]:
            bu['depletion_times'] += '{0:>8.4G}\n'.format(bt)

        bu['transmute_inventory'] = ''
        for nuc in self.env['core_transmute']:
            bu['transmute_inventory'] += '{0:>8}\n'.format(nuc)

        return bu


    def make_detector(self, nuc):
        """Generates a dictionary of values that fill the detector/cross-section portion of 
        the serpent template.  Requires the nuclide to be specified."""
        det = {}

        nuc_zz   = nucname.zzaaam(nuc)
        nuc_serp = nucname.serpent(nuc_zz)

        # Set the nuclide to calculate XS for
        det['xsnuc'] = "{0}.{1}".format(nuc_serp, self.env['temp_flag'])

        if nuc_zz in self.env['cool_form'].keys():
            det['detector_mat'] = 'coolant'
        elif nuc_zz in self.env['clad_form'].keys():
            det['detector_mat'] = 'cladding' 
        else:
            det['detector_mat'] = 'fuel' 

        # Setup detectors to calculate XS for
        det['xsdet'] = ''
        det_format = "det {tally_name} de energies dm {detector_mat} dr {tally_type} xsmat dt 3 phi\n"

        # Get tallies
        tallies = self.env['tallies']
        nuc_mts = self.env['nuc_mts'][nuc_zz]

        # Add tally line if MT number is valid
        for tally in tallies:
            if tallies[tally] in nuc_mts:
                det['xsdet'] += det_format.format(tally_name=tally, 
                                                  tally_type=tallies[tally], 
                                                  detector_mat=det['detector_mat'])

        # Add partial fission MTs, to sum later, if total fission tally is not avilable
        if (tallies['sigma_f'] not in nuc_mts) and (0 < len(partial_fission_mts & nuc_mts)):
            for mt in (partial_fission_mts & nuc_mts):
                det['xsdet'] += det_format.format(tally_name='sigma_f{0}'.format(mt), 
                                                  tally_type=mt, 
                                                  detector_mat=det['detector_mat'])
        return det


    def make_deltam(self, nuc, frac):
        """Generates a dictionary of values that fill the fuel mass stream portion of the 
        serpent template with this nuclide (zzaaam) pertubed to this mass value."""

        ihm_mat = 1.0 * self.ihm_mat

        pert_nuc = set([nuc])
        init_nuc_conc = {nuc: frac} 

        # generate a pertubed stream
        pert_stream = Material(init_nuc_conc)

        # generate a non-pertubed stream
        all_nucs = set(self.ihm_mat.comp.keys())
        non_pert_nucs = all_nucs - pert_nuc
        non_pert_stream = self.ihm_mat[non_pert_nucs]
        non_pert_stream.mass = 1.0 - pert_stream.mass

        # generate an initial heavy metal stream
        ihm_mat = pert_stream + non_pert_stream

        self.ihm_mat = ihm_mat

        # Convolve the streams
        atom_frac_fuel = {k: v for k, v in self.env['fuel_chemical_form'].items() if k != "IHM"}
        atom_frac_fuel[ihm_mat] = self.env['fuel_chemical_form'].get("IHM", 0.0)
        self.initial_fuel_stream = from_atom_frac(atom_frac_fuel)
        self.IHM_weight = self.ihm_mat.molecular_weight()
        self.fuel_weight = self.initial_fuel_stream.molecular_weight()
    

        # make burnup dictionary
        dm = {'fuel': self.make_input_fuel()}

        return dm


    def make_common_input(self, n):
        # Initial serpent fill dictionary
        serpent_fill = {
            'reactor': self.env['reactor'],
            'xsdata':  self.env['serpent_xsdata'],

            'fuel_density': '{0:.5G}'.format(self.pert_cols['fuel_density'][n]),
            'clad_density': '{0:.5G}'.format(self.pert_cols['clad_density'][n]),
            'cool_density': '{0:.5G}'.format(self.pert_cols['cool_density'][n]),

            'k_particles':   self.env['k_particles'],
            'k_cycles':      self.env['k_cycles'],
            'k_cycles_skip': self.env['k_cycles_skip'],
            }

        # Set the material lines
        self.update_initial_fuel(n)
        serpent_fill['fuel']     = self.make_input_fuel()
        serpent_fill['cladding'] = self.make_input_cladding()
        serpent_fill['coolant']  = self.make_input_coolant()

        # Add the geometry information
        serpent_fill.update(self.make_input_geometry(n))

        # Set the energy group structure
        serpent_fill.update(self.make_input_energy_groups())

        # Assign serpent_fill to the class
        self.serpent_fill = serpent_fill


    def make_burnup_input(self, n):
        self.serpent_fill.update(self.make_burnup(n))

        # Fill the burnup template
        with open(self.env['reactor'] + "_burnup_{0}".format(n), 'w') as f:
            f.write(self.env['burnup_template'].format(**self.serpent_fill))


    def make_xs_gen_input(self, nuc="U235", n=0):
        self.serpent_fill.update(self.make_detector(nuc))

        # Fill the XS template
        with open(self.env['reactor'] + "_xs_gen_{0}_{1}".format(nuc, n), 'w') as f:
            f.write(self.env['xs_gen_template'].format(**self.serpent_fill))


    def make_flux_g_input(self, n):
        # Make the flux only calculation with this energy group structure
        self.serpent_fill.update(self.make_input_energy_groups(self.hi_res_group_structure[::-1]))

        # Load the detectors to ensure valid entries,
        # Then wash out mutlipliers that are not just the flux
        self.serpent_fill.update(self.make_detector("U235"))
        self.serpent_fill['xsdet'] = "% Only calculating the flux here!"

        # Fill the XS template
        with open(self.env['reactor'] + "_flux_g_{0}".format(n), 'w') as f:
            f.write(self.env['xs_gen_template'].format(**self.serpent_fill))

        # Restore the detectors and energy groups to their default values
        self.serpent_fill.update(self.make_detector("U235"))
        self.serpent_fill.update(self.make_input_energy_groups())


    def make_deltam_input(self, nuc, n, s, frac):
        """While n indexs the perturbations, nuc is the nuclide (zzaaam) to perturb and 
        frac is the new mass fraction of this nuclide."""
        nuc_LL = nucname.name(nuc)
        self.serpent_fill.update(self.make_burnup(n))
        self.serpent_fill.update(self.make_deltam(nuc, frac[s]))

        # Fill the burnup template
        with open(self.env['reactor'] + "_deltam_{0}_{1}_{2}".format(nuc_LL, n, s), 'w') as f:
            f.write(self.env['burnup_template'].format(**self.serpent_fill))


    def get_mpi_flag(self):
        mpi_flag = ''

        if 'run_parallel' in self.env:
            run_flag = self.env['run_parallel'].upper()
        else:
            run_flag = ''

        if run_flag in ["MPI", "PBS"]:
            if 'number_cpus' in self.env:
                num_cpus = self.env['number_cpus']
            else:
                print(message("The number of cpus was not specified even though a multicore calculation was requested.\n"
                              "Setting the number of cpus to 1 for this calculation."))
                num_cpus = 1

            mpi_flag = '-mpi {0}'.format(self.env['number_cpus'])

        return mpi_flag


    #
    # Serpent run methods
    # 

    def run_script_walltime(self):
        # Set PBS_Walltime
        if 'walltime' in self.env:
            return self.env['walltime']
        else:
            # Get the number of perturbatiosn from the reactor file.
            rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'r')
            nperturbations = len(rx_h5.root.perturbations)
            rx_h5.close()

            # Take a guess as to the number of hours required to run
            return int(4.5 * nperturbations) + 1


    def run_script_fill_values(self):
        """Sets the fill values for running serpent."""

        rsfv = {}

        # Set Transport Job Context
        rsfv['transport_job_context'] = self.run_str + " -version"

        # Set Run_Commands 
        if self.env['options'].LOCAL:
            rsfv['run_commands'] = ''
        else:
            rsfv['run_commands'] = ''
            #rsfv['run_commands'] = 'lamboot\n'

        # Add burnup information
        if self.env['options'].RUN_BURNUP:
            #rsfv['run_commands'] += "{0} {1}_burnup {2}\n".format(self.run_str, self.env['reactor'], self.get_mpi_flag())
            rsfv['run_commands'] += "char --cwd -b defchar.py\n"

        # Add cross section information
        if self.env['options'].RUN_XS_GEN:
            rsfv['run_commands'] += "char --cwd -x defchar.py\n"

        # Add nuclide sensitivity analysis
        if self.env['options'].RUN_DELTAM:
            rsfv['run_commands'] += "char --cwd -m defchar.py\n"

        return rsfv


    def run_burnup_pert(self, n):
        """Runs a burnup perturbation step."""
        # Ensure that the burnup times are at t = 0
        if 0 != n%self.ntimes:
            raise IndexError("Burnups must be started at t = 0 perturbations.")

        self.env['logger'].info('Running burnup calculation at perturbation step {0}.'.format(n))

        # Make new input file
        self.make_common_input(n)
        self.make_burnup_input(n)

        # Run serpent on this input file as a subprocess
        if not self.env['options'].CACHE:
            argstr = "{0}_burnup_{1} {2}".format(self.env['reactor'], n, self.mpi_flag)
            rtn = self.run_serpent(argstr)

        # Parse & write this output to HDF5
        res, dep = self.parse_burnup(n)
        return res, dep



    def run_xs_gen_pert(self, nuc, n, ms_n, E_n=None, E_g=None, phi_n=None):
        """Runs the perturbation for an nuclide that is in serpent.

        nuc : nuclide identifier
        n : perterbation step number.
        ms_n : Mass stream of nuclides in serpent at this step.

        NOTE: This method adds filler fision products.
        If nuc is not zirconium, add Zr-90. If is zirconium, add Sr-90
        These two nuclides have almost the same mass and neutronic profile:
            http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc=Zr-90&n=2
            http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc=Sr-90&n=2
        We need to do this to preseve the atom density of the fuel, 
        while not inducing errors through self-shielding and strong absorbers.
        Basically, Zr-90 and Sr-90 become representative fision product pairs.
        
        WARNING: This is only suppossed to be a first order correction!
        Make sure that you include enough fission products in core_transmute.
        """
        nuc_zz = nucname.zzaaam(nuc)
        nuc_LL = nucname.name(nuc_zz)

        args_xs_gen = "{0}_xs_gen_{1}_{2} {3}".format(self.env['reactor'], nuc_LL, n, self.mpi_flag)

        info_str = 'Generating cross-sections for {0} at perturbation step {1} using serpent.'
        self.env['logger'].info(info_str.format(nuc_LL, n))

        # Add filler fision product
        top_up_mass = 1.0 - ms_n.mass
        if top_up_mass == 0.0:
            top_up = 0.0
        elif nucname.zz_name[nuc_zz//10000] == 'ZR':
            top_up = Material({380900: 90.0, 621480: 148.0}, top_up_mass)
        elif nucname.zz_name[nuc_zz//10000] == 'SM':
            top_up = Material({400900: 90.0, 601480: 148.0}, top_up_mass)
        else:
            top_up = Material({400900: 90.0, 621480: 148.0}, top_up_mass)

        ms = ms_n + top_up

        # Add a differential mass of this nuclide, if not otherwise present
        # Needed to capture sigma_s_gh effects via serpent.
        if (nuc_zz not in ms.comp) or (ms.comp[nuc_zz] == 0.0):
            # Grab the next (or prev) mass stream
            o = n + 1
            if (o == self.nperturbations) or (self.perturbations[o][-1] == 0.0):
                o = n - 1
            ms_o = Material()
            ms_o.from_hdf5(self.env['reactor'] + ".h5", "/Ti0", o, protocol=0)

            # make a new mass weight guess
            mw_nuc = 0.001 * ms_o.comp[nuc_zz]
            if mw_nuc == 0.0:
                mw_nuc = 0.001

            # Change the mass weight for this nuc
            ms_nuc = Material({nuc_zz: mw_nuc})
            ms = ms + ms_nuc

        # Set the final material
        atom_frac_fuel = {k: v for k, v in self.env['fuel_chemical_form'].items() if k != "IHM"}
        atom_frac_fuel[ms] = self.env['fuel_chemical_form'].get("IHM", 0.0)
        ms = from_atom_frac(atom_frac_fuel)

        # Make new input file
        self.make_common_input(n)
        self.serpent_fill['fuel'] = self.make_input_fuel(ms)
        self.make_xs_gen_input(nuc_LL, n)

        # Run serpent on this input file as a subprocess
        if not self.env['options'].CACHE:
            rtn = self.run_serpent(args_xs_gen)

        # Parse this run 
        res, det = self.parse_xs_gen(nuc_LL, n)

        # Prep for metastable tallies
        tallies = self.env['tallies']
        nuc_mts = self.env['nuc_mts'][nuc_zz]

        #
        # Use models, if serpent is not available
        #

        # Get absorption reaction tallies
        for tally in tallies:
            if (tallies[tally] not in nuc_mts) and (tallies[tally] is not None):
                if tally == 'sigma_f':
                    if (0 < len(partial_fission_mts & nuc_mts)):
                        sig_f  = np.zeros(len(E_g) - 1, dtype=float)
                        for mt in (partial_fission_mts & nuc_mts):
                            sig_f += det['DETsigma_f{0}'.format(mt)][:, 10]
                    else:
                        sig_f = pyne.xs.channels.sigma_f(nuc, E_n=E_n, E_g=E_g, phi_n=phi_n)
                    det['_sigma_f'] = sig_f
                elif tally == 'nubar_sigma_f':
                    # Do nubar later, after fission XS have been calculated
                    pass
                elif tally == 'sigma_a':
                    # Do absorption later, after other XS have been calculated
                    pass 
                elif tally == 'sigma_s_gh':
                    det['_sigma_s_gh'] = pyne.xs.channels.sigma_s_gh(nuc, 
                                                                     self.env['temperature'], 
                                                                     E_n=E_n, 
                                                                     E_g=E_g, 
                                                                     phi_n=phi_n)
                elif tally == 'sigma_s':
                    det['_sigma_s'] = pyne.xs.channels.sigma_s(nuc, self.env['temperature'], 
                                                               E_n=E_n, E_g=E_g, phi_n=phi_n)
                else:
                    tally_rx = tally.partition('_')[2]
                    try:
                        det['_'+tally] = pyne.xs.channels.sigma_a_reaction(nuc_zz, tally_rx, 
                                                                           E_n=E_n, E_g=E_g, 
                                                                           phi_n=phi_n)
                    except IndexError:
                        pass

        # Get (n, g *)
        if 'sigma_gamma_x' in tallies:
            if ('sigma_gamma' in tallies) and (tallies['sigma_gamma'] in nuc_mts):
                tot_sig_g = det['DETsigma_gamma'][:, 10]
                ms_rat = pyne.xs.channels.metastable_ratio(nuc_zz, 'gamma', E_g=E_g, 
                                                           E_n=E_n, phi_n=phi_n)
                sig_g = tot_sig_g / (1.0 + ms_rat)
                sig_g_x = ms_rat * sig_g
                det['_sigma_gamma_x'] = sig_g_x
                det['DETsigma_gamma'][:, 10] = sig_g
            else:
                det['_sigma_gamma_x'] = pyne.xs.channels.sigma_a_reaction(nuc_zz, 'gamma_x', 
                                                                          E_n=E_n, E_g=E_g, 
                                                                          phi_n=phi_n)

        # Get (n, 2n *)
        if 'sigma_2n_x' in tallies:
            if ('sigma_2n' in tallies) and (tallies['sigma_2n'] in nuc_mts):
                tot_sig_2n = det['DETsigma_2n'][:, 10]
                ms_rat = pyne.xs.channels.metastable_ratio(nuc_zz, '2n', E_g=E_g, E_n=E_n, 
                                                           phi_n=phi_n)
                sig_2n = tot_sig_2n / (1.0 + ms_rat)
                sig_2n_x = ms_rat * sig_2n
                det['_sigma_2n_x'] = sig_2n_x
                det['DETsigma_2n'][:, 10] = sig_2n
            else:
                det['_sigma_2n_x'] = pyne.xs.channels.sigma_a_reaction(nuc_zz, '2n_x', E_n=E_n, 
                                                                       E_g=E_g, phi_n=phi_n)

        # Get fission XS if MT not available
        if tallies['nubar_sigma_f'] not in nuc_mts:
            if (tallies['sigma_f'] in nuc_mts):
                sig_f = det['DETsigma_f'][:, 10]
            else:
                sig_f = det['_sigma_f']
            det['_nubar_sigma_f'] = 2.5 * sig_f

        # Get absoprtion XS if MT not available
        if tallies['sigma_a'] not in nuc_mts:
            det['_sigma_a'] = np.zeros(len(E_g) - 1, dtype=float)
            for tally in tally_types.sigma_a_tallies:
                if (tallies[tally] in nuc_mts):
                    det['_sigma_a'] += det['DET' + tally][:, 10]
                elif '_'+tally in det:
                    det['_sigma_a'] += det['_' + tally]

            #print(tally + " = " + repr(det['_sigma_a']))

        return res, det




    def run_flux_g_pert(self, n, ms_n):
        """Runs a perturbation of is high-resolution flux.

        n : perterbation step number.
        ms_n : Mass stream of nuclides in serpent at this step.
        """
        self.env['logger'].info("Generating high resolution flux for use with non-serpent models at at perturbation step {0}.".format(n))

        args_flux_g = "{0}_flux_g_{1} {2}".format(self.env['reactor'], n, self.mpi_flag)

        # Make mass stream 
        top_up_mass = 1.0 - ms_n.mass
        if top_up_mass == 0.0:
            top_up = 0.0
        else:
            top_up = Material({400900: 90.0, 621480: 148.0}, top_up_mass)

        ms = ms_n + top_up
        atom_frac_fuel = {k: v for k, v in self.env['fuel_chemical_form'].items() if k != "IHM"}
        atom_frac_fuel[ms] = self.env['fuel_chemical_form'].get("IHM", 0.0)
        ms = from_atom_frac(atom_frac_fuel)

        # Make new input file
        self.make_common_input(n)
        self.serpent_fill['fuel'] = self.make_input_fuel(ms)
        self.make_flux_g_input(n)

        # Run serpent on this input file as a subprocess
        if not self.env['options'].CACHE:
            rtn = self.run_serpent(args_flux_g)

        # Parse & write this output to HDF5
        res, det = self.parse_flux_g(n)

        return res, det



    def run_xs_mod_pert(self, nuc, n, E_n, E_g, phi_n):
        """Generates crosss sections for nuclides not in serpent.

        nuc : nuclide identifier
        n : perterbation step number.
        """
        info_str = 'Generating cross-sections for {0} at perturbation step {1} using models.'
        self.env['logger'].info(info_str.format(nucname.name(nuc), n))

        tallies = self.env['tallies']

        # Load cross-section cahce with proper values
        xs_cache['E_n'] = E_n
        xs_cache['E_g'] = E_g
        xs_cache['phi_n'] = phi_n

        xs_dict = {}

        # Add the cross-section data from models
        if 'sigma_f' in tallies:
            xs_dict['sigma_f'] = pyne.xs.channels.sigma_f(nuc)

        if 'sigma_a' in tallies:
            xs_dict['sigma_a'] = pyne.xs.channels.sigma_a(nuc)

        if 'sigma_s_gh' in tallies:
            xs_dict['sigma_s_gh'] = pyne.xs.channels.sigma_s_gh(nuc, self.env['temperature'])

        if 'sigma_s' in tallies:
            xs_dict['sigma_s'] = pyne.xs.channels.sigma_s(nuc, self.env['temperature'])

        if 'chi' in tallies:
            xs_dict['chi'] = pyne.xs.channels.chi(nuc)

        if 'sigma_gamma' in tallies:
            xs_dict['sigma_gamma'] = pyne.xs.channels.sigma_a_reaction(nuc, 'gamma')

        if 'sigma_2n' in tallies:
            xs_dict['sigma_2n'] = pyne.xs.channels.sigma_a_reaction(nuc, '2n')

        if 'sigma_3n' in tallies:
            xs_dict['sigma_3n'] = pyne.xs.channels.sigma_a_reaction(nuc, '3n')

        if 'sigma_alpha' in tallies:
            xs_dict['sigma_alpha'] = pyne.xs.channels.sigma_a_reaction(nuc, 'alpha')

        if 'sigma_proton' in tallies:
            xs_dict['sigma_proton'] = pyne.xs.channels.sigma_a_reaction(nuc, 'proton')

        if 'sigma_deut' in tallies:
            xs_dict['sigma_deut'] = pyne.xs.channels.sigma_a_reaction(nuc, 'deut')

        if 'sigma_trit' in tallies:
            xs_dict['sigma_trit'] = pyne.xs.channels.sigma_a_reaction(nuc, 'trit')

        if 'sigma_gamma_x' in tallies:
            xs_dict['sigma_gamma_x'] = pyne.xs.channels.sigma_a_reaction(nuc, 'gamma_x')

        if 'sigma_2n_x' in tallies:
            xs_dict['sigma_2n_x'] = pyne.xs.channels.sigma_a_reaction(nuc, '2n_x')

        if 'sigma_t' in tallies:
            xs_dict['sigma_t'] = pyne.xs.channels.sigma_t(nuc, self.env['temperature'])

        return xs_dict



    def run_deltam_pert(self, nuc, n, s, nuc_fracs):
        """Runs a sensitivity pertutbation."""
        # Ensure that pertubrations are at time t = 0
        if 0 != n%self.ntimes:
            raise IndexError("Sensitivities must be started at t = 0 perturbations.")

        nuc_LL = nucname.name(nuc)
        argstr = "{0}_deltam_{1}_{2}_{3} {4}".format(self.env['reactor'], nuc_LL, n, s, self.mpi_flag)

        info_str = 'Running {0} sensitivity study at mass fraction {1} at perturbation step {2}.'.format(nuc, nuc_fracs[s], n)
        self.env['logger'].info(info_str)

        # Make new input file
        self.make_common_input(n)
        self.make_deltam_input(nuc, n, s, nuc_fracs)

        # Run serpent on this input file as a subprocess
        if not self.env['options'].CACHE:
            rtn = self.run_serpent(argstr)

        # Parse & write this output to HDF5
        res, dep = self.parse_deltam(nuc, n, s)

        return res, dep


    #
    # Parsing functions
    #

    def parse_burnup(self, n):
        """Parse the burnup/depletion files into an equivelent python modules."""
        res_file = self.env['reactor'] + "_burnup_{0}_res".format(n)
        dep_file = self.env['reactor'] + "_burnup_{0}_dep".format(n)

        # Convert files
        res = pyne.serpent.parse_res(res_file + ".m")
        dep = pyne.serpent.parse_dep(dep_file + ".m")

        return res, dep


    def parse_xs_gen(self, nuc, n):
        """Parse the cross-section generation files into an equivelent python modules."""
        res_file = self.env['reactor'] + "_xs_gen_{0}_{1}_res".format(nuc, n)
        det_file = self.env['reactor'] + "_xs_gen_{0}_{1}_det0".format(nuc, n)

        # Convert files
        res = pyne.serpent.parse_res(res_file + ".m")
        det = pyne.serpent.parse_det(det_file + '.m')

        return res, det


    def parse_flux_g(self, n):
        """Parse the group flux only generation files into an equivelent python modules."""
        res_file = self.env['reactor'] + "_flux_g_{0}_res".format(n)
        det_file = self.env['reactor'] + "_flux_g_{0}_det0".format(n)

        # Convert files
        res = pyne.serpent.parse_res(res_file + ".m")
        det = pyne.serpent.parse_det(det_file + '.m')

        return res, det


    def parse_deltam(self, nuc, n, s):
        """Parse the sensitivity study files into an equivelent python modules."""
        res_file = self.env['reactor'] + "_deltam_{0}_{1}_{2}_res".format(nuc, n, s)
        dep_file = self.env['reactor'] + "_deltam_{0}_{1}_{2}_dep".format(nuc, n, s)

        # Convert files
        res = pyne.serpent.parse_res(res_file + ".m")
        dep = pyne.serpent.parse_dep(dep_file + ".m")

        return res, dep


    #
    # Init HDF5 groups and arrays
    #

    def init_h5(self):
        """Initialize a new HDF5 file in preparation for burnup and XS runs."""
        # Setup tables 
        desc = {}
        for n, param in enumerate(self.env['perturbation_params']):
            desc[param] = tb.Float64Col(pos=n)

        # Open HDF5 file
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'w')
        base_group = "/"

        # Add pertubation table index
        p = rx_h5.createTable(base_group, 'perturbations', desc)
        p.append(self.perturbations)

        # Add the nuclide tracking arrays.  
        rx_h5.createArray(base_group, 'nucstrack', np.array(self.env['core_transmute']), 
                          "Nuclides to track, copy of transmute_nucs_zz")

        rx_h5.createArray(base_group, 'load_nucs_zz', np.array(self.env['core_load']), 
                          "Core loading nuclides [zzaaam]")
        rx_h5.createArray(base_group, 'load_nucs_LL', np.array(\
                          [nucname.name(nuc) for nuc in self.env['core_load']]), 
                          "Core loading nuclides [LLAAAM]")

        rx_h5.createArray(base_group, 'transmute_nucs_zz', np.array(self.env['core_transmute']), 
                          "Core transmute nuclides [zzaaam]")
        rx_h5.createArray(base_group, 'transmute_nucs_LL', np.array(\
                          [nucname.name(nuc) for nuc in self.env['core_transmute']]), 
                          "Core transmute nuclides [LLAAAM]")

        # Close HDF5 file
        rx_h5.close()


    def init_array(self, rx_h5, base_group, array_name, init_array, array_string='Helpful Array is Helpful'):
        """Inits an array in an hdf5 file."""

        # Remove existing array
        if hasattr(rx_h5.getNode(base_group), array_name):
            rx_h5.removeNode(base_group, array_name, recursive=True)

        # Create new array
        rx_h5.createArray(base_group, array_name, init_array, array_string)


    def init_tally_group(self, rx_h5, base_group, tally, init_array, 
                         gstring='Group {tally}', astring='Array {tally} {nuc}'):
        """Inits a tally group in an hdf5 file."""

        # Remove existing group, create new group of same name.
        if hasattr(rx_h5.getNode(base_group), tally):
            rx_h5.removeNode(base_group, tally, recursive=True)
        tally_group = rx_h5.createGroup(base_group, tally, gstring.format(tally=tally))

        # Add nuclide arrays for this tally to this group.
        for nuc in self.env['core_transmute']:
            nuc_LL = nucname.name(nuc)
            rx_h5.createArray(tally_group, nuc_LL, init_array, astring.format(tally=tally, nuc=nuc_LL))


    def init_h5_burnup(self):
        """Initialize the hdf5 file for a set of burnup calculations based on the its input params."""
        # Open a new hdf5 file 
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'a')
        base_group = "/"

        # Number of energy groups
        G = self.G

        # Init the raw tally arrays
        neg1 = -1.0 * np.ones( (self.nperturbations, ) )
        negG = -1.0 * np.ones( (self.nperturbations, G) )
        negE = -1.0 * np.ones( (self.nperturbations, G+1) )


        # Add basic BU information
        self.init_array(rx_h5, base_group, 'BU0',   neg1, "Burnup of the initial core loading [MWd/kg]")
        self.init_array(rx_h5, base_group, 'time0', neg1, "Time after initial core loading [days]")

        # Add flux arrays
        self.init_array(rx_h5, base_group, 'phi',   neg1, "Total flux [n/cm2/s]")
        self.init_array(rx_h5, base_group, 'phi_g', negG, "Group fluxes [n/cm2/s]")

        # Create Fluence array
        self.init_array(rx_h5, base_group, 'Phi', neg1, "Fluence [n/kb]")

        # Energy Group bounds
        self.init_array(rx_h5, base_group, 'energy', negE, "Energy boundaries [MeV]")

        # Initialize transmutation matrix
        self.init_tally_group(rx_h5, base_group, 'Ti0', neg1, 
                              "Transmutation matrix from initial core loading [kg_i/kgIHM]", 
                              "Mass weight of {nuc} [kg/kgIHM]")

        self.init_array(rx_h5, base_group + '/Ti0', 'Mass', neg1, "Mass fraction of fuel [kg/kgIHM]")

        # close the file before returning
        rx_h5.close()


    def init_h5_xs_gen(self):
        """Initialize the hdf5 file for writing for the XS Gen stage.
        The shape of these arrays is dependent on the number of time steps."""
        # Open a new hdf5 file 
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'a')
        base_group = "/"

        # Number of energy groups
        G = self.G

        # Grab the tallies
        tallies = self.env['tallies']

        # Init the raw tally arrays
        neg1 = -1.0 * np.ones( (self.nperturbations, G) )
        negG = -1.0 * np.ones( (self.nperturbations, G, G) )

        for tally in tallies:
            if tally is not None:
                self.init_tally_group(rx_h5, base_group, tally, neg1, 
                                      "Microscopic Cross Section {tally} [barns]", 
                                      "Microscopic Cross Section {tally} for {nuc} [barns]")

        # Init aggregate tallies

        # nubar
        if ('sigma_f' in tallies) and ('nubar_sigma_f' in tallies) and ('nubar' not in tallies):
            self.init_tally_group(rx_h5, base_group, 'nubar', neg1, 
                                  "{tally} [unitless]", "{tally} for {nuc} [unitless]")

        # sigma_i
        if np.array(['sigma_i' in tally  for tally in tallies]).any() and ('sigma_i' not in tallies):
            self.init_tally_group(rx_h5, base_group, 'sigma_i', neg1, 
                                  "Microscopic Cross Section {tally} [barns]", 
                                  "Microscopic Cross Section {tally} for {nuc} [barns]")

        # sigma_s
        if ('sigma_s' not in tallies):
            self.init_tally_group(rx_h5, base_group, 'sigma_s', neg1, 
                                  "Microscopic Cross Section {tally} [barns]", 
                                  "Microscopic Cross Section {tally} for {nuc} [barns]")

        # Scattering kernel, sigma_s_gh
        self.init_tally_group(rx_h5, base_group, 'sigma_s_gh', negG, 
                              "Microscopic Scattering Kernel {tally} [barns]", 
                              "Microscopic Scattering Kernel {tally} for {nuc} [barns]")

        # Chi
        if 'chi' in tallies:
            self.init_tally_group(rx_h5, base_group, 'chi', neg1, 
                                  "Fission neutron energy spectrum {tally} [n/MeV]", 
                                  "Fission neutron energy spectrum {tally} for {nuc} [n/MeV]")

        # (n, gamma *)
        if 'sigma_gamma_x' in tallies:
            self.init_tally_group(rx_h5, base_group, 'sigma_gamma_x', neg1, 
                                  "Microscopic Cross Section (n, gamma *) [barns]", 
                                  "Microscopic Cross Section (n, gamma *) for {nuc} [barns]")

        # (n, 2n *)
        if 'sigma_2n_x' in tallies:
            self.init_tally_group(rx_h5, base_group, 'sigma_2n_x', neg1, 
                                  "Microscopic Cross Section (n, 2n *) [barns]", 
                                  "Microscopic Cross Section (n, 2n *) for {nuc} [barns]")


        # close the file before returning
        rx_h5.close()


    def init_h5_flux_g(self):
        """Initialize the hdf5 file for a set of high-resoultion flux calculations."""
        # Open a new hdf5 file 
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'a')
        base_group = "/"

        # Remove existing group, create new group of same name.
        if hasattr(rx_h5.getNode(base_group), "hi_res"):
            rx_h5.removeNode(base_group, "hi_res", recursive=True)
        hi_res_group = rx_h5.createGroup(base_group, "hi_res", "High Resolution Group Fluxes")

        # Get default group structure and number of energy groups
        group_structure = self.hi_res_group_structure
        G = len(group_structure) - 1

        # Init the raw tally arrays
        neg1 = -1.0 * np.ones( (self.nperturbations, ) )
        negG = -1.0 * np.ones( (self.nperturbations, G) )
        negE = -1.0 * np.ones( (self.nperturbations, G+1) )

        # Add flux arrays
        self.init_array(rx_h5, hi_res_group, 'phi',   neg1, "High-Resolution Total flux [n/cm2/s]")
        self.init_array(rx_h5, hi_res_group, 'phi_g', negG, "High-Resolution Group fluxes [n/cm2/s]")

        # Energy Group bounds
        self.init_array(rx_h5, hi_res_group, 'energy', group_structure, "High-Resolution Energy Boundaries [MeV]")

        # close the file before returning
        rx_h5.close()


    def init_h5_deltam(self):
        """Initialize the hdf5 file for nuclide sensitivity study."""
        deltam_desc = {
            'nuc_LL': tb.StringCol(6, pos=0),
            'nuc_zz': tb.Int32Col(pos=1),
                         
            'perturbation': tb.Int16Col(pos=2), 
            'ihm_mass_fraction': tb.Float64Col(pos=3),

            'reactivity': tb.Float64Col(shape=(self.ntimes, ), pos=4), 
            }

        # Open a new hdf5 file 
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'a')
        base_group = "/"

        # Remove existing group, create new group of same name.
        if hasattr(rx_h5.getNode(base_group), "nuclide_sensitivity"):
            rx_h5.removeNode(base_group, "nuclide_sensitivity", recursive=True)
        deltam_table = rx_h5.createTable(base_group, "nuclide_sensitivity", deltam_desc, "Nuclide Sensitivity")

        # close the file before returning
        rx_h5.close()


    #
    # Writing functions
    #

    def write_burnup(self, n, res, dep):
        """Writes the results of the burnup calculation to an hdf5 file.

        n : Perturbation index of first time step for this burnup calculation.
        res : A dictionary containing the results of the res file. 
        dep : A dictionary containing the results of the dep file. 
        """

        # Find end index 
        t = n + len(dep['DAYS'])

        # Open a new hdf5 file 
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'a')
        base_group = rx_h5.root

        # Grab pertubation columns
        pert_cols = rx_h5.root.perturbations.cols

        # Add basic BU information
        base_group.BU0[n:t] =  dep['BU']
        base_group.time0[n:t] = dep['DAYS']

        phi = res['TOT_FLUX'][:, ::2].flatten()
        base_group.phi[n:t] = phi
        base_group.phi_g[n:t] = res['FLUX'][:,::2][:, 1:]

        # Create Fluence array
        Phi = np.zeros(len(phi))
        Phi[1:] = cumtrapz(phi * (10.0**-21), dep['DAYS'] * (3600.0 * 24.0))
        base_group.Phi[n:t] = Phi

        # Energy Group bounds 
        base_group.energy[n:t] = res['GC_BOUNDS']

        # Calculate and store weight percents per kg fuel form (ie not per IHM)
        # Serepent masses somehow unnormalize themselves in all of these conversions, which is annoying.
        # This effect is of order 1E-5, which is large enough to be noticable.
        # Thus we have to go through two bouts of normalization here.
        mw_conversion = self.fuel_weight / (self.IHM_weight * dep['TOT_VOLUME'] * pert_cols.fuel_density[n])
        mw = dep['TOT_MASS'] * mw_conversion 

        nuc_LL = {}
        nuc_index = {}
        for nuc_zz in dep['ZAI']:
            # Find valid nuclide indeces
            try: 
                nuc_LL[nuc_zz] = nucname.name(int(nuc_zz))
            except:
                continue
            nuc_index[nuc_zz] = dep['i{0}'.format(nuc_zz)] - 1

        # Caclulate actual mass of nuclides present
        mass = mw[nuc_index.values()].sum(axis=0)   

        # Store normalized mass vector for each nuclide
        for nuc_zz in nuc_index:
            nuc_array = getattr(base_group.Ti0, nuc_LL[nuc_zz])

            mw_i =  mw[nuc_index[nuc_zz]] / mass[0]
            nuc_array[n:t] = mw_i

        # Renormalize mass
        mass = mass / mass[0]   
        base_group.Ti0.Mass[n:t] = mass

        # close the file before returning
        rx_h5.close()


    def write_xs_gen(self, nuc, n, res, det):
        # Convert nucname
        nuc_zz = nucname.zzaaam(nuc)
        nuc_LL = nucname.name(nuc_zz)

        nuc_is_fissionable = (86 <= nuc_zz/10000)

        # Add current working directory to path
        if sys.path[0] != os.getcwd():
            sys.path.insert(0, os.getcwd())

        # Open a new hdf5 file 
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'a')
        base_group = rx_h5.root

        # Grab the tallies
        tallies = self.env['tallies']

        # Write the raw tally arrays for this time and this nuc        
        for tally in tallies:
            tally_hdf5_group = getattr(base_group, tally)
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)

            # Make sure the detector was calculated, 
            # Or replace the tally with zeros
            if tallies[tally] in self.env['nuc_mts'][nuc_zz]:
                tally_serp_array = det['DET{0}'.format(tally)]
                tally_serp_array = tally_serp_array[::-1, 10]
            elif ('_'+tally in det):
                tally_serp_array = det['_'+tally][::-1]
            else:
                tally_serp_array = np.zeros(len(tally_hdf5_array[n]), dtype=float)

            # Zero out non-fissionable nuclide tallies
            if (tally in set(['sigma_f', 'nubar_sigma_f', 'nubar', 'chi'])) and (not nuc_is_fissionable):
                tally_serp_array = np.zeros(len(tally_hdf5_array[n]), dtype=float)

            # Make sure there are no NaNs
            mask = np.isnan(tally_serp_array)
            tally_serp_array[mask] = 0.0

            tally_hdf5_array[n] = tally_serp_array

        #
        # Write aggregate or otherwise special tallies
        #
        # nubar
        if ('sigma_f' in tallies) and ('nubar_sigma_f' in tallies) and ('nubar' not in tallies):
            tally_hdf5_group = getattr(base_group, 'nubar')
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)

            # Get fission XS
            if 'DETsigma_f' in det:
                sigma_f = det['DETsigma_f'][::-1, 10]
            elif '_sigma_f' in det:
                sigma_f = det['_sigma_f'][::-1]

            if 'DETnubar_sigma_f' in det:
                nubar_sigma_f = det['DETnubar_sigma_f'][::-1, 10]
            elif '_nubar_sigma_f' in det:
                nubar_sigma_f = det['_nubar_sigma_f'][::-1]

            if nuc_is_fissionable:
                nubar = nubar_sigma_f / sigma_f
            else:
                nubar = np.zeros(len(tally_hdf5_array[n]), dtype=float)            

            # Ensure nubar is well-formed            
            m = ((0.0 <= nubar) != True)
            nubar[m] = 0.0

            tally_hdf5_array[n] = nubar

        # sigma_i
        if ('sigma_i' not in tallies):
            tally_hdf5_group = getattr(base_group, 'sigma_i')
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)

            sigma_i = np.zeros(len(tally_hdf5_array[n]), dtype=float)

            # Sum all sigma_iN's together
            for tally in tallies:
                # Skip this tally when appropriate
                if 'sigma_i' not in tally:
                    continue

                if tallies[tally] not in self.env['nuc_mts'][nuc_zz]:
                    continue

                # Grab a sigma_iN array
                tally_serp_array = det['DET{0}'.format(tally)]

                # Add this array to the current sigma_i 
                sigma_i += tally_serp_array[::-1, 10]

            tally_hdf5_array[n] = sigma_i

        # sigma_s
        sigma_s = None
        if ('sigma_s' not in tallies):
            tally_hdf5_group = getattr(base_group, 'sigma_s')
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)

            if 'sigma_e' in tallies:
                sigma_e = det['DETsigma_e']
                sigma_e = sigma_e[::-1, 10]
            else:
                sigma_e = None

            if (sigma_i == None) and ('sigma_i' in tallies):
                sigma_i = det['DETsigma_i']
                sigma_i = sigma_e[::-1, 10]

            if (sigma_e == None) and (sigma_i == None):
                sigma_s = None
            elif (sigma_e == None) and (sigma_i != None):
                sigma_s = sigma_i
            elif (sigma_e != None) and (sigma_i == None):
                sigma_s = sigma_e
            else:
                sigma_s = sigma_e + sigma_i

            if sigma_s != None:
                tally_hdf5_array[n] = sigma_s

        # Scattering kernel, sigma_s_gh
        if sigma_s != None:
            tally_hdf5_group = getattr(base_group, 'sigma_s_gh')
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)

            gtp = res['GTRANSFP'][res['idx']][::2]
            G = len(sigma_s)
            gtp = gtp.reshape((G, G))

            # Fix mal-formed group transfer probs.
            cond = np.where(gtp.sum(axis=0) == 0.0)
            gtp[cond, cond] = 1.0

            sigma_s_gh = sigma_s * gtp

            tally_hdf5_array[n] = sigma_s_gh

        # chi
        if ('chi' in tallies):
            tally_hdf5_group = getattr(base_group, 'chi')
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)

            if nuc_is_fissionable:
                chi = res['CHI'][res['idx']][::2]
            else:
                chi = np.zeros(len(tally_hdf5_array[n]), dtype=float)                

            tally_hdf5_array[n] = chi

        # close the file before returning
        rx_h5.close()


    def write_flux_g(self, n, res, det):
        # Open a new hdf5 file 
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'a')
        base_group = rx_h5.root.hi_res

        # Grab the HDF5 arrays
        phi_hdf5_array = base_group.phi   
        phi_g_hdf5_array = base_group.phi_g

        # Grab the Serepent arrays
        phi_g_serp_array = det['DETphi']
        phi_g_serp_array = phi_g_serp_array[::-1, 10]

        phi_serp = phi_g_serp_array.sum()

        # Write the flux tallies
        phi_hdf5_array[n] = phi_serp
        phi_g_hdf5_array[n] = phi_g_serp_array

        # close the file before returning
        rx_h5.close()


    def write_xs_mod(self, nuc, n, xs_dict):
        """Writes cross sections from models."""

        # Convert nucname
        nuc_zz = nucname.zzaaam(nuc)
        nuc_LL = nucname.name(nuc_zz)

        nuc_is_fissionable = (86 <= nuc_zz/10000)

        # Open a new hdf5 file 
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'a')
        base_group = rx_h5.root

        # Grab the tallies
        tallies = self.env['tallies']

        # Write the tallies that were not calculated from models
        for tally in tallies:
            if tally in xs_dict:
                continue

            tally_hdf5_group = getattr(base_group, tally)
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)

            tally_hdf5_array[n] = np.zeros(len(tally_hdf5_array[n]), dtype=float)


        # Write the raw tally arrays for this time and this nuc        
        for tally in xs_dict:
            if tally not in tallies:
                continue

            tally_hdf5_group = getattr(base_group, tally)
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)

            tally_model_array = np.array(xs_dict[tally][::-1])

            # Make sure there are no NaNs
            mask = np.isnan(tally_model_array)
            tally_model_array[mask] = 0.0

            tally_hdf5_array[n] = tally_model_array

    
        #
        # Write special tallies
        #

        # nubar
        if ('sigma_f' in tallies) and ('sigma_f' in xs_dict) and  ('nubar_sigma_f' in tallies):
            # set the value of the number of neutrons per fission, based on 
            # whether or not the species actually fissions.
            if (xs_dict['sigma_f'] == 0.0).all():
                nubar = 0.0
            else:
                # Average value for thermal U-238 / Pu-239 systems
                # best guess for now.
                nubar = 2.871 

            # Write nubar
            nubar_array = nubar * np.ones(len(xs_dict['sigma_f']), dtype=float)

            tally_hdf5_group = getattr(base_group, 'nubar')
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)
            
            tally_hdf5_array[n] = nubar_array

            # Write nubar * sigma_f
            nubar_sigma_f_array = nubar * xs_dict['sigma_f'][::-1]

            tally_hdf5_group = getattr(base_group, 'nubar_sigma_f')
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)
            
            tally_hdf5_array[n] = nubar_sigma_f_array

        # sigma_s
        if ('sigma_s' in xs_dict):
            sigma_s = xs_dict['sigma_s'][::-1]
        elif ('sigma_s_gh' in xs_dict):
            sigma_s = xs_dict['sigma_s_gh'][::-1, ::-1].sum(axis=1)
        else:
            sigma_s = None

        tally_hdf5_group = getattr(base_group, 'sigma_s')
        tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)

        if sigma_s is not None:
            tally_hdf5_array[n] = sigma_s

        # Scattering kernel, sigma_s_gh
        if ('sigma_s_gh' in xs_dict):
            tally_hdf5_group = getattr(base_group, 'sigma_s_gh')
            tally_hdf5_array = getattr(tally_hdf5_group, nuc_LL)

            sigma_s_gh = xs_dict['sigma_s_gh'][::-1, ::-1].T

            tally_hdf5_array[n] = sigma_s_gh

        # close the file before returning
        rx_h5.close()


    def write_deltam(self, nuc, n, s, frac, res, dep):
        """Writes the results of a nuclide sensitivity study run to the hdf5 file.

        n : Perturbation index of first time step for this burnup calculation.
        nuc : The nuclide name in zzaaam form.
        frac : The mass fraction of the IHM of this nuctopr.
        """
        nuc_zz = nucname.zzaaam(nuc)
        nuc_LL = nucname.name(nuc_zz)

        # Open the hdf5 file 
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'a')
        base_group = rx_h5.root

        # Calculate the effectiv ereactivity
        keff = res['SIX_FF_KEFF'][:, ::2].flatten()
        rho = (keff - 1.0) / keff

        # Store this row
        nuc_sense_table = base_group.nuclide_sensitivity
        nuc_row = nuc_sense_table.row

        nuc_row['nuc_LL'] = nuc_LL
        nuc_row['nuc_zz'] = nuc_zz
                         
        nuc_row['perturbation'] = n
        nuc_row['ihm_mass_fraction'] = frac[s]

        nuc_row['reactivity'] = rho

        # Write this table out
        nuc_row.append()
        nuc_sense_table.flush()

        # close the file before returning
        rx_h5.close()


    #
    # Analysis functions
    # 

    def analyze_deltam(self):
        """Analyzes the results of the nuclide sensitivity study, producing a report to stdout."""

        # Open the hdf5 file 
        rx_h5 = tb.openFile(self.env['reactor'] + ".h5", 'a')
        base_group = rx_h5.root

        # Grab the appropriate table
        nuc_sense_table = base_group.nuclide_sensitivity

        # Get unique set of nuclides and perturbation numbers.
        nucs = np.unique(nuc_sense_table.cols.nuc_zz)
        ns = np.unique(nuc_sense_table.cols.perturbation)

        # Loop over all nucs and perturbations.
        rho_std_nuc = []
        for nuc in nucs:
            nuc_std = []
            for n in ns:
                where_cond = '(perturbation == {n}) & (nuc_zz == {nuc})'.format(n=n, nuc=nuc)
                rhos = np.array([row['reactivity'] for row in nuc_sense_table.where(where_cond)])

                # Take the standard deviation for each burn step, along this perturbations
                rho_std = np.std(rhos, axis=0)

                # Pick out the standard dev that is the highest among the burnstep
                nuc_std.append(np.max(rho_std))

            # Pick out the standard dev that is highest among all perturbations.
            # Couple that with this nuclide.
            rho_std_nuc.append((np.max(nuc_std), nuc))

        # Sort these reactivity standard deviations highest to lowest.
        rho_std_nuc.sort(reverse=True)

        print(message("Maximum standard deviation for nuclide sensitivity study:"))
        for sig, nuc_zz in rho_std_nuc:
            nuc_LL = nucname.name(int(nuc_zz))
            print("{0:<8}{1}".format(nuc_LL, sig))

        # close the file before returning
        rx_h5.close()
