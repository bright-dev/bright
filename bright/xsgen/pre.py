from __future__ import print_function
import re
import os
import sys
import shutil
import subprocess
from itertools import product
from collections import namedtuple

import numpy as np
from pyne import nucname
from pyne.material import Material, from_atom_frac
from pyne.utils import failure

import bright.xsgen.utils as utils
from bright.xsgen.utils import NotSpecified
from bright.xsgen.plugins import Plugin
from bright.xsgen.openmc_origen import OpenMCOrigen
from bright.xsgen.tally_types import restricted_tallies
from bright.xsgen.brightlite import BrightliteWriter

if sys.version_info[0] > 2:
    basestring = str

INITIAL_NUC_RE = re.compile('initial_([A-Za-z]{0,2}\d{1,7}[Mm]?)')

SOLVER_ENGINES = {'openmc+origen': OpenMCOrigen}

FORMAT_WRITERS = {
    'brightlite': BrightliteWriter,
    }

ensure_mat = lambda m: m is isinstance(m, Material) else Material(m)

def run_ui():
    """Runs the cross section user interface."""
    # Test to see if ui library is installed
    try:
        from bright.xsgen.ui import app
    except ImportError:
        sys.exit(failure("Please install the Enthought Tool Suite (ETS) for CHAR UI."))
    # Open UI
    application = app.Application()
    application.configure_traits()
    # Clean-up UI
    if application.rx_h5 is not None:
        application.rx_h5.close()
    sys.exit()

class XSGenPlugin(Plugin):

    requires = ('bright.xsgen.base',)

    defaultrc = {'solver': NotSpecified,
                 'formats': ('h5',),
                 'ui': False,
                 }

    rcdocs = { 
        'solver': ('The physics codes that are used to solve the '
                   'burnup-criticality problem and compute cross sections and '
                   'transmutation matrices.'),
        'formats': 'The output formats to write out.',
        }

    def update_argparser(self, parser):
        parser.add_argument("--ui", action="store_true", dest="ui",
            help="Launches the char ui.")
        parser.add_argument("-c", "--clean", action="store_true", dest="clean",
            help="Cleans the reactor directory of current files.")
        parser.add_argument('--solver', dest='solver', help=self.rcdocs['solver'])
        parser.add_argument('--formats', dest='formats', help=self.rcdocs['formats'], 
                            nargs='+')

    def setup(self, rc):
        if rc.ui:
            run_ui()

        self.ensure_rc(rc)
        self.make_states(rc)

        if rc.solver is NotSpecified:
            raise ValueError('a solver type must be specified')
        rc.engine = SOLVER_ENGINES[rc.solver](rc)

        rc.writers = [FORMAT_WRITERS[format](rc) for format in rc.formats]

    def ensure_rc(self, rc):
        self._ensure_bt(rc)
        self._ensure_nl(rc)
        self._ensure_temp(rc)
        self._ensure_smf(rc)
        self._ensure_av(rc)
        self._ensure_inp(rc)
        self._ensure_pp(rc)
        self._ensure_mats(rc)

    def _ensure_bt(self, rc):
        # Make Time Steps
        if 'burn_times' in rc:
            rc.burn_times = np.asarray(rc.burn_times, dtype=float)
        else:
            bt_upper_lim = rc.burn_time + rc.time_step/10.0
            rc.burn_times = np.arange(0, bt_upper_lim, rc.time_step)
        rc.burn_times_index = list(range(len(rc.burn_times)))

    def _ensure_nl(self, rc):
        # Make nuclide lists
        if isinstance(rc.core_load_nucs, basestring):
            core_load = load_nuc_file(rc.core_load_nucs)
        else:
            core_load = [nucname.zzaaam(nuc) for nuc in rc.core_load_nucs]
        rc.core_load = sorted(set(core_load))

        if isinstance(rc.core_transmute_nucs, basestring):
            core_transmute = load_nuc_file(rc.core_transmute_nucs)
        else:
            core_transmute = [nucname.zzaaam(nuc) for nuc in rc.core_transmute_nucs]
        rc.core_transmute = sorted(set(core_transmute))

    def _ensure_temp(self, rc):
        # Make temperature
        rc.temperature = rc.get('temperature', 600)

    def _ensure_smf(self, rc):
        # make sensitivity mass fractions
        if 'sensitivity_mass_fractions' in rc:
            rc.deltam = np.atleast_1d(rc.sensitivity_mass_fractions)
            rc.deltam.sort()

    def _ensure_av(self, rc):
        # Make arrays out of quatities that are allowed to vary.
        rc.fuel_density = np.atleast_1d(rc.fuel_density)
        rc.clad_density = np.atleast_1d(rc.clad_density)
        rc.cool_density = np.atleast_1d(rc.cool_density)
        rc.fuel_cell_radius = np.atleast_1d(rc.fuel_cell_radius)
        rc.void_cell_radius = np.atleast_1d(rc.void_cell_radius)
        rc.clad_cell_radius = np.atleast_1d(rc.clad_cell_radius)
        rc.unit_cell_pitch = np.atleast_1d(rc.unit_cell_pitch)
        rc.burn_regions = np.atleast_1d(rc.burn_regions)
        rc.fuel_specific_power = np.atleast_1d(rc.fuel_specific_power)

    def _ensure_inp(self, rc):
        # Grab initial nuc perturbation
        max_mass = 0.0
        initial_nuc_keys = []
        for key in rc:
            m = INITIAL_NUC_RE.match(key)
            if m is None:
                continue

            rc_initial_nuc = getattr(rc, key)
            rc_initial_nuc = np.atleast_1d(rc_initial_nuc)
            setattr(rc, key, rc_initial_nuc)

            initial_nuc_keys.append(key)
            max_mass += np.max(rc_initial_nuc)

        initial_nuc_keys.sort()
        rc.initial_nuc_keys = initial_nuc_keys

        if 1.0 < max_mass:
            msg = "The maxium mass of initial heavy metal perturbations exceeds 1.0 kg!"
            sys.exit(failure(msg))

    def _ensure_pp(self, rc):
        # Set up tuple of parameters to perform a burnup step for
        rc.perturbation_params = ['fuel_density', 'clad_density', 'cool_density',
            'fuel_cell_radius', 'void_cell_radius', 'clad_cell_radius',
            'unit_cell_pitch', 'burn_regions', 'fuel_specific_power',]

        rc.perturbation_params.extend(rc.initial_nuc_keys)
        # burn_times needs to be the last element
        rc.perturbation_params.append('burn_times')

    def _ensure_mats(self, rc):
        if 'fuel_material'in rc:
            rc.fuel_material = ensure_mat(rc.fuel_material)
        elif 'fuel_chemical_form' in rc and 'initial_heavy_metal' in rc:
            ihm_mat = Material(rc.initial_heavy_metal)
            atom_frac = {k: v for k, v in rc.fuel_chemical_form.items() if k != "IHM"}
            atom_frac[ihm_mat] = rc.fuel_chemical_form.get("IHM", 0.0)
            rc.fuel_material = from_atom_frac(atom_frac)
        else:
            raise ValueError("Please specify a fuel.")

        if 'clad_material' in rc: 
            rc.clad_material = ensure_mat(rc.clad_material)
        else: 
            rc.clad_material = Material({
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
                })

        if 'cool_material' in rc: 
            rc.cool_material = ensure_mat(rc.cool_material)
        else: 
            MW = (2 * 1.0) + (1 * 16.0) + (0.199 * 550 * 10.0**-6 * 10.0) + \
                                          (0.801 * 550 * 10.0**-6 * 11.0)
            rc.cool_material = Material({{
                10010: (2 * 1.0) / MW,
                80160: (1 * 16.0) / MW,
                50100: (0.199 * 550 * 10.0**-6 * 10.0) / MW,
                50110: (0.801 * 550 * 10.0**-6 * 11.0) / MW,
                })

    def make_states(self, rc):
        """Makes the reactor state table."""
        State = rc.State = namedtuple('State', rc.perturbation_params)
        data = [getattr(rc, a) for a in rc.perturbation_params]
        rc.states = [State(*p) for p in product(*data)]
        rc.nstates = len(rc.states)

