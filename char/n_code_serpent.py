"""A class to setup, run, and parse Serpent."""
from __future__ import print_function
import os
import sys
import subprocess
from itertools import product

import isoname
import numpy as np
import tables as tb
import metasci.nuke as msn
from metasci import SafeRemove
from MassStream import MassStream
from metasci.colortext import message, failure

from scipy.integrate import cumtrapz

import tally_types
from char import defchar
from n_code import NCode
from m2py import convert_res, convert_dep, convert_det


def zzaaam_2_serpent(iso):
    """Makes an isotope in serpent form."""
    if 0 == iso%10:
        iso_serp = iso/10
    else:
        iso_zz = iso/10000
        iso_LL = isoname.zzLL[iso_zz]
        iso_LL = iso_LL.capitalize()
        iso_aaa = (iso/10)%1000
        iso_serp = "{0}-{1}m".format(iso_LL, iso_aaa)

    return iso_serp


class NCodeSerpent(NCode):
    """A Serpent neutronics code wrapper class."""

    def __init__(self):
        # Reload defchar, just in case.
        global defchar
        from char import defchar

        self.name    = "Serpent"
        self.run_str = "sss-dev"

        if hasattr(defchar, 'ISO_FLAG'):
            self.iso_flag = defchar.ISO_FLAG
        else:
            self.iso_flag = ''

        # Remote file lists
        self.place_remote_files = ['.']
        self.fetch_remote_files = ['.']

    def make_input_material_weights(self, comp, mass_weighted=True):
        """This function takes an isotopic vector, comp, and returns a serpent string representation.
        Note that by default these are mass weights, rather than atomic weights."""
        comp_str = ''

        for iso in comp.keys():
            iso_serp = isoname.mixed_2_zzaaam(iso)
            iso_serp = zzaaam_2_serpent(iso_serp)

            if self.iso_flag == '':
                comp_str += "{0:>11}".format( "{0}".format(iso_serp) )
            else:
                comp_str += "{0:>11}".format( "{0}.{1}".format(iso_serp, self.iso_flag) )

            if mass_weighted:
                comp_str += "  -{0:.5G}\n".format(comp[iso])
            else:
                comp_str += "   {0:.5G}\n".format(comp[iso])
        return comp_str
            

    def make_input_fuel(self, ms=None):
        if hasattr(defchar, 'fuel_form_mass_weighted'):
            mass_weighted = defchar.fuel_form_mass_weighted
        else:
            mass_weighted = True

        if ms == None:
            comp_str = self.make_input_material_weights(defchar.initial_fuel_stream.comp, mass_weighted)
        else:
            comp_str = self.make_input_material_weights(ms.multByMass(), True)

        return comp_str

    def make_input_cladding(self):
        # Try to load cladding stream
        if hasattr(defchar, 'clad_form'):
            clad_stream = MassStream(defchar.clad_form)
        else:
            if 0 < defchar.verbosity:
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
            clad_stream = MassStream(clad_form)

        # Try to find if the cladding form is mass or atomic weighted
        if hasattr(defchar, 'clad_form_mass_weighted'):
            mass_weighted = defchar.clad_form_mass_weighted
        else:
            mass_weighted = True
        
        return self.make_input_material_weights(clad_stream.comp, mass_weighted)

    def make_input_coolant(self):
        # Try to load coolant stream
        if hasattr(defchar, 'cool_form'):
            cool_stream = MassStream(defchar.cool_form)
        else:
            if 0 < defchar.verbosity:
                print(message("Coolant not found.  Proceeding with borated light water."))
                print()
            MW = (2 * 1.0) + (1 * 16.0) + (0.199 * 550 * 10.0**-6 * 10.0) + (0.801 * 550 * 10.0**-6 * 11.0)
            cool_form = {
                10010: (2 * 1.0) / MW,
                80160: (1 * 16.0) / MW,
                50100: (0.199 * 550 * 10.0**-6 * 10.0) / MW,
                50110: (0.801 * 550 * 10.0**-6 * 11.0) / MW,
                }
            cool_stream = MassStream(cool_form)

        # Try to find if the coolant form is mass or atomic weighted
        if hasattr(defchar, 'cool_form_mass_weighted'):
            mass_weighted = defchar.cool_form_mass_weighted
        else:
            mass_weighted = True
        
        return self.make_input_material_weights(cool_stream.comp, mass_weighted)

    def make_input_geometry(self, n=0):
        # Open a new hdf5 file 
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'r')
        pert_cols = rx_h5.root.perturbations.cols

        # Require
        geom = {
            'fuel_radius': pert_cols.fuel_cell_radius[n],
            'clad_radius': pert_cols.clad_cell_radius[n],
            'void_radius': pert_cols.void_cell_radius[n],
            'cell_pitch':  pert_cols.unit_cell_pitch[n],
            }

        # Close hdf5 file
        rx_h5.close()

        # Tries to get the lattice specification
        # If it isn't present, use a default 17x17 PWR lattice
        if hasattr(defchar, 'lattice') and hasattr(defchar, 'lattice_xy'):
            geom['lattice']    = defchar.lattice
            geom['lattice_xy'] = defchar.lattice_xy
        else:
            if 0 < defchar.verbosity:
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

    def make_input_energy_groups(self):
        """Makes the energy group structure.

        CHAR and most other neutronics codes sepecify this using 
        upper energy bounds.  That way the number of bounds equals the number
        of groups G. Serpent, however, uses only the internal boundaries, so there
        are G-1 energies given for G groups.  Additionally, the highest energy group 
        has the range Bound[G-1] <= Group 1 < inifinity.  The lowest energy group 
        thus covers 0.0 MeV <= Group G < Bound[1].
        """
        e = {}
        gs = ["{0:.6G}".format(float(gb)) for gb in defchar.group_structure]

        # Set number of (serpent) groups
        e['n_groups'] = len(defchar.group_structure) - 1
        e['group_lower_bound'] = gs[0]
        e['group_upper_bound'] = gs[-1]
        e['group_inner_structure'] = "  " + "\n  ".join(gs[1:-1])
        e['group_structure'] = "  " + "\n  ".join(gs)

        return e        

    def make_burnup(self, n=0):
        """Generates a dictionary of values that fill the burnup portion of the serpent template."""
        # Open a new hdf5 file 
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'r')
        pert_cols = rx_h5.root.perturbations.cols

        # make burnup dictionary
        bu = {'decay_lib': defchar.serpent_decay_lib,
              'fission_yield_lib': defchar.serpent_fission_yield_lib,
              'num_burn_regions':  int(pert_cols.burn_regions[n]), 
              'fuel_specific_power': pert_cols.fuel_specific_power[n],
              }

        # Close hdf5 file
        rx_h5.close()

        bu['depletion_times'] = ''
        for ct in defchar.coarse_time[1:]:
            bu['depletion_times'] += '{0:>8.4G}\n'.format(ct)

        bu['transmute_inventory'] = ''
        for iso in defchar.core_transmute['zzaaam']:
            bu['transmute_inventory'] += '{0:>8}\n'.format(iso)

        return bu

    def make_detector(self, iso):
        """Generates a dictionary of values that fill the detector/cross-section portion of 
        the serpent template.  Requires the isotope to be specified."""
        det = {}

        iso_zz   = isoname.mixed_2_zzaaam(iso)
        iso_serp = zzaaam_2_serpent(iso_zz)

        # Set the isotope to calculate XS for
        if self.iso_flag == '':
            det['xsiso'] = "{0}".format(iso_serp) 
        else:
            det['xsiso'] = "{0}.{1}".format(iso_serp, self.iso_flag)

        # Setup detectors to calculate XS for
        det['xsdet'] = ''
        det_format = "det {tally_name} de energies dm fuel dr {tally_type} xsmat dt 3 phi\n"

        if hasattr(defchar, 'tallies'):
            tallies = defchar.tallies
        else:
            tallies = tally_types.serpent_default

        for tally in tallies:
            det['xsdet'] += det_format.format(tally_name=tally, tally_type=tallies[tally])

        return det

    def make_common_input(self, n=0):
        # Open a new hdf5 file 
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'r')
        pert_cols = rx_h5.root.perturbations.cols

        # Initial serpent fill dictionary
        serpent_fill = {
            'reactor': defchar.reactor,
            'xsdata':  defchar.serpent_xsdata,

            'fuel_density': '{0:.5G}'.format(pert_cols.fuel_density[n]),
            'clad_density': '{0:.5G}'.format(pert_cols.clad_density[n]),
            'cool_density': '{0:.5G}'.format(pert_cols.cool_density[n]),

            'k_particles':  defchar.k_particles,
            'k_cycles':     defchar.k_cycles,
            'k_cycles_skip': defchar.k_cycles_skip,
            }

        # Close hdf5 file.
        rx_h5.close()

        # Set the material lines
        serpent_fill['fuel']     = self.make_input_fuel()
        serpent_fill['cladding'] = self.make_input_cladding()
        serpent_fill['coolant']  = self.make_input_coolant()

        # Add the geometry information
        serpent_fill.update(self.make_input_geometry())

        # Set the energy group structure
        serpent_fill.update(self.make_input_energy_groups())

        # Assign serpent_fill to the class
        self.serpent_fill = serpent_fill


    def make_burnup_input(self, n=0):
        self.serpent_fill.update(self.make_burnup(n))

        # Fill the burnup template
        with open(defchar.reactor + "_burnup", 'w') as f:
            f.write(defchar.burnup_template.format(**self.serpent_fill))

        # Initilaize a new HDF5 file with this defchar info. 
        self.init_h5()


    def make_xs_gen_input(self, iso="U235"):
        self.serpent_fill.update( self.make_detector(iso) )

        # Fill the XS template
        with open(defchar.reactor + "_xs_gen", 'w') as f:
            f.write(defchar.xs_gen_template.format(**self.serpent_fill))


    def make_input(self):
        self.make_common_input()

        # Add burnup information
        if defchar.options.RUN_BURNUP:
            self.make_burnup_input()

        self.make_xs_gen_input()


    def get_mpi_flag(self):
        mpi_flag = ''

        if hasattr(defchar, 'run_parallel'):
            run_flag = defchar.run_parallel.upper()
        else:
            run_flag = ''

        #if defchar.scheduler == 'PBS':
        #    run_flag = 'PBS'

        if run_flag in ["MPI", "PBS"]:
            if hasattr(defchar, 'number_cpus'):
                num_cpus = defchar.number_cpus
            else:
                print(message("The number of cpus was not specified even though a multicore calculation was requested.\n"
                              "Setting the number of cpus to 1 for this calculation."))
                num_cpus = 1

            mpi_flag = '-mpi {0}'.format(defchar.number_cpus)

        return mpi_flag


    def run_script_walltime(self):
        # Set PBS_Walltime
        if defchar.scheduler in ["PBS"]:
            return 36
        else:
            return 36


    def run_script_fill_values(self):
        """Sets the fill values for running serpent."""

        rsfv = {}

        # Set Transport Job Context
        rsfv['transport_job_context'] = self.run_str + " -version"

        # Set Run_Commands 
        if defchar.options.Local:
            rsfv['run_commands'] = ''
        else:
            rsfv['run_commands'] = ''
            #rsfv['run_commands'] = 'lamboot\n'

        # Add burnup information
        if defchar.options.RunBU:
            #rsfv['run_commands'] += "{0} {1}_burnup {2}\n".format(self.run_str, defchar.reactor, self.get_mpi_flag())
            rsfv['run_commands'] += "char --cwd -bp defchar.py\n"

        #Add cross section information
        rsfv['run_commands'] += "char --cwd -x defchar.py\n"

        return rsfv


    def run_burnup(self):
        """Runs the burnup portion of CHAR."""
        # Initializa the common serpent_fill values
        self.make_common_input()
        run_command = "{0} {1}_burnup {2}".format(self.run_str, defchar.reactor, self.get_mpi_flag())

        # Open the hdf5 library
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'r')

        # Get the number of time points from the file
        nperturbations = len(rx_h5.root.perturbations)

        # close the file before returning
        rx_h5.close()

        # Initialize the hdf5 file to take XS data
        self.init_h5_burnup(nperturbations)

        # Loop over all non-burnup perturbations.
        data = [getattr(defchar, a) for a in defchar.perturbation_params[:-1]]
        pert = [r for r in product(*data)]
        for d


    def run_xs_gen(self):
        """Runs the cross-section generation portion of CHAR."""
        # Initializa the common serpent_fill values
        self.make_common_input()
        run_command = "{0} {1}_xs_gen {2}".format(self.run_str, defchar.reactor, self.get_mpi_flag())

        # Open the hdf5 library
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'r')

        # Get the number of time points from the file
        nperturbations = len(rx_h5.root.perturbations)

        # close the file before returning
        rx_h5.close()

        # Initialize the hdf5 file to take XS data
        self.init_h5_xs_gen(nperturbations)

        # Loop over all times
        for n in range(nperturbations):
            # Grab the MassStream at this time.
            ms_n = MassStream()
            ms_n.load_from_hdf5(defchar.reactor + ".h5", "/Ti0", n)

            # Loop over all output isotopes
            for iso in defchar.core_transmute['zzaaam']:
                defchar.logger.info('Generating cross-sections for {0} at perturbation step {1}.'.format(iso, n))

                # Add filler fision product
                # If iso is not zirconium, add Zr-90
                # If is zirconium, add Sr-90
                # These two isotopes have almost the same mass
                # and neutronic profile:
                #     http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc=Zr-90&n=2
                #     http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc=Sr-90&n=2
                # We need to do this to preseve the atom density of the fuel, 
                # while not inducing errors through self-shielding and strong absorbers.
                # Basically, Zr-90 and Sr-90 become representative fision product pairs.
                #
                # WARNING: This is only suppossed to be a first order correction!
                # Make sure that you include enough FP in core_transmute.
                top_up_mass = 1.0 - ms_n.mass
                if top_up_mass == 0.0:
                    top_up = 0.0
                elif isoname.zzLL[iso//10000] == 'ZR':
                    top_up = MassStream({380900: 90.0, 621480: 148.0}, top_up_mass)
                elif isoname.zzLL[iso//10000] == 'SM':
                    top_up = MassStream({400900: 90.0, 601480: 148.0}, top_up_mass)
                else:
                    top_up = MassStream({400900: 90.0, 621480: 148.0}, top_up_mass)

                ms = ms_n + top_up
                isovec, AW, MW = msn.convolve_initial_fuel_form(ms, defchar.fuel_chemical_form)
                ms = MassStream(isovec)

                # Update fuel in serpent_fill
                self.serpent_fill['fuel'] = self.make_input_fuel(ms)

                # Make new input file
                self.make_xs_gen_input(iso)

                # Run serpent on this input file as a subprocess
                rtn = subprocess.check_call(run_command, shell=True)

                # Parse & write this output to HDF5
                self.parse_xs_gen()
                self.write_xs_gen(iso, n)



    #
    # Parsing functions
    #

    def parse(self):
        """Convienence function to parse results."""
        self.parse_burnup()
        self.write_burnup()
        defchar.logger.info('Parsed burnup calculation.')


    def parse_burnup(self):
        """Parse the burnup/depletion files into an equivelent python modules.
        Writes the output to hdf5."""

        # Convert files
        convert_res(defchar.reactor + "_burnup_res.m")
        convert_dep(defchar.reactor + "_burnup_dep.m")


    def parse_xs_gen(self):
        """Parse the burnup/depletion files into an equivelent python modules.
        Writes the output to hdf5."""

        # Convert files
        convert_res(defchar.reactor + "_xs_gen_res.m")
        convert_det(defchar.reactor + "_xs_gen_det0.m")


    #
    # Init HDF5 groups and arrays
    #

    def init_h5(self):
        """Initialize a new HDF5 file in preparation for burnup and XS runs."""
        # Setup tables 
        desc = {}
        for n, param in enumerate(defchar.perturbation_params):
            desc[param] = tb.Float64Col(pos=n)

        # Open HDF5 file
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'w')
        base_group = "/"

        # Add pertubation table index
        p = rx_h5.createTable(base_group, 'perturbations', desc)

        # Add rows to table
        data = [getattr(defchar, a) for a in defchar.perturbation_params]
        rows = [r for r in product(*data)]
        p.append(rows)

        # Add the isotope tracking arrays.  
        rx_h5.createArray(base_group, 'isostrack', np.array(defchar.core_transmute['zzaaam']), 
                          "Isotopes to track, copy of transmute_isos_zz")

        rx_h5.createArray(base_group, 'load_isos_zz', np.array(defchar.core_load['zzaaam']), 
                          "Core loading isotopes [zzaaam]")
        rx_h5.createArray(base_group, 'load_isos_LL', np.array(defchar.core_load['LLAAAM']), 
                          "Core loading isotopes [LLAAAM]")

        rx_h5.createArray(base_group, 'transmute_isos_zz', np.array(defchar.core_transmute['zzaaam']), 
                          "Core transmute isotopes [zzaaam]")
        rx_h5.createArray(base_group, 'transmute_isos_LL', np.array(defchar.core_transmute['LLAAAM']), 
                          "Core transmute isotopes [LLAAAM]")

        # Close HDF5 file
        rx_h5.close()


    def init_tally_group(self, rx_h5, base_group, tally, init_array, 
                         gstring='Group {tally}', astring='Array {tally} {iso}'):
        """Inits a tally group in an hdf5 file."""

        # Remove existing group, create new group of same name.
        if hasattr(rx_h5.getNode(base_group), tally):
            rx_h5.removeNode(base_group, tally, recursive=True)
        tally_group = rx_h5.createGroup(base_group, tally, gstring.format(tally=tally))

        # Add isotopic arrays for this tally to this group.
        for iso_LL in defchar.core_transmute['LLAAAM']: 
            rx_h5.createArray(tally_group, iso_LL, init_array, astring.format(tally=tally, iso=iso_LL))


    def init_h5_burnup(self, nperturbations=1):
        """Initialize the hdf5 file for a set of burnup calculations based on the its input params."""
        # Open a new hdf5 file 
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'a')
        base_group = "/"

        # Number of energy groups
        G = self.serpent_fill['n_groups']

        # Init the raw tally arrays
        neg1 = -1.0 * np.ones( (nperturbations, ) )
        negG = -1.0 * np.ones( (nperturbations, G) )
        negE = -1.0 * np.ones( (nperturbations, G+1) )


        # Add basic BU information
        rx_h5.createArray(base_group, 'BU0',   neg1, "Burnup of the initial core loading [MWd/kg]")
        rx_h5.createArray(base_group, 'time0', neg1, "Time after initial core loading [days]")

        # Add flux arrays
        rx_h5.createArray(base_group, 'phi',   neg1, "Total flux [n/cm2/s]")
        rx_h5.createArray(base_group, 'phi_g', negG, "Group fluxes [n/cm2/s]")

        # Create Fluence array
        rx_h5.createArray(base_group, 'Phi', neg1, "Fluence [n/kb]")

        # Energy Group bounds
        rx_h5.createArray(base_group, 'energy', negE, "Energy boundaries [MeV]")

        # Initialize transmutation matrix
        self.init_tally_group(rx_h5, base_group, 'Ti0', neg1, 
                              "Transmutation matrix from initial core loading [kg_i/kgIHM]", 
                              "Mass weight of {iso} [kg/kgIHM]")

        rx_h5.createArray(base_group + '/Ti0', 'Mass', neg1, "Mass fraction of fuel [kg/kgIHM]")

        # close the file before returning
        rx_h5.close()


    def init_h5_xs_gen(self, nperturbations=1):
        """Initialize the hdf5 file for writing for the XS Gen stage.
        The shape of these arrays is dependent on the number of time steps."""
        # Open a new hdf5 file 
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'a')
        base_group = "/"

        # Number of energy groups
        G = self.serpent_fill['n_groups']

        # Grab the tallies
        if hasattr(defchar, 'tallies'):
            tallies = defchar.tallies
        else:
            tallies = tally_types.serpent_default
        # Init the raw tally arrays
        neg1 = -1.0 * np.ones( (nperturbations, G) )
        negG = -1.0 * np.ones( (nperturbations, G, G) )

        for tally in tallies:
            self.init_tally_group(rx_h5, base_group, tally, neg1, 
                                  "Microscopic Cross Section {tally} [barns]", 
                                  "Microscopic Cross Section {tally} for {iso} [barns]")

        # Init aggregate tallies

        # nubar
        if ('sigma_f' in tallies) and ('nubar_sigma_f' in tallies) and ('nubar' not in tallies):
            self.init_tally_group(rx_h5, base_group, 'nubar', neg1, 
                                  "{tally} [unitless]", "{tally} for {iso} [unitless]")

        # sigma_i
        if np.array(['sigma_i' in tally  for tally in tallies]).any() and ('sigma_i' not in tallies):
            self.init_tally_group(rx_h5, base_group, 'sigma_i', neg1, 
                                  "Microscopic Cross Section {tally} [barns]", 
                                  "Microscopic Cross Section {tally} for {iso} [barns]")

        # sigma_s
        if ('sigma_s' not in tallies):
            self.init_tally_group(rx_h5, base_group, 'sigma_s', neg1, 
                                  "Microscopic Cross Section {tally} [barns]", 
                                  "Microscopic Cross Section {tally} for {iso} [barns]")

        # Scattering kernel, sigma_s_gh
        self.init_tally_group(rx_h5, base_group, 'sigma_s_gh', negG, 
                              "Microscopic Scattering Kernel {tally} [barns]", 
                              "Microscopic Scattering Kernel {tally} for {iso} [barns]")

        # close the file before returning
        rx_h5.close()


    #
    # Writing functions
    #

    def write_burnup(self, n):
        """Writes the results of the burnup calculation to an hdf5 file.

        n : Perturbation index of first time step for this burnup calculation.
        """

        # Add current working directory to path
        if sys.path[0] != os.getcwd():
            sys.path.insert(0, os.getcwd())

        # Import data
        rx_res = __import__(defchar.reactor + "_burnup_res")
        rx_dep = __import__(defchar.reactor + "_burnup_dep")

        # Find end index 
        t = n + len(rx_dep.DAYS)

        # Open a new hdf5 file 
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'a')
        base_group = rx_h5.root

        # Add basic BU information
        base_group.BU0[n:t] =  rx_dep.BU
        base_group.time0[n:t] = rx_dep.DAYS

        phi = rx_res.TOT_FLUX[:, ::2].flatten()
        base_group.phi[n:t] = phi
        base_group.phi_g[n:t] = rx_res.FLUX[:,::2][:, 1:]

        # Create Fluence array
        Phi = np.zeros(len(phi))
        Phi[1:] = cumtrapz(phi * (10.0**-21), rx_dep.DAYS * (3600.0 * 24.0))
        base_group.Phi[n:t] = Phi

        # Energy Group bounds 
        base_group.energy[n:t] = rx_res.GC_BOUNDS

        # Calculate and store weight percents per IHM
        # Serepent masses somehow unnoprmalize themselves in all of these conversions, which is annoying.
        # This eefect is of order 1E-5, which is large enough to be noticable.
        # Thus we have to go through two bouts of normalization here.
        mw_conversion = defchar.fuel_weight / (defchar.IHM_weight * rx_dep.TOT_VOLUME * defchar.fuel_density)
        mw = rx_dep.TOT_MASS * mw_conversion 

        iso_LL = {}
        iso_index = {}
        for iso_zz in rx_dep.ZAI:
            # Find valid isotope indeces
            try: 
                iso_LL[iso_zz] = isoname.mixed_2_LLAAAM(int(iso_zz))
            except:
                continue
            iso_index[iso_zz] = getattr(rx_dep, 'i{0}'.format(iso_zz)) - 1

        mass = mw[iso_index.values()].sum(axis=0)   # Caclulate actual mass of isotopes present

        # Store normalized mass vector for each isotope
        for iso_zz in iso_index:
            iso_array = getattr(base_group.Ti0, iso_LL[iso_zz])

            mw_i =  mw[iso_index[iso_zz]] / mass[0]
            iso_array[n:t] = mw_i

        # Renormalize mass
        mass = mass / mass[0]   
        base_group.Ti0.Mass[n:t] = mass

        # close the file before returning
        rx_h5.close()


    def write_xs_gen(self, iso, n):
        # Convert isoname
        iso_zz = isoname.mixed_2_zzaaam(iso)
        iso_LL = isoname.zzaaam_2_LLAAAM(iso_zz)

        # Add current working directory to path
        if sys.path[0] != os.getcwd():
            sys.path.insert(0, os.getcwd())

        # Import data
        rx_res = __import__(defchar.reactor + "_xs_gen_res")
        rx_det = __import__(defchar.reactor + "_xs_gen_det0")

        # Open a new hdf5 file 
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'a')
        base_group = rx_h5.root

        # Grab the tallies
        if hasattr(defchar, 'tallies'):
            tallies = defchar.tallies
        else:
            tallies = tally_types.serpent_default

        # Write the raw tally arrays for this time and this iso        
        for tally in tallies:
            tally_hdf5_group = getattr(base_group, tally)
            tally_hdf5_array = getattr(tally_hdf5_group, iso_LL)

            tally_serp_array = getattr(rx_det, 'DET{0}'.format(tally))

            tally_hdf5_array[n] = tally_serp_array[::-1, 10]

        # Write aggregate tallies

        # nubar
        if ('sigma_f' in tallies) and ('nubar_sigma_f' in tallies) and ('nubar' not in tallies):
            tally_hdf5_group = getattr(base_group, 'nubar')
            tally_hdf5_array = getattr(tally_hdf5_group, iso_LL)

            sigma_f = getattr(rx_det, 'DETsigma_f')
            sigma_f = sigma_f[::-1, 10]

            nubar_sigma_f = getattr(rx_det, 'DETnubar_sigma_f')
            nubar_sigma_f = nubar_sigma_f[::-1, 10] 

            nubar = nubar_sigma_f / sigma_f

            tally_hdf5_array[n] = nubar

        # sigma_i
        sigma_i = None
        if ('sigma_i' not in tallies):
            tally_hdf5_group = getattr(base_group, 'sigma_i')
            tally_hdf5_array = getattr(tally_hdf5_group, iso_LL)

            # Sum all sigma_iN's together
            for tally in tallies:
                tally_serp_array = getattr(rx_det, 'DET{0}'.format(tally))

                if ('sigma_i' in tally) and (sigma_i == None):
                    sigma_i = tally_serp_array[::-1, 10]
                elif ('sigma_i' in tally) and (sigma_i != None):
                    sigma_i += tally_serp_array[::-1, 10]
                else:
                    pass

            tally_hdf5_array[n] = sigma_i

        # sigma_s
        sigma_s = None
        if ('sigma_s' not in tallies):
            tally_hdf5_group = getattr(base_group, 'sigma_s')
            tally_hdf5_array = getattr(tally_hdf5_group, iso_LL)

            if 'sigma_e' in tallies:
                sigma_e = getattr(rx_det, 'DETsigma_e')
                sigma_e = sigma_e[::-1, 10]
            else:
                sigma_e = None

            if (sigma_i == None) and ('sigma_i' in tallies):
                sigma_i = getattr(rx_det, 'DETsigma_i')
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
            tally_hdf5_array = getattr(tally_hdf5_group, iso_LL)

            gtp = rx_res.GTRANSFP[rx_res.idx][::2]
            G = len(sigma_s)
            gtp = gtp.reshape((G, G))

            sigma_s_gh = sigma_s * gtp

            tally_hdf5_array[n] = sigma_s_gh

        # close the file before returning
        rx_h5.close()
