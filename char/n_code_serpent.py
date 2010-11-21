"""A class to setup, run, and parse Serpent."""
from __future__ import print_function
import os
import sys

import isoname
import numpy as np
import tables as tb
import metasci.nuke as msn
from metasci import SafeRemove
from MassStream import MassStream
from metasci.colortext import message, failure

from scipy.integrate import cumtrapz

from char import defchar
from n_code import NCode
from m2py import convert_res, convert_dep

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
        self.place_remote_files = [defchar.reactor]
        self.fetch_remote_files = [defchar.reactor, 
                                   defchar.reactor + '.seed', 
                                   defchar.reactor + '_geom1.png', 
                                   defchar.reactor + '_mesh1.png', 
                                   defchar.reactor + '_lat_res.m', 
                                   ]

    def make_input_material_weights(self, comp, mass_weighted=True):
        """This function takes an isotopic vector, comp, and returns a serpent string representation.
        Note that by default these are mass weights, rather than atomic weights."""
        comp_str = ''

        for iso in comp.keys():
            iso_serp = isoname.mixed_2_zzaaam(iso)

            if 0 == iso_serp%10:
                iso_serp = iso_serp/10
            else:
                iso_zz = iso_serp/10000
                iso_LL = isoname.zzLL[iso_zz]
                iso_LL.capitalize()
                iso_aaa = (iso_serp/10)%1000
                iso_serp = "{0}-{1}m".format(iso_LL, iso_aaa)

            if self.iso_flag == '':
                comp_str += "{0:>11}".format( "{0}".format(iso_serp) )
            else:
                comp_str += "{0:>11}".format( "{0}.{1}".format(iso_serp, self.iso_flag) )

            if mass_weighted:
                comp_str += "  -{0:.5G}\n".format(comp[iso])
            else:
                comp_str += "   {0:.5G}\n".format(comp[iso])
        return comp_str
            

    def make_input_fuel(self):
        if hasattr(defchar, 'fuel_form_mass_weighted'):
            mass_weighted = defchar.fuel_form_mass_weighted
        else:
            mass_weighted = True
        
        return self.make_input_material_weights(defchar.initial_fuel_stream.comp, mass_weighted)

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

    def make_input_geometry(self):
        # Require
        geom = {
            'fuel_radius': defchar.fuel_cell_radius,
            'clad_radius': defchar.clad_cell_radius,
            'cell_pitch':  defchar.unit_cell_pitch,
            }

        # Tries to add a void region, 
        # If there isn't space, cladding is used instead.
        if hasattr(defchar, 'void_cell_radius'):
            geom['void_radius'] = defchar.void_cell_radius
        else:
            geom['void_radius'] = defchar.fuel_cell_radius + 0.0085
            if defchar.clad_cell_radius <= geom['void_radius']:
                geom['void_radius'] = defchar.fuel_cell_radius

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
        half_lat_pitch = (float(geom['lattice_xy']) * defchar.unit_cell_pitch) / 2.0
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

        # Set number of (serpent) groups
        e['n_groups'] = len(defchar.group_structure) + 1

        # Set serpent energy group bounds.
        gs = list(defchar.group_structure)
        gs = str(gs)
        gs = gs[1:-1]
        gs = gs.replace(',', '')
        e['group_structure'] = gs

        return e        

    def make_burnup(self):
        """Generates a dictionary of values that fill the burnup portion of the serpent template."""
        bu = {'decay_lib': defchar.serpent_decay_lib,
              'fission_yield_lib': defchar.serpent_fission_yield_lib,
              'num_burn_regions':  int(defchar.burn_regions), 
              'fuel_specific_power': defchar.fuel_specific_power,
              }

        bu['depletion_times'] = ''
        for ct in defchar.coarse_time:
            bu['depletion_times'] += '{0:>8.4G}\n'.format(ct)

        bu['transmute_inventory'] = ''
        for iso in defchar.core_transmute['zzaaam']:
            bu['transmute_inventory'] += '{0:>8}\n'.format(iso)

        return bu

    def make_input(self):
        serpent_fill = {
            'reactor': defchar.reactor,
            'xsdata':  defchar.serpent_xsdata,

            'fuel_density': '{0:.5G}'.format(defchar.fuel_density),
            'clad_density': '{0:.5G}'.format(defchar.clad_density),
            'cool_density': '{0:.5G}'.format(defchar.cool_density),

            'k_particles':  defchar.k_particles,
            'k_cycles':     defchar.k_cycles,
            'k_cycles_skip': defchar.k_cycles_skip,
            }

        # Set the material lines
        serpent_fill['fuel']     = self.make_input_fuel()
        serpent_fill['cladding'] = self.make_input_cladding()
        serpent_fill['coolant']  = self.make_input_coolant()

        # Add the geometry information
        serpent_fill.update(self.make_input_geometry())

        # Set the energy group structure
        serpent_fill.update(self.make_input_energy_groups())

        # Add burnup information
        if defchar.options.NoBurnBool:
            pass
        else:
            serpent_fill.update(self.make_burnup())

            # Fill the burnup template
            with open(defchar.reactor + "_burnup", 'w') as f:
                f.write(defchar.burnup_template.format(**serpent_fill))

        # Fill the XS template
        with open(defchar.reactor + "_xs_gen", 'w') as f:
            f.write(defchar.xs_gen_template.format(**serpent_fill))

        return


    def run_script_fill_values(self, runflag=''):
        """Sets the fill values for running serpent."""

        rsfv = {}

        # Set PBS_Walltime
        if runflag in ["PBS"]:
            rsfv['PBS_walltime'] = ",walltime={0:02G}:00:00\n".format(36)
        else:
            rsfv['PBS_walltime'] = '\n'
        
        # Set Schedular input an output files
        rsfv['PBS_stagein_settings']  = ''
        rsfv['PBS_stageout_settings'] = ''
        if runflag in ["PBS"]:
            for f in self.place_remote_files:
                rsfv['PBS_stagein_settings'] += "#PBS -W stagein=./{f}@{rg}:{rd}{f}\n".format(f=f, rd=RemoteDir, rg=RemoteGateway)

            for f in self.fetch_remote_files:
                if f not in self.place_remote_files:
                    rsfv['PBS_stageout_settings']  = "#PBS -W stageout=./{f}@{rg}:{rd}{f}\n".format(f=f, rd=RemoteDir, rg=RemoteGateway)

        # Set Transport Job Context
        rsfv['transport_job_context'] = self.run_str + " -version"

        # Set Run_Commands 
        mpi_flag = ''
        if runflag in ["MPI", "PBS"]:
            if hasattr(defchar, 'number_cpus'):
                num_cpus = defchar.number_cpus
            else:
                print(message("The number of cpus was not specified even though a multicore calculation was requested."))
                print(message("Setting the number of cpus to 1 for this calculation."))
                num_cpus = 1

            mpi_flag = '-mpi {0}'.format(NumberCPUs)

        rsfv['run_commands'] = "{0} {1}_burnup {2}\n".format(self.run_str, defchar.reactor, mpi_flag)

        return rsfv


    def parse(self):
        """Convienence function to parse results."""
        self.parse_burnup()


    def parse_burnup(self):
        """Parse the burnup/depletion files into an equivelent python modules.
        Writes the output to hdf5."""

        # Convert files
        convert_res(defchar.reactor + "_burnup_res.m")
        convert_dep(defchar.reactor + "_burnup_dep.m")

        self.write_burnup()


    def write_burnup(self):
        """Writes the results of the burnup calculation to an hdf5 file."""

        # Add current working directory to path
        sys.path.insert(0, os.getcwd())

        # Import data
        rx_res = __import__(defchar.reactor + "_burnup_res")
        rx_dep = __import__(defchar.reactor + "_burnup_dep")

        # Open a new hdf5 file 
        rx_h5 = tb.openFile(defchar.reactor + ".h5", 'w')
        base_group = "/"

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

        # Add basic BU information
        rx_h5.createArray(base_group, 'BU0', rx_dep.BU, "Burnup of the initial core loading [MWd/kg]")
        rx_h5.createArray(base_group, 'time0', rx_dep.DAYS, "Time after initial core loading [days]")

        phi = rx_res.TOT_FLUX[:, ::2].flatten()
        rx_h5.createArray(base_group, 'phi', phi, "Total flux [n/cm2/s]")
        rx_h5.createArray(base_group, 'phi_g', rx_res.FLUX[:,::2][:, 1:], "Group fluxes [n/cm2/s]")

        # Create Fluence array
        Phi = np.zeros(len(phi))
        Phi[1:] = cumtrapz(phi * (10.0**-21), rx_dep.DAYS * (3600.0 * 24.0))
        rx_h5.createArray(base_group, 'Phi', Phi, "Fluence [n/kb]")

        # Energy Group bounds
        rx_h5.createArray(base_group, 'energy', rx_res.GC_BOUNDS[0], "Energy boundaries [MeV]")

        # Calculate and store weight percents per IHM
        mw_conversion = defchar.fuel_weight / (defchar.IHM_weight * rx_dep.TOT_VOLUME * defchar.fuel_density)
        mw = rx_dep.TOT_MASS * mw_conversion 

        Ti0_group = rx_h5.createGroup(base_group, 'Ti0', "Transmutation matrix from initial core loading [kg_i/kgIHM]")

        for iso_zz in rx_dep.ZAI:
            try: 
                iso_LL = isoname.mixed_2_LLAAAM(int(iso_zz))
            except:
                continue

            i = getattr(rx_dep, 'i{0}'.format(iso_zz)) - 1

            rx_h5.createArray(Ti0_group, iso_LL, mw[i], "Mass weight of {0} [kg/kgIHM]".format(iso_LL))

        # close the file before returning
        rx_h5.close()
