"""A class to setup, run, and parse Serpent."""

import isoname
import numpy as np
import metasci.nuke as msn
from metasci import SafeRemove
from MassStream import MassStream

from char import reactor
from char import GroupStructure
from char import FMdic, dicFM

from char import FineTimeIndex, FineTime
from char import CoarseTimeIndex, CoarseTime
from char import CoreLoad_zzaaam, CoreLoad_LLAAAM, CoreLoad_MCNP
from char import CoreTran_zzaaam, CoreTran_LLAAAM, CoreTran_MCNP
from char import metastabletrak, metastableMCNP
from char import mat_number, number_mat
from char import InitialFuelStream
from char import kParticles, kCycles, kCyclesSkip

from n_code import NCode


class NCodeSerpent(NCode):
    """A Serpent neutronics code wrapper class."""

    def __init__(self):
        self.name    = "Serpent"
        self.run_str = "sss-dev"

        try:
            from defchar import ISO_FLAG
            self.iso_flag = ISO_FLAG
        except ImportError:
            self.iso_flag = ''

        # Remote file lists
        self.place_remote_files = [reactor]
        self.fetch_remote_files = [reactor]

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
                comp_str += " {0}".format(iso_serp)
            else:
                comp_str += " {0}.{1}".format(iso_serp, self.iso_flag)

            if mass_weighted:
                comp_str += "  -{0:.5G}\n".format(comp[iso])
            else:
                comp_str += "   {0:.5G}\n".format(comp[iso])

        return comp_str
            

    def make_input_fuel(self):
        try:
            from defchar import FuelForm_MassWeighted
            mass_weighted = FuelForm_MassWeighted
        except ImportError:
            mass_weighted = True
        
        return self.make_input_material_weights(InitialFuelStream.comp, mass_weighted)

    def make_input_cladding(self):
        # Try to load cladding stream
        try:
            from defchar import CladForm
            CladStream = MassStream(CladForm)
        except ImportError:
            print(message("Cladding not found.  Proceeding with standard zircaloy mixture."))
            CladForm = {
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
            CladStream = MassStream(CladForm)

        # Try to find if the cladding form is mass or atomic weighted
        try:
            from defchar import CladForm_MassWeighted
            mass_weighted = CladForm_MassWeighted
        except ImportError:
            mass_weighted = True
        
        return self.make_input_material_weights(CladStream.comp, mass_weighted)

    def make_input_coolant(self):
        # Try to load coolant stream
        try:
            from defchar import CoolForm
            CoolStream = MassStream(CoolForm)
        except ImportError:
            print(message("Coolant not found.  Proceeding with borated light water."))
            MW = (2 * 1.0) + (1 * 16.0) + (0.199 * 550 * 10.0**-6 * 10.0) + (0.801 * 550 * 10.0**-6 * 11.0)
            CoolForm = {
                10010: (2 * 1.0) / MW,
                80160: (1 * 16.0) / MW,
                50100: (0.199 * 550 * 10.0**-6 * 10.0) / MW,
                50110: (0.801 * 550 * 10.0**-6 * 11.0) / MW,
                }
            CoolStream = MassStream(CoolForm)

        # Try to find if the coolant form is mass or atomic weighted
        try:
            from defchar import CoolForm_MassWeighted
            mass_weighted = CoolForm_MassWeighted
        except ImportError:
            mass_weighted = True
        
        return self.make_input_material_weights(CoolStream.comp, mass_weighted)

    def make_input_geometry(self):
        # Require
        geom = {
            'FuelRadius': FuelCellRadius,
            'CladRadius': CladCellRadius,
            'CellPitch':  UnitCellPitch,
            }

        # Tries to add a void region, 
        # If there isn't space, cladding is used instead.
        try:
            from defchar import VoidCellRadius
             goem['VoidRadius'] = VoidCellRadius
        except ImportError:
            goem['VoidRadius'] = FuelCellRadius + 0.0085
            if CladCellRadius <= geom['VoidRadius']:
                geom['VoidRadius'] = FuelCellRadius

        # Tries to get the lattice specification
        # If it isn't present, use a default 17x17 PWR lattice
        try:
            from defchar import Lattice
            from defchar import LatticeXY
            goem['Lattice']   = Lattice
            goem['LatticeXY'] = LatticeXY
        except ImportError:
            if 0 < verbosity:
                print(message("Lattice specification not found."))
                print(message("Using a default 17x17 PWR lattice."))
            goem['LatticeXY'] = 17
            goem['Lattice']   = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" + \
                                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" + \
                                "1 1 1 1 1 2 1 1 2 1 1 2 1 1 1 1 1 \n" + \
                                "1 1 1 2 1 1 1 1 1 1 1 1 1 2 1 1 1 \n" + \
                                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" + \
                                "1 1 2 1 1 2 1 1 2 1 1 2 1 1 2 1 1 \n" + \
                                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" + \
                                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" + \
                                "1 1 2 1 1 2 1 1 2 1 1 2 1 1 2 1 1 \n" + \
                                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" + \
                                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" + \
                                "1 1 2 1 1 2 1 1 2 1 1 2 1 1 2 1 1 \n" + \
                                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" + \
                                "1 1 1 2 1 1 1 1 1 1 1 1 1 2 1 1 1 \n" + \
                                "1 1 1 1 1 2 1 1 2 1 1 2 1 1 1 1 1 \n" + \
                                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n" + \
                                "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n"

        # Determine if lattice is symetric
        lat = np.array(geom['Lattice'].split(), dtype='int')
        lat = lat.reshape((geom['LatticeXY'], geom['LatticeXY']))
        lat_trans = lat.transpose()
        if (lat == lat_trans).all():    # Condition for symmetry; A = A^T
            geom['SymFlag'] = ''
        else:
            geom['SymFlag'] = '% The lattice is not symmetric! Forced to use whole geometry.\n%'

        # Set half of the lattice pitch.
        half_lat_pitch = (float(geom['LatticeXY']) * CellPitch) / 2.0
        geom['HalfLatticePitch'] = "{0:.5G}".format(half_lat_pitch)

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
        e['NGroups'] = len(GroupStructure) + 1

        # Set serpent energy group bounds.
        gs = list(GroupStructure)
        gs = str(gs)
        gs = gs[1:-1]
        gs = gs.replace(',', '')
        e['GroupStructure'] = gs

        return e        

    def make_input(self):
        serpent_fill = {
            'reactor':     reactor,

            'FuelDensity': '{0:.5G}'.fromat(FuelDensity),
            'CladDensity': '{0:.5G}'.fromat(CladDensity),
            'CoolDensity': '{0:.5G}'.fromat(CoolDensity),

            'kParticles':  kParticles,
            'kCycles':     kCycles,
            'kCyclesSkip': kCyclesSkip,
            }

        # Set the material lines
        serpent_fill['fuel']     = self.make_input_fuel()
        serpent_fill['cladding'] = self.make_input_cladding()
        serpent_fill['coolant']  = self.make_input_coolant()

        # Add the geometry information
        serpent_fill.update(self.make_input_geometry())

        # Set the energy group structure
        serpent_fill.update(self.make_input_energy_groups())

        # Fill the template
        with open('templates/{0}.serpent.template'.format(reactor), 'r') as f:
            template_file = f.read()

        with open(reactor, 'r') as f:
            f.write(template_file.format(**serpent_fill))

        return
