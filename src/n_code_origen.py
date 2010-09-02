"""A class to setup, run, and parse ORIGEN."""
import shutil

import isoname
import metasci.nuke as msn
import metasci.nuke.Origen as msno
from metasci import SafeRemove
from metasci.colortext import *

from char import reactor

from char import FineTimeIndex, FineTime
from char import CoarseTimeIndex, CoarseTime
from char import CoreLoad_zzaaam, CoreLoad_LLAAAM, CoreLoad_MCNP
from char import CoreTran_zzaaam, CoreTran_LLAAAM, CoreTran_MCNP

from n_code import NCode


class NCodeORIGEN(NCode):
    """An ORIGEN neutronics code wrapper class."""

    def __init__(self):
        self.name    = "ORIGEN"

        # Try to find & set the origen fast or therm flag.
        try:
            from defchar import ORIGEN_FASTorTHERM
        except NameError:
            ORIGEN_FASTorTHERM = ''

        self.fast_or_therm = ORIGEN_FASTorTHERM.lower()
        if not (ORIGEN_FASTorTHERM in ['fast', 'therm']):
            print(message("ORIGEN_FASTorTHERM not properly set: {0}", ORIGEN_FASTorTHERM))
            print(message("Setting to 'therm'."))
            self.fast_or_therm = 'therm'

        self.run_str = "o2"

        # Try to find and set the origenc conectration cut-off
        try:
            from defchar import ORIGEN_Concentration_Cut_Off
        except NameError:
            print(message("ORIGEN_Concentration_Cut_Off not set."))
            print(message("Using default value of 1.0E-10."))
            ORIGEN_Concentration_Cut_Off = 10.0**-10

        try:
            self.cut_off = "{0:.3E}".format(ORIGEN_Concentration_Cut_Off)
        except:
            print(message("ORIGEN_Concentration_Cut_Off not properly set: {0}", ORIGEN_Concentration_Cut_Off))
            print(message("Using default value of 1.0E-10."))
            self.cut_off = '1.-10'

        # Remote file lists
        self.place_remote_files = []
        self.fetch_remote_files = []

        return

    def make_input_tape4(self, isovec, n = "TAPE4.INP"):
        """Writes a TAPE4.INP file."""
        msno.writeTAPE4(isovec, n)
        return

    def make_input_tape5(self, t, p = ""):
        """Writes a TAPE5 files to the current directory + p.
        Only use this function from within the libs/ORIGEN/ directory."""

        # Grab important data from library file.
        libfile = tb.openFile('../{0}.h5'.format(reactor), 'r')
        lfr = libfile.root
        n_t = FineTime.index(t)
        IRFtime = FineTime[n_t] - FineTime[n_t-1]
        IRFflux = lfr.Fine.flux[n_t]
        libfile.close()

        # Grab the library numbers used in the original TAPE9 file...
        NLB = []
        tape9_orig = open("../../{0}.tape9".format(reactor), 'r')
        for line in tape9_orig:
            ls = line.split()
            if ls == []:
                continue
            elif (3 < int(ls[0])) and not (ls[0] in NLB):
                NLB.append(ls[0])
        tape9_orig.close()

        # Make template file fill-value dictionary
        tape5_kw = {
            'CUT_OFF': self.cut_off,
            'IRFtime': '{0:.10E}'.format(IRFtime),
            'IRFflux': '{0:.10E}'.format(IRFflux),
            'NLB1':    NLB[0],
            'NLB2':    NLB[1],
            'NLB3':    NLB[2],
            }

        # Fill the template
        with open("../../../templates/{0}.tape5.origen.template".format(reactor), 'r') as f:
            template_file = f.read()

        with open("{0}{1}_T{2}.tape5".format(p, reactor, t), 'w') as f:
            f.write(template_file.format(**tape5_kw))

        return

    def make_input_tape9(self, t, p = ""):
        """Makes a new ORIGEN TAPE9.INP file based on original TAPE9 file and stored XS.
        Places the new TAPE9.INP library in the current directory + path p.

        Args:
            * `t` (float): Time at which to evaluate the new library.

        Keyword Args:
    	    * `p` (str):  Path to append to the new tape9 file.

        Only use this function from within the libs/ORIGEN/ directory.
        """

        libfile = tb.openFile('../{0}.h5'.format(reactor), 'r')
        lfr = libfile.root
        lfrf = lfr.Fine

        n_t    = FineTime.index(t)
        phi_t  = lfrf.flux[n_t]
        phi_gt = lfrf.flux_g[n_t]

        tape9_kw = {
            'name_new': "{0}{1}_T{2}.tape9".format(p, reactor, t),
            'name_org': "../../../templates/{0}.tape9.origen.template".format(reactor),
            'iso_list': list(lfr.CoreLoad_zzaaam),
            'SNG':   {},
            'SN2N':  {},
            'SN3N':  {},
            'SNF':   {},
            'SNA':   {},
            'SNP':   {},
            'SNGX':  {},
            'SN2NX': {},
            }

        # Get Nomral XS
        for iso in tape9_kw['iso_list']:
            iso_LLAAAM = isoname.zzaaam_2_LLAAAM(iso)
            tape9_kw['SNG'][iso]  = msn.GroupCollapse(getattr(lfrf.sigma_gamma,  iso_LLAAAM)[n_t], phi_gt, phi_t)
            tape9_kw['SN2N'][iso] = msn.GroupCollapse(getattr(lfrf.sigma_2n,     iso_LLAAAM)[n_t], phi_gt, phi_t)
            tape9_kw['SN3N'][iso] = msn.GroupCollapse(getattr(lfrf.sigma_3n,     iso_LLAAAM)[n_t], phi_gt, phi_t)
            tape9_kw['SNF'][iso]  = msn.GroupCollapse(getattr(lfrf.sigma_f,      iso_LLAAAM)[n_t], phi_gt, phi_t)
            tape9_kw['SNA'][iso]  = msn.GroupCollapse(getattr(lfrf.sigma_alpha,  iso_LLAAAM)[n_t], phi_gt, phi_t)
            tape9_kw['SNP'][iso]  = msn.GroupCollapse(getattr(lfrf.sigma_proton, iso_LLAAAM)[n_t], phi_gt, phi_t)

        # (n, g*) XS
        for metaiso_table in lfrf.BranchRatio.NG:
            branch_ratio = metaiso_table[n_t]
            metaiso = isoname.mixed_2_zzaaam(metaiso_table.name)
            initiso = int(metaiso/10)*10 - 10

            if initiso in tape9_kw['iso_list']:
                tape9_kw['SNGX'][initiso] = tape9_kw['SNG'][initiso] * branch_ratio
        
        # (n, 2n*) XS
        for metaiso_table in lfrf.BranchRatio.N2N:
            branch_ratio = metaiso_table[n_t]
            metaiso = isoname.mixed_2_zzaaam(metaiso_table.name)
            initiso = int(metaiso/10)*10 + 10

            if initiso in tape9_kw['iso_list']:
                tape9_kw['SN2NX'][initiso] = tape9_kw['SN2N'][initiso] * branch_ratio

        # Close out the library file
        libfile.close()

        # Finally, make the tape9 file
        msno.writeTAPE9(**tape9_kw)

        return

    def run(self):
        """Runs the ORIGEN Burnup Calculations."""
        os.chdir('libs/ORIGEN/')

        # Grab General Data from the HDF5 File
        libfile = tb.openFile("../{0}.h5".format(reactor), 'r')
        CoreLoadIsos = list(libfile.root.CoreLoad_zzaaam)
        libfile.close()

        if 0 < verbosity:
            print(message("Preping the ORIGEN Directories..."))
        t1 = time.time()
        for t in FineTime[1:]:
            self.make_input_tape5(t)
            self.make_input_tape9(t)
        for iso in CoreLoadIsos: 
            os.mkdir("{0}".format(iso))
        t2 = time.time()
        if 0 < verbosity:
            print(message("...Done!  That only took {0:time} min.\n", "{0:.3G}".format((t2-t1)/60.0) ))

        if 0 < verbosity:
            print(message("  ~~~~~  Starting ORIGEN Runs  ~~~~~  "))
        orit1 = time.time()

        # Initialize the data structures
        self.BU  = {}
        self.k   = {}
        self.Pro = {}
        self.Des = {}
        self.Tij = {}

        for iso in CoreLoadIsos:
            isoLL = isoname.zzaaam_2_LLAAAM(iso)
            if 0 < verbosity:
                print(message("  ~~~~~  Now on {0:iso}  ~~~~~  \n", "Isotope {0}".format(isoLL)))
            isot1 = time.time()

            # Initilize iso data, for t = 0
            self.BU[iso]  = [0.0]
            self.k[iso]   = [0.0]
            self.Pro[iso] = [0.0]
            self.Des[iso] = [0.0]
            self.Tij[iso] = [{iso: 1000.0}]

            for t in FineTime[1:]:
                if 0 < verbosity:
                   print(message("Starting ORIGEN run for {0:iso} at {1:time}...", isoLL, "Time {0}".format(t)))
                t1 = time.time()

                os.chdir("{0}".format(iso))

                # Make/Get Input Decks
                self.make_input_tape4(Tij[iso][-1])
                shutil.copy("../{0}_T{1}.tape5".format(reactor, t), "TAPE5.INP")
                shutil.copy("../{0}_T{1}.tape9".format(reactor, t), "TAPE9.INP")

                # Run ORIGEN
                subprocess.call("o2_{0}_linux.exe".format(ORIGEN_FASTorTHERM), shell=True)

                # Parse Output
                parsed = parsechar.Parse_TAPE6()
                self.BU[iso].append(  BU[iso][-1] + parsed[0] )
                self.k[iso].append(   parsed[1] )
                self.Pro[iso].append( parsed[2] )
                self.Des[iso].append( parsed[3] )
                self.Tij[iso].append( parsed[4] )

                # Clean up the directory
                for f in os.listdir('.'):
                    if f[-4:] in ['.INP', '.OUT']:
                        metasci.SafeRemove(f)
                os.chdir('../') # Back to ORIGEN Directory

                t2 = time.time()
                if 0 < verbosity:
                    print(message("ORIGEN run completed in {0:time} min!", "{0:.3G} min".format((t2-t1)/60.0) ))
    
            isot2 = time.time()
            if 0 < verbosity:
                print(message("  ~~~~~  Isotope {0:iso} took {1:time} min!  ~~~~~  \n", isoLL, "{0:.3G} min".format((isot2-isot1)/60.0) ))


        # Kludge to put Tij in the right units and form
        allORIGENisoList = []
        for iso in CoreLoadIsos:
            for t in Tij[iso]:
                for j in t.keys():
                    if (j not in allORIGENisoList):
                        allORIGENisoList.append(j)
        for iso in CoreLoadIsos:
            for n_t in range(len(Tij[iso])):
                for j in allORIGENisoList:
                    if j in Tij[iso][n_t].keys():
                        Tij[iso][n_t][j] = Tij[iso][n_t][j] / (10.0**3)
                    else:
                        Tij[iso][n_t][j] = 0.0
    
        orit2 = time.time()
        if 0 < verbosity:
            print(message("  ~~~~~  ORIGEN took {0:time} to run!  ~~~~~  ", "{0:.3G} min".format((orit2-orit1)/60.0) ))

        os.chdir('../../') #Back to 'reactor' root
        return 
