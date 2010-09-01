"""A class to setup, run, and parse MCNP."""

from math import pi

import isoname
import metasci.nuke as msn
from metasci import SafeRemove

from char import reactor
from char import UnitCellHeight, FuelCellRadius, CladCellRadius, UnitCellPitch, \
    FuelDensity, CladDensity, CoolDensity, FuelSpecificPower
from char import ENDF_FLAG
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


class NCodeMCNP(NCode):
    """An MCNP neutronics code wrapper class."""

    def __init__(self):
        self.name    = "MCNP"
        self.run_str = "mcnp"

        ### Calculate Cell Volumes ###
        # Division by 4 because template only looks at a quarter of the unit cell
        self.FuelCellVolume = 0.25 * UnitCellHeight * pi * (FuelCellRadius**2) 
        self.CladCellVolume = (0.25 * UnitCellHeight * pi * (CladCellRadius**2)) - self.FuelCellVolume
        self.UnitCellVolume = 0.25 * UnitCellHeight * (UnitCellPitch**2)
        self.CoolCellVolume = self.UnitCellVolume - self.FuelCellVolume - self.CladCellVolume

        self.UnitCellHalfPitch = UnitCellPitch / 2.0 

        # Remote file lists
        self.place_remote_files = [reactor + '.i']
        self.fetch_remote_files = [reactor + '.i', reactor + '.o', 'runlog.txt', "run_{0}.*".format(reactor)]

    def make_input_burn(self):
        """Generates the Burn card string for an MCNP inputfile."""

        # make list of time steps
        TimeSteps = []
        for n in CoarseTimeIndex:
            if (0 == n%2) and not (n == 0):
                	TimeSteps.append(CoarseTime[n] - CoarseTime[n-2])

        # Loads the burn card.
        burnline  = "BURN TIME=" + str(TimeSteps)[1:-1].replace(',', '') 
        burnline += " MAT=1 MATVOL=" + "{0:G}".format(self.FuelCellVolume) 
        burnline += " POWER={0:G}".format(msn.CellPower(InitialFuelStream.comp, FuelSpecificPower, self.FuelCellVolume, FuelDensity)) 
        burnline += " PFRAC="
        for t in TimeSteps:
            burnline += "1.0 " 
        burnline += "OMIT=1 27 6014 7016 8018 9018 87223 88227 89228 90234 91229 91230 91232 92229 92230 92231 94245 "
        burnline += "95240 95642 95644 97245 97246 97247 97248 99249 99250 99251 99252 99253 "
        burnline += "BOPT=1.0 -04 -1"  

        burnline = msn.Line2MCNP(burnline)

        return burnline

    def make_input_mat1(self):
        """Generates the Material 1 card string for an MCNP inputfile."""

        # writes m1, the initial fuel composition
        m1line = "m1 "
        for key in InitialFuelStream.comp.keys():
            m1line = m1line + '{0}{1} -{2:G} '.format(isoname.zzaaam_2_MCNP( key ), ENDF_FLAG, InitialFuelStream.comp[key])
        for iso in CoreLoad_MCNP:
            if not (iso in isoname.zzaaam_2_MCNP_List(InitialFuelStream.comp.keys()) ) and not (iso in metastableMCNP):
                m1line = m1line + '{0}{1} -{2:G} '.format(iso, ENDF_FLAG, 1.0E-36)
        m1line = msn.Line2MCNP(m1line)

        return m1line

    def make_input_mat_other(self):
        """Generates the Other Material card string for an MCNP inputfile."""

        # writes remaining materials
        moline = ""
        for mnum in range(moffset+1, moffset+1+len(mat_number)):
            for iso in mat_number.keys():
                if mat_number[iso] == mnum:
                    moline = moline + "m{0} {1}{2} 1.0\n".format(mat_number[iso], iso, ENDF_FLAG) 
                    break
        return moline[:-1]

    def make_input_pert(self):
        """Generates the perturbation cards string for an MCNP inputfile."""

        # writes remaining materials
        pline  = ""
        pline += "PERT1:N CELL=1 RHO=-{0:G}\n".format(FuelDensity * 0.0)
        pline += "PERT2:N CELL=1 RHO=-{0:G}\n".format(FuelDensity * 0.5)
        pline += "PERT3:N CELL=1 RHO=-{0:G}\n".format(FuelDensity * 1.5)
        pline += "PERT4:N CELL=1 RHO=-{0:G}\n".format(FuelDensity * 2.0)
        return pline[:-1]

    def make_input_tally(self):
        """Generates the Tally card string for an MCNP inputfile."""

        tline  = ""

        # Write Group Flux Tally
        tline += "f14:n 1\n"
        tline += "FC14 flux\n"
        tline += "E14 {0} NT\n".format(GroupStructure)
        tline += "sd14 {0:G}\n".format(self.FuelCellVolume)
        tline += "c \n" 

        # Write Multiplier Tallies
        tline += "f24:n 1\n"
        tline += "FC24 XS\n"
        l2w = "fm24 "
        for mnum in sorted(mat_number.values()):
            for fm in FMdic.values():
                l2w = l2w + "(1 {0:d} {1:d}) ".format(mnum, fm)
        tline += msn.Line2MCNP(l2w) 
        tline += "\n" 
        tline += "E24 {0} NT\n".format(GroupStructure) 
        tline += "sd24 {0:G}\n".format(self.FuelCellVolume) 
        tline += "c \n"

        #Write Group Flux Tally for CINDER data
        #This is to fix meta-stable read-in flags for ORIGEN
        tline += "f34:n 1\n" 
        tline += "FC34 cinderflux\n" 
        CinderGroupLine = "E34 1.00000E-11 5.00000E-09 1.00000E-08 1.50000E-08 2.00000E-08 2.50000E-08 " + \
                        "3.00000E-08 3.50000E-08 4.20000E-08 5.00000E-08 5.80000E-08 6.70000E-08 " + \
                        "8.00000E-08 1.00000E-07 1.52000E-07 2.51000E-07 4.14000E-07 6.83000E-07 " + \
                        "1.12500E-06 1.85500E-06 3.05900E-06 5.04300E-06 8.31500E-06 1.37100E-05 " + \
                        "2.26000E-05 3.72700E-05 6.14400E-05 1.01300E-04 1.67000E-04 2.75400E-04 " + \
                        "4.54000E-04 7.48500E-04 1.23400E-03 2.03500E-03 2.40400E-03 2.84000E-03 " + \
                        "3.35500E-03 5.53100E-03 9.11900E-03 1.50300E-02 1.98900E-02 2.55400E-02 " + \
                        "4.08700E-02 6.73800E-02 1.11100E-01 1.83200E-01 3.02000E-01 3.88700E-01 " + \
                        "4.97900E-01 6.39279E-01 8.20850E-01 1.10803E+00 1.35335E+00 1.73774E+00 " + \
                        "2.23130E+00 2.86505E+00 3.67879E+00 4.96585E+00 6.06500E+00 1.00000E+01 " + \
                        "1.49182E+01 1.69046E+01 2.00000E+01 2.50000E+01 NT"
        tline += msn.Line2MCNP(CinderGroupLine)
        tline += "\n" 
        tline += "sd34 {0:G}\n".format(self.FuelCellVolume)
        tline += "c " 

        return tline


    def make_input_scat_tally(self):
        """Generates Group-to-Group Scattering XS Tallies."""

        tline  = ""

        tline += "F44:n 1\n"
        tline += "FC44 Group-to-Group Scattering XS\n"
        l2w = "FM44 "
        for mnum in sorted(mat_number.values()):
            l2w = l2w + "(1 {0:d} 2:4) ".format(mnum)
        tline += msn.Line2MCNP(l2w) 
        tline += "\n" 
        tline += "E44 {0} NT\n".format(GroupStructure) 
        tline += "SD44 {0:G}\n".format(self.FuelCellVolume) 
        tline += "FT44 SCX 1 \n"
        tline += "c \n"

        tline += "c " 

        return tline

    def make_input(self, NoBurnBool=False, NoPertBool=False):
        """Make the MCNP input file."""

        SafeRemove(reactor + ".i")
        mcnp_fill = {
            'FuelDensity': '{0:G}'.format(FuelDensity), 
            'CladDensity': '{0:G}'.format(CladDensity), 
            'CoolDensity': '{0:G}'.format(CoolDensity), 

            'FuelCellVolume': '{0:G}'.format(self.FuelCellVolume), 

            'FuelCellRadius':    '{0:G}'.format(FuelCellRadius), 
            'CladCellRadius':    '{0:G}'.format(CladCellRadius), 
            'UnitCellHalfPitch': '{0:G}'.format(self.UnitCellHalfPitch), 
            'UnitCellHeight':    '{0:G}'.format(UnitCellHeight), 

            'GroupStructure': GroupStructure, 

            'kParticles':  '{0:G}'.format(kParticles),    
            'kCycles':     '{0:G}'.format(kCycles),    
            'kCyclesSkip': '{0:G}'.format(kCyclesSkip),    
            }

        # Make the MCNP input deck fill values
        if NoBurnBool:
            mcnp_fill['Burn'] = 'c No Burn Card'
        else:
            mcnp_fill['Burn'] = self.make_input_burn()

        mcnp_fill['Mat1']     = self.make_input_mat1()
        mcnp_fill['MatOther'] = self.make_input_mat_other()

        if NoPertBool:
            mcnp_fill['Pert'] = 'c No Pert Cards'
        else:
            mcnp_fill['Pert'] = self.make_input_pert()

        mcnp_fill['Tally'] = self.make_input_tally()

        mcnp_fill['TallyScat'] = self.make_input_scat_tally()

        # Fill the template
        with open(reactor + '.i.template', 'r') as f:
            template_file = f.read()

        with open(reactor + '.i', 'w') as f:
            f.write(template_file.format(**mcnp_fill))

    return


    def run_script_fill_values(self):
        """Fills the runscript with values appropriate to MCNP."""

        rsfv = {}

        # Set PBS_Walltime
        if runflag in ["PBS"]:
            #Walltime = 4 hr/burn-step * Num burn-steps (1 + Numb Particle / (3000 part/hr/cpu) / NumCpu )
            #rsfv['PBS_Walltime'] = ",walltime={0:02G}:00:00\n".format(4*len(CoarseTime)*(1 + kParticles*kCycles/3000/NumberCPUs))
            rsfv['PBS_Walltime'] = ",walltime={0:02G}:00:00\n".format(36)
        else:
            rsfv['PBS_Walltime'] = '\n'

        # Set PBS_Stagein and PBS_Stageout
        rdict = {
            'RDir': RemoteDir,
            'RGateway': RemoteGateway, 
            'reactor': reactor, 
            }        

        if runflag in ["PBS"]:
            rsfv['PBS_Stagein'] = "#PBS -W stagein=./{reactor}.i@{RGateway}:{RDir}{reactor}.i\n".format(**rdict)

            rsfv['PBS_Stageout']  = "#PBS -W stageout=./{reactor}.o@{RGateway}:{RDir}{reactor}.o\n".format(**rdict)
            rsfv['PBS_Stageout'] += "#PBS -W stageout=./{reactor}.s@{RGateway}:{RDir}{reactor}.s\n".format(**rdict)
            rsfv['PBS_Stageout'] += "#PBS -W stageout=./{reactor}.m@{RGateway}:{RDir}{reactor}.m\n".format(**rdict)
            rsfv['PBS_Stageout'] += "#PBS -W stageout=./{reactor}.r@{RGateway}:{RDir}{reactor}.r\n".format(**rdict)
            rsfv['PBS_Stageout'] += "#PBS -W stageout=./CHAR_{reactor}.o*@{RGateway}:{RDir}CHAR_{reactor}.o*\n".format(**rdict)
        else:
            rsfv['PBS_Stagein']  = ''
            rsfv['PBS_Stageout'] = ''
    
        # Set Transport Job Context
        rsfv['Transport_Job_Context'] = "echo \"DATAPATH is ${DATAPATH}\""

        # Set Run_Commands 
        rsfv['Run_Commands']  = ''

        if runflag in ["PBS"]:
            rsfv['Run_Commands'] += "### Set MCNPX datapath variable\n"
            if localflag:
                PathDATAPATH = os.getenv("DATAPATH")
            else:
                PathDATAPATH = RemoteDATAPATH
            rsfv['Run_Commands'] += "export DATAPATH={0}\n".format(PathDATAPATH)
            rsfv['Run_Commands'] += "\n"

        if localflag:
            PathMPI  = LocalPathMPI
            PathMCNP = LocalPathMCNP
        else:
            PathMPI  = RemotePathMPI
            PathMCNP = RemotePathMCNP

        if runflag in ["MPI", "PBS"]:
            rsfv['Run_Commands'] += "### Run MCNP with MPI\n"
            rsfv['Run_Commands'] += "{0} \\\n".format(PathMPI)
            rsfv['Run_Commands'] += "-machinefile $PBS_NODEFILE \\\n"
            rsfv['Run_Commands'] += "{0} \\\n".format(PathMCNP)
            rsfv['Run_Commands'] += "i={0}.i \\\n".format(reactor)
            rsfv['Run_Commands'] += "o={0}.o \\\n".format(reactor)
            rsfv['Run_Commands'] += "s={0}.s \\\n".format(reactor)
            rsfv['Run_Commands'] += "m={0}.m \\\n".format(reactor)
            rsfv['Run_Commands'] += "r={0}.r   \n".format(reactor)
        else:
            rsfv['Run_Commands'] += "{0} inp={1}.i n={1}. ".format(PathMCNP, reactor)
            rsfv['Run_Commands'] += "\n"

        return rsfv


    def parse_init(self):
        """Initializes variables for MCNP parsing"""

        self.InFlux = False
        self.InXS = False
        self.InBurnSum = False
        self.AtMult = False
        self.AtBin = False
        self.LibFileOpen = False
        self.libfile = None
        self.mat = None
        self.fm = None
        self.n = 0
        self.flux_g = []
        self.E_up = []
        self.mult = {}
        self.burn = {
            "step":     [], 
            "duration": [], 
            "time":     [], 
            "power":    [], 
            "keff":     [], 
            "flux":     [], 
            "ave. nu":  [], 
            "ave. q":   [], 
            "burnup":   [], 
            "source":   [],
            }

        self.InFluxCinder = False
        self.E_up_cinder = []
        self.flux_g_cinder = []
        self.flux_cinder = []

        self.BR_NG  = {}
        self.BR_N2N = {}

        return

    def parse(self):
        """Parse MCNPX Output"""

        if not (self.mult == {}):
            self.parse_init()

        outfile = open(reactor + '.o', 'r')

        for line in outfile:
            ls = line.split()
            lslen = len(ls)

            if lslen == 0:
                if (self.InFlux or self.InFluxCinder) and self.AtBin:
                    self.InFlux = False
                    self.InFluxCinder = False
                    self.AtBin = False
                    continue
                elif self.AtMult and self.AtBin:
                    self.AtMult = False
                    self.AtBin = False
                    self.mat = None
                    self.fm = None
                    self.n = 0
                    continue
                elif self.InBurnSum:
                    self.InBurnSum = False
                else:
                    continue

            # Get Flux Data
            elif ls[0] == 'energy' and (self.InFlux or self.InFluxCinder):
                self.AtBin = True
                continue
            elif (1 < lslen) and (ls[1] == 'flux'):
                self.InFlux = True
                self.flux_g.append([])
                continue
            elif (1 < lslen) and (ls[1] == 'cinderflux'):
                self.InFluxCinder = True
                self.flux_g_cinder.append([])
                continue
            elif self.InFlux and self.AtBin:
                if len(self.flux_g) == 1:
                    self.E_up.append(float(ls[0]))
                self.flux_g[-1].append(float(ls[1]))
                continue
            elif self.InFluxCinder and self.AtBin:
                if len(self.flux_g_cinder) == 1:
                    self.E_up_cinder.append(float(ls[0]))
                self.flux_g_cinder[-1].append(float(ls[1]))
                continue

            # Get Multiplier Data
            elif (1 < lslen) and (ls[1] == 'XS'):
                self.InXS = True
                continue
            elif InXS and (ls[0] == 'multiplier'):
                self.AtMult = True
                self.mat = isoname.MCNP_2_LLAAAM( number_mat[ ls[-2] ] )
                self.fm = dicFM[ ls[-1] ]
                lib = '{0}_{1}.lib'.format(mat, fm)
                if lib in self.mult.keys():
                    self.mult[lib].append([])
                else:
                    self.mult[lib] = [[]]
                continue
            elif ls[0] == 'energy' and self.AtMult:
                self.AtBin = True
                continue
            elif self.AtMult and self.AtBin:
                self.mult[lib][-1].append( float(ls[1]) / self.flux_g[-1][n] )
                n = n + 1

            # Get Burnup Summary Table
            elif ls == ["(days)", "(days)", "(MW)", "(GWd/MTU)", "(nts/sec)"]:
                self.InBurnSum = True
                continue
            elif self.InBurnSum and (0 == len(burn['step'])):
                self.burn["step"].append( float(ls[0]) )
                self.burn["duration"].append( float(ls[1]) )
                self.burn["time"].append( float(ls[2]) )
                self.burn["power"].append( float(ls[3]) )
                self.burn["keff"].append( float(ls[4]) )
                self.burn["flux"].append( float(ls[5]) )
                self.burn["ave. nu"].append( float(ls[6]) )
                self.burn["ave. q"].append( float(ls[7]) )
                self.burn["burnup"].append( float(ls[8]) )
                self.burn["source"].append( float(ls[9]) )
            elif self.InBurnSum and (float(ls[0]) <= self.burn['step'][-1]):
                continue
            elif self.InBurnSum and (self.burn['step'][-1] < float(ls[0])):
                # Half Step Data Interpolation
                self.burn["step"].append( (self.burn["step"][-1] + float(ls[0])) / 2.0)
                self.burn["duration"].append( (self.burn["duration"][-1] + float(ls[1])) / 2.0)
                self.burn["time"].append( (self.burn["time"][-1] + float(ls[2])) / 2.0)
                self.burn["power"].append( (self.burn["power"][-1] + float(ls[3])) /2.0)
                self.burn["keff"].append( (self.burn["keff"][-1] + float(ls[4])) / 2.0)
                self.burn["flux"].append( (self.burn["flux"][-1] + float(ls[5])) / 2.0)
                self.burn["ave. nu"].append( (self.burn["ave. nu"][-1] + float(ls[6])) / 2.0)
                self.burn["ave. q"].append( (self.burn["ave. q"][-1] + float(ls[7])) / 2.0)
                self.burn["burnup"].append( (self.burn["burnup"][-1] + float(ls[8])) / 2.0)
                self.burn["source"].append( (self.burn["source"][-1] + float(ls[9])) / 2.0)

                # Full Step Data
                self.burn["step"].append( float(ls[0]) )
                self.burn["duration"].append( float(ls[1]) )
                self.burn["time"].append( float(ls[2]) )
                self.burn["power"].append( float(ls[3]) )
                self.burn["keff"].append( float(ls[4]) )
                self.burn["flux"].append( float(ls[5]) )
                self.burn["ave. nu"].append( float(ls[6]) )
                self.burn["ave. q"].append( float(ls[7]) )
                self.burn["burnup"].append( float(ls[8]) )
                self.burn["source"].append( float(ls[9]) )
                continue

        outfile.close()

        # Calculate Total Fluxes	
        # User-specified groups
        for gth_flux_spectrum in self.flux_g:
            f = 0.0
            for f_g in gth_flux_spectrum:
                f = f + f_g
            self.flux.append(f)
        # Cinder Energy Groups
        for gth_flux_spectrum in self.flux_g_cinder:
            f = 0.0
            for f_g in gth_flux_spectrum:
                f = f + f_g
            self.flux_cinder.append(f)

        # Renormalize Fluxes if burned...
        if 0 < len( self.burn["step"] ):
            # User-specified groups
            for t in range(len(self.flux_g)):
                for g in range(len(self.flux_g[t])):
                    self.flux_g[t][g] = self.flux_g[t][g] * self.burn["flux"][t] / self.flux[t]
            self.flux = self.burn["flux"]
            # Cinder Energy Groups
            for t in range(len(self.flux_g_cinder)):
                for g in range(len(self.flux_g_cinder[t])):
                    self.flux_g_cinder[t][g] = self.flux_g_cinder[t][g] * self.burn["flux"][t] / self.flux_cinder[t]
            self.flux_cinder = self.burn["flux"]

        # Calculate Meta-Stable to Ground State branch ratios, if metastables are suppossed to be tracked...
        for iso in metastabletrak:
            anum = (iso%10000)/10
            znum = iso/10000
            NG_grd_str  = "{0}{1:03d}0 produced by the following C-X  (g   )".format(anum, znum)
            NG_mss_str  = "{0}{1:03d}1 produced by the following C-X  (g  *)".format(anum, znum)
            N2N_grd_str = "{0}{1:03d}0 produced by the following C-X  (2n  )".format(anum, znum)
            N2N_mss_str = "{0}{1:03d}1 produced by the following C-X  (2n *)".format(anum, znum)

            InGrndNG  = False
            InMetaNG  = False
            InGrndN2N = False
            InMetaN2N = False

            GrndXS_NG  = []
            MetaXS_NG  = []
            GrndXS_N2N = []
            MetaXS_N2N = []

            cinderfile = open(CINDER_DAT, 'r')
            for line in cinderfile:
                if len(line) <= 4:
                    continue
                elif NG_grd_str in line:
                    InGrndNG = True
                    InMetaNG = False
                    continue
                elif NG_mss_str in line:
                    InGrndNG = False
                    InMetaNG = True
                    continue
                elif N2N_grd_str in line:
                    InGrndN2N = True
                    InMetaN2N = False
                    continue
                elif N2N_mss_str in line:
                    InGrndN2N = False
                    InMetaN2N = True
                    continue
                elif "_______________________" in line:
                    InGrndNG  = False
                    InMetaNG  = False
                    InGrndN2N = False
                    InMetaN2N = False
                    continue
                elif "   #" == line[:4]:
                    InGrndNG  = False
                    InMetaNG  = False
                    InGrndN2N = False
                    InMetaN2N = False
                    continue
                elif InGrndNG:
                    for xs in line.split():
                        GrndXS_NG.append( float(xs) )
                elif InMetaNG:
                    for xs in line.split():
                        MetaXS_NG.append( float(xs) )
                elif InGrndN2N:
                    for xs in line.split():
                        GrndXS_N2N.append( float(xs) )
                elif InMetaN2N:
                    for xs in line.split():
                        MetaXS_N2N.append( float(xs) )
            cinderfile.close()

            if not (GrndXS_NG == []) and not (MetaXS_NG == []):
                self.BR_NG[iso] = []
                for t in range( len(flux_g_cinder) ):
                    gXS_NG = msn.GroupCollapse(GrndXS_NG, flux_g_cinder[t], flux_cinder[t])
                    mXS_NG = msn.GroupCollapse(MetaXS_NG, flux_g_cinder[t], flux_cinder[t])
                    self.BR_NG[iso].append(mXS_NG / gXS_NG)

            if not (GrndXS_N2N == []) and not (MetaXS_N2N == []):
                self.BR_N2N[iso] = []
                for t in range( len(flux_g_cinder) ):
                    gXS_N2N = msn.GroupCollapse(GrndXS_N2N, flux_g_cinder[t], flux_cinder[t])
                    mXS_N2N = msn.GroupCollapse(MetaXS_N2N, flux_g_cinder[t], flux_cinder[t])
                    self.BR_N2N[iso].append(mXS_N2N / gXS_N2N)
        return

    def write_text_lib(self):
        """Writes MCNP output to text libraries."""
    
        if self.mult == {}:
            self.parse()

        if not ( 'libs' in os.listdir('.') ):
            os.mkdir('libs/')
        os.chdir('libs/')
        for lib in os.listdir('.'):
            if lib[-4:] == ".lib":
                metasci.SafeRemove(lib)

        headline = "E_up"
        for t in CoarseTime:
            headline = headline + "\t\t" + str(t)
        headline = headline + "\n"

        ### Write Flux Files
        # Group Fluxes
        libfile = open('flux_g.lib', 'w')
        libfile.write(headline)
        for e in range(len(self.E_up)):
            line = '%.5E'%E_up[e]
            for f in range(len(self.flux_g)):
                line = line + "\t{0:.6E}".format(self.flux_g[f][e])
            libfile.write('{0}\n'.format(line))
        libfile.close()	 

        # Total Flux
        libfile = open('flux.lib', 'w')
        libfile.write('Time\t' + headline[5:])
        line = "Flux"
        for f in range(len(self.flux)):
            line = line + "\t{0:.6E}".format(self.flux[f])
        libfile.write('{0}\n'.format(line))
        libfile.close()	 

        # Cinder Group Fluxes
        libfile = open('flux_g_cinder.lib', 'w')
        libfile.write(headline)
        for e in range(len(self.E_up_cinder)):
            line = '{0:.5E}'.format(self.E_up_cinder[e])
            for f in range(len(self.flux_g_cinder)):
                line = line + "\t{0:.6E}".format(self.flux_g_cinder[f][e])
            libfile.write('{0}\n'.format(line))
        libfile.close()	 

        # Cinder Total Flux
        libfile = open('flux_cinder.lib', 'w')
        libfile.write('Time\t' + headline[5:])
        line = "Flux"
        for f in range(len(self.flux_cinder)):
            line = line + "\t{0:.6E}".format(self.flux_cinder[f])
        libfile.write('{0}\n'.format(line))
        libfile.close()	 

        # Write Multiplier Files
        for lib in mult.keys():
            libfile = open(lib, 'w')
            libfile.write(headline)
            for e in range(len(self.E_up)):
                line = '{0:.5E)'.format(self.E_up[e])
                for ml in range(len(self.mult[lib])):
                    line = line + "\t{0:.6E}".format(self.mult[lib][ml][e])
                libfile.write('{0}\n'.format(line))
            libfile.close()

        # Meta-Stable Branch Ratio Files
        for iso in BR_NG.keys():
            libfile = open(isoname.zzaaam_2_LLAAAM(iso) + '_NGamma_Branch_Ratio.lib', 'w')
            libfile.write('Time\t' + headline[5:])
            line = "Ratio"
            for br in range(len(self.BR_NG[iso])):
                line = line + "\t{0:.6E}".format(self.BR_NG[iso][br])
            libfile.write('{0}\n'.format(line))
            libfile.close()	 
        for iso in BR_N2N.keys():
            libfile = open(isoname.zzaaam_2_LLAAAM(iso) + '_N2N_Branch_Ratio.lib', 'w')
            libfile.write('Time\t' + headline[5:])
            line = "Ratio"
            for br in range(len(self.BR_N2N[iso])):
                line = line + "\t{0:.6E}".format(self.BR_N2N[iso][br])
            libfile.write('{0}\n'.format(line))
            libfile.close()	 

        os.chdir('..')
        return 

    def write_hdf5_lib(self):
        """Writes MCNP output to an HDF5 library."""
    
        if self.mult == {}:
            self.parse()

        if not ( 'libs' in os.listdir('.') ):
            os.mkdir('libs/')
        os.chdir('libs/')
        for h5 in os.listdir('.'):
            if h5[-3:] == ".h5":
                metasci.SafeRemove(h5)

        # Initialize the HDF5 file
        libfile = tb.openFile(reactor + ".h5", mode = "w", title = '[CHAR] {0}'.format(reactor))
        root = libfile.root

        ##########################################
        ### Make arrays that apply to all data ###
        ##########################################
        # Add Energy Group Array
        libfile.createArray(root, "E_up", E_up, "Upper Energy Limit [MeV]")                                     #Upper Energy Limit
        libfile.createArray(root, "CoreLoad_zzaaam", CoreLoad_zzaaam, "Core Loading Isotopes [zzaaam]")         #Core load isotopes in zzaaam form
        libfile.createArray(root, "CoreLoad_LLAAAM", CoreLoad_LLAAAM, "Core Loading Isotopes [LLAAAM]")         #Core load isotopes in LLAAAM form
        libfile.createArray(root, "CoreLoad_MCNP",   CoreLoad_MCNP,   "Core Loading Isotopes [MCNP]")           #Core load isotopes in MCNP form
        libfile.createArray(root, "CoreTran_zzaaam", CoreTran_zzaaam, "Core Transformation Isotopes [zzaaam]")	#Core transformation isotopes in zzaaam form
        libfile.createArray(root, "CoreTran_LLAAAM", CoreTran_LLAAAM, "Core Transformation Isotopes [LLAAAM]")	#Core transformation isotopes in LLAAAM form
        libfile.createArray(root, "CoreTran_MCNP",   CoreTran_MCNP,   "Core Transformation Isotopes [MCNP]")	#Core transformation isotopes in MCNP form

        # Add Coarsely binned time and flux data
        libfile.createGroup(root, "Coarse")
        libfile.createArray("/Coarse", "time",   CoarseTime,  "Burn Time, Coarse [days]")		        #Burn up times
        libfile.createArray("/Coarse", "flux",   self.flux,   "Neutron Flux, Coarse [n/s/cm2]")	 	    #Flux
        libfile.createArray("/Coarse", "flux_g", self.flux_g, "Neutron Group Flux, Coarse [n/s/cm2]")	#Group fluxes

        # Add Finely Binned time and flux data
        libfile.createGroup(root, "Fine")
        flux_fine = []
        flux_g_fine = []
        for t in FineTime:
            n = 0
            while (CoarseTime[n] <= t) and not (t <= CoarseTime[n+1]):
                n = n + 1
            flux_fine.append( metasci.SolveLine(t, CoarseTime[n+1], self.flux[n+1], CoarseTime[n], self.flux[n]) )
            flux_g_fine.append([])
            for e in range( len(self.E_up) ):
                flux_g_fine[-1].append(metasci.SolveLine(t, CoarseTime[n+1], self.flux_g[n+1][e], CoarseTime[n], self.flux_g[n][e]))
        libfile.createArray("/Fine", "time",   FineTime,    "Burn Time, Fine [days]")			    #Burn up times
        libfile.createArray("/Fine", "flux",   flux_fine,   "Neutron Flux, Fine [n/s/cm2]") 		#Flux
        libfile.createArray("/Fine", "flux_g", flux_g_fine, "Neutron Group Flux, Fine [n/s/cm2]")	#Group fluxes


        #############################################
        ### Make arrays that apply to Cinder data ###
        #############################################
        libfile.createGroup(root, "CINDER")
        libfile.createArray("/CINDER", "time",   CoarseTime,         "Cinder Burn Time [days]")			        #Burn up times
        libfile.createArray("/CINDER", "E_up",   self.E_up_cinder,   "Upper Energy Limit - Cinder [MeV]")       #Upper Energy Limit
        libfile.createArray("/CINDER", "flux_g", self.flux_g_cinder, "Neutron Group Flux - Cinder [n/s/cm2]")   #Group fluxes
        libfile.createArray("/CINDER", "flux",   self.flux_cinder,   "Neutron Flux - Cinder [n/s/cm2]")		    #Flux

        ##########################################
        ### Make Meta-Stable Branch Ratio Data ###
        ##########################################
        libfile.createGroup("/Coarse", "BranchRatio")
        libfile.createGroup("/Coarse/BranchRatio", "NG")
        libfile.createGroup("/Coarse/BranchRatio", "N2N")
        libfile.createGroup("/Fine", "BranchRatio")
        libfile.createGroup("/Fine/BranchRatio", "NG")
        libfile.createGroup("/Fine/BranchRatio", "N2N")

        # Coarse Data
        for iso in BR_NG.keys():
            iname = isoname.zzaaam_2_LLAAAM(iso)
            libfile.createArray("/Coarse/BranchRatio/NG", iname, self.BR_NG[iso], "(n, gamma) Meta-Stable to Ground {0} Branch Ratio".format(iname))	#Meta-Stable Branch Ratio Data
        for iso in BR_N2N.keys():
            iname = isoname.zzaaam_2_LLAAAM(iso)
            libfile.createArray("/Coarse/BranchRatio/N2N", iname, self.BR_N2N[iso], "(n, 2n) Meta-Stable to Ground {0} Branch Ratio".format(iname))	#Meta-Stable Branch Ratio Data

        # Fine Data
        BR_NG_Fine  = {}
        for iso in BR_NG.keys():
            finetemp = []
            for t in FineTime:
                n = 0
                while (CoarseTime[n] <= t) and not (t <= CoarseTime[n+1]):
                    n = n + 1
                    finetemp.append( metasci.SolveLine(t, CoarseTime[n+1], self.BR_NG[iso][n+1], CoarseTime[n], self.BR_NG[iso][n]) )
            iname = isoname.zzaaam_2_LLAAAM(iso)
            libfile.createArray("/Fine/BranchRatio/NG", iname, finetemp, "(n, gamma) Meta-Stable to Ground {0} Branch Ratio".format(iname))	#Meta-Stable Branch Ratio Data
        BR_N2N_Fine = {}
        for iso in BR_N2N.keys():
            finetemp = []
            for t in FineTime:
                n = 0
                while (CoarseTime[n] <= t) and not (t <= CoarseTime[n+1]):
                    n = n + 1
                    finetemp.append( metasci.SolveLine(t, CoarseTime[n+1], self.BR_N2N[iso][n+1], CoarseTime[n], self.BR_N2N[iso][n]) )
            iname = isoname.zzaaam_2_LLAAAM(iso)
            libfile.createArray("/Fine/BranchRatio/N2N", iname, finetemp, "(n, 2n) Meta-Stable to Ground {0} Branch Ratio".format(iname))	#Meta-Stable Branch Ratio Data

        ##############################
        ### Make Multiplier Arrays ###
        ##############################
        for key in FMdic.keys():
            libfile.createGroup("/Coarse", key)
        for key in FMdic.keys():
            libfile.createGroup("/Fine", key)

        # Coarse Data
        for key in mult.keys():
            iso, nada, mult_type = key[:-4].partition("_")
            if mult_type[:5] == "sigma":
                mult_str = mult_type + " [barns]"
            else:
                mult_str = mult_type + " [unitless]"
            libfile.createArray("/Coarse/"+mult_type, iso, self.mult[key], mult_str)

            # Fine Data
            finetemp = []
            for t in FineTime:
                n = 0
                while (CoarseTime[n] <= t) and not (t <= CoarseTime[n+1]):
                     n = n + 1
                    finetemp.append([])
                    for e in range( len(E_up) ):
                        finetemp[-1].append( metasci.SolveLine(t, CoarseTime[n+1], self.mult[key][n+1][e], CoarseTime[n], self.mult[key][n][e]) )
            libfile.createArray("/Fine/"+mult_type, iso, finetemp, mult_str)
    
        libfile.close()
        os.chdir('..')
        return

