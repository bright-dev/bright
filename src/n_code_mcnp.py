"""An class to setup, run, and parse MCNP."""

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

    def make_input(NoBurnBool=False, NoPertBool=False):
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


