#############################
### General specifcations ###
#############################
reactor = "lwr"
coreloadtrackfile = "loadtrack.txt"				#Initial core loading nuclide list file
coretrantrackfile = "transformtrack.txt"		#Transmutation tracking nuclide list file
BurnTime   = 4200								#Number of days to burn the material [days]	
CoarseStep = 840								#Coarse Time step by which to increment the burn, for MCNP [days]
FineStep   = 30         						#Shorter time step for ORIGEN runs [days]
#CoarseStep = 4200								#Coarse Time step by which to increment the burn, for MCNP [days]
#FineStep   = 150								#Shorter time step for ORIGEN runs [days]
email      = "scopatz@gmail.com"		    	#E-mail address to send job information to.

NumberCPUs  = 4									#Number of CPUs to run transport code on.
CPUsPerNode = 4                                 #Processors per node

verbosity = 100
verbosity = 0

################################
### Unit Cell Sepcifications ###
################################
FuelCellRadius = 0.410
VoidCellRadius = 0.4185
CladCellRadius = 0.475
UnitCellPitch  = 0.65635 * 2.0 
UnitCellHeight = 10.0

FuelDensity = 10.7					#Denisty of Fuel
CladDensity = 5.87                  #Cladding Density
CoolDensity = 0.73                  #Coolant Density

###########################
### MCNPX Specification ###
###########################
moffset = 3 						#Number of materials used by default template.  
InitialFuelForm = {                 #Dictionary of initial fuel loading in weight fraction. 
    80160: 2.0 * 16.0 / 270.0, 
    922380: 0.96*238.0/270.0, 
    922350: 0.04*238.0/270.0
    }	
FuelSpecificPower = 40.0 / 1000.0   #Power garnered from fuel [MW / kgIHM]

kParticles  = 5000      #Number of particles to run per kcode cycle
#kParticles  = 1000      #Number of particles to run per kcode cycle
#kParticles  = 100       #Number of particles to run per kcode cycle
kCycles     = 130       #Number of kcode cycles to run
kCyclesSkip = 30        #Number of kcode cycles to run but not tally at the begining.

CINDER_DAT = "/usr/share/MCNPX/v260/Data/cinder.dat" 	#path to cinder.dat file, needed for metastables...

#GroupStructure = "1.0e-9 98log 10.0"						#Any Valid MCNP Group structure.  Remember, these are upper energy bounds.
GroupStructure = "1.0e-9 4log 10.0"						#Any Valid MCNP Group structure.  Remember, these are upper energy bounds.

FMdic = {"sigma_t": -1}
#FMdic = {"sigma_t": -1, "sigma_f": -2, "nubar": -3, "chi": -4, "sigma_a": -5}
#FMdic = {"sigma_t": -1, "sigma_f": -2, "nubar": -3, "chi": -4, "sigma_a": -5, "sigma_e": 2, "sigma_i": 4}
#FMdic = {"sigma_t": -1, "sigma_f": -2, "nubar": -3, "chi": -4, "sigma_a": -5, "sigma_e": 2, "sigma_i": 4, "sigma_2n": 16, "sigma_3n": 17, "sigma_gamma": 102, "sigma_proton": 103, "sigma_alpha": 107}

#ISO_FLAG = ".60c"
#ISO_FLAG = ".70c"
ISO_FLAG = "" 

############################
### ORIGEN Specification ###
############################
ORIGEN_FASTorTHERM = "THERM"							#Specifies whether the reactor is thermal or fast, and thus what origen executable to use...
ORIGEN_Concentration_Cut_Off = 10.0**-10                #Default is 1E-10

##################################
### Local Server Specification ###
##################################
LocalPathMPI  = "mpirun"							#Local Path to 'mpirun'
LocalPathMCNP = "mcnpx"								#Local Path to MCNP

###################################
### Remote Server Specification ###
###################################
RemoteURL  = "nukestar.me.utexas.edu"						#Remote server address
RemoteUser = "scopatz"								#Remoter username
RemoteDir  = "/home/scopatz/"							#Remote directory to store files.

RemoteGateway = 'nukestar01'

RemotePathMPI  = "/usr/local/bin/mpiexec"					#Remote Path to 'mpirun'
RemotePathMCNP = "/usr/share/mcnpxv260/bin/mcnpx260"				#Remote Path to mcnp
RemoteDATAPATH = "/usr/share/mcnpxv260/lib/"					#Remote DATAPATH enviromental variable
