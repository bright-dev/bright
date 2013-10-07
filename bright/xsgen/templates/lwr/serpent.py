xs_gen = """\
set title "[CHAR] {reactor} Cross Section Generator"

set acelib "{xsdata}"

% --- Matrial Definitions ---

% Initial Fuel Stream
mat fuel -{fuel_density}
{fuel}

% Cladding Stream
mat cladding   -{clad_density}
{cladding}

% Coolant Stream
mat coolant  -{cool_density} moder lwtr 1001
{coolant}

therm lwtr lwj3.20t

% --- Run Specification ---

% Periodic boundary conditions
set bc 3

% Fuel universe
set gcu 100

% 1/8 square symmetry
{sym_flag}set sym 8

% Group Stucture
set egrid 5E-05 {group_lower_bound} {group_upper_bound}
set nfg {n_groups} 
{group_inner_structure}


% Criticality calc
set pop {k_particles} {k_cycles} {k_cycles_skip}

% --- Geometry ---
pin 1
fill 100 {fuel_radius}
void     {void_radius}
cladding {clad_radius}
coolant  

pin 2
coolant

surf 100 inf
cell 110 100 fuel   -100

lat 10 1 0.0 0.0 {lattice_xy} {lattice_xy} {cell_pitch}

{lattice}

surf 3000  sqc  0.0 0.0 {half_lattice_pitch}
cell 300   0  fill 10   -3000
cell 301   0  outside    3000

% --- Graphs ---
%plot 3 800 800
%mesh 3 800 800

% --- Group Constant Generation ---

% Energy group structure
ene energies 1
{group_structure}

% Total flux in {detector_mat}
det phi de energies dm {detector_mat}

% Group constant material
mat xsmat 1.0 {xsnuc} 1.0

% Set group transfer probability to this material
set gtpmat xsmat

% Specify the detectors
{xsdet}
"""

burnup = """\
set title "[CHAR] {reactor} Burnup Calculation"

set acelib "{xsdata}"

% --- Matrial Definitions ---

% Initial Fuel Stream
mat fuel -{fuel_density} burn {num_burn_regions}
{fuel}

% Cladding Stream
mat cladding   -{clad_density}
{cladding}

% Coolant Stream
mat coolant  -{cool_density} moder lwtr 1001
{coolant}

therm lwtr lwj3.20t

% --- Run Specification ---

% Periodic boundary conditions
set bc 3

% 1/8 square symmetry
{sym_flag}set sym 8

% Group Stucture
set egrid 5E-05 {group_lower_bound} {group_upper_bound}
set nfg {n_groups}
{group_inner_structure}


% Criticality calc
set pop {k_particles} {k_cycles} {k_cycles_skip}


% --- Geometry ---
pin 1
fuel     {fuel_radius}
void     {void_radius}
cladding {clad_radius}
coolant  

pin 2
coolant

lat 10 1 0.0 0.0 {lattice_xy} {lattice_xy} {cell_pitch}

{lattice}

surf 3000  sqc  0.0 0.0 {half_lattice_pitch}
cell 300   0  fill 10   -3000
cell 301   0  outside    3000

% --- Graphs ---
%plot 3 800 800
%mesh 3 800 800

% Decay and fission yield libraries
 
set declib "{decay_lib}"
set nfylib "{fission_yield_lib}"
 
% Burnup calculation options
 
set bumode  2  % CRAM method
set pcc     1  % Predictor-corrector calculation on
set xscalc  2  % Calc cross sections from spectrum (fast)

set powdens {fuel_specific_power} % Fuel specific power [W/g]

% Depletion cycle
dep daytot

{depletion_times}

% Nuclide inventory
 
set inventory

{transmute_inventory}
"""


