template = """\
set title "[CHAR] {reactor}"

set acelib "{xsdata}"

% --- Matrial Definitions ---

% Initial Fuel Stream
mat fuel -{fuel_density} {num_burn_regions}
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
set nfg {n_groups} {group_structure}


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
plot 3 800 800
mesh 3 800 800

% --- Burnup ---

{burnup}
"""
