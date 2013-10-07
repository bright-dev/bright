==================
Char Configuration
==================
Char is highly configurable, which enables easy small and large changes to the reactor of interest.  
All of char's configuration parameters are stored in a file, which is by default called ``defchar.py``.
However, the path to this definition file is given as a runtime argument, so the file may have any name.

As its name implies, ``defchar.py`` is a python module that serves as a cookie jar for all of the
parameters used.  Because it is a python file and is executed at runtime, the underlying variables 
that char uses may be calculated in as complex or as simple a way as the user desires.  Here, 
only simple examples will be demonstrated.

Char also comes with a variety of default parameters.  In some cases, these defaults may be imported
directly into ``defchar.py``.  In other instances, defaults are selected automatically if a value is 
not provided.

What follows is an explanation of all of the input paramters that may be specified in a ``defchar.py``
file.  

.. note::
    Any parameter whose type informatation includes the optional ``[sequence of]`` flag may be 
    confugured using a scalar value of this type or by a list or tuple, whose elements are all
    of this same type.

    Parameters flagged in this way are included in the perturbation table.  The perturbation table
    itself is calculated at the outproduct of all of these parameters.  This is the fundemental 
    method by which char calculates different cross section sets in the same library.


---------------------
General Specification
---------------------
* **reactor** (str): A string name for the reactor.  This is used as an identified throughout
  much of char.  Only letters and underscores are valid characters.  Exmples include
  ``"lwr"``, ``"fast_reactor"``, and ``"LWR_more_perts"``.
* **burn_regions** ([sequence of] int): This specifies the number of annular burn regions in a fuel pin.
  The absolut minimum number is 1, however three to four are probably required for decent statistics.
  The serpent manual recomends ten. 
* **burn_time** (int or float): This is the total time, in days, that all burnup calculations should 
  be run for. As such, it is a surrogate for the fluence.  For a light water reactor, a value of 
  ``4200`` is a decent guess.  Generally, you want to overshoot rather than undershoot this value.
* **time_step** (int or float): Time step by which to increment ``burn_time``.  The code generally
  works best when ``(0 == burn_time % time_step)``.  It is usually better to guess smaller values
  for this parameter rather than larger ones.
* **burn_times** (sequence): This is a list, array, or other sequence of burnup times, in days. 
  If present, this overrides the ``burn_time`` and ``time_step`` variables.  Zero should 
  be the first entry and it should be strictly monotonicly increasing.
* **email** (str): An email address to send runtime updates to, "char@zeon.gov".
* **scheduler** (str): If present, this value specifices which scheduler to use.  Accepted values include 
  an empty string ``''`` to indicate no scheduler and ``"PBS"`` to use the torque scheduler.
* **number_cpus** (int): This is the number of CPUs to run the transport code on if parallel processing 
  is available.
* **cpus_per_node** (int): If being run on a cluster, this value indicates the number of processors 
  available per node.
* **verbosity** (int): Indicates how much information should be printed to stdout.  The default level 
  is zero. You may set this arbitrarily high.

For example, one may set up the general portion as follows::

    reactor = "lwr"         # Run identifier
    burn_regions = 1        # Number of burnup annular regions.
    burn_time = 4200        # Number of days to burn the material [days]    
    time_step = 840         # Time step by which to increment the burn [days]
    email = "char@zeon.gov" # E-mail address to send job information to.
    verbosity = 100         # Print absolutely everything

    # A commented out line
    #scheduler = "PBS"

    number_cpus   = 3   # Number of CPUs to run transport code on.
    cpus_per_node = 4   # Processors per node


.. _nuclide-tracking:

----------------
Nuclide Tracking
----------------
* **core_load_nucs** (list or path):  This parameter specifies all of the nuclides that may 
  be present in the initial core loading of the reactor.  As such, it is these nuclides that 
  Char worries about transmuting from.
* **core_transmute_nucs** (list or path):  This parameter gives all of the nulcides that char
  tracks throughout a burnup-depeletion calcultion.  As such, it is this set that determines 
  which nuclides are transmuted to.  Cross sections are calculated for all nuclides in this
  list.

There are a few ways to instantiate these two variables.  
The preffered way is to import some of the predefined lists that have already been specified
in char.  This prevents the user from having to build undo specification::

    # Set isotopes to track
    from char.nuc_track import load, transmute
    core_load_nucs      = load
    core_transmute_nucs = transmute

Valid lists that may be imported are ``load``, ``transmute``, ``actinides``, ``uranium``, 
and ``u235``.  The contents of these lists may be seen at 
`the github website <https://github.com/scopatz/char/blob/master/char/nuc_track.py>`_.

Moreover, the user can simply pass in python lists:: 

    # Set isotopes to track
    core_load_nucs = ['U235', 922380, 'O16', 10010]
    core_transmute_nucs = ['Am245', 'PU240', 'tc99', ...]

Additionally, suppose that you have file that is a whitespace separated list of isotope names.
A string-valued path to this file may also be passed into these parameters.

    ``nucs_to_track.txt``::

        U235 
        Th232 80160
        AM242m

    Relevant part of ``defchar.py``::

        # Set isotopes to track
        core_load_nucs = "/path/to/nucs_to_track.txt"
        core_transmute_nucs = "/path/to/nucs_to_track.txt"


.. _calc_mode_templates:

--------------------------
Calculation Mode Templates
--------------------------
* **burnup_template** (str): An unformatted string template that is used for burnup-depletion
  calculations.
* **xs_gen_template** (str): An unformatted string template that is used for cross-section generation
  calculations.

Like in :ref:`nuclide-tracking`, the preffered method of supplying templates is from pre-defined
versions in char itself.  Current values may be seen 
`at github <https://github.com/scopatz/char/blob/master/char/templates/lwr/serpent.py>`_::

    # Load stock template string from char
    from char.templates.lwr import serpent
    burnup_template = serpent.burnup
    xs_gen_template = serpent.xs_gen

Of course, the user could generate their own template strings and place them here.
If one wehere to do this, the avaiable fill variables are listed below.  All of the
values of these variables are strings, or like integers or floats, easily coerced 
to strings.  For more information please refer to the
`Python manual <http://docs.python.org/library/string.html#format-specification-mini-language>`_.

    **Burnup**:
        * ``reactor``
        * ``xsdata``
        * ``decay_lib``
        * ``fission_yield_lib``
        * ``fuel``
        * ``fuel_radius``
        * ``fuel_density``
        * ``fuel_specific_power``
        * ``num_burn_regions``
        * ``cladding``
        * ``clad_radius``
        * ``clad_density``
        * ``coolant``
        * ``cool_radius``
        * ``cool_density``
        * ``sym_flag``
        * ``n_groups``
        * ``group_lower_bound``
        * ``group_upper_bound``
        * ``group_inner_structure``
        * ``k_cycles``
        * ``k_cycles_skip``
        * ``k_particles``
        * ``lattice``
        * ``lattice_xy``
        * ``cell_pitch``
        * ``half_lattice_pitch``
        * ``depletion_times``
        * ``transmute_inventory``

    **Cross Section Generation**:
        * ``reactor``
        * ``xsdata``
        * ``fuel``
        * ``fuel_radius``
        * ``fuel_density``
        * ``fuel_specific_power``
        * ``num_burn_regions``
        * ``cladding``
        * ``clad_radius``
        * ``clad_density``
        * ``coolant``
        * ``cool_radius``
        * ``cool_density``
        * ``sym_flag``
        * ``n_groups``
        * ``group_structure``
        * ``group_lower_bound``
        * ``group_upper_bound``
        * ``group_inner_structure``
        * ``k_cycles``
        * ``k_cycles_skip``
        * ``k_particles``
        * ``lattice``
        * ``lattice_xy``
        * ``cell_pitch``
        * ``half_lattice_pitch``
        * ``xsiso``
        * ``xsdet``


------------------------
Unit Cell Sepcifications
------------------------
* **fuel_cell_radius** ([sequence of] float): Fuel cell radius [cm].
* **void_cell_radius** ([sequence of] float): Void cell radius [cm].  Must be greater than or 
  equal to the ``fuel_cell_radius``. 
* **clad_cell_radius** ([sequence of] float): Cladding cell radius [cm].  Must be greater than or 
  equal to the ``void_cell_radius``. 
* **unit_cell_pitch** ([sequence of] float): The length of the unit cell box [cm].
* **unit_cell_height** ([sequence of] float): The length of the z-direction of the lattice of 
  interest [cm].
* **fuel_density** ([sequence of] float): Denisty of fuel region [g/cm^3].  
* **clad_density** ([sequence of] float): Denisty of cladding region [g/cm^3].  
* **cool_density** ([sequence of] float): Denisty of coolant region [g/cm^3].  
* **fuel_specific_power** ([sequence of] float): Mass-normalized power from a unit of fuel [W/g].
  Required for burnup and sensitivity calculations.

An example that is representative of a light water reactor is as follows::

    fuel_cell_radius = 0.410
    void_cell_radius = 0.4185
    clad_cell_radius = 0.475
    unit_cell_pitch  = 0.65635 * 2.0
    unit_cell_height = 10.0

    fuel_density = [10.7, 10.7*0.9, 10.7*1.1]   # Run 3 different fuel densities
    clad_density = 5.87
    cool_density = 0.73

    fuel_specific_power = 40.0 / 1000.0


---------------------
Lattice Specification
---------------------
* **lattice** (str):  This is a string that represents the lattice that is used by serpent.
  While interally char does some analysis of this string (to set the symmetric lattice flag ``sym_flag``), 
  this string is passed directly into serpent.  When using the 
  :ref:`default serpent templates <calc_mode_templates>`, the number ``1`` represents a fuel pin, 
  while the number ``2`` represents a coolant pin.  These are material numbers defined in the templates.
  New rows must be separated by newline characters.
  Optional, if not provided, a 17x17 pressurized water reactor assembly is subsitutied.
* **lattice_xy** (int):  The number of rows and columns in a fuel assembly.  
  For instance, this number would be 17 for a 17x17 assembly or 9 for 9x9 assembly.
  Optional, if not provided, this value is given as 17 to match ``lattice``.

More information on how to set up lattices is available in the serpent manual.
The default values are as follows::

    lattice_xy = 17
    lattice    = ("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n"
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


-----------------------
Mass Stream Information
-----------------------
* **initial_heavy_metal** (dict): A ditionary that specifies the initial heavy metal 
  concentrations of each isotope in a pure fuel stream.  The keys of this isotope are 
  integers in ``zzaaam``-form.  The values are floats on the range ``[0, 1]``.
  The sum of all values here should equal 1.  Therefore, this represents 1 [kgIHM].
* **fuel_chemical_form** (dict): This is another python dictionary that gives the 
  chemical composition of the fuel region.  Keys are either isotopes in ``zzaaam``-form
  or the string ``"IHM"``, which is a placeholder for the ``initial_heavy_metal`` stream.
  Values are floats that repesent the number of atoms in this chemical form.
* **clad_form** (dict): A python dictionary that represents the cladding material. The keys of 
  this isotope are integers in ``zzaaam``-form.  The values are floats on the range ``[0, 1]``.
  The sum of all values here should equal 1.
  Optional, if not present in ``defchar.py`` a default zircaloy will be substituted.
* **cool_form** (dict): A python dictionary that represents the cladding material. The keys of 
  this isotope are integers in ``zzaaam``-form.  The values are floats on the range ``[0, 1]``.
  The sum of all values here should equal 1.
  Optional, if not present in ``defchar.py`` borated light water will be substituted.

The following values represent a light water reactor::

    # Low enriched uranium
    initial_heavy_metal = {
        922350: 0.04,
        922380: 0.96,
        }

    # Uranium oxide
    fuel_chemical_form = {
        80160: 2.0,
        "IHM": 1.0,
        }

    # Default zircaloy
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

    # Default borated light water
    MW = (2 * 1.0) + (1 * 16.0) + (0.199 * 550 * 10.0**-6 * 10.0) + (0.801 * 550 * 10.0**-6 * 11.0)
    cool_form = {
        10010: (2 * 1.0) / MW,
        80160: (1 * 16.0) / MW,
        50100: (0.199 * 550 * 10.0**-6 * 10.0) / MW,
        50110: (0.801 * 550 * 10.0**-6 * 11.0) / MW,
        }


Additionally, if for some reason the user does not wish to supply a mass-weighted fuel stream
and instead spcifies atom fractions, the following parameters may be set away from their default 
``True`` value.  If any of these are ``False``, atom fractions wil be used.  This is not
recomended in general.

* **fuel_form_mass_weighted** (bool): Fuel stream mass-weighted / atom-fraction flag.
* **clad_form_mass_weighted** (bool): Cladding stream mass-weighted / atom-fraction flag.
* **cool_form_mass_weighted** (bool): Coolant stream mass-weighted / atom-fraction flag.


-------------------------------------------
Initial Nuclide Concentration Perturbations
-------------------------------------------
There are two main (optional) ways to pertub nuclides.  The first is such that 
a pertubed isotopic vector shows up in the outer product perturbation set.

* **initial_{nuc}** ([sequence of] float):  The mass fraction value(s) of this isotope that 
  should be replaced in the ``initial_heavy_metal`` material.  The ``{nuc}`` term specifices 
  the name of the isotope to be pertubed in ``zzaaam``- or ``name``-form (i.e. 'U235' or 922350).

For example, take the following hypothetical mass stream::

    # Initial heavy metal mass fraction distribution
    initial_heavy_metal = {
        922340: 0.01,
        922350: 0.04,
        "U238": 0.95,
        }

    # Pertub some of these nuclides
    initial_U234 = [0.01, 0.015]
    initial_922350 = [0.02, 0.04, 0.06]

The above would produce 6 initial heavy metal streams, one for each (U234, U235) combination, and
generate burnups or cross sections for each of these mass streams.  The remaining isotopes 
(U238 here) would have their mass fractions altered to accomdated the perturbed material.


The second way that initial nuclides are perturbed is during the isotopic ssensitivity study 
(``-m``, calculation).  Every nuclide present in the ``initial_heavy_metal`` stream is pertubed by a set 
amounts.  These relative amounts are given via the following parameter.

* **sensitivity_mass_fractions** ([sequence of] float): The relative amount by which to 
  perturb each nuclide in an isotopic sensitivity study, ``-m``.

An example that will perturb each nuclide's initial amount by +/-10% from the initial 
value for every other pertubation is as follows::

    sensitivity_mass_fractions = [1.1, 0.9]


--------------------
Crtiticality Control
--------------------
* **k_cycles** (int): The total number of criticality cycles to run.
* **k_cycles_skip** (int): The number of cycles to run but not tally at the begining.
  This number must be strictly less than ``k_cycles``.
* **k_particles** (int): The number of source particles to run per cycle.

These are some decent values for the criticality calculation::

    k_particles   = 1000
    k_cycles      = 130
    k_cycles_skip = 30


-------------------------
Other Physical Parameters
-------------------------
* **group_structure** (sequence of floats): Defines the energy group structure but setting bounds 
  on each group.  Energies are given in [MeV].  For ``G`` energy groups, this list must be of length 
  ``G+1``.  Note that this list is not allowed to change.  
* **temperature** (int or float): This is the value of the temperature [K] of the fuel region.  This 
  parameter determines which cross-section set to use and is therefor not a continuous variable.
  Though char does not check, ``temperature`` should be a positive multiple of 300 K (ie 300, 600, 900, etc).
  If a value other than this is supplied, char will likely run, but use physical models everywhere instead 
  of serpent generated values.

Examples of the above parameters are as follows::

    # A log-spaced 10-group structure
    group_structure = [1.0E-09, 1.0E-08, 1.0E-07, 1.0E-06, 1.0E-05, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0]

    # Temperature at 600 K
    temperature = 600



---------------------------
Remote Server Specification
---------------------------
* **remote_url** (str): Remote server web- or ip-address. Optional.
* **remote_user** (str): Login name to the remote server. Optional.
* **remote_dir** (str): Directory on remote filesystem to store char files. Optional.
* **remote_gateway** (str): Local hostname of the remote machine that the user logs in to. Optional.
  This is especially useful when running char on a remote cluster where the gateway machine
  might not have the same name locally that it provides to the outside world. Optional.

Here is a Mr. Aznable might set up his remote configuration::

    remote_url  = "mycluster.zeon.gov"
    remote_user = "char"
    remote_dir  = "/home/char/"
    remote_gateway = 'mycluster01'


---------------------------------------
Serpent Data Library Path Specification
---------------------------------------
* **serpent_xsdata** (str):  Path to the serpent cross section library to use.
* **serpent_decay_lib** (str): Path to the serpent decay library to use.  Optional, only 
  required for burnup calculations.
* **serpent_fission_yield_lib** (str): Path to the serpent fission yeild library to use.  Optional, only
  required for burnup calculations.

Mr. Aznable's friend, Casval, might want to use a common set of libraries. These parameters are therefore
set as follows::

    # Use ENDF libraries.
    serpent_xsdata = "/usr/share/serpent/xsdata/endf7.xsdata"

    # The following two are only needed for burnup runs
    serpent_decay_lib = "/usr/share/serpent/xsdata/sss_endfb7.dec"
    serpent_fission_yield_lib = "/usr/share/serpent/xsdata/sss_endfb7.nfy"


========================
Bringing it all together
========================
Using the above specification information, a typical ``defchar.py`` might look like the 
following.  This configuration can be veiwed a a work in progress, since the user has
chosen to comment some values out, but still wants to keep them around::

    #############################
    ### General specifcations ###
    #############################
    reactor = "lwr"
    burn_regions = 1
    burn_time   = 4200
    #time_step = 840
    email      = "char@zeon.gov"

    #scheduler = "PBS"

    number_cpus   = 3   # Number of CPUs to run transport code on.
    cpus_per_node = 4   # Processors per node

    verbosity = 100

    # Set isotopes to track
    from char.iso_track import load, transmute
    core_load_isos      = load
    core_transmute_isos = transmute


    # Load stock template string from char
    # Having this allows users to specify other templates
    from char.templates.lwr import serpent
    xs_gen_template = serpent.xs_gen
    burnup_template = serpent.burnup


    ################################
    ### Unit Cell Sepcifications ###
    ################################
    fuel_cell_radius = 0.410
    void_cell_radius = 0.4185
    clad_cell_radius = 0.475
    unit_cell_pitch  = 0.65635 * 2.0 
    unit_cell_height = 10.0

    fuel_density = [10.7, 10.7*0.9, 10.7*1.1]   # Denisty of Fuel
    clad_density = 5.87                         # Cladding Density
    cool_density = 0.73                         # Coolant Density

    fuel_specific_power = 40.0 / 1000.0   # Power garnered from fuel [W / g]


    #################################
    ### Mass Stream Specification ###
    #################################
    # LEU
    initial_heavy_metal = {
        922350: 0.04, 
        922380: 0.96, 
        }

    initial_U235 = [0.02, 0.04, 0.06]

    # UOX
    fuel_chemical_form = {
        80160: 2.0, 
        "IHM": 1.0, 
        }	


    fuel_form_mass_weighted = True


    ###############################
    ### Criticality Information ###
    ###############################
    #k_particles  = 500
    #k_particles  = 5000
    k_particles  = 1000
    k_cycles      = 130
    k_cycles_skip = 30


    # A log-spaced 10-group structure
    group_structure = [1.0E-09, 1.0E-08, 1.0E-07, 1.0E-06, 1.0E-05, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0]

    # Temperature at 600 K
    temperature = 600


    ###################################
    ### Remote Server Specification ###
    ###################################
    remote_url  = "mycluster.zeon.gov"
    remote_user = "char"
    remote_dir  = "/home/char/"
    remote_gateway = 'mycluster01'


    #############################
    ### Serpent Specification ###
    #############################
    serpent_xsdata = "/usr/share/serpent/xsdata/endf7.xsdata"
    #serpent_xsdata = "/usr/share/serpent/xsdata/jeff311.xsdata"

    # The following two are only needed for burnup runs
    serpent_decay_lib = "/usr/share/serpent/xsdata/sss_endfb7.dec"
    serpent_fission_yield_lib = "/usr/share/serpent/xsdata/sss_endfb7.nfy"
    #serpent_decay_lib = "/usr/share/serpent/xsdata/sss_jeff311.dec"
    #serpent_fission_yield_lib = "/usr/share/serpent/xsdata/sss_jeff311.nfy"

