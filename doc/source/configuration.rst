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


.. _isotope-tracking:

----------------
Isotope Tracking
----------------
* **core_load_isos** (list or path):  This parameter specifies all of the isotopes that may 
  be present in the initial core loading of the reactor.  As such, it is these isotopes that 
  char worries about transmuting from.
* **core_transmute_isos** (list or path):  This parameter gives all of the nulcides that char
  tracks throughout a burnup-depeletion calcultion.  As such, it is this set that determines 
  which isotopes are transmuted to.  Cross sections are calculated for all isotopes in this
  list.

There are a few ways to instantiate these two variables.  
The preffered way is to import some of the predefined lists that have already been specified
in char.  This prevents the user from having to build undo specification::

    # Set isotopes to track
    from char.iso_track import load, transmute
    core_load_isos      = load
    core_transmute_isos = transmute

Valid lists that may be imported are ``load``, ``transmute``, ``actinides``, ``uranium``, and ``u235``.
The contents of these lists may be seen at 
`the github website <https://github.com/scopatz/char/blob/master/char/iso_track.py>`_.

Moreover, the user can simply pass in python lists:: 

    # Set isotopes to track
    core_load_isos = ['U235', 922380, 'O16', 10010]
    core_transmute_isos = ['Am245', 'PU240', 'tc99', ...]

Additionally, suppose that you have file that is a whitespace separated list of isotope names.
A string-valued path to this file may also be passed into these parameters.

    ``isos_to_track.txt``::

        U235 
        Th232 80160
        AM242m

    Relevant part of ``defchar.py``::

        # Set isotopes to track
        core_load_isos = "/path/to/isos_to_track.txt"
        core_transmute_isos = "/path/to/isos_to_track.txt"


--------------------------
Calculation Mode Templates
--------------------------
* **burnup_template** (str): An unformatted string template that is used for burnup-depletion
  calculations.
* **xs_gen_template** (str): An unformatted string template that is used for cross-section generation
  calculations.

Like in :ref:`isotope-tracking`, the preffered method of supplying templates is from pre-defined
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

An example that is representative of a light water reactor is as follows::

    fuel_cell_radius = 0.410
    void_cell_radius = 0.4185
    clad_cell_radius = 0.475
    unit_cell_pitch  = 0.65635 * 2.0
    unit_cell_height = 10.0

    fuel_density = [10.7, 10.7*0.9, 10.7*1.1]   # Run 3 different fuel densities
    clad_density = 5.87
    cool_density = 0.73

