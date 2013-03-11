class_ds = \
"""Fuel Fabrication Fuel Cycle Component Class.  Daughter of FCComp class.

Parameters
----------
mats : dict or map or None, optional 
    A dictionary whose keys are string labels (eg, "U-235", "TRU", "My Fuel") 
    and whose values are mass streams.  For example::

        materials = {
            "U235": Material({922350: 1.0}, 1.0, "U-235"),
            "U236": Material({922360: 1.0}, 1.0, "U-236"),
            "U238": Material({922380: 1.0}, 1.0, "U-238"),
            }

    would be valid for a light water reactor.
mws_in : dict or map or None, optional 
    A dictionary whose keys are the same as for materials and whose values are the 
    associated weight (float) for that stream.  If a material should be allowed to 
    vary (ie optimized over), specify its weight as a negative number. For instance::

        mass_weights_in = {
            "U235": -1.0,
            "U236": 0.005,
            "U238": -1.0,
            }

    would be valid for a light water reactor with half a percent of U-236 always 
    present.
r : Reactor1G or None, optional 
    An instance of a Reactor1G class to fabricate fuel for.
paramtrack : list of str or None, optional 
    Additional parameters to track, if any.    
n : str, optional 
    The name of the fuel fabrication fuel cycle component instance.

"""

desc = {
    'docstrings': {
        'module': "Python wrapper for fuel fabrication.",
        'class': class_ds,
        'attrs': {},
        'methods': {},
        },
    'attrs': {},
    'extra': {},
    }

desc['docstrings']['attrs']['materials'] = \
"A mapping of materials which are mixed to create a valid fuel form for reactor."

desc['docstrings']['attrs']['mass_weights_in'] = \
"""A mapping representing the initial specification of the mass weights.  
Exactly two of these should have negative values."""

desc['docstrings']['attrs']['mass_weights_out'] = \
"""A mapping representing the mass weights that are calculated to generate 
a valid fuel from the materials  provided."""

desc['docstrings']['attrs']['deltaRs'] = \
"""A mapping representing the deltaR values of each of the materials."""

desc['docstrings']['attrs']['reactor'] = \
"""An instance of a Reactor1G class that the mat_prod is valid as a fuel for."""

desc['docstrings']['methods']['initialize'] = \
"""The initialize() method takes the appropriate materials, input mass weights,
and a reactor and resets the current FuelFabrication instance.

Parameters
----------
mats : dict
    A dictionary whose keys are string labels and whose values are mass streams.  
mws_in : dict 
    A dictionary whose keys are the same as for materials and whose values are the 
    associated weight (float) for that stream.
r : Reactor1G 
    An instance of a Reactor1G class to fabricate fuel for.

"""

desc['docstrings']['methods']['calc_params'] = \
"""Here the parameters for FuelFabrication are set.  For example::

    self.params_prior_calc["Weight_U235"] = self.mass_weights_in["U235"]
    self.params_after_calc["Weight_U235"] = self.mass_weights_out["U235"]

    self.params_prior_calc["deltaR_U235"] = self.deltaRs["U235"]
    self.params_after_calc["deltaR_U235"] = self.deltaRs["U235"]

    self.params_prior_calc["Weight_U238"] = self.mass_weights_in["U238"]
    self.params_after_calc["Weight_U238"] = self.mass_weights_out["U238"]

    self.params_prior_calc["deltaR_U238"] = self.deltaRs["U238"]
    self.params_after_calc["deltaR_U238"] = self.deltaRs["U238"]

"""

desc['docstrings']['methods']['calc_deltaRs'] = \
"""Computes deltaRs for each mass stream."""


desc['docstrings']['methods']['calc_core_input'] = \
"""Computes the core input material that becomes mat_prod based on materials and 
mass_weights_out.

Returns
-------
core_input : Material 
    mat_prod.
"""

desc['docstrings']['methods']['calc_mass_ratios'] = \
"""Calculates mass_weights_out by varying the values of the two parameter which had 
negative values in mass_weights_in.  Therefore, this is the portion of the code that 
performs the optimization calculation.
"""

desc['docstrings']['methods']['calc'] = \
"""This method performs an optimization calculation on all input materials to 
determine the mass ratios that generate the correct fuel form for the reactor.
It then compiles the fuel and returns the resultant material. 

Parameters
----------
mats : dict or None, optional
    A dictionary whose keys are string labels and whose values are materials.  
mws_in : dict or None, optional
    A dictionary whose keys are the same as for materials and whose values are 
    the associated weight (float) for that material.
r : Reactor1G or None, optional 
    An instance of a Reactor1G class to fabricate fuel for.

Returns
-------
core_input : Material 
    mat_prod

"""
