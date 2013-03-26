class_ds = \
"""Storage Fuel Cycle Component Class.  Daughter of FCComp class.

Parameters
----------
n : str, optional
    The name of the storage fuel cycle component instance.

"""

desc = {
    'docstrings': {
        'module': "Python wrapper for Storage.",
        'class': class_ds,
        'attrs': {},
        'methods': {},
        },
    'methods': {
        ('calc', ('t', 'f8', '0.0')): 'Material', 
        },
    }

mod = {'Storage': desc,
       'docstring': "Python wrapper for Storage.",}

desc['docstrings']['attrs']['decay_time'] = \
"""This the float (double) attribute that represents how long an input fuel 
mass should be stored for.  This time is represented in seconds, so be sure 
to convert to the proper units before using.  Consider using pyne.utils.to_sec() 
for second conversions.
"""

desc['docstrings']['methods']['calc_params'] = \
"""Here the parameters for Storage are set.  For storage, this amounts to just
a "Mass" parameter::

    self.params_prior_calc["Mass"]  = self.mat_feed.mass
    self.params_after_calc["Mass"] = self.mat_prod.mass

"""

desc['docstrings']['methods']['calc'] = \
"""As usual, calc sets up the Storage component's input stream and calculates 
the corresponding output Material.  Here, this amounts to calling bateman() 
for every nuclide in mat_feed, for each chain that ends with a nuclide in track_nucs.

Parameters
----------
input : dict or Material or None, optional 
    If input is present, it set as the component's mat_feed.  If input is a 
    isotopic dictionary (zzaaam keys, float values), this dictionary is first 
    converted into a Material before being set as mat_feed.
decay_time : float or None, optional 
    decay_time is set to the time value here prior to any other calculations.  
    This time has units of seconds.

Returns
-------
output : Material
    mat_prod

"""


