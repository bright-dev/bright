class_ds = \
"""A One-Group Light Water Reactor Fuel Cycle Component.  This is a daughter 
class of Reactor1G and a granddaughter of FCComp.

Parameters
----------
lib : str, optional
    The path the the LWR HDF5 data library.  This value is set to 
    Reactor1G.libfile and used by Reactor1G.loadlib().
rp : ReactorParameters, optional
    The physical reactor parameter data to initialize this LWR instance with.  
    If this argument is not provided, default values are taken.
n : str, optional
    The name of this LWR instance.

"""

desc = {
    'docstrings': {
        'class': class_ds,
        'attrs': {},
        'methods': {},
        },
    'attrs': {},
    'extra': {},
    }

mod = {'LightWaterReactor1G': desc,
       'docstring': "Python wrapper for LWR1G.",}


desc['docstrings']['methods']['calc_params'] = \
"""Along with its own parameter set to track, the LWR model implements its own 
function to set these parameters.  This function is equivalent to the following::

    self.params_prior_calc["BUd"]  = 0.0
    self.params_after_calc["BUd"] = self.BUd

    self.params_prior_calc["U"]  = self.mat_feed_u.mass
    self.params_after_calc["U"] = self.mat_prod_u.mass

    self.params_prior_calc["TRU"]  = self.mat_feed_tru.mass
    self.params_after_calc["TRU"] = self.mat_prod_tru.mass

    self.params_prior_calc["ACT"]  = self.mat_feed_act.mass
    self.params_after_calc["ACT"] = self.mat_prod_act.mass

    self.params_prior_calc["LAN"]  = self.mat_feed_lan.mass
    self.params_after_calc["LAN"] = self.mat_prod_lan.mass

    self.params_prior_calc["FP"]  = 1.0 - self.mat_feed_act.mass  - self.mat_feed_lan.mass

"""
