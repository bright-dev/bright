class_ds = \
"""A One-Group Fast Reactor Fuel Cycle Component.  This is a daughter class 
of Reactor1G and a granddaughter of FCComp.

Parameters
----------
lib : str, optional
    The path the the FR HDF5 data library.  This value is set to Reactor1G.libfile 
    and used by Reactor1G.loadlib().
rp : ReactorParameters, optional
    The physical reactor parameter data to initialize this FR instance with.  
    If this argument is not provided, default values are taken.
n : str, optional
    The name of this FR instance.

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

mod = {'FastReactor1G': desc,
       'docstring': "Python wrapper for FR1G.",}

desc['docstrings']['methods']['calc_params'] = \
"""Along with its own parameter set to track, the FR model implements its own 
function to set these parameters.  This function is equivalent to the following::

    self.params_prior_calc["BUd"]  = 0.0
    self.params_after_calc["BUd"] = self.BUd

    self.params_prior_calc["TRUCR"]  = 0.0
    self.params_after_calc["TRUCR"] = self.calc_tru_cr()

    self.params_prior_calc["P_NL"]  = 0.0
    self.params_after_calc["P_NL"] = self.P_NL

    self.params_prior_calc["U"]  = self.mat_feed_u.mass
    self.params_after_calc["U"] = self.mat_prod_u.mass

    self.params_prior_calc["TRU"]  = self.mat_feed_tru.mass
    self.params_after_calc["TRU"] = self.mat_prod_tru.mass

    self.params_prior_calc["ACT"]  = self.mat_feed_act.mass
    self.params_after_calc["ACT"] = self.mat_prod_act.mass

    self.params_prior_calc["LAN"]  = self.mat_feed_lan.mass
    self.params_after_calc["LAN"] = self.mat_prod_lan.mass

    self.params_prior_calc["FP"]  = 1.0 - self.mat_feed_act.mass  - self.mat_feed_lan.mass
    self.params_after_calc["FP"] = 1.0 - self.mat_prod_act.mass - self.mat_prod_lan.mass

"""


